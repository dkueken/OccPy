import numpy as np
import pandas as pd
from scipy import interpolate
import rasterio
from rasterio.fill import fillnodata

import laspy
from occpy.TerrainModel import TerrainModel

# import plotting functions
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from matplotlib.transforms import Affine2D
from matplotlib.collections import PathCollection
import matplotlib.lines as mlines
import seaborn as sns

import time
import os
import glob

# we need to add the path to the src folder so raytr module can be found
import sys

sys.path.append(r"./src")  # TODO: check if src folder should be better placed into occpy package
from occpy.src.raytr import PyRaytracer
from occpy.TerrainModel import TerrainModel

is_sorted = lambda a: np.all(a[:-1] <= a[1:])

def last_nonzero(arr, axis, invalid_val=-1):
    """_summary_

    Args:
        arr (_type_): _description_
        axis (_type_): _description_
        invalid_val (int, optional): _description_. Defaults to -1.

    Returns:
        _type_: _description_
    """
    mask = arr!=0
    val = arr.shape[axis] - np.flip(mask, axis=axis).argmax(axis=axis) - 1
    return np.where(mask.any(axis=axis), val, invalid_val)

def lineplot_plusplus(orientation = "horizontal", **kwargs):
    line = sns.lineplot(**kwargs)

    r = Affine2D().scale(sx=1, sy=-1).rotate_deg(90)
    for x in line.images + line.lines + line.collections:
        trans = x.get_transform()
        x.set_transform(r+trans)
        if isinstance(x, PathCollection):
            transoff = x.get_offset_transform()
            x._transOffset = r+transoff

    old = line.axis()
    line.axis(old[2:4] + old[0:2])
    xlabel = line.get_xlabel()
    line.set_xlabel(line.get_ylabel())
    line.set_ylabel(xlabel)

    return line


class OccPy:
    def __init__(self, laz_in, out_dir, vox_dim=0.1, lower_threshold=1, points_per_iter=10000000, plot_dim=None):
        self.laz_in_f = laz_in
        self.out_dir = out_dir
        os.makedirs(out_dir, exist_ok=True)
        self.vox_dim = vox_dim  # voxel dimension (cubic) in meters TODO: maybe implement non-cubic voxels?
        self.lower_threshold = lower_threshold  # lower threshold above ground to exclude TODO: check if necessary
        self.points_per_iter = points_per_iter  # number of points read in from laz file in each iteration
        self.is_mobile = False
        self.traj_f = None
        self.traj = None
        self.senspos_f = None
        self.senspos = None
        self.single_return = None
        # TODO: find a better way to link scan position id between laz file and scan pos file!
        self.scan_pos_id_stridx = 0

        # some parameters that will be filled during function calls
        self.TotalVolume = 0
        self.Volume0_3 = 0
        self.Volume3_10 = 0
        self.Volume10_max = 0
        self.TotalOcclusion = 0
        self.Occlusion0_3 = 0
        self.Occlusion3_10 = 0
        self.Occlusion10_max = 0

        self.OcclFrac2D = None


        # 3DGrids
        self.Nhit = None
        self.Nmiss = None
        self.Nocc = None
        self.Classification = None

        # 2D Grids
        self.dtm = None
        self.dsm = None
        self.chm = None

        if plot_dim is None:
            # TODO: Test if this works!
            with laspy.open(laz_in) as file:
                hdr = file.header
                self.PlotDim = dict(minX=hdr.x_min,
                                    maxX=hdr.x_max,
                                    minY=hdr.y_min,
                                    maxY=hdr.y_max,
                                    minZ=hdr.z_min,
                                    maxZ=hdr.z_max)

        else:
            # we expect the format of plot_dim to be [minX, minY, minZ, maxX, maxY, maxZ]
            self.PlotDim = dict(minX=plot_dim[0],
                                maxX=plot_dim[3],
                                minY=plot_dim[1],
                                maxY=plot_dim[4],
                                minZ=plot_dim[2],
                                maxZ=plot_dim[5])

        self.grid_dim = dict(nx=int((self.PlotDim['maxX'] - self.PlotDim['minX']) / self.vox_dim),
                             ny=int((self.PlotDim['maxY'] - self.PlotDim['minY']) / self.vox_dim),
                             nz=int((self.PlotDim['maxZ'] - self.PlotDim['minZ']) / self.vox_dim))

        # initialize RayTr Object
        self.RayTr = PyRaytracer()

        # Define Grid
        minBound = np.array([self.PlotDim['minY'], self.PlotDim['minX'], self.PlotDim['minZ']])
        maxBound = np.array([self.PlotDim['maxY'], self.PlotDim['maxX'], self.PlotDim['maxZ']])
        self.RayTr.defineGrid(minBound, maxBound, self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz'],
                              self.vox_dim)

    def read_trajectory_file(self, path2traj, delimiter=" ", hdr_time='%time', hdr_x='x', hdr_y='y', hdr_z='z'):
        self.is_mobile = True
        traj_in = pd.read_csv(path2traj, sep=delimiter)

        d = {'time': traj_in[hdr_time], 'sensor_x': traj_in[hdr_x], 'sensor_y': traj_in[hdr_y],
             'sensor_z': traj_in[hdr_z]}

        self.traj = pd.DataFrame(data=d)

    def read_sensorpos_file(self, path2senspos, delimiter=" ", hdr_scanpos_id='', hdr_x='', hdr_y='', hdr_z='', sens_pos_id_offset=0, str_idx_ScanPosID=0, single_return=False):
        self.is_mobile = False
        sens_pos_in = pd.read_csv(path2senspos, sep=delimiter)

        d = {'ScanPos': sens_pos_in[hdr_scanpos_id]+sens_pos_id_offset,
             'sensor_x': sens_pos_in[hdr_x], 'sensor_y': sens_pos_in[hdr_y], 'sensor_z': sens_pos_in[hdr_z]}

        self.senspos = pd.DataFrame(data=d)
        self.scan_pos_id_stridx = str_idx_ScanPosID
        self.single_return = single_return

    def interpolate_traj(self, traj_time, traj_x, traj_y, traj_z, pts_gpstime):
        f_x = interpolate.interp1d(traj_time, traj_x, kind='linear', fill_value="extrapolate")
        sensor_x = f_x(pts_gpstime)
        f_y = interpolate.interp1d(traj_time, traj_y, kind='linear', fill_value="extrapolate")
        sensor_y = f_y(pts_gpstime)
        f_z = interpolate.interp1d(traj_time, traj_z, kind='linear', fill_value="extrapolate")
        sensor_z = f_z(pts_gpstime)

        d = {'time': pts_gpstime, 'sensor_x': sensor_x, 'sensor_y': sensor_y, 'sensor_z': sensor_z}

        df = pd.DataFrame(data=d)

        return df

    def do_raytracing(self):
        run_raytraycing_after_loading = False
        if os.path.isdir(self.laz_in_f):
            ## get list of laz files in input directory
            fCont = glob.glob(f"{self.laz_in_f}\*.laz")

            for f in fCont:
                # get scan position
                lastdir_ind = f.rfind("\\")
                scan_name = f[lastdir_ind + 1:]
                # TODO: This is very specific to the test data and needs to be made generic!
                end_of_scanID_idx = scan_name.rfind("_")
                scan_id = int(scan_name[self.scan_pos_id_stridx:end_of_scanID_idx])

                print(f"###############################")
                print(f"##### Processing {scan_name}...")
                print(f"###############################")

                scanpos_X = self.senspos.loc[self.senspos['ScanPos'] == scan_id, 'sensor_x'].values[0]
                scanpos_Y = self.senspos.loc[self.senspos['ScanPos'] == scan_id, 'sensor_y'].values[0]
                scanpos_Z = self.senspos.loc[self.senspos['ScanPos'] == scan_id, 'sensor_z'].values[0]

                # read in laz file
                tic = time.time()
                with laspy.open(f) as file:
                    count = 0
                    for points in file.chunk_iterator(self.points_per_iter):

                        if self.single_return:
                            sorted = True

                            x = points.x.copy()
                            y = points.y.copy()
                            z = points.z.copy()

                            sensor_x = np.ones(x.shape) * scanpos_X
                            sensor_y = np.ones(x.shape) * scanpos_Y
                            sensor_z = np.ones(x.shape) * scanpos_Z

                            gps_time = np.linspace(start=count + 1, stop=count + len(x), num=len(x), endpoint=True)

                            # run raytracing algorithme using singleReturnPulses version
                            #print("Do raytracing with all pulses in batch")
                            tic_r = time.time()
                            self.RayTr.doRaytracing_singleReturnPulses(x, y, z, sensor_x, sensor_y,
                                                                  sensor_z, gps_time)
                            toc_r = time.time()
                            #print("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))

                        else:
                            x = points.x.copy()
                            y = points.y.copy()
                            z = points.z.copy()
                            gps_time = points.gps_time.copy()
                            return_number = points.return_number.copy()
                            number_of_returns = points.number_of_returns.copy()

                            # check if gps_time is sorted
                            if count == 0:  # only check sort state in the first iteration of the for loop
                                if not is_sorted(gps_time):
                                    print(
                                        f"!!!!! input laz file is not sorted along gps_time. The algorithm will still run. However, the "
                                        f"performance will be greatly decreased as the entire content of the laz file has to be read into "
                                        f"the system memory. If you have multi return data, consider sorting your laz data first, e.g. using "
                                        f"LASTools lassort: lassort -i laz_in -gps_time -return_number -odix _sort -olaz -v !!!!")
                                    sorted = False
                                else:
                                    sorted = True

                            sensor_x = np.ones(gps_time.shape) * scanpos_X
                            sensor_y = np.ones(gps_time.shape) * scanpos_Y
                            sensor_z = np.ones(gps_time.shape) * scanpos_Z

                            self.RayTr.addPointData(x, y, z, sensor_x, sensor_y, sensor_z, gps_time, return_number,
                                               number_of_returns)

                            if sorted:  # only if pulses are sorted run raytracing now. Otherwise we have to read in the entire dataset first!
                                # Get report on pulse dataset - comment this out once everythin is working or TODO: add a verbose flag!
                                #self.RayTr.getPulseDatasetReport()

                                # run raytracing on added points
                                #print("Do raytracing with stored pulses")
                                tic_r = time.time()
                                self.RayTr.doRaytracing()
                                toc_r = time.time()
                                #print("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))

                                self.RayTr.clearPulseDataset()

                                # Check if traversed pulses have been deleted from map - comment this out once everything is working or TODO: add a verbose flag!
                                # RayTr.getPulseDatasetReport()

                        count = count + len(gps_time)

            toc = time.time()
            if sorted and not self.single_return:
                # optional: incomplete pulses can occur if the data has been filtered (either actively or during black box processing
                # of the processing software. We could actively turn the incomplete pulses into complete ones and do the raytracing
                # for them!
                print("convert incomplete pulses to complete ones - be cautious with that!")
                # RayTr.getPulseDatasetReport()
                self.RayTr.cleanUpPulseDataset()
                # RayTr.getPulseDatasetReport()
                print("Run raytracing for incomplete pulses")
                tic_r = time.time()
                self.RayTr.doRaytracing()
                toc_r = time.time()
                print("Time elapsed for raytracing incomplete pulses: {:.2f} seconds".format(toc_r - tic_r))
                print("Time elapsed for reading and raytracing entire data: {:.2f} seconds".format(toc_r - tic))
            elif not sorted and not self.single_return:
                print("Time elapsed for reading in data: {:.2f} seconds".format(toc - tic))

                # RayTr.getPulseDatasetReport()

                print("Clean up pulse dataset in order to handle incomplete pulses")
                self.RayTr.cleanUpPulseDataset()

                self.RayTr.getPulseDatasetReport()

                print("Do actual raytracing with all pulses")
                tic = time.time()
                self.RayTr.doRaytracing()
                toc = time.time()
                print("Time elapsed for raytracing: {:.2f} seconds".format(toc - tic))

        else: # if input is a single laz file
            with laspy.open(self.laz_in_f) as file:
                count = 0
                for points in file.chunk_iterator(self.points_per_iter):
                    print("{:.2f}%".format(count / file.header.point_count * 100))

                    # For performance we need to use copy
                    # so that the underlying arrays are contiguous
                    x = points.x.copy()
                    y = points.y.copy()
                    z = points.z.copy()
                    gps_time = points.gps_time.copy()
                    return_number = points.return_number.copy()
                    number_of_returns = points.number_of_returns.copy()

                    if np.max(
                            return_number) == 0:  # a not very nice hack for the special case where return_number and number_of_returns are all 0 for Horizon measurements - TODO: figure out why!
                        return_number[:] = 1
                        number_of_returns[:] = 1

                    # for the case of mobile acquisitions, inerpolate trajectory for gps_time
                    if self.is_mobile:
                        SensorPos = self.interpolate_traj(self.traj['time'], self.traj['sensor_x'], self.traj['sensor_y'],
                                                          self.traj['sensor_z'], gps_time)
                    # call interpolate function for trajectory to extract sensor position for each gps_time



                    if np.max(number_of_returns) == 1 or np.max(return_number) == 1:
                        run_raytraycing_after_loading = False
                        self.RayTr.doRaytracing_singleReturnPulses(x, y, z, SensorPos['sensor_x'], SensorPos['sensor_y'],
                                                                   SensorPos['sensor_z'], gps_time)
                    else:
                        # check if gps_time is sorted
                        if count == 0:  # only check sort state in the first iteration of the for loop
                            if not is_sorted(gps_time):
                                print(
                                    f"!!!!! input laz file is not sorted along gps_time. The algorithm will still run. However, the "
                                    f"performance will be greatly decreased as the entire content of the laz file has to be read into "
                                    f"the system memory. If you have multi return data, consider sorting your laz data first, e.g. using "
                                    f"LASTools lassort: lassort -i laz_in -gps_time -return_number -odix _sort -olaz -v !!!!")
                                sorted = False
                                run_raytraycing_after_loading=True


                            else:
                                sorted = True
                                run_raytraycing_after_loading=False

                        self.RayTr.addPointData(x, y, z, SensorPos['sensor_x'], SensorPos['sensor_y'],
                                                SensorPos['sensor_z'],
                                                gps_time, return_number, number_of_returns)

                        if sorted:  # only if pulses are sorted run raytracing now. Otherwise we have to read in the entire dataset first!
                            # Get report on pulse dataset - comment this out once everythin is working or TODO: add a verbose flag!
                            # self.RayTr.getPulseDatasetReport()

                            # run raytracing on added points
                            # print("Do raytracing with stored pulses")
                            tic_r = time.time()
                            self.RayTr.doRaytracing()
                            toc_r = time.time()
                            # print("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))

                            self.RayTr.clearPulseDataset() # clear out data that have been traced.

                            # Check if traversed pulses have been deleted from map - comment this out once everything is working or TODO: add a verbose flag!
                            # self.RayTr.getPulseDatasetReport()

                    count = count + len(gps_time)

        if run_raytraycing_after_loading:
            self.RayTr.doRaytracing()

        self.get_raytracing_report()
        self.save_raytracing_output()

    def get_raytracing_report(self):
        # Get report on traversal
        self.RayTr.reportOnTraversal()

    def save_raytracing_output(self):
        print("Extracting Nhit")
        tic = time.time()
        self.Nhit = self.RayTr.getNhit()
        self.Nhit = np.array(self.Nhit, dtype=np.int32)

        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        print("Extracting Nocc")
        tic = time.time()
        self.Nocc = self.RayTr.getNocc()
        self.Nocc = np.array(self.Nocc, dtype=np.int32)

        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        print("Extracting Nmiss")
        tic = time.time()
        self.Nmiss = self.RayTr.getNmiss()
        self.Nmiss = np.array(self.Nmiss, dtype=np.int32)

        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        print("Saving Occlusion Outputs")
        tic = time.time()
        np.save(f"{self.out_dir}/Nhit.npy", self.Nhit)
        np.save(f"{self.out_dir}/Nmiss.npy", self.Nmiss)
        np.save(f"{self.out_dir}/Nocc.npy", self.Nocc)
        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        # Create Classification grid
        print("Classify Grid")
        tic = time.time()
        self.Classification = np.zeros((self.grid_dim['ny'], self.grid_dim['nx'], self.grid_dim['nz']), dtype=int)

        self.Classification[np.logical_and.reduce((self.Nhit > 0, self.Nmiss >= 0, self.Nocc >= 0))] = 1  # voxels that were observed
        self.Classification[np.logical_and.reduce((self.Nhit == 0, self.Nmiss > 0, self.Nocc >= 0))] = 2  # voxels that are empty
        self.Classification[
            np.logical_and.reduce((self.Nhit == 0, self.Nmiss == 0, self.Nocc > 0))] = 3  # voxels that are hidden (occluded)
        self.Classification[np.logical_and.reduce((self.Nhit == 0, self.Nmiss == 0,
                                              self.Nocc == 0))] = 4  # voxels that were not observed # TODO: Figure out, why this overwrites voxels that are classified as occluded! -> this was because np.logical_and only takes in 2 arrays as input, not 3! use np.logical_and.reduce() for that!

        np.save(f"{self.out_dir}/Classification.npy", self.Classification)
        toc = time.time()
        print("Elapsed Time: " + str(toc - tic) + " seconds")

    def normalize_occlusion_output(self, input_folder, dtm_file, dsm_file=None):
        """
        normalize_occlusion_output normalizes all occlusion output grids (Nhit, Nmiss, Nocc, Classification) with the specified DTM
        This function also calculates occlusion statistics for the total canopy volume (defined by the volume between DTM
        and DSM). Currently only binary occlusion is analysed at the moment (TODO: implement also fractional occlusion),
        i.e. only voxels that are completely occluded (Nhit==0 and Nmiss==0 and Nocc >0)
        :param input_folder: directory to the output of the raytracing algorithm
        :param dtm_file: DTM file (.tif) of the area of interest. Currently, both dimensions and pixel size should match the output grids
        :param dsm_file: DSM file (.tif) of the area of interest. Currently, both dimensions and pixel size should match the output grids
        :return:
        """

        self.Nhit = np.load(f"{input_folder}\\Nhit.npy")
        self.Nmiss = np.load(f"{input_folder}\\Nmiss.npy")
        self.Nocc = np.load(f"{input_folder}\\Nocc.npy")
        self.Classification = np.load(f"{input_folder}\\Classification.npy")

        # This is a bit of a quick and dirty solution to check on the compatibility of voxel size and pixel size of terrain models. TODO: improve that!
        dtm = TerrainModel(dtm_file)
        gt = dtm.dtm.res
        pix_size = gt[0]

        if pix_size != self.vox_dim:
            dtm.crop2extent(extent=(self.PlotDim['minX'], self.PlotDim['maxY'], self.PlotDim['maxX'], self.PlotDim['minY']),
                                   out_file=f"{dtm_file[:-4]}_resc_{self.vox_dim}.tif",
                                   res=self.vox_dim)

        ext = dtm.get_extent()
        extent_dtm = (ext.left, ext.top, ext.right, ext.bottom)
        extent_voxgrid = (self.PlotDim['minX'], self.PlotDim['maxY'], self.PlotDim['maxX'], self.PlotDim['minY'])
        if extent_dtm != extent_voxgrid:
            dtm.crop2extent(extent=extent_voxgrid,
                            out_file = f"{dtm.get_terrainmodel_path()[:-4]}_clipped.tif",
                            res=self.vox_dim)

        dtm_file = dtm.get_terrainmodel_path()

        dtm_src = rasterio.open(dtm_file)
        self.dtm = np.flipud(dtm_src.read(
            1))  # we need to flip the terrain models in order to make them compatible with the Occlusion output
        # fill in data gaps in dtm
        self.dtm = fillnodata(self.dtm, mask=self.dtm != dtm_src.get_nodatavals()[0])


        if dsm_file is not None:

            dsm = TerrainModel(dsm_file)
            # check on pixel size
            gt = dsm.dtm.res
            pix_size = gt[0]

            if pix_size != self.vox_dim:
                dsm.crop2extent(extent=(self.PlotDim['minX'], self.PlotDim['maxY'], self.PlotDim['maxX'], self.PlotDim['minY']),
                                   out_file=f"{dtm_file[:-4]}_resc_{self.vox_dim}.tif",
                                   res=self.vox_dim)

            ext = dsm.get_extent()
            extent_dsm = (ext.left, ext.top, ext.right, ext.bottom)
            if extent_dsm != extent_voxgrid:
                dsm.crop2extent(extent=extent_voxgrid,
                                out_file=f"{dsm.get_terrainmodel_path()[:-4]}_clipped.tif",
                                res=self.vox_dim)

            dsm_file = dsm.get_terrainmodel_path()
            dsm_src = rasterio.open(dsm_file)
            self.dsm = np.flipud(dsm_src.read(
                1))  # we need to flip the terrain models in order to make them compatible with the Occlusion output
            self.dsm = fillnodata(self.dsm, mask=self.dsm != dsm_src.get_nodatavals()[0])

            self.chm = self.dsm - self.dtm

            self.Nhit_norm = np.zeros((self.dtm.shape[0], self.dtm.shape[1], int(np.ceil(np.amax(self.chm) / self.vox_dim))), dtype=int)
            self.Nmiss_norm = np.zeros((self.dtm.shape[0], self.dtm.shape[1], int(np.ceil(np.amax(self.chm) / self.vox_dim))), dtype=int)
            self.Nocc_norm = np.zeros((self.dtm.shape[0], self.dtm.shape[1], int(np.ceil(np.amax(self.chm) / self.vox_dim))), dtype=int)
            self.Classification_norm = np.zeros((self.dtm.shape[0], self.dtm.shape[1], int(np.ceil(np.amax(self.chm) / self.vox_dim))),
                                           dtype=int)

            self.OcclFrac2D = np.zeros(self.dtm.shape)

            for y in range(0, self.dsm.shape[0], 1):
                for x in range(0, self.dsm.shape[1], 1):
                    # get zind where DTM is located in grid at x,y
                    zind_dtm = int(np.floor((self.dtm[y, x] - self.PlotDim['minZ']) / self.vox_dim))
                    zind_dsm = int(np.floor((self.dsm[y, x] - self.PlotDim['minZ']) / self.vox_dim))

                    # extract profile from grids
                    prof_class = self.Classification[y, x, zind_dtm:zind_dsm]
                    prof_class_buf = self.Classification[y, x, zind_dtm+int(np.ceil(self.lower_threshold/self.vox_dim)):zind_dsm]

                    self.Classification_norm[y, x, 0:len(prof_class)] = prof_class
                    # Calculate occlusion fraction for z profile
                    num_occl =sum(prof_class_buf == 3)

                    if len(prof_class_buf)==0:
                        self.OcclFrac2D[y, x] = 0
                    else:
                        self.OcclFrac2D[y, x] = num_occl/len(prof_class_buf)

                    self.Nhit_norm[y, x, 0:len(prof_class)] = self.Nhit[y, x, zind_dtm:zind_dsm]
                    self.Nmiss_norm[y, x, 0:len(prof_class)] = self.Nmiss[y, x, zind_dtm:zind_dsm]
                    self.Nocc_norm[y, x, 0:len(prof_class)] = self.Nocc[y, x, zind_dtm:zind_dsm]

                    self.__updateOccl_Volumes(prof_class)


        else:
            # as we do not know the height of the scene a priori, we will initialize a 3 D grid with the same dimensions
            # as the unnormalized grids, introducing quite some overhead...
            self.Nhit_norm = np.zeros(self.Nhit.shape, dtype=int)
            self.Nmiss_norm = np.zeros(self.Nmiss.shape, dtype=int)
            self.Nocc_norm = np.zeros(self.Nocc.shape, dtype=int)
            self.Classification_norm = np.zeros(self.Classification.shape, dtype=int)

            self.OcclFrac2D = np.zeros(self.dtm.shape)
            self.chm = np.zeros(self.dtm.shape)

            max_len_prof = 0
            for y in range(0, self.dtm.shape[0], 1):
                for x in range(0, self.dtm.shape[1], 1):
                    # get zind where DTM is located in grid at x,y
                    zind_dtm = int(np.floor((self.dtm[y, x] - self.PlotDim['minZ']) / self.vox_dim))
                    zind_dsm = last_nonzero(self.Nhit[y, x, :], axis=0) # If no dsm is provided, we take the DSM from the same
                    # acquisition. This will introduce an under estimation of occlusion for ground based acquisitions
                    # as occlusion on top of canopy is not counted.
                    self.chm[y,x] = (zind_dsm - zind_dtm) * self.vox_dim

                    # extract profile from grids
                    prof_class = self.Classification[y, x, zind_dtm:zind_dsm]
                    prof_class_buf = self.Classification[y, x,
                                     zind_dtm + int(np.ceil(self.lower_threshold / self.vox_dim)):zind_dsm]

                    self.Classification_norm[y, x, 0:len(prof_class)] = prof_class
                    # Calculate occlusion fraction for z profile
                    num_occl = sum(prof_class_buf == 3)

                    if len(prof_class_buf) == 0:
                        self.OcclFrac2D[y, x] = 0
                    else:
                        self.OcclFrac2D[y, x] = num_occl / len(prof_class_buf)

                    if len(prof_class) > max_len_prof:
                        max_len_prof = len(prof_class)

                    self.Classification_norm[y, x, 0:len(prof_class)] = self.Classification[y, x, zind_dtm:zind_dsm]
                    self.Nhit_norm[y, x, 0:len(prof_class)] = self.Nhit[y, x, zind_dtm:zind_dsm]
                    self.Nmiss_norm[y, x, 0:len(prof_class)] = self.Nmiss[y, x, zind_dtm:zind_dsm]
                    self.Nocc_norm[y, x, 0:len(prof_class)] = self.Nocc[y, x, zind_dtm:zind_dsm]

                    self.__updateOccl_Volumes(prof_class)

            # get rid of the excessive height of the grid
            self.Classification_norm = self.Classification_norm[:, :, 0:max_len_prof]
            self.Nhit_norm = self.Nhit_norm[:, :, 0:max_len_prof]
            self.Nmiss_norm = self.Nmiss_norm[:, :, 0:max_len_prof]
            self.Nocc_norm = self.Nocc_norm[:, :, 0:max_len_prof]

        print(f"Saving normalized output files into input directory...")
        np.save(f"{input_folder}/Nhit_norm.npy", self.Nhit_norm)
        np.save(f"{input_folder}/Nmiss_norm.npy", self.Nmiss_norm)
        np.save(f"{input_folder}/Nocc_norm.npy", self.Nocc_norm)
        np.save(f"{input_folder}/Classification_norm.npy", self.Classification_norm)

    #def visualize_2d_occlusion_map(self, out_fig):

    def __updateOccl_Volumes(self, prof_class):
        # update canopy volume and occlusion statistics
        self.TotalVolume = self.TotalVolume + len(prof_class)
        self.TotalOcclusion = self.TotalOcclusion + sum(prof_class==3)

        # Update Volume per strata
        if len(prof_class) >= 10 / self.vox_dim:
            self.Volume10_max = self.Volume10_max + len(prof_class[int(10 / self.vox_dim):])
            self.Occlusion10_max = self.Occlusion10_max + sum(prof_class[int(10 / self.vox_dim):] == 3)
            self.Volume3_10 = self.Volume3_10 + 7 / self.vox_dim  # The entire volume from 3 to 10 m (7m) is within the canopy
            self.Occlusion3_10 = self.Occlusion3_10 + sum(
                prof_class[int(3 / self.vox_dim):int(10 / self.vox_dim)] == 3)
            self.Volume0_3 = self.Volume0_3 + 3 / self.vox_dim  # The entire volume from 0 to 3 m (3m) is within the canopy
            self.Occlusion0_3 = self.Occlusion0_3 + sum(prof_class[0:int(3 / self.vox_dim)] == 3)
        elif len(prof_class) < 10 / self.vox_dim and len(prof_class) >= 3 / self.vox_dim:
            self.Volume3_10 = self.Volume3_10 + len(prof_class[int(3 / self.vox_dim):])
            self.Occlusion3_10 = self.Occlusion3_10 + sum(prof_class[int(3 / self.vox_dim):] == 3)
            self.Volume0_3 = self.Volume0_3 + 3 / self.vox_dim
            self.Occlusion0_3 = self.Occlusion0_3 + sum(prof_class[0:int(3 / self.vox_dim)] == 3)
        else:
            self.Volume0_3 = self.Volume0_3 + len(prof_class)
            self.Occlusion0_3 = self.Occlusion0_3 + sum(prof_class == 3)

    def get_Occl_TransectFigure(self, start_ind, end_ind, axis=0, format="png", show_plots=False):

        chm_slice_ref = None
        if axis==0: # get a slice of Y-Axis
            Nhit_Slice = np.sum(self.Nhit_norm[start_ind:end_ind, :, :], axis=axis)
            OcclFrac_Slice = np.sum(self.Classification_norm[start_ind:end_ind, :, :]==3, axis=axis) / (end_ind - start_ind)
            if self.chm is not None:
                chm_slice_ref = np.max(self.chm[start_ind:end_ind, :], axis=0)
        elif axis==1:
            Nhit_Slice = np.sum(self.Nhit_norm[:,start_ind:end_ind,:], axis=axis)
            OcclFrac_Slice = np.sum(self.Classification_norm[:, start_ind:end_ind, :] == 3, axis=axis) / (
                        end_ind - start_ind)
            if self.chm is not None:
                chm_slice_ref = np.max(self.chm[:, start_ind:end_ind], axis=1)
        else:
            Nhit_Slice = np.sum(self.Nhit_norm[:, :, start_ind:end_ind], axis=axis)
            OcclFrac_Slice = np.sum(self.Classification_norm[:, :, start_ind:end_ind] == 3, axis=axis) / (
                    end_ind - start_ind)

        NHits_Slice_log = np.log10(Nhit_Slice, where=(Nhit_Slice != 0))

        # we need to rotate the slice for visualization purposes
        OcclFrac_Slice = np.rot90(OcclFrac_Slice)
        NHits_Slice_log = np.rot90(NHits_Slice_log)

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1,1,1)
        x_axis_vect = None
        if axis==0:
            ax.set_xlabel(f"X [m]")
            ax.set_ylabel(f"Z [m]")
            extent = [self.PlotDim['minX'], self.PlotDim['maxX'], 0, OcclFrac_Slice.shape[0]*self.vox_dim]
            ax.axis(extent)
            x_axis_vect = np.linspace(start=self.PlotDim['minX'], stop=self.PlotDim['maxX'], num=self.grid_dim['nx'])

        elif axis==1:
            ax.set_xlabel(f"Y [m]")
            ax.set_ylabel(f"Z [m]")
            extent = [self.PlotDim['minY'], self.PlotDim['maxY'], 0, OcclFrac_Slice.shape[0] * self.vox_dim]
            ax.axis(extent)
            x_axis_vect = np.linspace(start=self.PlotDim['minY'], stop=self.PlotDim['maxY'], num=self.grid_dim['ny'])
        else:
            ax.set_xlabel(f"X [m]")
            ax.set_ylabel(f"Y [m]")
            extent = [self.PlotDim['minX'], self.PlotDim['maxX'], self.PlotDim['minY'], self.PlotDim['maxY']]
            ax.axis(extent)

        reds_cmap = plt.get_cmap(name='inferno_r')
        reds_cmap.set_under('k', alpha=0)
        greens_cmap = plt.get_cmap(name='Greens_r')
        greens_cmap.set_under('k', alpha=0)
        #plot raster data
        im1 = ax.imshow(NHits_Slice_log, cmap=greens_cmap, clim=[0.1, np.amax(NHits_Slice_log)], interpolation='none',
                        extent=extent)
        im2 = ax.imshow(OcclFrac_Slice * 100, cmap=reds_cmap, vmin=1, vmax=50, clim=[1, 100], interpolation='none',
                        alpha=0.8,
                        extent=extent)
        if x_axis_vect is not None:
            chm_ref_plot = ax.plot(x_axis_vect, chm_slice_ref, label="Ref CHM")
            # chm_comp_plot = ax.plot(x_axis_vect, chm_slice_comp, label="Comp CHM", linestyle='--') #TODO: implement that!
            ax.legend(handles=[chm_ref_plot[0]], loc='upper left', bbox_to_anchor=(0.375, 1))

        # define colorbars with position and dimension
        axins1 = inset_axes(
            ax,
            width="35%",
            height="5%",
            loc="upper right",
        )
        axins1.xaxis.set_ticks_position("bottom")
        fig.colorbar(im1, cax=axins1, orientation='horizontal', label="Log Nr.Hits")
        axins2 = inset_axes(
            ax,
            width="35%",
            height="5%",
            loc="upper left",
        )
        axins2.xaxis.set_ticks_position("bottom")
        fig.colorbar(im2, cax=axins2, orientation='horizontal', label="Occlusion [%]")




        # save figure
        if axis==0:
            plt.savefig(
                f"{self.out_dir}/Occlusion_Slice_XZ_{start_ind}_{end_ind}_voxels.{format}",
                dpi=300, format=format)
        elif axis==1:
            plt.savefig(
                f"{self.out_dir}/Occlusion_Slice_YZ_{start_ind}_{end_ind}_voxels.{format}",
                dpi=300, format=format)
        else:
            plt.savefig(
                f"{self.out_dir}/Occlusion_Slice_XY_{start_ind}_{end_ind}_voxels.{format}",
                dpi=300, format=format)

        if show_plots:
            plt.show()
        else:
            plt.close()

    def get_Occlusion_Profile(self, format="png", show_plots=False):

        OcclVertProf = np.sum(self.Classification_norm == 3, axis=0)
        OcclVertProf = np.sum(OcclVertProf, axis=0)
        OcclVertProf_Rel = OcclVertProf / ( (self.grid_dim['nx']) *
                                            (self.grid_dim['ny']))


        sns.set_theme(style="darkgrid")

        vert_vect = np.arange(start=0, stop=self.Classification_norm.shape[2] * self.vox_dim, step=self.vox_dim)
        # a hack to make sure that vert_vect is of the same length as OcclVertProf TODO: this has to be checked if it is generic!
        vert_vect = vert_vect[0:len(OcclVertProf)]

        # get rid of voxels below self.lower_threshold
        vert_vect = vert_vect[int(self.lower_threshold/self.vox_dim):]
        OcclVertProf_Rel = OcclVertProf_Rel[int(self.lower_threshold/self.vox_dim):]
        OcclVertProf = OcclVertProf[int(self.lower_threshold/self.vox_dim):]

        # extract height of max occlusion - we exclude the lowest 2 m to exclude occlusion from the ground
        ind_max_occl = np.argmax(OcclVertProf_Rel)
        h_max_occl = vert_vect[ind_max_occl]

        occl_vert_prof = pd.DataFrame(
            data={'vert_vect': vert_vect, 'Occl_Sum': OcclVertProf, 'OcclRel': OcclVertProf_Rel * 100})

        # Plot the vertical occlusion profile
        fig2 = plt.figure(figsize=(4.5, 7))
        ax2 = fig2.add_subplot(1, 1, 1)
        ax2.set_ylabel(f"Occlusion [%]")
        ax2.set_xlabel(f"Height above ground [m]")
        line = lineplot_plusplus(x="vert_vect", y="OcclRel", data=occl_vert_prof, orientation="vertical", color='blue')
        # add horizontal line at mean canopy height
        mean_canopy_h = self.chm.mean()
        mcl = line.axhline(mean_canopy_h, color='r', linestyle='--', label='Mean canopy height')
        line_proxy = mlines.Line2D([], [], color='blue', label="MLS Occlusion")
        ax2.legend(handles=[line_proxy, mcl], loc='upper right', labels=["MLS Occlusion", "Mean ALS canopy height"])
        plt.savefig(f"{self.out_dir}/OcclusionVertProf.{format}", dpi=300, format=format)
        if show_plots:
            plt.show()
        else:
            plt.close()
       

        return occl_vert_prof

    def clean_up_RayTr(self):
        del self.RayTr
