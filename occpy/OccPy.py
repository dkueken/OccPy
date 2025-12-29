import time
import os
import glob

import numpy as np
import pandas as pd
import laspy
import OSToolBox as ost
from tqdm import tqdm

# plotting functions
import matplotlib
try:
    matplotlib.use('TkAgg')
except ImportError:
    print("couldn't change matplotlib backend")
    
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import to_rgb
import seaborn as sns

from raytr import PyRaytracer
from occpy.TerrainModel import TerrainModel
from occpy.util import prepare_ply, read_trajectory_file, read_sensorpos_file, interpolate_traj, last_nonzero

# TODO: change print statements to logging like in occpyRIEGL

is_sorted = lambda a: np.all(a[:-1] <= a[1:])

# Function to darken an RGB color
def darken_color(color, amount=0.6):
    r, g, b = to_rgb(color)
    return (r * amount, g * amount, b * amount)

""" Moved to visualization module. delete once everything is working.
def get_Occlusion_ProfileFigure(Classification, plot_dim, vox_dim, out_dir, low_thresh=0, vertBuffer=0, max_percentage=100, fig_prop=None, show_plots=False):

    grid_dim = (int((plot_dim[3] - plot_dim[0]) / vox_dim), int((plot_dim[4] - plot_dim[1]) / vox_dim),
                int((plot_dim[5] - plot_dim[2]) / vox_dim))

    vert_vect = np.arange(start=low_thresh, stop=Classification.shape[2] * vox_dim, step=vox_dim)
    Classification = Classification[:,:,int(low_thresh / vox_dim):]
    # a hack to make sure that vert_vect is of the same length as OcclVertProf TODO: this has to be checked if it is generic!


    OcclVertProf = np.sum(Classification == 3, axis=0)
    OcclVertProf = np.sum(OcclVertProf, axis=0)
    OcclVertProf_Rel = OcclVertProf / ((grid_dim[0]) * (grid_dim[1]))

    FilledVertProf = np.sum(Classification == 1, axis=0)
    FilledVertProf = np.sum(FilledVertProf, axis=0)
    FilledVertProf_Rel = FilledVertProf / (grid_dim[0] * grid_dim[1])

    EmptyVertProf = np.sum(np.logical_or(Classification == 2, Classification==0), axis=0)
    EmptyVertProf = np.sum(EmptyVertProf, axis=0)
    EmptyVertProf_Rel = EmptyVertProf / (grid_dim[0] * grid_dim[1])

    heights = vert_vect[0:len(OcclVertProf)]

    percentages = np.column_stack([FilledVertProf_Rel*100, OcclVertProf_Rel*100, EmptyVertProf_Rel*100])
    categories = ['Filled', 'Occluded', 'Empty']
    colors = ['skyblue', 'salmon', 'lightgreen']

    # Compute cumulative percentages for stacking
    cumulative = np.cumsum(percentages, axis=1)

    palette = sns.color_palette('colorblind', n_colors=len(categories))
    palette[2] = (1.0, 1.0, 1.0) # white for empty

    fig, ax = plt.subplots(figsize=fig_prop['fig_size'])

    for i, cat in enumerate(categories):
        left = cumulative[:, i - 1] if i > 0 else np.zeros_like(heights)
        face_color = palette[i]
        edge_color = darken_color(face_color, 0.8)  # slightly darker for lines

        # Fill area
        ax.fill_betweenx(
            heights, left, cumulative[:, i],
            color=face_color, alpha=0.6
        )
        # Outline
        ax.plot(cumulative[:, i], heights, color=edge_color, linewidth=1.5, label="_nolegend_")

    ax.set_xlabel('Percentage of voxels [%]', fontsize=fig_prop['label_size'])
    ax.set_ylabel('Height above ground [m]', fontsize=fig_prop['label_size'])
    ax.set_xlim(0.1,max_percentage)
    ax.set_ylim(0,np.max(heights) + vertBuffer)
    plt.xticks(fontsize=fig_prop['label_size_ticks'])
    plt.yticks(fontsize=fig_prop['label_size_ticks'])
    ax.legend(categories[0:2], fontsize=fig_prop['label_size_ticks'])
    plt.tight_layout()

    plt.savefig(f"{out_dir}/OcclusionVertProf.{fig_prop['out_format']}", dpi=300, format=fig_prop['out_format'])
    if show_plots:
        plt.show(block=True)
    else:
        plt.close()
"""




class OccPy:
    def __init__(self, laz_in, out_dir, vox_dim=0.1, lower_threshold=1, points_per_iter=10000000, plot_dim=None, output_voxels=False):
        """
        initialize OccPy object

        Parameters
        ----------
        laz_in: str
            path to laz file or directory with multiple laz files
        out_dir: str
            path to output directory for occlusion maps
        vox_dim: float [default: 0.1]
            voxel dimension in meters. We currently assume cubic voxels.
        lower_threshold: float [default: 1]
            lower threshold above ground to exclude from occlusion mapping in voxels. TODO: figure out, if this is needed at this place!
        points_per_iter: int [default: 10000000
            number of points read in from laz file in each iteration. Warning: this only applies to LAZ files with either
            only single returns or when points have been sorted by GPSTime and return number.
        plot_dim: list [default: None]
            corner coordinates of the plot in question. Expected format is: [minX, minY, minZ, maxX, maxY, maxZ]
        output_voxels: bool [default: False]
            whether voxel grid should be outputed as ply file. WARNING: this requires significant resources and takes a long time. Has not been properly tested!
        """
        self.laz_in_f = laz_in
        self.out_dir = out_dir
        os.makedirs(out_dir, exist_ok=True)
        self.vox_dim = vox_dim  # voxel dimension (cubic) in meters TODO: maybe implement non-cubic voxels?
        self.lower_threshold = lower_threshold  # lower threshold above ground to exclude TODO: check if necessary
        self.points_per_iter = points_per_iter  # number of points read in from laz file in each iteration
        self.output_voxels = output_voxels
        self.is_mobile = False
        self.traj_f = None
        self.traj = None
        self.senspos_f = None
        self.senspos = None
        self.single_return = None
        # TODO: find a better way to link scan position id between laz file and scan pos file!
        self.scan_pos_id_stridx = 0
        self.scan_pos_id_endstridx = 0

        # some parameters that will be filled during function calls
        self.TotalVolume = 0
        self.TotalOcclusion = 0

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
        minBound = np.array([self.PlotDim['minX'], self.PlotDim['minY'], self.PlotDim['minZ']])
        maxBound = np.array([self.PlotDim['maxX'], self.PlotDim['maxY'], self.PlotDim['maxZ']])
        self.RayTr.defineGrid(minBound, maxBound, self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz'],
                              self.vox_dim)

    def define_sensor_pos(self, path2file, is_mobile, single_return=None, delimiter=" ", hdr_time='%time', hdr_scanpos_id='', hdr_x='x', hdr_y='y', hdr_z='z', sens_pos_id_offset=0, str_idx_ScanPosID=0, str_end_idx_ScanPosID=0):
        """
        defines sensor position based on the provided csv file. CSV file needs to include
        Parameters
        ----------
        path2file: str [mandatory]
            path to csv file with sensor position information
        is_mobile: bool [mandatory]
            True or False whether platform is mobile (MLS, ULS) or static (TLS)
        single_return: bool [adviced]
            True or False whether data is single return or multi return data
        delimiter: str [default: " "]
            csv delimiter
        hdr_time: str
            column header for time (only needed for mobile acquisitions)
        hdr_scanpos_id str
            column header for scan pos id (only needed for static acquisitions -> equivalent to hdr_time in mobile acquisitions
        hdr_x: str [default 'x']
            column header for x coordinates
        hdr_y: str [default 'y']
            column header for y coordinates
        hdr_z: str [default 'z']
            column header for z coordinates
        sens_pos_id_offset: int
            Very specific use case where Scan Pos ID in position file does not correspond with Scan Pos ID in LAZ files and we need to add an offset
        str_idx_ScanPosID: int
            string index of where the scan position identifier is written in the laz file name TODO: find a better way to handle this!
        str_end_idx_ScanPosID: int
            string index of where the scan position identifier ends in the laz file name TODO: find a better way to handle this!

        Returns
        -------

        """
        self.is_mobile = is_mobile
        self.single_return = single_return

        if is_mobile: # case of mobile acquisitions (MLS, ULS)
            self.traj = read_trajectory_file(path2traj=path2file, delimiter=delimiter, hdr_time=hdr_time, hdr_x=hdr_x, hdr_y=hdr_y, hdr_z = hdr_z)
        else: # case of static acquisition (TLS)
            self.senspos = read_sensorpos_file(path2senspos=path2file, delimiter=delimiter, hdr_scanpos_id=hdr_scanpos_id, hdr_x=hdr_x, hdr_y=hdr_y, hdr_z=hdr_z, sens_pos_id_offset=sens_pos_id_offset)
            self.scan_pos_id_stridx = str_idx_ScanPosID
            self.scan_pos_id_endstridx = str_end_idx_ScanPosID

    # TODO: change to structure of occpyRIEGL: read input las files and link to scan positions/ trajectory first (how to do this for MLS/ULS?) -> allows us to error out/log early when position link is not properly found

    # TODO: ugly workaround for the case where a single laz file from a single TLS position should be run
    def define_sensor_pos_singlePos(self, scan_pos_id, x, y, z):
        """
        defines the scanner position of a single TLS scan position. This is currently just a work-around where we have
        the case of a single laz file and a single position without a text file defining e.g. multiple scan positions.
        Writes scanner position into self.senspos

        Parameters
        ----------
        scan_pos_id: int
            Scan Position Identificaiton number
        x: float
            X-Coordinates of scanner position
        y: float
            Y-Coordinates of scanner position
        z: float
            Z-Coordinates of scanner position

        Returns
        -------

        """
        d = {'ScanPos': scan_pos_id,
             'sensor_x': x, 'sensor_y': y, 'sensor_z': z}

        senspos = pd.DataFrame(data=d, index=[0])

        self.senspos = senspos


    def do_raytracing(self):
        """
        Perform ray tracing.

        This method processes either a directory of LAZ files (for TLS with known scan positions) or a single LAZ file
        (Single TLS position or MLS/ULS with a trajectory, depends on self.is_mobile).
        In the case of TLS, for each LAZ file, it extracts point positions and sensor positions, and
        then performs ray tracing, accounting for single or multi-return pulse information. Multi-return handling
        supports on-the-fly processing if the data is sorted by GPS time; otherwise, the full dataset must be loaded first.

        Raises
        ------
        FileNotFoundError
            TODO: Needs to be implemented: If `self.laz_in_f` is not a valid file or directory.

        RuntimeWarning
            TODO: Needs to be implemented: If multi-return data is detected but the LAZ file is not sorted by GPS time.

        """

        run_raytraycing_after_loading = False
        if os.path.isdir(self.laz_in_f):
            ## get list of laz files in input directory
            fCont = glob.glob(f"{self.laz_in_f}/*.laz")

            for f in fCont:
                # get scan position
                # lastdir_ind = f.rfind("/")
                scan_name = os.path.basename(f)
                # TODO: This is very specific to the test data and needs to be made generic!
                #end_of_scanID_idx = scan_name.rfind("_")
                scan_id = int(scan_name[self.scan_pos_id_stridx:self.scan_pos_id_endstridx])

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

                        # if self.single_return is not set, check data to find out if there are multiple returns stored per pulse
                        if self.single_return == None:
                            print(
                                f"TIPP: it was not specified if the point cloud stores single return data. to avoid any "
                                f"issues where there are only single return pulses in the first chunk, put there are "
                                f"multi return pulses in the whole dataset, please set the single_return variable within "
                                f"the define_sensor_pos function")
                            max_ret_num = np.max(points.return_number)
                            if max_ret_num > 1:
                                self.single_return = False
                            else:
                                self.single_return = True

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
                with tqdm(total=file.header.point_count, desc="Tracing Pulses...", unit="pulses") as pbar:
                    for points in file.chunk_iterator(points_per_iteration=self.points_per_iter):

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
                            # call interpolate function for trajectory to extract sensor position for each gps_time
                            SensorPos = interpolate_traj(self.traj['time'], self.traj['sensor_x'], self.traj['sensor_y'],
                                                              self.traj['sensor_z'], gps_time)

                        else:
                            SensorPos = self.senspos

                            # TODO: figure out, why we have to repeat scan position to the shape of input point cloud and try to implement it, so that we only have to pass one position
                            SensorPos = pd.DataFrame(data={'ScanPos': np.ones(gps_time.shape) * self.senspos['ScanPos'].values[0],
                                                           'sensor_x': np.ones(gps_time.shape) * self.senspos['sensor_x'].values[0],
                                                           'sensor_y': np.ones(gps_time.shape) * self.senspos['sensor_y'].values[0],
                                                           'sensor_z': np.ones(gps_time.shape) * self.senspos['sensor_z'].values[0]})




                        if np.max(number_of_returns) == 1 or np.max(return_number) == 1:
                            run_raytraycing_after_loading = False
                            self.RayTr.doRaytracing_singleReturnPulses(x, y, z,  SensorPos['sensor_x'], SensorPos['sensor_y'],
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
                        pbar.update(len(points))

        if run_raytraycing_after_loading:
            self.RayTr.doRaytracing()

        self.get_raytracing_report()
        self.save_raytracing_output()

    def get_raytracing_report(self):
        """
        prints a report on the voxel traversal to the console

        Returns
        -------

        """
        # Get report on traversal
        self.RayTr.reportOnTraversal()

    def save_raytracing_output(self):
        """
        Extract and save the outputs of the ray-tracing process.
        This method performs the following steps:
        1. Extracts the voxel-wise hit (`Nhit`), miss (`Nmiss`), and occlusion (`Nocc`) voxelgrids and saves as .npy in self.out_dir
        2. Creates a voxel classification grid based on the `Nhit`, `Nmiss`, and `Nocc` values:
        - 1 = observed (hit > 0)
        - 2 = empty (miss > 0, hit == 0)
        - 3 = occluded (occlusion > 0, hit == 0, miss == 0)
        - 4 = unobserved (all three == 0)
        3. Writes `.ply` files for all voxel outputs if `self.output_voxels` is True (takes long and creates large files)

        Returns
        -------

        """
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

        print("Saving Occlusion Outputs As .npy")
        tic = time.time()
        np.save(f"{self.out_dir}/Nhit.npy", self.Nhit)
        np.save(f"{self.out_dir}/Nmiss.npy", self.Nmiss)
        np.save(f"{self.out_dir}/Nocc.npy", self.Nocc)
        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        # Create Classification grid
        print("Classify Grid")
        tic = time.time()
        self.Classification = np.zeros((self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz']), dtype=int)

        self.Classification[np.logical_and.reduce((self.Nhit > 0, self.Nmiss >= 0, self.Nocc >= 0))] = 1  # voxels that were observed
        self.Classification[np.logical_and.reduce((self.Nhit == 0, self.Nmiss > 0, self.Nocc >= 0))] = 2  # voxels that are empty
        self.Classification[
            np.logical_and.reduce((self.Nhit == 0, self.Nmiss == 0, self.Nocc > 0))] = 3  # voxels that are hidden (occluded)
        self.Classification[np.logical_and.reduce((self.Nhit == 0, self.Nmiss == 0,
                                              self.Nocc == 0))] = 4  # voxels that were not observed # TODO: Figure out, why this overwrites voxels that are classified as occluded! -> this was because np.logical_and only takes in 2 arrays as input, not 3! use np.logical_and.reduce() for that!

        np.save(f"{self.out_dir}/Classification.npy", self.Classification)
        toc = time.time()
        print("Elapsed Time: " + str(toc - tic) + " seconds")

        # write ply file
        if self.output_voxels:
            print("Saving Occlusion Outputs As .ply")
            tic = time.time()
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nhit)
            ost.write_ply(f"{self.out_dir}/Nhit.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nmiss)
            ost.write_ply(f"{self.out_dir}/Nmiss.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nocc)
            ost.write_ply(f"{self.out_dir}/Nocc.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Classification)
            ost.write_ply(f"{self.out_dir}/Classification.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            self.occl = np.zeros(shape=self.Classification.shape)
            x4, y4, z4 = np.where(self.Classification == 4)
            self.occl[x4, y4, z4] = self.Classification[x4, y4, z4]
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.occl)
            ost.write_ply(f"{self.out_dir}/Occl.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            toc = time.time()
            print("Elapsed Time: " + str(toc - tic) + " seconds")


    def get_chm(self):
        """
        returns chm.
        TODO: check if this is still needed!
        Returns
        -------
        reference to canopy height model

        """
        if self.chm is None:
            print("No CHM was defined. To define CHM ")
        return self.chm

    def clean_up_RayTr(self):
        """
        Free up memory after raytracing.

        Returns
        -------

        """
        del self.RayTr

    def getGridDimensions(self):
        """
        Get the grid dimensions
        Returns
        -------
        vector with grid diemsnions [minx, maxx, miny, maxy]

        """
        gridDim = self.RayTr.getGridDimensions()
        gridDim = np.asarray(gridDim, dtype=np.int32)

        return gridDim
