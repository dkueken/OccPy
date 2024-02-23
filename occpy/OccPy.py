import numpy as np
import pandas as pd
from scipy import interpolate

import laspy
from occpy.TerrainModel import TerrainModel

import time
import os

# we need to add the path to the src folder so raytr module can be found
import sys
sys.path.append(r"./src") # TODO: check if src folder should be better placed into occpy package
from occpy.src.raytr import PyRaytracer


class OccPy:
    def __init__(self, laz_in, out_dir, vox_dim=0.1, lower_threshold=1, points_per_iter=10000000, plot_dim=None):
        self.laz_in_f = laz_in
        self.out_dir = out_dir
        os.makedirs(out_dir, exist_ok=True)
        self.vox_dim = vox_dim    # voxel dimension (cubic) in meters TODO: maybe implement non-cubic voxels?
        self.lower_threshold = lower_threshold    # lower threshold above ground to exclude TODO: check if necessary
        self.points_per_iter = points_per_iter    # number of points read in from laz file in each iteration
        self.is_mobile = False
        self.traj_f = None
        self.traj = None
        self.senspos_f = None
        self.senspos = None

        if plot_dim is None:
            #TODO: Test if this works!
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

        self.grid_dim = dict(nx=int((self.PlotDim['maxX'] - self.PlotDim['minX'])/self.vox_dim),
                             ny=int((self.PlotDim['maxY'] - self.PlotDim['minY'])/self.vox_dim),
                             nz=int((self.PlotDim['maxZ'] - self.PlotDim['minZ'])/self.vox_dim))

        # initialize RayTr Object
        self.RayTr = PyRaytracer()

        # Define Grid
        minBound = np.array([self.PlotDim['minY'], self.PlotDim['minX'], self.PlotDim['minZ']])
        maxBound = np.array([self.PlotDim['maxY'], self.PlotDim['maxX'], self.PlotDim['maxZ']])
        self.RayTr.defineGrid(minBound, maxBound, self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz'], self.vox_dim)


    def read_trajectory_file(self, path2traj, delimiter=" ", hdr_time='%time', hdr_x='x', hdr_y='y', hdr_z='z'):
        self.is_mobile = True
        traj_in = pd.read_csv(path2traj, sep=delimiter)

        d = {'time': traj_in[hdr_time], 'sensor_x': traj_in[hdr_x], 'sensor_y': traj_in[hdr_y], 'sensor_z': traj_in[hdr_z]}

        self.traj = pd.DataFrame(data=d)

    def read_sensorpos_file(self, path2senspos, delimiter=" ", hdr_time='', hdr_x='', hdr_y='', hdr_z=''):
        print('TODO!')

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
                    SensorPos = self.interpolate_traj(self.traj['time'], self.traj['sensor_x'], self.traj['sensor_y'], self.traj['sensor_z'], gps_time)
                # call interpolate function for trajectory to extract sensor position for each gps_time


                if np.max(number_of_returns)==1 or np.max(return_number)==1:
                    self.RayTr.doRaytracing_singleReturnPulses(x, y, z, SensorPos['sensor_x'], SensorPos['sensor_y'],
                                                      SensorPos['sensor_z'], gps_time)
                else:
                    self.RayTr.addPointData(x, y, z, SensorPos['sensor_x'], SensorPos['sensor_y'], SensorPos['sensor_z'],
                                       gps_time, return_number, number_of_returns)

                count = count + len(gps_time)

        self.get_raytracing_report()
        self.save_raytracing_output()

    def get_raytracing_report(self):
        # Get report on traversal
        self.RayTr.reportOnTraversal()

    def save_raytracing_output(self):
        print("Extracting Nhit")
        tic = time.time()
        Nhit= self.RayTr.getNhit()
        Nhit = np.array(Nhit, dtype=np.int32)

        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        print("Extracting Nocc")
        tic = time.time()
        Nocc = self.RayTr.getNocc()
        Nocc = np.array(Nocc, dtype=np.int32)

        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        print("Extracting Nmiss")
        tic = time.time()
        Nmiss = self.RayTr.getNmiss()
        Nmiss = np.array(Nmiss, dtype=np.int32)

        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        print("Saving Occlusion Outputs")
        tic = time.time()
        np.save(f"{self.out_dir}/Nhit.npy", Nhit)
        np.save(f"{self.out_dir}/Nmiss.npy", Nmiss)
        np.save(f"{self.out_dir}/Nocc.npy", Nocc)
        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        # Create Classification grid
        print("Classify Grid")
        tic = time.time()
        Classification = np.zeros((self.grid_dim['ny'], self.grid_dim['nx'], self.grid_dim['nz']), dtype=int)

        Classification[np.logical_and.reduce((Nhit > 0, Nmiss >= 0, Nocc >= 0))] = 1  # voxels that were observed
        Classification[np.logical_and.reduce((Nhit == 0, Nmiss > 0, Nocc >= 0))] = 2  # voxels that are empty
        Classification[
            np.logical_and.reduce((Nhit == 0, Nmiss == 0, Nocc > 0))] = 3  # voxels that are hidden (occluded)
        Classification[np.logical_and.reduce((Nhit == 0, Nmiss == 0,
                                              Nocc == 0))] = 4  # voxels that were not observed # TODO: Figure out, why this overwrites voxels that are classified as occluded! -> this was because np.logical_and only takes in 2 arrays as input, not 3! use np.logical_and.reduce() for that!

        np.save(f"{self.out_dir}/Classification.npy", Classification)
        toc = time.time()
        print("Elapsed Time: " + str(toc - tic) + " seconds")

    def clean_up_RayTr(self):
        del self.RayTr



