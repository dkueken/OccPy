import numpy as np
import pickle
import pandas as pd
import laspy
import time
import zarr

from raytr import PyRaytracer

is_sorted = lambda a: np.all(a[:-1] <= a[1:])

# Define parameters
laz_in = r"D:\_tmp_wdir\OcclusionMapping_VZ400i_Rameren\ScanPos_001.laz"
out_dir = r"D:\_tmp_wdir\OcclusionMapping_VZ400i_Rameren\OcclusionMapping_ScanPos001"

points_per_iter = 10000000
ScanPos = dict(X=2676485.407984,
               Y=1246058.337989,
               Z=551.136297)

PlotDim = dict(minX=2676480,
               maxX=2676550,
               minY=1246050,
               maxY=1246120,
               minZ=540,
               maxZ=585)

vox_dim = 0.1

# Prepare Raytracer
gridDim = dict(nx=int((PlotDim['maxX'] - PlotDim['minX']) / vox_dim),
               ny=int((PlotDim['maxY'] - PlotDim['minY']) / vox_dim),
               nz=int((PlotDim['maxZ'] - PlotDim['minZ']) / vox_dim)
               )

# initialize output grids
Nhit = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)
Nmiss = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)
Nocc = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)


# read in laz file and create points dataset to be stored as a npy file with sorted gpstime/returnnumber
with laspy.open(laz_in) as file:
    count = 0
    for points in file.chunk_iterator(points_per_iteration=points_per_iter):
        print("{:.2f}%".format(count / file.header.point_count * 100))

        x = points.x.copy()
        y = points.y.copy()
        z = points.z.copy()
        gps_time = points.gps_time.copy()
        return_number = points.return_number.copy()
        number_of_returns = points.number_of_returns.copy()

        sensor_x = np.ones(gps_time.shape) * ScanPos['X']
        sensor_y = np.ones(gps_time.shape) * ScanPos['Y']
        sensor_z = np.ones(gps_time.shape) * ScanPos['Z']

        # Performe raytracing
        tic_tot = time.time()

        # init Raytracer object
        print("Initializing Raytracer")
        RayTr = PyRaytracer()

        tic = time.time()
        print("Define Grid")
        minBound = np.array([PlotDim['minY'], PlotDim['minX'], PlotDim['minZ']])
        maxBound = np.array([PlotDim['maxY'], PlotDim['maxX'], PlotDim['maxZ']])
        RayTr.defineGrid(minBound, maxBound, gridDim['nx'], gridDim['ny'], gridDim['nz'], vox_dim)
        toc = time.time()
        print("Time elapsed: {:.2f} seconds".format(toc - tic))

        RayTr.addPointData(x, y, z, sensor_x, sensor_y, sensor_z, gps_time, return_number, number_of_returns)

        RayTr.getPulseDatasetReport()


