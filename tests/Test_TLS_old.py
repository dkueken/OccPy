import os

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

# read in laz file and create points dataset to be stored as a npy file with sorted gpstime/returnnumber
with laspy.open(laz_in) as file:
    count = 0
    for points in file.chunk_iterator(points_per_iteration=points_per_iter):
        print("{:.2f}%".format(count / file.header.point_count * 100))

        if count==0:
            x = points.x.copy()
            y = points.y.copy()
            z = points.z.copy()
            gps_time = points.gps_time.copy()
            return_number = points.return_number.copy()
            number_of_returns = points.number_of_returns.copy()

        else:
            x = np.vstack((x, points.x.copy()))
            y = np.vstack((y, points.y.copy()))
            z = np.vstack((z, points.z.copy()))
            gps_time = np.vstack((gps_time, points.gps_time.copy()))
            return_number = np.vstack((return_number, points.return_number.copy()))
            number_of_returns = np.vstack((number_of_returns, points.number_of_returns.copy()))

        count += len(points)

# sort stacked pdata according to gps_time and return_number
if ~is_sorted(gps_time):
    sort_ind = np.lexsort((return_number, gps_time)) # sort first by gps_time, then by return_number ... somehow the key for ordering first is metnioned last ...
    x = x[sort_ind]
    y = y[sort_ind]
    z = z[sort_ind]
    gps_time = gps_time[sort_ind]
    return_number = return_number[sort_ind]
    number_of_returns = number_of_returns[sort_ind]


# write sorted pdata into zarr file
zarr.save(f"{out_dir}/PulseData.zip", gps_time=gps_time, return_number=return_number,
          number_of_returns=number_of_returns, x=x, y=y, z=z)

del gps_time, return_number, number_of_returns, x, y, z

# Prepare Raytracer
gridDim = dict(nx=int((PlotDim['maxX'] - PlotDim['minX'])/vox_dim),
               ny=int((PlotDim['maxY'] - PlotDim['minY'])/vox_dim),
               nz=int((PlotDim['maxZ'] - PlotDim['minZ'])/vox_dim)
               )

# initialize output grids
Nhit = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)
Nmiss = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)
Nocc = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)

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

# Load the previously generated pulse dataset with sorted pulses
pdata = zarr.load(f"{out_dir}/PulseData.zip")

if points_per_iter > len(pdata['gps_time']):
    gps_time = pdata['gps_time']
    return_number = pdata['return_number']
    number_of_returns = pdata['number_of_returns']
    x = pdata['x']
    y = pdata['y']
    z = pdata['z']

    SensorPosX = np.ones(gps_time.shape)*ScanPos['X']
    SensorPosY = np.ones(gps_time.shape)*ScanPos['Y']
    SensorPosZ = np.ones(gps_time.shape)*ScanPos['Z']

    print("Do actual raytracing with all pulses")
    tic = time.time()
    RayTr.doRaytracing(x, y, z, SensorPosX, SensorPosY, SensorPosZ, gps_time)
    toc = time.time()

    print("Elapsed Time: {:.2f} seconds".format(toc - tic))

print("Extracting Nhit")
tic = time.time()
Nhit_tmp = RayTr.getNhit()
Nhit_tmp = np.array(Nhit_tmp, dtype=np.int32)
Nhit = Nhit + Nhit_tmp
del Nhit_tmp

toc = time.time()
print("Elapsed Time: {:.2f} seconds".format(toc - tic))

print("Extracting Nocc")
tic = time.time()
Nocc_tmp = RayTr.getNocc()
Nocc_tmp = np.array(Nocc_tmp, dtype=np.int32)
Nocc = Nocc + Nocc_tmp
del Nocc_tmp

toc = time.time()
print("Elapsed Time: {:.2f} seconds".format(toc - tic))

print("Extracting Nmiss")
tic = time.time()
Nmiss_tmp = RayTr.getNmiss()
Nmiss_tmp = np.array(Nmiss_tmp, dtype=np.int32)
Nmiss = Nmiss + Nmiss_tmp
del Nmiss_tmp

toc = time.time()
print("Elapsed Time: {:.2f} seconds".format(toc - tic))

del RayTr

toc_tot = time.time()
print("Total Elapsed Time for Voxel Traversal: {:.2f} seconds".format(toc_tot - tic_tot))






