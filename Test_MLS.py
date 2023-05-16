import os

import numpy as np
import pickle
import pandas as pd
import laspy
import time
from prepare_trajectory import interpolate_traj

from raytr import PyRaytracer

# Input parameters
CLNR = 7379
OP = 2
laz_in = r"\\speedy11-12-fs\data_15\_PLS\20220503_LFI_CH_Campaign\CLNR_7379\Operator_2\LAZ\CLNR_7379_TRAN_OP02_2022-05-03_13-17-04_100pct_height_world_rot2LV95.laz"
laz_clip = r"\\speedy11-12-fs\data_15\_PLS\20220503_LFI_CH_Campaign\CLNR_7379\Operator_2\LAZ\CLNR_7379_TRAN_OP02_2022-05-03_13-17-04_100pct_height_world_rot2LV95_clip2plot.laz"
traj_in = r"\\speedy11-12-fs\data_15\_PLS\20220503_LFI_CH_Campaign\CLNR_7379\Operator_2\Trajectories\CLNR_7379_TRAN_OP02_2022-05-03_13-17-04_results_traj_rot2LV95.txt"
plot_dim_pkl = r"\\speedy11-12-fs\data_15\_PLS\20220503_LFI_CH_Campaign\CLNR_7379\Operator_2\Occlusion_Mapping\PlotDimensions.pkl"
dtm_file = r"\\speedy11-12-fs\data_15\_PLS\20220503_LFI_CH_Campaign\CLNR_7379\Operator_2\DEM\DTM_7379_swissAlti3D_10cm.tif"
dsm_file = r"\\speedy11-12-fs\data_15\_PLS\20220503_LFI_CH_Campaign\CLNR_7379\Operator_2\DEM\DSM_7379_swissSurface3D_10cm.tif"

out_dir = r"D:\Projects\OcclusionMapping_Timelapse\CLNR_7379\Operator_2\Frames\\"
os.makedirs(os.path.dirname(out_dir), exist_ok=True)

parameters = dict(
    sensorType=3,
    voxDim=0.1,
    lower_threshold=1,
    points_per_iter=10000000,
)

# read in trajectory file
traj = pd.read_csv(traj_in, sep=" ")

# Read in Plot dimensions from pickle
with open(plot_dim_pkl, 'rb') as f:
    PlotDim = pickle.load(f)

gridDim = dict(nx=int((PlotDim['maxX'] - PlotDim['minX'])/parameters['voxDim']),
               ny=int((PlotDim['maxY'] - PlotDim['minY'])/parameters['voxDim']),
               nz=int((PlotDim['maxZ'] - PlotDim['minZ'])/parameters['voxDim'])
               )


# initialize output grids
Nhit = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)
Nmiss = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)
Nocc = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)

# read in laz iteratively and run raytracing tool
tic_tot = time.time()
with laspy.open(laz_in) as file:

    # init Raytracer object
    print("Initializing Raytracer")
    RayTr = PyRaytracer()

    tic = time.time()
    print("Define Grid")
    minBound = np.array([PlotDim['minY'], PlotDim['minX'], PlotDim['minZ']])
    maxBound = np.array([PlotDim['maxY'], PlotDim['maxX'], PlotDim['maxZ']])
    RayTr.defineGrid(minBound, maxBound, gridDim['nx'], gridDim['ny'], gridDim['nz'], parameters['voxDim'])
    toc = time.time()
    print("Time elapsed: {:.2f} seconds".format(toc - tic))

    count = 0
    for points in file.chunk_iterator(parameters['points_per_iter']):
        print("{:.2f}%".format(count / file.header.point_count * 100))

        # For performance we need to use copy
        # so that the underlying arrays are contiguous
        x, y, z, gps = points.x.copy(), points.y.copy(), points.z.copy(), points.gps_time.copy()

        # call interpolate function for trajectory
        SensorPos = interpolate_traj(traj, gps)

        print("Do actual raytracing with all pulses")
        tic = time.time()
        RayTr.doRaytracing(x, y, z, SensorPos['sensor_x'], SensorPos['sensor_y'], SensorPos['sensor_z'], SensorPos['time'])
        toc = time.time()

        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        count += len(points)

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


