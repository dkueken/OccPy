import os

import numpy as np
import pickle
import pandas as pd
import laspy
import time
from prepare_trajectory import interpolate_traj

from BOVwriter import writeBOV

# we need to add the path to the src folder so raytr module can be found
import sys
sys.path.append(r".\src")

from raytr import PyRaytracer



# Input parameters
laz_in = r"D:\data\occpy\Rameren_FP05_2022-04-27_14-32-22_100pct_height_world_rot2LV95.laz"
traj_in = r"D:\data\occpy\Rameren_FP05_2022-04-27_14-32-22_results_traj_rot_LV95.txt"

out_dir = r"D:\data\occpy\output"
os.makedirs(os.path.dirname(out_dir), exist_ok=True)

parameters = dict(
    voxDim=0.1,
    lower_threshold=1,
    points_per_iter=10000000,
)

# read in trajectory file
traj = pd.read_csv(traj_in, sep=";")

# Plot Dimensions - Plot Dim FP05
PlotDim = dict(minX=2676541,
               maxX=2676591,
               minY=1246160,
               maxY=1246210,
               minZ=540,
               maxZ=615)

gridDim = dict(nx=int((PlotDim['maxX'] - PlotDim['minX'])/parameters['voxDim']),
               ny=int((PlotDim['maxY'] - PlotDim['minY'])/parameters['voxDim']),
               nz=int((PlotDim['maxZ'] - PlotDim['minZ'])/parameters['voxDim'])
               )

# initialize output grids
Nhit = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)
Nmiss = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)
Nocc = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)

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

starttime = time.time()

# read in laz iteratively and run raytracing tool
print("Reading in LAZ data")
tic_tot = time.time()
with laspy.open(laz_in) as file:
    count = 0
    for points in file.chunk_iterator(parameters['points_per_iter']):
        print("{:.2f}%".format(count / file.header.point_count * 100))

        # For performance we need to use copy
        # so that the underlying arrays are contiguous
        x = points.x.copy()
        y = points.y.copy()
        z = points.z.copy()
        gps_time = points.gps_time.copy()
        return_number = points.return_number.copy()
        number_of_returns = points.number_of_returns.copy()

        if np.max(return_number)==0:  # a not very nice hack for the special case where return_number and number_of_returns are all 0 for Horizon measurements - TODO: figure out why!
            return_number[:] = 1
            number_of_returns[:] = 1

        # call interpolate function for trajectory to extract sensor position for each gps_time
        SensorPos = interpolate_traj(traj['X.time'], traj['x'], traj['y'], traj['z'], gps_time)

        RayTr.doRaytracing_singleReturnPulses(x, y, z, SensorPos['sensor_x'], SensorPos['sensor_y'], SensorPos['sensor_z'], gps_time)

        count = count + len(gps_time)

print("RAYTRACING Elapsed Time: " + str(totaltime) + " seconds")

# Get report on traversal
RayTr.reportOnTraversal()

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

print("Saving Occlusion Outputs")
tic = time.time()
np.save(f"{out_dir}/Nhit.npy", Nhit)
np.save(f"{out_dir}/Nmiss.npy", Nmiss)
np.save(f"{out_dir}/Nocc.npy", Nocc)
toc = time.time()
print("Elapsed Time: {:.2f} seconds".format(toc - tic))

# Create Classification grid
print("Classify Grid")
tic = time.time()
Classification = np.zeros((gridDim['ny'], gridDim['nx'], gridDim['nz']), dtype=int)

Classification[np.logical_and.reduce((Nhit>0, Nmiss>=0, Nocc>=0))] = 1 # voxels that were observed
Classification[np.logical_and.reduce((Nhit==0, Nmiss>0, Nocc>=0))] = 2 # voxels that are empty
Classification[np.logical_and.reduce((Nhit==0, Nmiss==0, Nocc>0))] = 3 # voxels that are hidden (occluded)
Classification[np.logical_and.reduce((Nhit==0, Nmiss==0, Nocc==0))] = 4 # voxels that were not observed # TODO: Figure out, why this overwrites voxels that are classified as occluded! -> this was because np.logical_and only takes in 2 arrays as input, not 3! use np.logical_and.reduce() for that!

np.save(f"{out_dir}/Classification.npy", Classification)
toc = time.time()
print("Elapsed Time: " + str(toc - tic) + " seconds")


# write BOV file -> Actually VISIT can read in npy files. However, axis definition is wrong.
print("Writing BOV File")
tic = time.time()
writeBOV(out_dir + '\\', "Classification", "Classification", 'i', Classification)
writeBOV(out_dir + '\\', "Nhit", "Nhit", 'i', Nhit)
writeBOV(out_dir + '\\', "Nmiss", "Nmiss", 'i', Nmiss)
writeBOV(out_dir + '\\', "Nocc", "Nocc", 'i', Nocc)
toc = time.time()
print("Elapsed Time: " + str(toc - tic) + " seconds")

totaltime = time.time() - starttime

