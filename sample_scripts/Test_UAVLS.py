import os

import numpy as np
import pickle
import pandas as pd
import laspy
import time
from occpy.OccPy import interpolate_traj

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

from raytr import PyRaytracer

is_sorted = lambda a: np.all(a[:-1] <= a[1:])

# convenience functions for testing purposes. not actually needed!
# function to check if points are collinear
def are_points_collinear(points, threshold):
    # Ensure at least 3 points are provided
    if len(points) < 3:
        return False

    # Calculate vecotrs between consecutive points
    vectors = []
    for i in range(len(points) - 1):
        vector = [points[i+1][j] - points[i][j] for j in range(3)]
        vectors.append(vector)

    # Calculate cross products between consecutive vectors
    cross_products = []
    for i in range(len(vectors) - 1):
        cross_product = [
            vectors[i][1] * vectors[i+1][2] - vectors[i][2] * vectors[i+1][1],
            vectors[i][2] * vectors[i+1][0] - vectors[i][0] * vectors[i+1][2],
            vectors[i][0] * vectors[i+1][1] - vectors[i][1] * vectors[i+1][0]
        ]
        cross_products.append(cross_product)

    # Check if all cross products are zero
    return all(
        abs(cross_product[0]) <= threshold and
        abs(cross_product[1]) <= threshold and
        abs(cross_product[2]) <= threshold
        for cross_product in cross_products
    )

def move_point_to_collinear(points, selected_index, threshold=1e-6, max_iterations=100):
    # Ensure at least 3 points are provided
    if len(points) < 3:
        return points

    iterations = 0
    while True:
        # Calculate the average position of the other points in the first and second dimensions
        avg_position = np.mean(points[:selected_index] + points[selected_index+1:], axis=0)[:2]


        # Calculate the direction vector from the selected point to the average position in the first and second dimensions
        direction = avg_position - points[selected_index][:2]

        # Normalize the direction vector
        direction /= np.linalg.norm(direction)
        print(direction)

        # Define the desired distance to move the point in the first and second dimensions (adjust this as needed)
        # distance = 100

        # Calculate the new position of the selected point in the first and second dimensions
        new_position = points[selected_index][:2] + direction

        # Create a new updated points list with the selected point's position modified in the first and second dimensions
        updated_points = points.copy()
        updated_points[selected_index] = (new_position[0], new_position[1], points[selected_index][2])

        # check if the updated points are collinear
        if are_points_collinear(updated_points, threshold):
            print("Collinearity reached!")
            return updated_points

        # Exit the loop if the maximum number of iterations is reached
        iterations += 1
        if iterations >= max_iterations:
            print("Maximum iterations has been reached!")
            break

        points[selected_index] = (new_position[0], new_position[1], points[selected_index][2])

    return points

def move_point_to_line(point, line_points):
    # Convert the points to Numpy arrays for easier calculations
    point = np.array(point)
    line_points = np.array(line_points)

    # Calculate the line direction vector
    line_direction = line_points[-1] - line_points[0]

    # Calculate the point's projection onto the line
    t = np.dot(point - line_points[0], line_direction) / np.dot(line_direction, line_direction)
    projected_point = line_points[0] + t * line_direction

    return projected_point


# Define parameters
laz_in = r"path_to_input_laz_file.laz"

# If you have multiple trajectory files, define their paths here and combine them using pandas as seen below.
traj_file_1 = r"path_to_trajectory_file.txt"
traj_file_2 = r""

out_dir = r"path_to_output_directory"
os.makedirs(os.path.dirname(out_dir), exist_ok=True)

parameters = dict(
    voxDim=0.1,
    lower_threshold=1,
    points_per_iter=10000000,
)

# read in trajectory file
traj1 = pd.read_csv(traj_file_1, sep=",")
if len(traj_file_2)!=0:
    traj2 = pd.read_csv(traj_file_2, sep=",")

    # combine the two trajectories
    traj = pd.concat([traj1, traj2], ignore_index=True)
else:
    traj = traj1

# Define plot dimensions
# Plot Dim FP05
"""
PlotDim = dict(minX=2676541,
               maxX=2676591,
               minY=1246160,
               maxY=1246210,
               minZ=540,
               maxZ=615)
"""
# Plot Dim FP08
PlotDim = dict(minX=2676489,
               maxX=2676539,
               minY=1246115,
               maxY=1246165,
               minZ=540,
               maxZ=590)

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

# read in laz file
tic = time.time()
print("Reading in point cloud data")
with laspy.open(laz_in) as file:
    count = 0
    for points in file.chunk_iterator(points_per_iteration=parameters['points_per_iter']):
        print("{:.2f}%".format(count / file.header.point_count * 100))

        x = points.x.copy()
        y = points.y.copy()
        z = points.z.copy()
        gps_time = points.gps_time.copy()
        return_number = points.return_number.copy()
        number_of_returns = points.number_of_returns.copy()

        # call interpolate function for trajectory to extract sensor position for each gps_time
        SensorPos = interpolate_traj(traj['Time[s]'], traj['Easting[m]'], traj['Northing[m]'], traj['Height[m]'], gps_time)

        # check if gps_time is sorted
        if count==0:  # only check sort state in the first iteration of the for loop
            if not is_sorted(gps_time):
                print(f"!!!!! input laz file is not sorted along gps_time. The algorithm will still run. However, the "
                      f"performance will be greatly decreased as the entire content of the laz file has to be read into "
                      f"the system memory. If you have multi return data, consider sorting your laz data first, e.g. using "
                      f"LASTools lassort: lassort -i laz_in -gps_time -return_number -odix _sort -olaz -v !!!!")
                sorted = False
            else:
                sorted = True

        RayTr.addPointData(x, y, z, SensorPos['sensor_x'], SensorPos['sensor_y'], SensorPos['sensor_z'], gps_time, return_number, number_of_returns)

        if sorted: # only if pulses are sorted run raytracing now. Otherwise we have to read in the entire dataset first!
            # Get report on pulse dataset - comment this out once everythin is working or TODO: add a verbose flag!
            # RayTr.getPulseDatasetReport()

            # For UAVLS data, interpolated Sensor Positions can lead to sensor positions, which are not collinear with the
            # laser returns. This can lead to not-traversed returns during voxel traversal (return lies in another voxel
            # than the traversal vector). For laser pulses with two or more returns we can move the Sensor position so it is
            # collinear with its returns. TODO: Check on the implications for doing that! What are the added uncertainties?
            # TODO: Check on how far we move the Sensor positions!

            RayTr.moveSensorPos2Collinearity()
            # get report on sensor shifts
            sensor_shifts_tmp = RayTr.reportSensorShifts()
            sensor_shifts_tmp = np.array(sensor_shifts_tmp, float)
            if count==0:
                SensorShifts = sensor_shifts_tmp
            else:
                SensorShifts = np.vstack((SensorShifts, sensor_shifts_tmp))

            # run raytracing on added points
            print("Do raytracing with stored pulses")
            tic_r = time.time()
            RayTr.doRaytracing()
            toc_r = time.time()
            print("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))

            # clear up pulsedataset in RayTr object
            RayTr.clearPulseDataset()

            # Check if traversed pulses have been deleted from map - comment this out once everything is working or TODO: add a verbose flag!
            # RayTr.getPulseDatasetReport()


        count = count + len(gps_time)


toc = time.time()
if sorted:
    # optional: incomplete pulses can occur if the data has been filtered (either actively or during black box processing
    # of the processing software. We could actively turn the incomplete pulses into complete ones and do the raytracing
    # for them!
    print("convert incomplete pulses to complete ones - be cautious with that!")
    RayTr.getPulseDatasetReport()
    RayTr.cleanUpPulseDataset()
    RayTr.getPulseDatasetReport()
    print("Run raytracing for incomplete pulses")
    tic_r = time.time()
    RayTr.doRaytracing()
    toc_r = time.time()
    print("Time elapsed for raytracing incomplete pulses: {:.2f} seconds".format(toc_r - tic_r))
    print("Time elapsed for reading and raytracing entire data: {:.2f} seconds".format(toc-tic))
else:
    print("Time elapsed for reading in data: {:.2f} seconds".format(toc-tic))

    RayTr.getPulseDatasetReport()

    print("Clean up pulse dataset in order to handle incomplete pulses")
    RayTr.cleanUpPulseDataset()

    RayTr.getPulseDatasetReport()

    RayTr.moveSensorPos2Collinearity()
    # get report on sensor shifts
    sensor_shifts_tmp = RayTr.reportSensorShifts()
    SensorShifts = np.array(sensor_shifts_tmp, float)

    print("Do actual raytracing with all pulses")
    tic = time.time()
    RayTr.doRaytracing()
    toc = time.time()
    print("Time elapsed for raytracing: {:.2f} seconds".format(toc-tic))


# Get report on traversal
RayTr.reportOnTraversal()

# save data to sensor shifts
np.save(f"{out_dir}/SensorShifts.npy", SensorShifts)


print(f"### Raytracing is complete - Extracting 3D voxel grid")
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



