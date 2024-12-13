import numpy as np
import os
import pandas as pd
import laspy
import time
from BOVwriter import writeBOV

# adding the path to the src folder in case raytr module cannot be found
import sys
sys.path.append(r"src")
from raytr import PyRaytracer

is_sorted = lambda a: np.all(a[:-1] <= a[1:])

# Define parameters
# for performance reasons it would be best that multi return datasets are sorted based on their GPS time and return number
# you could use LASTools lassort function for that: lassort -i laz_in -gps_time -return_number -odix _sort -olaz -v
laz_in = r"path_to_input_laz_file.laz"
out_dir = r"path_to_output_directory"

single_return = True    # set this to true, if your TLS data is single return data (e.g. FARO -> they will also not have
                        # the attributes gps_time, return_number and number_of_returns stored in the laz file...

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

points_per_iter = 10000000 # define how many points should be read from the laz file simultaneously

# Define the coordinates of the Scanner position
ScanPos = dict(X=2676566.32042296,
               Y=1246174.808677466,
               Z=553.1482458720556)

# Define the coordinates of the plot
PlotDim = dict(minX=ScanPos['X'] - 25,
               maxX=ScanPos['X'] + 25,
               minY=ScanPos['Y'] - 25,
               maxY=ScanPos['Y'] + 25,
               minZ=540,
               maxZ=585)

vox_dim = 0.1  # voxel dimension in meters

# Prepare Raytracer
gridDim = dict(nx=int((PlotDim['maxX'] - PlotDim['minX']) / vox_dim),
               ny=int((PlotDim['maxY'] - PlotDim['minY']) / vox_dim),
               nz=int((PlotDim['maxZ'] - PlotDim['minZ']) / vox_dim)
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
RayTr.defineGrid(minBound, maxBound, gridDim['nx'], gridDim['ny'], gridDim['nz'], vox_dim)
toc = time.time()
print("Time elapsed: {:.2f} seconds".format(toc - tic))


# read in laz file
tic = time.time()
print("Reading in point cloud data")
with laspy.open(laz_in) as file:
    count = 0
    for points in file.chunk_iterator(points_per_iteration=points_per_iter):
        print("{:.2f}%".format(count / file.header.point_count * 100))

        if single_return:
            sorted = True

            x = points.x.copy()
            y = points.y.copy()
            z = points.z.copy()

            sensor_x = np.ones(x.shape) * ScanPos['X']
            sensor_y = np.ones(x.shape) * ScanPos['Y']
            sensor_z = np.ones(x.shape) * ScanPos['Z']

            gps_time = np.linspace(start=count+1, stop=count+len(x), num=len(x), endpoint=True)

            # run raytracing algorithme using singleReturnPulses version
            print("Do raytracing with all pulses in batch")
            tic_r = time.time()
            RayTr.doRaytracing_singleReturnPulses(x, y, z, sensor_x, sensor_y,
                                                  sensor_z, gps_time)
            toc_r = time.time()
            print("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))


        else:
            x = points.x.copy()
            y = points.y.copy()
            z = points.z.copy()
            gps_time = points.gps_time.copy()
            return_number = points.return_number.copy()
            number_of_returns = points.number_of_returns.copy()

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


            sensor_x = np.ones(gps_time.shape) * ScanPos['X']
            sensor_y = np.ones(gps_time.shape) * ScanPos['Y']
            sensor_z = np.ones(gps_time.shape) * ScanPos['Z']

            RayTr.addPointData(x, y, z, sensor_x, sensor_y, sensor_z, gps_time, return_number, number_of_returns)

            if sorted: # only if pulses are sorted run raytracing now. Otherwise we have to read in the entire dataset first!
                # Get report on pulse dataset - comment this out once everythin is working or TODO: add a verbose flag!
                RayTr.getPulseDatasetReport()

                # run raytracing on added points
                print("Do raytracing with stored pulses")
                tic_r = time.time()
                RayTr.doRaytracing()
                toc_r = time.time()
                print("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))

                # Check if traversed pulses have been deleted from map - comment this out once everything is working or TODO: add a verbose flag!
                RayTr.getPulseDatasetReport()

        count = count + len(gps_time)


toc = time.time()
if sorted and not single_return:
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
elif not sorted and not single_return:
    print("Time elapsed for reading in data: {:.2f} seconds".format(toc-tic))

    RayTr.getPulseDatasetReport()

    print("Clean up pulse dataset in order to handle incomplete pulses")
    RayTr.cleanUpPulseDataset()

    RayTr.getPulseDatasetReport()

    print("Do actual raytracing with all pulses")
    tic = time.time()
    RayTr.doRaytracing()
    toc = time.time()
    print("Time elapsed for raytracing: {:.2f} seconds".format(toc-tic))

# Get report on traversal
RayTr.reportOnTraversal()


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






