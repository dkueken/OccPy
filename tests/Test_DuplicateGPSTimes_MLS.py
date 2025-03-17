import numpy as np
import pandas as pd
import laspy
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('TkAgg')

laz_in = r"D:\_tmp_wdir\OcclusionMappingTests\Rameren\Horizon\FP05_2022\Rameren_FP05_2022-06-01_07-32-05_100pct_height_world_rot2LV95_new.laz"

with laspy.open(laz_in) as file:
    count = 0
    for points in file.chunk_iterator(points_per_iteration=10000000):
        print("{:.2f}%".format(count / file.header.point_count * 100))

        if count==0:
            x = points.x.copy()
            y = points.y.copy()
            z = points.z.copy()
            gps_time = points.gps_time.copy()
            return_number = points.return_number.copy()
            number_of_returns = points.number_of_returns.copy()

        else:
            x = np.hstack((x, points.x.copy()))
            y = np.hstack((y, points.y.copy()))
            z = np.hstack((z, points.z.copy()))
            gps_time = np.hstack((gps_time, points.gps_time.copy()))
            return_number = np.hstack((return_number, points.return_number.copy()))
            number_of_returns = np.hstack((number_of_returns, points.number_of_returns.copy()))

        count += len(points)


# get unique gps_times in sorted array
gps_time_un, un_indices, un_count = np.unique(gps_time, return_index=True, return_counts=True)

max_count = np.max(un_count)

gps_time_maxret = gps_time_un[un_count==max_count]

# create a (apparently) pulse subset
x_sub = x[gps_time==gps_time_maxret[0]]
y_sub = y[gps_time==gps_time_maxret[0]]
z_sub = z[gps_time==gps_time_maxret[0]]

# plot this point subset with the same GPS Time
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(x_sub, y_sub, z_sub, marker='x')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()

print('test')


