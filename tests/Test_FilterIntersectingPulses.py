from occpy.util import filterPointsIntersectingBox, read_trajectory_file

import numpy as np


laz_in = rf".\data\MLS\MLS_TestData_20perc_FP10_2025.laz"
traj_in = r".\data\MLS\MLS_TestData_traj_FP10_2025.txt"

traj = read_trajectory_file(path2traj=traj_in, delimiter=" ", hdr_time="//world_time", hdr_x="x", hdr_y="y", hdr_z="z")


minbounds = np.array([1246170, 2676560, 550])
maxbounds = np.array([1246180, 2676570, 610])


pulses_intersecting = filterPointsIntersectingBox(laz_in=laz_in, laz_out=f"{laz_in[:-4]}_filteredBox.laz", min_bound=minbounds, max_bound=maxbounds, traj_in=traj, points_per_iter=100000)

print(pulses_intersecting)
