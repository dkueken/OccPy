from occpy.OccPy import OccPy

import numpy as np


laz_in = r".\data\MLS\OccPy_MLS_test_data.laz"
traj_in = r".\data\MLS\Rameren_FP05_2022-06-01_07-32-05_results_traj_rot2LV95.txt"

traj = OccPy.read_trajectory_file(path2traj=traj_in)

minbounds = np.array([1246170, 2676560, 550])
maxbounds = np.array([1246180, 2676570, 610])


pulses_intersecting = OccPy.filterPointsIntersectingBox(laz_in=laz_in, laz_out=f"{laz_in[:-4]}_filteredBox.laz", min_bound=minbounds, max_bound=maxbounds, traj_in=traj, points_per_iter=100000)


