from occpy.OccPy import OccPy
from occpy import TerrainModel


test = OccPy(laz_in=r'D:\_tmp_wdir\TestData_occpy\MLS\Rameren_FP05_2022-04-06_07-39-44_100pct_height_world_rot2LV95.laz',
             out_dir=r'D:\_tmp_wdir\TestData_occpy\MLS\OccPy_Test',
             vox_dim=0.1,
             lower_threshold=1,
             points_per_iter=1000000,
             plot_dim=[2676541,
                       1246160,
                       540,
                       2676591,
                       1246210,
                       615])

test.read_trajectory_file(path2traj=r'D:\_tmp_wdir\TestData_occpy\MLS\Rameren_FP05_2022-04-06_07-39-44_results_traj_rot2LV95.txt')

test.do_raytracing()

