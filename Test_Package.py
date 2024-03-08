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

#test.read_trajectory_file(path2traj=r'D:\_tmp_wdir\TestData_occpy\MLS\Rameren_FP05_2022-04-06_07-39-44_results_traj_rot2LV95.txt')

#test.do_raytracing()

#test.save_raytracing_output()

test.normalize_occlusion_output(input_folder=r'D:\_tmp_wdir\TestData_occpy\MLS\OccPy_Test', dtm_file=r'D:\_tmp_wdir\TestData_occpy\MLS\DTM_FP05_swissAlti3D_10cm.tif', dsm_file=r'D:\_tmp_wdir\TestData_occpy\MLS\DSM_FP05_swissSurface3D_10cm.tif')

# Get some occlusion statistics
print(f"Total canopy volume of the plot: {test.TotalVolume * (test.vox_dim**3)} m3")
print(f"Total occluded volume of the plot: {test.TotalOcclusion * (test.vox_dim**3)} m3")
print(f"Average occlusion fraction: {test.OcclFrac2D.mean()}")
print(f"Max occlusion fraction: {test.OcclFrac2D.max()}")

# Test occlusion visualization
test.get_Occl_TransectFigure(start_ind=200, end_ind=300, axis=0)
test.get_Occl_TransectFigure(start_ind=200, end_ind=300, axis=1)
test.get_Occl_TransectFigure(start_ind=50, end_ind=100, axis=2)

occl_prof = test.get_Occlusion_Profile()

print(occl_prof)

