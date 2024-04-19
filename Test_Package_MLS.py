from occpy.OccPy import OccPy
from occpy import TerrainModel


test = OccPy(laz_in=r'Z:\Data\Occlusion_TestData\OccPy_TestData_RamerenWald_FP05\Data\RamerenWald\MLS\ZebHorizon\FP05\LAZ\Rameren_FP05_2022-06-01_07-32-05_100pct_height_world_rot2LV95.laz',
             out_dir=r'Z:\Data\Occlusion_TestData\OccPy_TestData_RamerenWald_FP05\Data\RamerenWald\MLS\ZebHorizon\FP05\OccPy_Test',
             vox_dim=0.1,
             lower_threshold=1,
             points_per_iter=1000000,
             plot_dim=[2676541,
                       1246160,
                       540,
                       2676591,
                       1246210,
                       615])

test.read_trajectory_file(path2traj=r'Z:\Data\Occlusion_TestData\OccPy_TestData_RamerenWald_FP05\Data\RamerenWald\MLS\ZebHorizon\FP05\Trajectories\Rameren_FP05_2022-06-01_07-32-05_results_traj_rot2LV95.txt')

test.do_raytracing()

test.save_raytracing_output()

test.normalize_occlusion_output(input_folder=test.out_dir,
                                dtm_file=r'Z:\Data\Occlusion_TestData\OccPy_TestData_RamerenWald_FP05\Data\RamerenWald\DTM_FP05_swissAlti3D_10cm.tif',
                                dsm_file=r'Z:\Data\Occlusion_TestData\OccPy_TestData_RamerenWald_FP05\Data\RamerenWald\DSM_FP05_swissSurface3D_10cm.tif')

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

