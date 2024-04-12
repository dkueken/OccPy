from occpy.OccPy import OccPy
from occpy import TerrainModel


test = OccPy(laz_in=r'Z:\Data\Occlusion_TestData\OccPy_TestData_RamerenWald_FP05\Data\RamerenWald\UAVLS\20240410\20240410_1040_Ramerenwald_LFI_FP05_07_08_10_0_45_3.1_44.2lps_200pts_60m_90_clip_gS1_heb.laz',
             out_dir=r'Z:\Data\Occlusion_TestData\OccPy_TestData_RamerenWald_FP05\Data\RamerenWald\UAVLS\20240410\OccPyTest',
             vox_dim=0.1,
             lower_threshold=1,
             points_per_iter=1000000,
             plot_dim=[2676541,
                       1246160,
                       540,
                       2676591,
                       1246210,
                       615])

test.read_trajectory_file(path2traj=r'Z:\Data\Occlusion_TestData\OccPy_TestData_RamerenWald_FP05\Data\RamerenWald\UAVLS\20240410\trajectory_20240410_Ramerenwald_LFI_FP05_07_08_10_0°_combined.txt',
                          delimiter=',',
                          hdr_time='Time[s]', hdr_x='Easting[m]', hdr_y='Northing[m]', hdr_z='Height[m]')

test.do_raytracing()

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

