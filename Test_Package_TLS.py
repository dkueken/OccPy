from occpy.OccPy import OccPy
from occpy import TerrainModel


test = OccPy(laz_in=r'Z:\Data\Occlusion_TestData\OccPy_TestData_RamerenWald_FP05_shared\Data\RamerenWald\TLS\Leica_BLK360_Oct2020\LAZ',
             out_dir=r'Z:\Data\Occlusion_TestData\OccPy_TestData_RamerenWald_FP05_shared\Data\RamerenWald\TLS\Leica_BLK360_Oct2020\OccPy_Test',
             vox_dim=0.5,
             lower_threshold=1,
             points_per_iter=100000000,
             plot_dim=[2676541,
                       1246160,
                       540,
                       2676591,
                       1246210,
                       615])

test.read_trajectory_file(path2traj=r'D:\_tmp_wdir\TestData_occpy\MLS\Rameren_FP05_2022-04-06_07-39-44_results_traj_rot2LV95.txt')
test.read_sensorpos_file(path2senspos=r'Z:\Data\Occlusion_TestData\OccPy_TestData_RamerenWald_FP05\Data\RamerenWald\TLS\Leica_BLK360_Oct2020\ScanPositions.txt',
                         delimiter=" ",
                         hdr_scanpos_id="ScanPosID",
                         hdr_x="X",
                         hdr_y="Y",
                         hdr_z="Z",
                         str_idx_ScanPosID=23,
                         sens_pos_id_offset=1,
                         single_return=True) # This is a very specific use case, where we had an offset of scan pos id between laz files and scan pos id in the scanner position file.

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

