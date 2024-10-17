from occpy.OccPy import OccPy
from occpy import TerrainModel
import argparse
import os
from pathlib import Path

parser = argparse.ArgumentParser(description="OccPy")
parser.add_argument('root_folder', type=str, help="Path to the folder containing relevant files")
args = parser.parse_args()

root_folder = Path(args.root_folder)

test = OccPy(laz_in= os.path.join(root_folder, "Occpy_MLS_test_data.laz"),
             out_dir=os.path.join(str(root_folder.parent.parent), "output"),
             vox_dim=0.1,
             lower_threshold=1,
             points_per_iter=1000000,
             plot_dim=[2676541,
                       1246160,
                       540,
                       2676591,
                       1246210,
                       615])

test.read_trajectory_file(path2traj=os.path.join(root_folder, "Rameren_FP05_2022-06-01_07-32-05_results_traj_rot2LV95.txt"))

test.do_raytracing()

test.save_raytracing_output()

test.normalize_occlusion_output(input_folder=test.out_dir,
                                dtm_file=os.path.join(str(root_folder.parent), "DTM_FP05_swissAlti3D_10cm.tif"),
                                dsm_file=os.path.join(str(root_folder.parent), "DSM_FP05_swissSurface3D_10cm.tif"))

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

