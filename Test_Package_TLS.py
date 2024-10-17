from occpy.OccPy import OccPy
from occpy import TerrainModel
import os
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description="OccPy")
parser.add_argument('root_folder', type=str, help="Path to the folder containing relevant files")
args = parser.parse_args()

root_folder = Path(args.root_folder)

test = OccPy(laz_in=os.path.join(root_folder, "LAZ"),
             out_dir=os.path.join(str(root_folder.parent.parent), "output"),
             vox_dim=0.1,
             lower_threshold=1,
             points_per_iter=100000000,
             plot_dim=[2676541,
                       1246160,
                       540,
                       2676591,
                       1246210,
                       615])


test.read_sensorpos_file(path2senspos=os.path.join(root_folder, "ScanPositions.txt"),
                         delimiter=" ",
                         hdr_scanpos_id="ScanPosID",
                         hdr_x="X",
                         hdr_y="Y",
                         hdr_z="Z",
                         str_idx_ScanPosID=23,
                         sens_pos_id_offset=1,  # This is a very specific use case, where we had an offset of scan pos id between laz files and scan pos id in the scanner position file.
                         single_return=True)

test.do_raytracing()

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

