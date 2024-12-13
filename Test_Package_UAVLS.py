from occpy.OccPy import OccPy
from occpy import TerrainModel
import os
import argparse
from pathlib import Path
import json

parser = argparse.ArgumentParser(description="OccPy")
parser.add_argument('config_file', help="JSON configuration file of the OccPy command")
args = parser.parse_args()

def load_config(config_file):
    with open(config_file, 'r') as file:
        return json.load(file)

config = load_config(args.config_file)
root_folder = config['root_folder']

test = OccPy(laz_in=f"{root_folder}{config['laz_in']}",
             out_dir=f"{root_folder}{config['out_dir']}",
             vox_dim=config['vox_dim'],
             lower_threshold=config['lower_threshold'],
             points_per_iter=config['points_per_iter'],
             plot_dim=config['plot_dim'],
             output_voxels=config['output_voxels'])

test.read_trajectory_file(path2traj=f"{root_folder}{config['ScanPos']}",
                          delimiter=',',
                          hdr_time='Time[s]', hdr_x='Easting[m]', hdr_y='Northing[m]', hdr_z='Height[m]')

test.do_raytracing()

test.normalize_occlusion_output(input_folder=test.out_dir,
                                dtm_file=f"{root_folder}{config['tif_in']['DTM']}",
                                dsm_file=f"{root_folder}{config['tif_in']['DSM']}")

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

