from occpy.OccPy import OccPy
from occpy import TerrainModel
import os
import argparse
from pathlib import Path
import json

import time

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
             )

test.define_sensor_pos(path2file=os.path.join(root_folder, "OccPy_MLS_test_data_FP05_trajectory.txt"),
                       is_mobile=True,
                       single_return=True,
                       delimiter=" ",
                       hdr_time='%time',
                       hdr_x='x',
                       hdr_y='y',
                       hdr_z = 'z')

tic = time.time()
test.do_raytracing()
toc = time.time()
print(f"Elapsed time: {toc-tic} seconds")

#test.save_raytracing_output()

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

