from occpy.OccPy import OccPy
from occpy.OccPy import normalize_occlusion_output
from occpy.OccPy import get_Occl_TransectFigure
from occpy.OccPy import get_Occl_TransectFigure_BinaryOcclusion
from occpy.OccPy import get_Occlusion_ProfileFigure

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

test = OccPy(laz_in=f"{config['laz_in']}",
             out_dir=config['out_dir'],
             vox_dim=config['vox_dim'],
             lower_threshold=config['lower_threshold'],
             points_per_iter=config['points_per_iter'],
             plot_dim=config['plot_dim'],)

test.define_sensor_pos(path2file=config['ScanPos'],         # Path to the trajectory file
                       is_mobile=True,                      # whether acquisition is mobile. Always true for MLS or ULS
                       single_return=False,                 # wheter the data is single or multi return
                       delimiter=",",                       # delimiter used in the trajectory file 
                       hdr_time='Time[s]',                  # column header for the time information in the trajectory file
                       hdr_x='Easting[m]',                  # column header for the x coordinate in the trajectory file
                       hdr_y='Northing[m]',                 # column header for the y coordinate in the trajectory file
                       hdr_z='Height[m]',)                  # column header for the z coordinate in the trajectory file

test.do_raytracing()

Nhit_norm, Nmiss_norm, Nocc_norm, Classification_norm, chm = normalize_occlusion_output(input_folder=config['out_dir'],
                                                                                        PlotDim=config['plot_dim'],
                                                                                        vox_dim=config['vox_dim'],
                                                                                        dtm_file=config['tif_in']['DTM'],
                                                                                        dsm_file=config['tif_in']['DSM'],
                                                                                        lower_threshold=config['lower_threshold'],
                                                                                        output_voxels=config['output_voxels'])

# If instead of a binary voxel grid classification, you are intereded in a fraction of occluded pulses per voxel, you can use the following function:
OcclFrac_norm = Nocc_norm.astype(float) / (Nhit_norm.astype(float) + Nmiss_norm.astype(float) + Nocc_norm.astype(float))

# Test occlusion visualization
# first define figure properties
fig_prop = dict(fig_size=(3.5, 3.2),        # figure size in inches
                label_size=8,               # label size in pts. for e.g. axis labels
                label_size_ticks=6,         # label size in pts. for e.g. axis ticks
                label_size_tiny=5,          # label size in pts. for e.g. ticks in colorbar
                out_format='png',)          # output format of the figure, can be 'png', 'pdf', 'svg', etc.


get_Occl_TransectFigure_BinaryOcclusion(Nhit_norm, Classification_norm, plot_dim=config['plot_dim'], vox_dim=config['vox_dim'],
                                        out_dir=config['out_dir'], axis=0, start_ind=0, end_ind=100, chm=chm, vertBuffer=10, fig_prop=fig_prop, show_plots=True)

# Also get a profile figure of the transect
from occpy.OccPy import get_Occlusion_ProfileFigure

# define other figure properties for the smaller figure
fig_prop = dict(fig_size=(1.75, 3.2),
                label_size=8,
                label_size_ticks=6,
                label_size_tiny=5,
                out_format='png', )
get_Occlusion_ProfileFigure(Classification_norm, plot_dim=config['plot_dim'], vox_dim=config['vox_dim'], out_dir=config['out_dir'], low_thresh=0, vertBuffer=10, max_percentage=100, fig_prop=fig_prop, show_plots=True)

