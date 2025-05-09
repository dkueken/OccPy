from occpy.OccPy import OccPy
from occpy import TerrainModel

import os.path

import numpy as np
import pandas as pd

import json
import glob
import rasterio

import subprocess

import time

### Input parameters
in_dir = r'E:\Data\BEF_Lagern\MLS\2024_02_21_LeafOff'
out_dir = r'E:\Data\BEF_Lagern\MLS\2024_02_21_LeafOff'
uls_dir = r'E:\Data\BEF_Lagern\UAVLS'
occpy_dir = r"C:\Users\Kueken\dev\occPy"
focal_tree_coords_file = r'E:\Data\BEF_Lagern\FocalTreeCoords\BEFLagern_FocalTree_Coords.csv'
plot_corner_coords_file = r'E:\Data\BEF_Lagern\PlotCoordiantes\CornerCoordinates_FieldAssessment.csv'

json_template = r"C:\Users\Kueken\dev\occPy\settings_MLS.JSON"
# load config template
with open(json_template, 'r') as file:
    config_template = json.load(file)

PLOT_CENTERS = pd.read_csv(focal_tree_coords_file)
plot_corner_coords = pd.read_csv(plot_corner_coords_file)

parameters = dict(vox_size = 0.1,
                  plot_buffer_if = 0)

# get all plots
plots = glob.glob(f"{in_dir}\\BEF*")

for plot in plots:
    last_dir_ind = plot.rfind('\\')
    plot_id = plot[last_dir_ind + 1:]
    print(plot_id)

    #if plot_id == "BEF_31A":
    #    print(f"!!!! There is an issue with processing {plot_id} for LeafOn... skipping for the moment. TODO: check input laz file!")
    #    continue

    if not os.path.exists(f"{plot}\\OcclusionMap"):
        os.makedirs(f"{plot}\\OcclusionMap")
        os.makedirs(f"{plot}\\OcclusionMap\\config")

    if os.path.exists(f"{plot}\\OcclusionMap\\Classification.npy"):
        print(f"{plot_id} has already been processd... skipping")
        continue

    minx = np.floor(PLOT_CENTERS.loc[PLOT_CENTERS['ID'] == plot_id, 'X'].values[0] - (
            25 + parameters['plot_buffer_if']))
    maxx = np.floor(PLOT_CENTERS.loc[PLOT_CENTERS['ID'] == plot_id, 'X'].values[0] + (
            25 + parameters['plot_buffer_if']))
    miny = np.floor(PLOT_CENTERS.loc[PLOT_CENTERS['ID'] == plot_id, 'Y'].values[0] - (
            25 + parameters['plot_buffer_if']))
    maxy = np.floor(PLOT_CENTERS.loc[PLOT_CENTERS['ID'] == plot_id, 'Y'].values[0] + (
            25 + parameters['plot_buffer_if']))

    # read in DTM and DSM to get min and max Z coordinates
    dtm_file = f"{uls_dir}\\{plot_id}\\Grids\\{plot_id}_DTM_ULS_clip2IF.tif"
    if not os.path.exists(dtm_file):
        dtm_file = f"{plot}\\Grids\\{plot_id}_DTM_MLS_clip2IF.tif"

    with rasterio.open(dtm_file) as dtm_f:
        dtm = dtm_f.read(1)
        dtm_mask_na = np.ma.masked_where(dtm == -9999,
                                         dtm,
                                         copy=True)

        minz = np.floor(dtm_mask_na.min()-1)
    dsm_file = f"{uls_dir}\\{plot_id}\\Grids\\{plot_id}_DSM_ULS_clip2IF.tif"
    if not os.path.exists(dsm_file):
        dsm_file = f"{plot}\\Grids\\{plot_id}_DSM_MLS_clip2IF.tif"

    with rasterio.open(dsm_file) as dsm_f:
        dsm = dsm_f.read(1)

        maxz = np.floor(dsm.max() + 5)

    # write config file
    cfg = config_template
    cfg['root_folder'] = r"C:\Users\Kueken\dev\occPy"
    cfg['laz_in'] = glob.glob(f"{plot}\\LAZ\\{plot_id}_rot2LV95.laz")[0] # we expect only one file here
    cfg['tif_in']['DTM'] = dtm_file
    cfg['tif_in']['DSM'] = dsm_file
    cfg['out_dir'] = f"{plot}\\OcclusionMap"
    cfg['vox_dim'] = float(parameters['vox_size'])
    cfg['lower_threshold'] = int(1)
    cfg['points_per_iter'] = int(1000000)
    cfg['plot_dim'] = [float(minx), float(miny), float(minz), float(maxx), float(maxy), float(maxz)]
    cfg['ScanPos'] = glob.glob(f"{plot}\\Trajectories\\{plot_id}_traj_rot2LV95.txt")[0]

    # write cfg to json file
    cfg_out = f"{plot}\\OcclusionMap\\config\\{plot_id}_occpy_cfg.json"
    with open(cfg_out, 'w') as json_out:
        json.dump(cfg, json_out, indent=4)

    occpy = OccPy(laz_in=f"{cfg['laz_in']}",
                 out_dir=f"{cfg['out_dir']}",
                 vox_dim=cfg['vox_dim'],
                 lower_threshold=cfg['lower_threshold'],
                 points_per_iter=cfg['points_per_iter'],
                 plot_dim=cfg['plot_dim'],
                 )

    # check whether we have %time or //world_time as time header in ScanPos file
    traj = pd.read_csv(cfg['ScanPos'], delimiter=" ")

    occpy.define_sensor_pos(path2file=cfg['ScanPos'],
                           is_mobile=True,
                           single_return=True,
                           delimiter=" ",
                           hdr_time=traj.columns[0],
                           hdr_x='x',
                           hdr_y='y',
                           hdr_z='z')

    tic = time.time()
    occpy.do_raytracing()
    toc = time.time()
    print(f"Elapsed time: {toc - tic} seconds")

    # test.save_raytracing_output()

    occpy.normalize_occlusion_output(input_folder=occpy.out_dir,
                                    dtm_file=f"{cfg['tif_in']['DTM']}",
                                    dsm_file=f"{cfg['tif_in']['DSM']}")

    # Get some occlusion statistics
    print(f"Total canopy volume of the plot: {occpy.TotalVolume * (occpy.vox_dim ** 3)} m3")
    print(f"Total occluded volume of the plot: {occpy.TotalOcclusion * (occpy.vox_dim ** 3)} m3")
    print(f"Average occlusion fraction: {occpy.OcclFrac2D.mean()}")
    print(f"Max occlusion fraction: {occpy.OcclFrac2D.max()}")















