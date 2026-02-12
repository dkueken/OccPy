from occpy.OccPy import OccPy
from occpy import TerrainModel
from occpy.util import normalize_occlusion_output

import os.path

import numpy as np
import pandas as pd

import json
import glob
import rasterio

import shutil

import subprocess

import time

def copy_file(src, dst):
    shutil.copyfile(src, dst)


def get_data(local_dir, main_dir, plot_id):
    # Create directory for local directory
    print(f"### Getting data for {plot_id}")
    os.makedirs(f"{local_dir}\\{plot_id}", exist_ok=True)

    # get path to input laz file
    laz_in = glob.glob(f"{main_dir}\\{plot_id}\\LAZ\\*_heb_sort.laz")[
        0]  # we assume there is only one file that fits this condition

    # get path to input traj file
    traj_in = glob.glob(f"{main_dir}\\{plot_id}\\Trajectories\\*.txt")[
        0]  # we assume there is only one file that fits this condition

    # copy laz file to local
    copy_file(src=laz_in, dst=f"{local_dir}\\{plot_id}\\{os.path.basename(laz_in)}")
    copy_file(src=traj_in, dst=f"{local_dir}\\{plot_id}\\{os.path.basename(traj_in)}")

def cleanup_output_dir(dir):
    """
    function to delete all the contents in the defined directory
    """
    for item in os.listdir(dir):
        item_path = os.path.join(dir, item)
        if os.path.isfile(item_path) or os.path.islink(item_path):
            os.unlink(item_path)
        elif os.path.isdir(item_path):
            shutil.rmtree(item_path)


if __name__ == "__main__":

    ### Input parameters
    in_dir = r'\\speedy16-36\data_15\_PROJEKTE\20231211_1058_BEF_Laegern\_analysis\UAVLS'
    out_dir = r'\\speedy16-36\data_15\_PROJEKTE\20231211_1058_BEF_Laegern\_analysis\UAVLS'
    local_dir = r'C:\Users\kueken\dev\occpy_data'
    uls_dir = r'\\speedy16-36\data_15\_PROJEKTE\20231211_1058_BEF_Laegern\_analysis\UAVLS'
    #occpy_dir = r"C:\Users\Kueken\dev\occPy"
    focal_tree_coords_file = r'\\speedy16-36\data_15\_PROJEKTE\20231211_1058_BEF_Laegern\_analysis\FocalTreeCoords\BEFLagern_FocalTree_Coords.csv'
    plot_corner_coords_file = r'\\speedy16-36\data_15\_PROJEKTE\20231211_1058_BEF_Laegern\_analysis\PlotCoordiantes\CornerCoordinates_FieldAssessment.csv'

    overwrite = False
    cleanup = False

    json_template = r"C:\Users\kueken\dev\OccPy\config\settings_ULS.JSON"
    # load config template
    with open(json_template, 'r') as file:
        config_template = json.load(file)

    PLOT_CENTERS = pd.read_csv(focal_tree_coords_file)
    plot_corner_coords = pd.read_csv(plot_corner_coords_file)

    parameters = dict(vox_size = 0.05,
                      plot_buffer_if = 0)

    # get all plots
    plots = glob.glob(f"{in_dir}\\BEF*")

    for plot in plots:
        last_dir_ind = plot.rfind('\\')
        plot_id = plot[last_dir_ind + 1:]
        print(plot_id)

        if os.path.exists(f"{plot}\\OcclusionMap_{int(parameters['vox_size']*100)}cm\\Classification.npy") and cleanup:
            cleanup_output_dir(f"{plot}\\OcclusionMap_{int(parameters['vox_size']*100)}cm")
            os.makedirs(f"{plot}\\OcclusionMap_{int(parameters['vox_size']*100)}cm\\config")

        if not os.path.exists(f"{plot}\\OcclusionMap_{int(parameters['vox_size']*100)}cm"):
            os.makedirs(f"{plot}\\OcclusionMap_{int(parameters['vox_size']*100)}cm")
            os.makedirs(f"{plot}\\OcclusionMap_{int(parameters['vox_size']*100)}cm\\config")

        if os.path.exists(f"{plot}\\OcclusionMap_{int(parameters['vox_size']*100)}cm\\Classification_{int(parameters['vox_size']*100)}cm.npy") and not overwrite:
            print(f"{plot_id} has already been processd... skipping")
            continue

        # create local directory to work in and copy necessary data
        if not os.path.exists(os.path.join(local_dir, plot_id)):
            os.makedirs(os.path.join(local_dir, plot_id, "LAZ"))
            os.makedirs(os.path.join(local_dir, plot_id, "Trajectories"))
            os.makedirs(os.path.join(local_dir, plot_id, f"OcclusionMap_{int(parameters['vox_size'] * 100)}cm"))
            os.makedirs(
                os.path.join(local_dir, plot_id, f"OcclusionMap_{int(parameters['vox_size'] * 100)}cm", "config"))

        get_data(local_dir, in_dir, plot_id)


        """# first sort laz file according to gps time - not needed anymore!
        if not glob.glob(f"{plot}\\LAZ\\*_clip_gs1_heb_sort.laz"):
            laz_in = glob.glob(f"{plot}\\LAZ\\*_clip_gS1_heb.laz")[0]  # we expect only one file here
            cmd = f"lassort -i {laz_in} -gps_time -return_number -odix _sort -olaz -cpu64 -v"
            subprocess.run(cmd)

        """
        minx = np.floor(PLOT_CENTERS.loc[PLOT_CENTERS['ID'] == plot_id, 'X'].values[0] - (
                25 + parameters['plot_buffer_if']))
        maxx = np.floor(PLOT_CENTERS.loc[PLOT_CENTERS['ID'] == plot_id, 'X'].values[0] + (
                25 + parameters['plot_buffer_if']))
        miny = np.floor(PLOT_CENTERS.loc[PLOT_CENTERS['ID'] == plot_id, 'Y'].values[0] - (
                25 + parameters['plot_buffer_if']))
        maxy = np.floor(PLOT_CENTERS.loc[PLOT_CENTERS['ID'] == plot_id, 'Y'].values[0] + (
                25 + parameters['plot_buffer_if']))

        # read in DTM and DSM to get min and max Z coordinates
        dtm_file = f"{uls_dir}\\{plot_id}\\Grids\\{plot_id}_DTM_miniVUX.tif"
        if not os.path.exists(dtm_file):
            dtm_file = f"{plot}\\Grids\\{plot_id}_DTM_MLS.tif"

        with rasterio.open(dtm_file) as dtm_f:
            dtm = dtm_f.read(1)
            dtm_mask_na = np.ma.masked_where(dtm == -9999,
                                             dtm,
                                             copy=True)

            minz = np.floor(dtm_mask_na.min() - 1)
        dsm_file = f"{uls_dir}\\{plot_id}\\Grids\\{plot_id}_DSM_miniVUX.tif"
        if not os.path.exists(dsm_file):
            dsm_file = f"{plot}\\Grids\\{plot_id}_DSM_MLS.tif"

        with rasterio.open(dsm_file) as dsm_f:
            dsm = dsm_f.read(1)

            maxz = np.floor(dsm.max() + 5)

        # write config file
        cfg = config_template
        cfg['root_folder'] = r"C:\Users\Kueken\dev\occpy_dev\OccPy"
        cfg['laz_in'] = glob.glob(os.path.join(local_dir, plot_id), "*_heb_sort.laz")[0]  # we expect only one file here
        cfg['tif_in']['DTM'] = dtm_file
        cfg['tif_in']['DSM'] = dsm_file
        cfg['out_dir'] = os.path.join(local_dir, plot_id, f"OcclusionMap_{int(parameters['vox_size']*100)}cm")
        cfg['vox_dim'] = float(parameters['vox_size'])
        cfg['lower_threshold'] = int(1)
        cfg['points_per_iter'] = int(1000000)
        cfg['plot_dim'] = [float(minx), float(miny), float(minz), float(maxx), float(maxy), float(maxz)]
        cfg['ScanPos'] = glob.glob(os.path.join(local_dir, plot_id), f"*.txt")[0]

        # write cfg to json file
        cfg_out = os.path.join(local_dir, plot_id, f"OcclusionMap_{int(parameters['vox_size']*100)}cm", "config", f"{plot_id}_occpy_cfg.json")

        with open(cfg_out, 'w') as json_out:
            json.dump(cfg, json_out, indent=4)

        occpy = OccPy(laz_in=f"{cfg['laz_in']}",
                      out_dir=f"{cfg['out_dir']}",
                      vox_dim=cfg['vox_dim'],
                      lower_threshold=cfg['lower_threshold'],
                      points_per_iter=cfg['points_per_iter'],
                      plot_dim=cfg['plot_dim'],
                      output_voxels=False
                      )

        traj = pd.read_csv(cfg['ScanPos'], delimiter=",")

        occpy.define_sensor_pos(path2file=cfg['ScanPos'],
                               is_mobile=True,
                               single_return=False,
                               delimiter=",",
                               hdr_time=traj.columns[0],
                               hdr_x=traj.columns[1],
                               hdr_y=traj.columns[2],
                               hdr_z=traj.columns[3])

        tic = time.time()
        occpy.do_raytracing()
        toc = time.time()
        print(f"Elapsed time: {toc - tic} seconds")

        # test.save_raytracing_output()

        normalize_occlusion_output(input_folder=occpy.out_dir,
                                   PlotDim=cfg['plot_dim'],
                                   vox_dim=cfg['vox_dim'],
                                   dtm_file=f"{cfg['tif_in']['DTM']}",
                                   dsm_file=f"{cfg['tif_in']['DSM']}")

        # move output folder to remote
        shutil.copytree(src=os.path.join(local_dir, plot_id, f"OcclusionMap_{int(parameters['vox_size'] * 100)}cm"),
                        dst=os.path.join(in_dir, plot_id, f"OcclusionMap_{int(parameters['vox_size'] * 100)}cm"),
                        dirs_exist_ok=True)

        # delete local folder
        shutil.rmtree(os.path.join(local_dir, plot_id))
