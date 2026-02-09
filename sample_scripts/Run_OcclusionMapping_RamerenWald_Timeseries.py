from occpy.OccPy import OccPy
from occpy.OccPy import normalize_occlusion_output
from occpy.OccPy import get_Occl_TransectFigure_BinaryOcclusion
from occpy.OccPy import get_Occlusion_ProfileFigure

import numpy as np
import pandas as pd
import json
import laspy
import rasterio
from rasterio.windows import from_bounds

import os
import shutil
import glob

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import time


def copy_file(src, dst):
    shutil.copyfile(src, dst)


def get_data(local_dir, main_dir, date, plot_id):
    # Create directory for local directory
    print(f"### Getting data for {plot_id} at {date}")
    os.makedirs(f"{local_dir}\\{date}\\{plot_id}", exist_ok=True)

    # get path to input laz file
    laz_in = glob.glob(f"{main_dir}\\{date}\\{plot_id}\\LAZ\\*_rot2LV95.laz")[
        0]  # we assume there is only one file that fits this condition

    # get path to input traj file
    traj_in = glob.glob(f"{main_dir}\\{date}\\{plot_id}\\Trajectories\\*_rot2LV95.txt")[
        0]  # we assume there is only one file that fits this condition

    # copy laz file to local
    copy_file(src=laz_in, dst=f"{local_dir}\\{date}\\{plot_id}\\{os.path.basename(laz_in)}")
    copy_file(src=traj_in, dst=f"{local_dir}\\{date}\\{plot_id}\\{os.path.basename(traj_in)}")


if __name__ == "__main__":

    years2process = [2023, 2024, 2025]

    local_wdir = r"C:\Users\kueken\dev\pai_estimation_amapvox\data\OcclusionMapping"

    plot_corner_coords_file = r'\\speedy16-36\data_15\_PLS\20250217_RamerenWald_TimeSeries\PlotCoordinates.csv'

    plot_coords = pd.read_csv(plot_corner_coords_file)

    parameters = dict(vox_dim=0.1,
                      lower_threshold = 0,
                      points_per_iter = 1000000,
                      plot_buffer=0)

    for year in years2process:
        print(f"##### Processing all acquisitions of year {year}")

        if year==2023:
            main_dir = fr"\\speedy16-36\data_15\_PLS\20230329_RamerenWald_TimeSeries"
            local_wdir = r"C:\Users\kueken\dev\pai_estimation_amapvox\data\OcclusionMapping\2023"
            # we always use the DTM from 2025, to be consistent. also, 2025 DTM has the least holes in it.
            uls_ref_dtm = fr"\\speedy16-36\data_15\_PLS\20250217_RamerenWald_TimeSeries\_references\UAVLS\Grids\20250305_1040_Ramerenwald_LFI_FP05_07_08_10_0_45_3.1_44.2lps_200pts_60m_90_clip_gS1_heb_dtm_0.50.tif"
            uls_ref_dsm = fr"\\speedy16-36\data_15\_PLS\20230329_RamerenWald_TimeSeries\REF_DEMs\20230216_Ramerenwald_LFI_FP05_07_08_10_0_45_doubleGrid_3.1_44.2lps_200pts_60m_90_clip_gS1_dsm_0.50.tif"
            hdr_time = "%time" # header of trajectory files are different between years
        elif year==2024:
            main_dir = fr"\\speedy16-36\data_15\_PLS\20240319_RamerenWald_TimeSeries"
            local_wdir = r"C:\Users\kueken\dev\pai_estimation_amapvox\data\OcclusionMapping\2024"
            # we always use the DTM from 2025, to be consistent. also, 2025 DTM has the least holes in it.
            uls_ref_dtm = fr"\\speedy16-36\data_15\_PLS\20250217_RamerenWald_TimeSeries\_references\UAVLS\Grids\20250305_1040_Ramerenwald_LFI_FP05_07_08_10_0_45_3.1_44.2lps_200pts_60m_90_clip_gS1_heb_dtm_0.50.tif"
            uls_ref_dsm = fr"\\speedy16-36\data_15\_PLS\20240319_RamerenWald_TimeSeries\_references\UAVLS\Grids\20240215_1040_Ramerenwald_LFI_FP05_07_08_10_0_45_3.1_44.2lps_200pts_60m_90_clip_gS1_heb_dsm_0.50.tif"
            hdr_time = "//world_time"  # header of trajectory files are different between years
        elif year==2025:
            main_dir = fr"\\speedy16-36\data_15\_PLS\20250217_RamerenWald_TimeSeries"
            local_wdir = r"C:\Users\kueken\dev\pai_estimation_amapvox\data\OcclusionMapping\2025"
            # we always use the DTM from 2025, to be consistent. also, 2025 DTM has the least holes in it.
            uls_ref_dtm = fr"\\speedy16-36\data_15\_PLS\20250217_RamerenWald_TimeSeries\_references\UAVLS\Grids\20250305_1040_Ramerenwald_LFI_FP05_07_08_10_0_45_3.1_44.2lps_200pts_60m_90_clip_gS1_heb_dtm_0.50.tif"
            uls_ref_dsm = fr"\\speedy16-36\data_15\_PLS\20250217_RamerenWald_TimeSeries\_references\UAVLS\Grids\20250305_1040_Ramerenwald_LFI_FP05_07_08_10_0_45_3.1_44.2lps_200pts_60m_90_clip_gS1_heb_dsm_0.50.tif"
            hdr_time = "//world_time"  # header of trajectory files are different between years

        # get acquisition dates for year
        dates = glob.glob(f"{main_dir}\\{year}*")

        json_template = r"C:\Users\kueken\dev\OccPy\config\settings_MLS_tutorial.JSON"

        for date_dir in dates:
            date = os.path.basename(date_dir)


            # get plots for date
            plots = glob.glob(f"{date_dir}\\FP*")

            for plot_dir in plots:
                plot_id = os.path.basename(plot_dir)

                if os.path.exists(f"{date_dir}\\{plot_id}\\OcclusionMapping\\Classification.npy"):
                    print(f"{plot_id} for acquisition date {date} has already been processed ... skipping")
                    continue

                print(f"##### Processing plot {plot_id} for acquisition date {date}")

                # get necessary data from speedy and store them locally
                get_data(local_dir=local_wdir, main_dir=main_dir, date=date, plot_id=plot_id)

                if not os.path.exists(f"{local_wdir}\\{date}\\{plot_id}\\OcclusionMapping"):
                    os.makedirs(f"{local_wdir}\\{date}\\{plot_id}\\OcclusionMapping")
                    os.makedirs(f"{local_wdir}\\{date}\\{plot_id}\\OcclusionMapping\\config")


                Plot_Center_X = plot_coords.loc[plot_coords['Plot'] == plot_id, 'cenX'].values[0]
                Plot_Center_Y = plot_coords.loc[plot_coords['Plot'] == plot_id, 'cenY'].values[0]


                minx = np.floor(Plot_Center_X - (
                        25 + parameters['plot_buffer']))
                maxx = np.floor(Plot_Center_X + (
                        25 + parameters['plot_buffer']))
                miny = np.floor(Plot_Center_Y - (
                        25 + parameters['plot_buffer']))
                maxy = np.floor(Plot_Center_Y + (
                        25 + parameters['plot_buffer']))

                # read in DTM and DSM to get min and max Z coordinates
                dtm_file = uls_ref_dtm

                with rasterio.open(dtm_file) as dtm_f:
                    # read and clip data to plot extent
                    bbox = (minx, miny, maxx, maxy)
                    window = from_bounds(*bbox, transform=dtm_f.transform)

                    dtm = dtm_f.read(1, window=window, masked=True)

                    minz = np.floor(dtm.min() - 1)

                dsm_file = uls_ref_dsm

                with rasterio.open(dsm_file) as dsm_f:
                    # read and clip data to plot extent
                    bbox = (minx, miny, maxx, maxy)
                    dsm = dsm_f.read(1, window=window, masked=True)

                    # get rid of

                    maxz = np.floor(dsm.max() + 5)

                # load config template and adapt accordingly
                with open(json_template, 'r') as file:
                    config = json.load(file)

                config['root_folder'] = os.path.dirname(os.getcwd())
                config['laz_in'] = glob.glob(f"{local_wdir}\\{date}\\{plot_id}\\*_rot2LV95.laz")[0] # There should only be one file
                config['ScanPos'] = glob.glob(f"{local_wdir}\\{date}\\{plot_id}\\*_rot2LV95.txt")[0] # There should only be one file
                config['tif_in']['DTM'] = uls_ref_dtm
                config['tif_in']['DSM'] = uls_ref_dsm
                config['out_dir'] = f"{local_wdir}\\{date}\\{plot_id}\\OcclusionMapping"
                config['vox_dim'] = parameters['vox_dim']
                config['lower_threshold'] = parameters['lower_threshold']
                config['points_per_iter'] = parameters['points_per_iter']
                config['plot_dim'] = [int(minx), int(miny), int(minz), int(maxx), int(maxy), int(maxz)]
                config['output_voxels'] = False

                # write json to file
                with open(rf"{local_wdir}\\{date}\\{plot_id}\\OcclusionMapping\\config\\settings_occpy_MLS_{date}_{plot_id}.JSON", 'w') as file:
                    json.dump(config, file, indent=4)

                # create OccPy Instance to initiate OccPy
                occpy = OccPy(laz_in = config['laz_in'],
                              out_dir = config['out_dir'],
                              vox_dim = config['vox_dim'],
                              lower_threshold = config['lower_threshold'],
                              points_per_iter = config['points_per_iter'],
                              plot_dim = config['plot_dim'])

                # Read in trajectory file
                occpy.define_sensor_pos(path2file=config['ScanPos'],
                                        is_mobile=True,
                                        single_return=True,
                                        delimiter=" ",
                                        hdr_time=hdr_time,
                                        hdr_x='x',
                                        hdr_y='y',
                                        hdr_z='z')

                # run occpy with specified parameters
                tic = time.time()
                occpy.do_raytracing()
                toc = time.time()
                print(f"Raytracing took {toc-tic:.2f} seconds.")


                # normalize outputs
                Nhit_norm, Nmiss_norm, Nocc_norm, Classification_norm, chm = normalize_occlusion_output(input_folder=config['out_dir'],
                                                                                                        PlotDim=config['plot_dim'],
                                                                                                        vox_dim=config['vox_dim'],
                                                                                                        dtm_file=config['tif_in']['DTM'],
                                                                                                        dsm_file=config['tif_in']['DSM'],
                                                                                                        lower_threshold=config['lower_threshold'],
                                                                                                        output_voxels=config['output_voxels'])

                # Create figures
                fig_prop = dict(fig_size=(3.5, 3.2),
                                label_size=8,
                                label_size_ticks=6,
                                label_size_tiny=5,
                                out_format='png')

                # Get a transect figure for Y axis (North-South) around plot center of 10 m depth
                get_Occl_TransectFigure_BinaryOcclusion(Nhit_norm, Classification_norm, plot_dim=config['plot_dim'], vox_dim=config['vox_dim'],
                                                        out_dir=config['out_dir'], axis=0, start_ind=200, end_ind=300, chm=chm, vertBuffer=10, fig_prop=fig_prop, show_plots=False)

                # Also get a profile figure of the transect
                fig_prop = dict(fig_size=(1.75, 3.2),
                                label_size=8,
                                label_size_ticks=6,
                                label_size_tiny=5,
                                out_format='png')

                get_Occlusion_ProfileFigure(Classification_norm, plot_dim=config['plot_dim'], vox_dim=config['vox_dim'], out_dir=config['out_dir'],
                                            low_thresh=0, vertBuffer=10, max_percentage=60, fig_prop=fig_prop, show_plots=False)


                # Copy/Move everything back to main_dir
                shutil.copytree(src=f"{local_wdir}\\{date}\\{plot_id}\\OcclusionMapping",
                                dst=f"{main_dir}\\{date}\\{plot_id}\\OcclusionMapping",
                                dirs_exist_ok=True)

                # delete local folder
                shutil.rmtree(f"{local_wdir}\\{date}\\{plot_id}")

















