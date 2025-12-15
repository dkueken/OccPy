from occpy.OccPy import OccPy
from occpy.OccPy import normalize_occlusion_output
from occpy.OccPy import get_Occl_TransectFigure
from occpy.OccPy import get_Occl_TransectFigure_BinaryOcclusion
from occpy.OccPy import get_Occlusion_ProfileFigure

import numpy as np
from occpy import TerrainModel
import os
import glob
import shutil
from datetime import datetime
import argparse
from pathlib import Path
import json
import laspy

import rasterio
from rasterio.fill import fillnodata

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


if __name__=="__main__":
    years2analyse = [2024]
    plots2analyse = ['FP10']
    overwrite = False

    vox_dim = 0.1

    lower_threshold = 1  # buffer above DTM to exclude the lowest layer which could be influenced by coregistration issues. buffer in meters

    height_prof = np.arange(lower_threshold,45,vox_dim)


    for year in years2analyse:

        if year == 2023:
            main_dir_remote = rf"\\speedy16-36\data_15\_PLS\20230329_RamerenWald_TimeSeries"
            main_dir_local = rf"E:\Data\_PLS\20230329_RamerenWald_TimeSeries"
        elif year == 2024:
            main_dir_remote = rf"\\speedy16-36\data_15\_PLS\20240319_RamerenWald_TimeSeries"
            main_dir_local = rf"E:\Data\_PLS\20240319_RamerenWald_TimeSeries"
        elif year == 2025:
            main_dir_remote = rf"\\speedy16-36\data_15\_PLS\20250217_RamerenWald_TimeSeries"
            main_dir_local = rf"E:\Data\_PLS\20250217_RamerenWald_TimeSeries"

        # get dates of the year
        dates = glob.glob(f"{main_dir_local}\\{year}*")

        for date_dir in dates:
            date = os.path.basename(date_dir)

            d2 = datetime(year, int(date[4:6]), int(date[-2:]))
            d1 = datetime(year, 1, 1)
            doy = (d2 - d1).days + 1

            # get plots for date
            plots = plots2analyse

            for plot_dir in plots:
                plot_id = os.path.basename(plot_dir)

                print(f"##### Processing {plot_id} for acquisition date {date}")
                if os.path.exists(f"{main_dir_local}\\{date}\\{plot_id}\\OcclusionMapping\\PointCloud_Slice_YZ_2676515_2676525.png"):
                    print(f"{plot_id} for {date} has already been process... skipping...")
                    continue


                config_file = f"{main_dir_remote}\\{date}\\{plot_id}\\OcclusionMapping\\config\\settings_occpy_MLS_{date}_{plot_id}.JSON"
                def load_config(config_file):
                    with open(config_file, 'r') as file:
                        return json.load(file)


                config = load_config(config_file)
                root_folder = config['root_folder']
                out_dir_local = f"{main_dir_local}\\{date}\\{plot_id}\\OcclusionMapping"
                # we have a problem with the 2024 dsm file, so we take the 2025 instead for this visualization purpose
                # dsm_file = glob.glob(f"{main_dir_local}\\{date}\\{plot_id}\\OcclusionMapping\\*_dsm_0.50_resc_0.1.tif")[0]
                dsm_file = glob.glob(
                    f"E:\\Data\\_PLS\\20250217_RamerenWald_TimeSeries\\20250217\\{plot_id}\\OcclusionMapping\\*_dsm_0.50_resc_0.1.tif")[
                    0]


                Nhit_norm, Nmiss_norm, Nocc_norm, Classification_norm, chm = normalize_occlusion_output(
                    input_folder=f"{main_dir_remote}\\{date}\\{plot_id}\\OcclusionMapping",
                    PlotDim=config['plot_dim'],
                    vox_dim=config['vox_dim'],
                    dtm_file=config['tif_in']['DTM'],
                    dsm_file=dsm_file,
                    lower_threshold=1,
                    output_folder=out_dir_local,
                    output_voxels=False)
                config['out_dir'] = out_dir_local


                # Get Occlusion Fraction
                OcclFrac = Nocc_norm.astype(float) / (
                            Nhit_norm.astype(float) + Nmiss_norm.astype(float) + Nocc_norm.astype(float))

                # Test occlusion visualization

                fig_prop = dict(fig_size=(3.5, 3.2),
                                label_size=8,
                                label_size_ticks=6,
                                label_size_tiny=5,
                                out_format='png', )

                get_Occl_TransectFigure_BinaryOcclusion(Nhit_norm, Classification_norm, plot_dim=config['plot_dim'],
                                                        vox_dim=config['vox_dim'], out_dir=config['out_dir'], axis=0,
                                                        start_ind=250, end_ind=350, chm=chm, vertBuffer=0,
                                                        fig_prop=fig_prop, show_plots=False)

                fig_prop = dict(fig_size=(1.75, 3.2),
                                label_size=8,
                                label_size_ticks=6,
                                label_size_tiny=5,
                                out_format='png', )
                get_Occlusion_ProfileFigure(Classification_norm, plot_dim=config['plot_dim'], vox_dim=config['vox_dim'],
                                            out_dir=config['out_dir'], low_thresh=0, vertBuffer=0, max_percentage=35,
                                            fig_prop=fig_prop, show_plots=False)

                # Create a point cloud slice of the shown transect
                laz_file = glob.glob(f"{main_dir_local}\\{date}\\{plot_id}\\LAZ\\*_ground_heb.laz")[0]


                xmin = 2676515 # equivalent to figure in perspective paper
                xmax = 2676525
                ymin = config['plot_dim'][1] + 0.81
                ymax = config['plot_dim'][4]

                laz = laspy.read(laz_file)
                # filter points based on corner coordinates
                mask = (
                    (laz.x >= xmin) & (laz.x <= xmax) &
                    (laz.y >= ymin) & (laz.y <= ymax)
                )
                subset = laz.points[mask]

                points = np.vstack((subset.x - xmin, subset.y - ymin, subset["height above ground"])).T
                intensity = subset.intensity.astype(np.float32)

                q1, q25, q75, q99 = np.percentile(intensity, [1, 25, 75, 99])
                iqr = q75 - q25
                p_low = q25 - 1.5 * iqr
                if p_low < q1:
                    p_low = q1
                p_high = q75 + 1.5 * iqr
                if p_high > q99:
                    p_high = q99
                intensity_clipped = np.clip(intensity, p_low, p_high)

                # Normalize intensity for coloring into range [0 1]
                gamma = 0.5
                intensity_norm = (intensity_clipped - np.min(intensity_clipped)) / (
                            np.max(intensity_clipped) - np.min(intensity_clipped))
                intensity_norm = intensity_norm ** gamma

                # reduce amoutn of points
                fraction = 1
                num_points = points.shape[0]
                keep_idx = np.random.choice(num_points, int(num_points * fraction), replace=False)

                points_sub = points[keep_idx]
                intensity_sub = intensity_norm[keep_idx]

                # we need to adapt plot_dimensions as we used a smaller slice to show transects

                fig_prop = dict(fig_size=(3.5, 3.2),
                                label_size=8,
                                label_size_ticks=6,
                                label_size_tiny=5,
                                out_format='png', )
                grey_cmap = plt.get_cmap(name='Grays_r')

                # define plot extent
                extent = [ymin - ymin, ymax - ymin, 0, Nhit_norm.shape[2] * config['vox_dim']]

                fig, ax = plt.subplots(figsize=fig_prop['fig_size'])
                sc = plt.scatter(points[:, 1], points[:, 2],
                                 c=intensity_norm, cmap=grey_cmap, marker=',', s=0.1, edgecolors='none')

                plt.xlabel("Y [m]", fontsize=fig_prop['label_size'])
                plt.ylabel("Height a.g. [m]", fontsize=fig_prop['label_size'])
                # define tick label size
                plt.yticks(fontsize=fig_prop['label_size_ticks'])
                plt.xticks(fontsize=fig_prop['label_size_ticks'])
                ax.axis(extent)

                # Add a label stating sensor and season similar to legend on occlusion figure
                slot_width = 0.28
                margin = 0.05
                slot_height = 0.05
                y_pos_axes = 0.98
                n_slots = 3
                gap = (1 - 2 * margin - n_slots * slot_width) / (n_slots - 1)

                # create dummy plot to get a legend like label
                dummy_handle = Line2D([], [], color='none', marker=None, linewidth=0)
                legend_ax = ax.inset_axes([margin, y_pos_axes, slot_width, slot_height])
                legend_ax.axis("off")
                legend = legend_ax.legend(
                    handles=[dummy_handle],
                    labels=[f"DOY: {doy}"],
                    loc='center',
                    frameon=True,
                    ncol=1,
                    fontsize=fig_prop['label_size_ticks'],
                    handlelength=0,  # no reserved space for the symbol
                    handletextpad=0  # adjust space between handle and text (0 = none)
                )
                legend.get_frame().set_alpha(1)

                # Colorbar for Intensity
                start_pos = 2 * (slot_width + gap)
                cax2 = ax.inset_axes([2 * (slot_width + gap) + margin, y_pos_axes, slot_width, slot_height])
                cb2 = plt.colorbar(sc, cax=cax2, orientation='horizontal')
                cb2.set_label("Intensity [-]", size=fig_prop['label_size_ticks'])
                cb2.ax.tick_params(labelsize=fig_prop['label_size_tiny'])

                # tight layout
                plt.tight_layout()

                # plt.colorbar(sc, label="Normalized Intensity")
                plt.savefig(f"{config['out_dir']}/PointCloud_Slice_YZ_{xmin}_{xmax}.{fig_prop['out_format']}",
                            dpi=600, format=fig_prop['out_format'])

                print('test')
                del laz, points, Nhit_norm, Nmiss_norm, Nocc_norm, Classification_norm, OcclFrac





