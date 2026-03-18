import os

import numpy as np
import pandas as pd
import open3d as o3d
import pyvista as pv
import laspy
import json
from functools import partial

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.colors import LogNorm
from matplotlib.colors import to_rgb
import seaborn as sns

from mpl_toolkits.axes_grid1.inset_locator import inset_axes



# TODO: temp, give as parameter
VOX_DIM = 0.1

def get_Occl_TransectFigure(Nhit, Classification, OcclFrac, plot_dim, vox_dim, out_dir, start_ind=None, end_ind=None, axis=0, chm=None, vertBuffer=0, fig_prop=None, show_plots=False):
    """
    get_Occl_TransectFigure creates a matplotlib figure of a defined transect through the occlusion mapping output grid
    TODO: this function should be implemented in a more generic way!


    Parameters
    ----------
    Nhit: np.ndarray
        3D numpy array with number of hits in each grid cell (voxel)
    Classification: np.ndarray
        3D numpy array with voxel Classification (Observed with hit = 1, Observed & empty = 2, Occluded = 3, Unobserved = 4)
    OcclFrac: np.ndarray
        3D numpy array with Occlusion fraction
    plot_dim: np.ndarray
        plot dimension of the input grid, as in [minX, minY, minZ, maxX, maxY, maxZ]
    vox_dim: float
        voxel dimensions in meters (cubic voxel are assumed)
    out_dir: str
        path to output directory
    start_ind: int [default: None]
        voxel index of where the transect should start. If None [default] start_ind = 0
    end_ind: int [default: None]
        voxel index of where the transect should end. If None [defaulte] end_ind = Nhit.shape[axis]
    axis: int [0, 1, 2]
        axis index, either 0 (X-Axis), 1 (Y-Axis) or 2 (Z-Axis)
    chm: np.ndarray [default=None]
        2D canopy height model raster. if chm is None no CHM line will be plotted.
    vertBuffer: float [default=0]
        optional vertical buffer added to the figure, if axis=0 or axis=1. This adds a padding above the canopy, so legend
        entries are not overlapping the actual transect.
    fig_prop: dict [default=None]
        python dictionary with figure properties. If fig_prop = None [default], the following settings will be defined:
        fig_prop = dict(fig_size=(3.14, 2.25),  # figure size in inch
                        lable_size=8,           # font size for labels (e.g. x, y, z-axis labels=
                        label_size_ticks=6,     # font size for tick-labels
                        label_size_tiny=4,      # font size for other labels (e.g. legend labels)
                        out_format='png')       # output format of figure file

    show_plots: bool [default=False]
        Whether output figures should be shown [will pause the execution until figure is closed] or not.

    Returns
    -------

    """

    if fig_prop is None:
        fig_prop = dict(fig_size=(3.14, 2.25),
                        label_size=8,
                        label_size_ticks=6,
                        label_size_tiny=4,
                        out_format='png', )


    if start_ind is None:
        start_ind = 0
    if end_ind is None:
        end_ind = Nhit.shape[axis]

    grid_dim = (int((plot_dim[3] - plot_dim[0]) / vox_dim), int((plot_dim[4] - plot_dim[1]) / vox_dim), int((plot_dim[5] - plot_dim[2]) / vox_dim))

    chm_slice_ref = None
    if axis == 0:  # get YZ, project axis X
        Nhit_Slice = np.sum(Nhit[start_ind:end_ind, :, :], axis=axis)
        OcclFrac_Slice = np.sum(Classification[start_ind:end_ind, :, :] == 3, axis=axis) / (
                    end_ind - start_ind)
        if chm is not None:
            # chm is [ny, nx] so to get YZ we project axis 1
            chm_slice_ref = np.max(chm[:, start_ind:end_ind], axis=1)
    elif axis == 1:  # get XZ, project axis Y
        Nhit_Slice = np.sum(Nhit[:, start_ind:end_ind, :], axis=axis)
        # OcclFrac_Slice = np.sum(Classification[:, start_ind:end_ind, :] == 3, axis=axis) / (
        #         end_ind - start_ind)
        OcclFrac = OcclFrac[:, start_ind:end_ind, :]
        mask = (OcclFrac >= 0.8)

        # sum only where mask is True
        sum_vals = np.sum(np.where(mask, OcclFrac, 0), axis=axis)

        # count matching values along the axis
        count_vals = np.sum(mask, axis=axis)

        # Safe division: avoid divide by zero and assign default where count == 0
        with np.errstate(divide='ignore', invalid='ignore'):
            OcclFrac_Slice = np.divide(sum_vals, count_vals)
            OcclFrac_Slice[count_vals == 0] = 0

        if chm is not None:
            # chm is [ny, nx] so to get XZ we project axis 0
            chm_slice_ref = np.max(chm[start_ind:end_ind, :], axis=0)
    else:  # get a slice of Z-Axis
        Nhit_Slice = np.sum(Nhit[:, :, start_ind:end_ind], axis=axis)
        OcclFrac_Slice = np.sum(Classification[:, :, start_ind:end_ind] == 3, axis=axis) / (
                end_ind - start_ind)

    #NHits_Slice_log = np.log10(Nhit_Slice, where=(Nhit_Slice != 0))

    # we need to rotate the slice for visualization purposes
    OcclFrac_Slice = np.rot90(OcclFrac_Slice)
    NHit_Slice = np.rot90(Nhit_Slice)

    fig = plt.figure(figsize=fig_prop['fig_size'])
    ax = fig.add_subplot(1, 1, 1)
    x_axis_vect = None
    if axis == 0:
        ax.set_xlabel(f"Y [m]", fontsize=fig_prop['label_size'])
        ax.set_ylabel(f"Height a.g. [m]", fontsize=fig_prop['label_size'])
        extent = [plot_dim[1]-plot_dim[1], plot_dim[4]-plot_dim[1], 0, OcclFrac_Slice.shape[0] * vox_dim]
        if vertBuffer != 0:
            extent_buf = extent.copy()
            extent_buf[3] = extent_buf[3] + vertBuffer
            ax.axis(extent_buf)
        else:
            ax.axis(extent)
        x_axis_vect = np.linspace(start=plot_dim[1]-plot_dim[1], stop=plot_dim[4]-plot_dim[1], num=grid_dim[1])

    elif axis == 1:
        ax.set_xlabel(f"X [m]", fontsize=fig_prop['label_size'])
        ax.set_ylabel(f"Height a.g. [m]", fontsize=fig_prop['label_size'])
        extent = [plot_dim[0]-plot_dim[0], plot_dim[3]-plot_dim[0], 0, OcclFrac_Slice.shape[0] * vox_dim]
        if vertBuffer!=0:
            extent_buf = extent.copy()
            extent_buf[3] = extent_buf[3] + vertBuffer
            ax.axis(extent_buf)
        else:
            ax.axis(extent)
        x_axis_vect = np.linspace(start=plot_dim[0]-plot_dim[0], stop=plot_dim[3]-plot_dim[0], num=grid_dim[0])
    else:
        ax.set_xlabel(f"X [m]", fontsize=fig_prop['label_size'])
        ax.set_ylabel(f"Y [m]", fontsize=fig_prop['label_size'])
        extent = [plot_dim[0]-plot_dim[0], plot_dim[3]-plot_dim[0], plot_dim[1]-plot_dim[1], plot_dim[4]-plot_dim[1]]
        ax.axis(extent)

    # define tick label size
    plt.yticks(fontsize=fig_prop['label_size_ticks'])
    plt.xticks(fontsize=fig_prop['label_size_ticks'])

    reds_cmap = plt.get_cmap(name='inferno_r')
    reds_cmap.set_under('k', alpha=0)
    grey_cmap = plt.get_cmap(name='Grays_r')
    grey_cmap.set_under('k', alpha=0)
    # plot raster data

    p50, p99 = np.percentile(OcclFrac_Slice*100, [50, 99])

    im1 = ax.imshow(NHit_Slice, cmap=grey_cmap, norm=LogNorm(vmin=1, vmax=np.amax(NHit_Slice)), interpolation='none',
                    extent=extent, alpha=1, aspect='auto')
    im2 = ax.imshow(OcclFrac_Slice * 100, cmap=reds_cmap, vmin=p50, vmax=p99, clim=[p50, p99], interpolation='none',
                    alpha=0.75, aspect='auto',
                    extent=extent)

    # Define equally spaced horizontal slots for two colorbars and one legend
    n_slots = 3
    slot_width = 0.28
    margin = 0.05
    gap = (1 - 2*margin - n_slots * slot_width) / (n_slots - 1)
    slot_height = 0.05
    y_pos_axes = 0.98


    if x_axis_vect is not None:
        if chm is not None:
            chm_ref_plot = ax.plot(x_axis_vect, chm_slice_ref, label="ULS CHM")
            # chm_comp_plot = ax.plot(x_axis_vect, chm_slice_comp, label="Comp CHM", linestyle='--') #TODO: implement that!

            legend_ax = ax.inset_axes([slot_width + gap + margin, y_pos_axes, slot_width, slot_height])
            legend_ax.axis("off")
            legend = legend_ax.legend(handles=[chm_ref_plot[0]], loc='center', frameon=True, ncol=1, fontsize=fig_prop['label_size_ticks'])
            legend.get_frame().set_alpha(1)


    # define colorbars with position and dimension
    start_pos = 0
    cax1 = ax.inset_axes([margin, y_pos_axes, slot_width, slot_height])
    cb1 = plt.colorbar(im1, cax=cax1, orientation='horizontal')
    cb1.ax.tick_params(labelsize=fig_prop['label_size_tiny'])
    cb1.set_label("Nr. Hits", size=fig_prop['label_size_ticks'])

    # Change ticks to actual values
    ticks = [10, 100, 1000]
    cb1.set_ticks(ticks)
    cb1.set_ticklabels([str(t) for t in ticks])


    # Second colorbar for Occlusion
    start_pos = 2 * (slot_width + gap)
    cax2 = ax.inset_axes([2 * (slot_width + gap) + margin, y_pos_axes, slot_width, slot_height])
    cb2 = plt.colorbar(im2, cax=cax2, orientation='horizontal')
    cb2.set_label("Occlusion [%]", size=fig_prop['label_size_ticks'])
    cb2.ax.tick_params(labelsize=fig_prop['label_size_tiny'])


    # tight layout
    plt.tight_layout()

    # save figure
    if axis == 0:
        plt.savefig(
            os.path.join(out_dir, f"Occlusion_Slice_YZ_{start_ind}_{end_ind}_voxels.{fig_prop['out_format']}"),
            dpi=300, format=fig_prop['out_format'])
    elif axis == 1:
        plt.savefig(
            os.path.join(out_dir, f"Occlusion_Slice_XZ_{start_ind}_{end_ind}_voxels.{fig_prop['out_format']}"),
            dpi=300, format=fig_prop['out_format'])
    else:
        plt.savefig(
            os.path.join(out_dir, f"Occlusion_Slice_XY_{start_ind}_{end_ind}_voxels.{fig_prop['out_format']}"),
            dpi=300, format=fig_prop['out_format'])

    if show_plots:
        plt.show(block=True)
    else:
        plt.close()

def get_Occl_TransectFigure_BinaryOcclusion(Nhit, Classification, plot_dim, vox_dim, out_dir, start_ind=None, end_ind=None, axis=0, chm=None, vertBuffer=0, nhit_max=100000, nhit_min=1, fig_prop=None, show_plots=False):
    """
        get_Occl_TransectFigure_BinaryOcclusion creates a matplotlib figure of a defined transect through the occlusion mapping output grid
        TODO: this function should be implemented in a more generic way and potentially be integrated into get_Occl_TransectFigure()


        Parameters
        ----------
        Nhit: np.ndarray
            3D numpy array with number of hits in each grid cell (voxel)
        Classification: np.ndarray
            3D numpy array with voxel Classification (Observed with hit = 1, Observed & empty = 2, Occluded = 3, Unobserved = 4)
        plot_dim: np.ndarray
            plot dimension of the input grid, as in [minX, minY, minZ, maxX, maxY, maxZ]
        vox_dim: float
            voxel dimensions in meters (cubic voxel are assumed)
        out_dir: str
            path to output directory
        start_ind: int [default: None]
            voxel index of where the transect should start. If None [default] start_ind = 0
        end_ind: int [default: None]
            voxel index of where the transect should end. If None [defaulte] end_ind = Nhit.shape[axis]
        axis: int [0, 1, 2]
            axis index, either 0 (X-Axis), 1 (Y-Axis) or 2 (Z-Axis)
        chm: np.ndarray [default=None]
            2D canopy height model raster. TODO: check and implement behavior if chm is not provided.
        vertBuffer: float [default=0]
            optional vertical buffer added to the figure, if axis=0 or axis=1. This adds a padding above the canopy, so legend
            entries are not overlapping the actual transect.
        fig_prop: dict [default=None]
            python dictionary with figure properties. If fig_prop = None [default], the following settings will be defined:
            fig_prop = dict(fig_size=(3.14, 2.25),  # figure size in inch
                            lable_size=8,           # font size for labels (e.g. x, y, z-axis labels=
                            label_size_ticks=6,     # font size for tick-labels
                            label_size_tiny=4,      # font size for other labels (e.g. legend labels)
                            out_format='png')       # output format of figure file
        show_plots: bool [default=False]
            Whether output figures should be shown [will pause the execution until figure is closed] or not.

        Returns
        -------

        """
    if fig_prop is None:
        fig_prop = dict(fig_size=(3.14, 2.25),
                        label_size=8,
                        label_size_ticks=6,
                        label_size_tiny=4,
                        out_format='png', )

    if start_ind is None:
        start_ind = 0
    if end_ind is None:
        end_ind = Nhit.shape[axis]

    grid_dim = (int((plot_dim[3] - plot_dim[0]) / vox_dim), int((plot_dim[4] - plot_dim[1]) / vox_dim), int((plot_dim[5] - plot_dim[2]) / vox_dim))

    chm_slice_ref = None
    if axis == 0:  # get YZ, project axis X
        Nhit_Slice = np.sum(Nhit[start_ind:end_ind, :, :], axis=axis)
        OcclFrac_Slice = np.sum(Classification[start_ind:end_ind, :, :] == 3, axis=axis) / (
                    end_ind - start_ind)
        if chm is not None:
            # chm is [ny, nx] so to get YZ we project axis 1
            chm_slice_ref = np.max(chm[:, start_ind:end_ind], axis=1)
    elif axis == 1:  # get XZ, project axis Y
        Nhit_Slice = np.sum(Nhit[:, start_ind:end_ind, :], axis=axis)
        OcclFrac_Slice = np.sum(Classification[:, start_ind:end_ind, :] == 3, axis=axis) / (end_ind - start_ind)
        if chm is not None:
            # chm is [ny, nx] so to get XZ we project axis 0
            chm_slice_ref = np.max(chm[start_ind:end_ind, :], axis=0)
    else:  # get a slice of Z-Axis
        Nhit_Slice = np.sum(Nhit[:, :, start_ind:end_ind], axis=axis)
        OcclFrac_Slice = np.sum(Classification[:, :, start_ind:end_ind] == 3, axis=axis) / (
                end_ind - start_ind)

    #NHits_Slice_log = np.log10(Nhit_Slice, where=(Nhit_Slice != 0))

    # we need to rotate the slice for visualization purposes
    OcclFrac_Slice = np.rot90(OcclFrac_Slice)
    NHit_Slice = np.rot90(Nhit_Slice)

    fig = plt.figure(figsize=fig_prop['fig_size'])
    ax = fig.add_subplot(1, 1, 1)
    x_axis_vect = None
    if axis == 0:
        ax.set_xlabel(f"Y [m]", fontsize=fig_prop['label_size'])
        ax.set_ylabel(f"Height a.g. [m]", fontsize=fig_prop['label_size'])
        extent = [plot_dim[1]-plot_dim[1], plot_dim[4]-plot_dim[1], 0, OcclFrac_Slice.shape[0] * vox_dim]
        if vertBuffer != 0:
            extent_buf = extent.copy()
            extent_buf[3] = extent_buf[3] + vertBuffer
            ax.axis(extent_buf)
        else:
            ax.axis(extent)
        x_axis_vect = np.linspace(start=plot_dim[1]-plot_dim[1], stop=plot_dim[4]-plot_dim[1], num=grid_dim[1])

    elif axis == 1:
        ax.set_xlabel(f"X [m]", fontsize=fig_prop['label_size'])
        ax.set_ylabel(f"Height a.g. [m]", fontsize=fig_prop['label_size'])
        extent = [plot_dim[0]-plot_dim[0], plot_dim[3]-plot_dim[0], 0, OcclFrac_Slice.shape[0] * vox_dim]
        if vertBuffer!=0:
            extent_buf = extent.copy()
            extent_buf[3] = extent_buf[3] + vertBuffer
            ax.axis(extent_buf)
        else:
            ax.axis(extent)
        x_axis_vect = np.linspace(start=plot_dim[0]-plot_dim[0], stop=plot_dim[3]-plot_dim[0], num=grid_dim[0])
    else:
        ax.set_xlabel(f"X [m]", fontsize=fig_prop['label_size'])
        ax.set_ylabel(f"Y [m]", fontsize=fig_prop['label_size'])
        extent = [plot_dim[0]-plot_dim[0], plot_dim[3]-plot_dim[0], plot_dim[1]-plot_dim[1], plot_dim[4]-plot_dim[1]]
        ax.axis(extent)

    # define tick label size
    plt.yticks(fontsize=fig_prop['label_size_ticks'])
    plt.xticks(fontsize=fig_prop['label_size_ticks'])

    reds_cmap = plt.get_cmap(name='plasma_r')
    reds_cmap.set_under('k', alpha=0)
    grey_cmap = plt.get_cmap(name='Grays_r')
    grey_cmap.set_under('k', alpha=0)
    # plot raster data


    im1 = ax.imshow(NHit_Slice, cmap=grey_cmap, norm=LogNorm(vmin=nhit_min, vmax=nhit_max), interpolation='none',
                    extent=extent, alpha=1, aspect='auto')
    im2 = ax.imshow(OcclFrac_Slice * 100, cmap=reds_cmap, vmin=1, vmax=50, clim=[1, 50], interpolation='none',
                    alpha=0.75, aspect='auto',
                    extent=extent)

    # Define equally spaced horizontal slots for two colorbars and one legend
    n_slots = 3
    slot_width = 0.28
    margin = 0.05
    gap = (1 - 2*margin - n_slots * slot_width) / (n_slots - 1)
    slot_height = 0.05
    y_pos_axes = 0.98


    if x_axis_vect is not None:
        if chm is not None:
            chm_ref_plot = ax.plot(x_axis_vect, chm_slice_ref, label="ULS CHM")
            # chm_comp_plot = ax.plot(x_axis_vect, chm_slice_comp, label="Comp CHM", linestyle='--') #TODO: implement that!

            legend_ax = ax.inset_axes([slot_width + gap + margin, y_pos_axes, slot_width, slot_height])
            legend_ax.axis("off")
            legend = legend_ax.legend(handles=[chm_ref_plot[0]], loc='center', frameon=True, ncol=1, fontsize=fig_prop['label_size_ticks'])
            legend.get_frame().set_alpha(1)


    # define colorbars with position and dimension
    start_pos = 0
    cax1 = ax.inset_axes([margin, y_pos_axes, slot_width, slot_height])
    cb1 = plt.colorbar(im1, cax=cax1, orientation='horizontal')
    cb1.ax.tick_params(labelsize=fig_prop['label_size_tiny'])
    cb1.set_label("Nr. Hits", size=fig_prop['label_size_ticks'])

    # Change ticks to actual values
    ticks = [10, 100, 1000]
    cb1.set_ticks(ticks)
    cb1.set_ticklabels([str(t) for t in ticks])


    # Second colorbar for Occlusion
    start_pos = 2 * (slot_width + gap)
    cax2 = ax.inset_axes([2 * (slot_width + gap) + margin, y_pos_axes, slot_width, slot_height])
    cb2 = plt.colorbar(im2, cax=cax2, orientation='horizontal')
    cb2.set_label("Occlusion [%]", size=fig_prop['label_size_ticks'])
    cb2.ax.tick_params(labelsize=fig_prop['label_size_tiny'])

    # tight layout
    plt.tight_layout()

    # save figure
    if axis == 0:
        plt.savefig(
            os.path.join(out_dir, f"Occlusion_Slice_YZ_{start_ind}_{end_ind}_voxels_binary.{fig_prop['out_format']}"),
            dpi=300, format=fig_prop['out_format'])
    elif axis == 1:
        plt.savefig(
            os.path.join(out_dir, f"Occlusion_Slice_XZ_{start_ind}_{end_ind}_voxels_binary.{fig_prop['out_format']}"),
            dpi=300, format=fig_prop['out_format'])
    else:
        plt.savefig(
            os.path.join(out_dir, f"Occlusion_Slice_XY_{start_ind}_{end_ind}_voxels_binary.{fig_prop['out_format']}"),
            dpi=300, format=fig_prop['out_format'])

    if show_plots:
        plt.show(block=True)
    else:
        plt.close()

# Function to darken an RGB color
def darken_color(color, amount=0.6):
    """
    helper function to make provided color darker by amount

    Parameters
    ----------
    color:
    amount

    Returns
    -------
    list of RGB colors (r, g, b)

    """
    r, g, b = to_rgb(color)
    return (r * amount, g * amount, b * amount)

def get_Occlusion_ProfileFigure(Classification, plot_dim, vox_dim, out_dir, low_thresh=0, vertBuffer=0, max_percentage=100, fig_prop=None, show_plots=False):
    """get_Occlusion_ProfileFigure produces a profile figure of Occluded, filled, and empty voxels

    Parameters
    ----------
    Classification: np.ndarray
        3D voxel grid with voxel classification (1=Observed with hit, 2 = observed empty, 3 = occlusion, 4 = unobserved)
    plot_dim: np.ndarray
        plot dimension of the input grid, as in [minX, minY, minZ, maxX, maxY, maxZ]
    vox_dim: float
        voxel dimension in meters (assuming cubic voxels)
    out_dir: str
        directory for figure output
    low_thresh: float, default 0
        cut-off for lower part of the grid to exclude e.g. high occlusion towards the ground
    vertBuffer: float, default 0
        vertical buffer to add ontop of highest canopy point. e.g. needed to align Y-Axis with transect figure made with
        get_Occl_TransectFigure
    max_percentage: float, default 100
        maximum volume percentage to be shown on x-Axis
    fig_prop: dict, default None
            python dictionary with figure properties. If fig_prop = None [default], the following settings will be defined:
            fig_prop = dict(fig_size=(3.14, 2.25),  # figure size in inch
                            lable_size=8,           # font size for labels (e.g. x, y, z-axis labels=
                            label_size_ticks=6,     # font size for tick-labels
                            label_size_tiny=4,      # font size for other labels (e.g. legend labels)
                            out_format='png')       # output format of figure file
    show_plots: bool, default False
            Whether output figures should be shown [will pause the execution until figure is closed] or not.

    Returns
    -------

    """

    grid_dim = (int((plot_dim[3] - plot_dim[0]) / vox_dim), int((plot_dim[4] - plot_dim[1]) / vox_dim),
                int((plot_dim[5] - plot_dim[2]) / vox_dim))

    vert_vect = np.arange(start=low_thresh, stop=Classification.shape[2] * vox_dim, step=vox_dim)
    Classification = Classification[:,:,int(low_thresh / vox_dim):]
    # a hack to make sure that vert_vect is of the same length as OcclVertProf TODO: this has to be checked if it is generic!


    OcclVertProf = np.sum(Classification == 3, axis=0)
    OcclVertProf = np.sum(OcclVertProf, axis=0)
    OcclVertProf_Rel = OcclVertProf / ((grid_dim[0]) * (grid_dim[1]))

    FilledVertProf = np.sum(Classification == 1, axis=0)
    FilledVertProf = np.sum(FilledVertProf, axis=0)
    FilledVertProf_Rel = FilledVertProf / (grid_dim[0] * grid_dim[1])

    EmptyVertProf = np.sum(np.logical_or(Classification == 2, Classification==0), axis=0)
    EmptyVertProf = np.sum(EmptyVertProf, axis=0)
    EmptyVertProf_Rel = EmptyVertProf / (grid_dim[0] * grid_dim[1])

    heights = vert_vect[0:len(OcclVertProf)]

    percentages = np.column_stack([FilledVertProf_Rel*100, OcclVertProf_Rel*100, EmptyVertProf_Rel*100])
    categories = ['Filled', 'Occluded', 'Empty']
    colors = ['skyblue', 'salmon', 'lightgreen']

    # Compute cumulative percentages for stacking
    cumulative = np.cumsum(percentages, axis=1)

    palette = sns.color_palette('colorblind', n_colors=len(categories))
    palette[2] = (1.0, 1.0, 1.0) # white for empty

    fig, ax = plt.subplots(figsize=fig_prop['fig_size'])

    for i, cat in enumerate(categories):
        left = cumulative[:, i - 1] if i > 0 else np.zeros_like(heights)
        face_color = palette[i]
        edge_color = darken_color(face_color, 0.8)  # slightly darker for lines

        # Fill area
        ax.fill_betweenx(
            heights, left, cumulative[:, i],
            color=face_color, alpha=0.6
        )
        # Outline
        ax.plot(cumulative[:, i], heights, color=edge_color, linewidth=1.5, label="_nolegend_")

    ax.set_xlabel('Percentage of voxels [%]', fontsize=fig_prop['label_size'])
    ax.set_ylabel('Height above ground [m]', fontsize=fig_prop['label_size'])
    ax.set_xlim(0.1,max_percentage)
    ax.set_ylim(0,np.max(heights) + vertBuffer)
    plt.xticks(fontsize=fig_prop['label_size_ticks'])
    plt.yticks(fontsize=fig_prop['label_size_ticks'])
    ax.legend(categories[0:2], fontsize=fig_prop['label_size_ticks'])
    plt.tight_layout()

    plt.savefig(os.path.join(out_dir, f"OcclusionVertProf.{fig_prop['out_format']}"), dpi=300, format=fig_prop['out_format'])
    if show_plots:
        plt.show(block=True)
    else:
        plt.close()

def interactive_figure(output_dir, axis=0):
    """
    Create an interactive slice viewer for voxel grids with occlusion overlay.

    Loads Classification.npy and Nhit.npy, builds interactive
    figure with sliders to select the slice center and projection depth, and
    visualizes:
      - log10(Nhit) as a grayscale heatmap, and
      - fraction of voxels classified as occluded as a colored overlay.

    Axis selects the slicing plane:
      - 0: YZ slice 
      - 1: XZ slice
      - 2: XY slice

    Parameters
    ----------
    output_dir : str
        Directory containing Classification.npy and Nhit.npy arrays.
    axis : int, default 0
        Axis orthogonal to the slicing plane (0, 1, or 2).
    """

    classification_arr = np.load(os.path.join(output_dir, "Classification.npy"))
    nhit_arr = np.load(os.path.join(output_dir, "Nhit.npy"))

    print("Shapes: (X,Y,Z)")
    print(f"Classification: {classification_arr.shape}")
    print(f"NHIT: {nhit_arr.shape}")
    
    # -- function to generate plot for given parameters
    
    def generate_image(center, depth, axis):
        start_ind = int(max(0, center-depth/2))
        max_ind = classification_arr.shape[axis] - 1 
        end_ind = int(min(max_ind, start_ind+depth))

        if axis==0: # get a slice of X-Axis, YZ image
            Nhit_Slice = np.sum(nhit_arr[start_ind:end_ind, :, :], axis=axis)
            OcclFrac_Slice = np.sum(classification_arr[start_ind:end_ind, :, :]==3, axis=axis) / (end_ind - start_ind)
        elif axis==1: # XZ image
            Nhit_Slice = np.sum(nhit_arr[:,start_ind:end_ind,:], axis=axis)
            OcclFrac_Slice = np.sum(classification_arr[:, start_ind:end_ind, :] == 3, axis=axis) / (end_ind - start_ind)
        else: # XY image
            Nhit_Slice = np.sum(nhit_arr[:, :, start_ind:end_ind], axis=axis)
            OcclFrac_Slice = np.sum(classification_arr[:, :, start_ind:end_ind] == 3, axis=axis) / (end_ind - start_ind)
        
        NHits_Slice_log = np.log10(Nhit_Slice, where=(Nhit_Slice != 0))
        OcclFrac_Slice = np.rot90(OcclFrac_Slice)
        NHits_Slice_log = np.rot90(NHits_Slice_log)
        vlim_occl = np.ceil(np.amax(OcclFrac_Slice*100) / 10.0) * 10
        MAX_value_nhits = np.amax(NHits_Slice_log)

        return NHits_Slice_log, OcclFrac_Slice, vlim_occl, MAX_value_nhits

    # -- init

    init_center = int(round(classification_arr.shape[axis]/2))
    init_depth = 50

    # ---- define figure thingies

    fig = plt.figure(figsize=(6, 4))
    # fig = plt.figure(figsize=(4, 8))
    # figure properties
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal', adjustable='box')
    # adjust the main plot to make room for the sliders
    # fig.subplots_adjust(left=0.1)
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.xaxis.label.set_size(12)
    ax.yaxis.label.set_size(12)
    ax.spines[['right', 'top']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(.8)
    ax.spines[['left', 'bottom']].set_color('k')
    reds_cmap = plt.get_cmap(name='viridis_r')
    reds_cmap.set_under('k', alpha=0)
    greens_cmap = plt.get_cmap(name='Greys_r')
    greens_cmap.set_under('k', alpha=0)

    # --- define extent
    if axis==0:
        ax.set_xlabel(f"Y (m)")
        ax.set_ylabel(f"Z (m)")
        extent = [0, classification_arr.shape[1]*VOX_DIM, 0, classification_arr.shape[2]*VOX_DIM]
        ax.axis(extent)
    elif axis==1:
        ax.set_xlabel(f"X (m)")
        ax.set_ylabel(f"Z (m)")
        extent = [0, classification_arr.shape[0]*VOX_DIM, 0, classification_arr.shape[2]*VOX_DIM]
        ax.axis(extent)
    else:
        ax.set_xlabel(f"X (m)")
        ax.set_ylabel(f"Y (m)")
        extent = [0, classification_arr.shape[0]*VOX_DIM, 0, classification_arr.shape[1]*VOX_DIM]
        ax.axis(extent)
    
    # --- generate initial image and show
    
    nhit_img, occl_img, lim_ocll, lim_nhit = generate_image(init_center, init_depth, axis=axis)
    im1 = ax.imshow(nhit_img, cmap=greens_cmap, clim=[0.1, lim_nhit], interpolation='none',
                    alpha=1, extent=extent)
    im2 = ax.imshow(occl_img * 100, cmap=reds_cmap, clim=[1, lim_ocll], interpolation='none',
                    alpha=1, extent=extent)
    
    # --- add colorbars
    # TODO: update these based on the slice?
    axins1 = inset_axes(
        ax,
        width="35%",
        height="5%",
        loc="upper right",
    )
    axins1.xaxis.set_ticks_position("bottom")
    fig.colorbar(im1, cax=axins1, orientation='horizontal', label="Log Nr.Hits")
    axins2 = inset_axes(
        ax,
        width="35%",
        height="5%",
        loc="upper left",
    )
    axins2.xaxis.set_ticks_position("bottom")
    fig.colorbar(im2, cax=axins2, orientation='horizontal', label="Occluded voxels (%)")

    # --- add sliders for center and depth

    ax_center = fig.add_axes([0.15, 0.1, 0.0225, 0.63])
    center_slider = Slider(
        ax=ax_center,
        label='Center \nvoxel',
        valmin=0,
        valmax=classification_arr.shape[axis]-1,
        valinit=init_center,
        orientation="vertical"
    )
    center_slider.label.set_size(8)
    ax_depth = fig.add_axes([0.05, 0.1, 0.0225, 0.63])
    depth_slider = Slider(
        ax=ax_depth,
        label="Depth of \nprojection \nin #voxels",
        valmin=1,
        valmax=100,
        valinit=init_depth,
        orientation="vertical"
    )
    depth_slider.label.set_size(8)

    # -- define update function

    def update(val):
        nhit_img, occl_img, lim_ocll, lim_nhit = generate_image(center_slider.val, depth_slider.val, axis)
        im1 = ax.imshow(nhit_img, cmap=greens_cmap, clim=[0.1, lim_nhit], interpolation='none',
                    alpha=1, extent=extent)
        im2 = ax.imshow(occl_img * 100, cmap=reds_cmap, clim=[1, lim_ocll], interpolation='none',
                        alpha=1, extent=extent)
        
    # --- show

    # plt.show()

    out_file = "test_out/TEMP_slice_fig.png"
    plt.savefig(out_file, dpi=300, format="png", bbox_inches="tight")

def plot_riegl_grid(data : pd.DataFrame, max_scanline, max_scanline_idx, image2=None, out_path=None):
    """
    Plot a scanline-by-index occupancy grid from RIEGL data.

    Builds a boolean image with shape (max_scanline_idx+1, max_scanline+1) marking
    where (scanline, scanline_idx) pairs exist in `data`. Optionally overlays a
    second boolean image and saves the figure to `out_path`.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame containing columns 'scanline' and 'scanline_idx'.
    max_scanline : int
        Maximum scanline index on the horizontal axis.
    max_scanline_idx : int
        Maximum scanline_idx on the vertical axis.
    image2 : array-like of bool, optional
        Secondary image to overlay (same shape as the grid).
    out_path : str, optional
        Path to save the resulting figure.

    """
    scanline_np = data[["scanline"]].to_numpy()
    scanline_idx_np = data[["scanline_idx"]].to_numpy()
    # scanline_idx_np = np.where(scanline_idx_np > max_scanline_idx, max_scanline_idx, scanline_idx_np)
    extent = [0, max_scanline+1, 0, max_scanline_idx+1]
    img = np.zeros(shape=(max_scanline_idx+1, max_scanline+1), dtype=bool)
    img[scanline_idx_np, scanline_np] = True
    figsize=(12,5)
    cmap = matplotlib.colors.ListedColormap(['white', 'red'])
    fig, ax = plt.subplots(ncols=1, nrows=1, squeeze=False, 
                           sharex=False, sharey=False, figsize=figsize)

    with plt.style.context('seaborn-v0_8-notebook'):
        ax[0,0].imshow(img, interpolation='nearest', extent=extent, 
                    clim=[0,1], cmap=plt.get_cmap(cmap, 2), vmin=0, vmax=1, alpha=1)
        ax[0,0].set(adjustable='box', aspect='equal')
        ax[0,0].set(xlabel="Scanline", ylabel="Scanline index")
        ax[0,0].set_facecolor('white')
        cmap_blue = matplotlib.colors.ListedColormap(['white', 'blue'])
        if image2 is not None:
            ax[0,0].imshow(image2, interpolation='nearest', extent=extent, clim=[0,1], cmap=plt.get_cmap(cmap_blue,2), vmin=0,vmax=1, alpha=0.5)
    fig.tight_layout()
    plt.show()
    if out_path is not None:
        fig.savefig(out_path)

def vis_pv_static_bounds(occmap_file, 
                       min_bound_voxel, 
                       max_bound_voxel, 
                       config_file,
                       pointcloud_file=None,
                       opacity_occluded=0.2,
                       opacity_hit=0.1,
                       opacity_unobserved=0.2,
                       point_size=2,
                       point_color=(0,0,0),
                       return_plotter=True):
    
    """
    Visualize an occupancy map with PyVista.

    Loads a voxel occupancy grid (.npy) and optionally point cloud (.las), crops a region based on max_bound and min_bound,
    and displays occluded (red), hit (green), and unobserved (blue) voxels as
    semi-transparent meshes, with optionally the point cloud overlaid.

    Ensure min_bound and max_bound are within the occupancy grid dimensions. Large regions may be quite slow to render.

    Parameters
    ----------
    occmap_file : str
        Path to the .npy file containing the occupancy map (3D array).
    min_bound_voxel : array-like of int, shape (3,)
        Minimum XYZ voxel coordinates to visualize
    max_bound_voxel : array-like of int, shape (3,)
        Maximum XYZ voxel coordinates to visualize
    config_file : str
        Path to config file containing occpy run parameters.
    pointcloud_file : str, default None
        If provided, visualize the point cloud from this .las file.
    opacity_occluded : float, default 0.2
        Opacity for occluded voxels (red).
    opacity_hit : float, default 0.1
        Opacity for hit voxels (green).
    opacity_unobserved : float, default 0.2
        Opacity for unobserved voxels (blue).
    point_size : float, default 2
        Point size for rendering the point cloud.
    point_color : tuple of float, default (0,0,0)
        RGB color for the point cloud points.
    return_plotter : bool, default False
        If True, return the PyVista plotter object for further manipulation instead of showing the plot
    """

    # check paths
    if not os.path.exists(occmap_file):
        raise FileNotFoundError(f"Occupancy map file not found: {occmap_file}")
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Config file not found: {config_file}")
    if pointcloud_file is not None and not os.path.exists(pointcloud_file):
        raise FileNotFoundError(f"Point cloud file not found: {pointcloud_file}")
    
    # read json config file
    with open(config_file) as file:
        settings = json.load(file)

    # check required keys
    if "vox_dim" not in settings or "plot_dim" not in settings:
        raise ValueError(f"Config file must contain 'vox_dim' and 'plot_dim' keys. Found keys: {list(settings.keys())}")
    
    # check if min_bound_voxel and max_bound_voxel are valid
    if not all(isinstance(x, int) for x in min_bound_voxel) or not all(isinstance(x, int) for x in max_bound_voxel):
        raise ValueError(f"min_bound_voxel and max_bound_voxel must be array-like of int. Got min_bound_voxel={min_bound_voxel}, max_bound_voxel={max_bound_voxel}")
    min_bound_voxel = np.array(min_bound_voxel)
    max_bound_voxel = np.array(max_bound_voxel)

    vox_dim = settings["vox_dim"]
    plot_dim = settings["plot_dim"]
    # get min and max bounds
    min_bound = np.array(plot_dim[:3])
    max_bound = np.array(plot_dim[3:6])

    occmap = np.load(occmap_file)
    dims = occmap.shape

    print(f"Occlusion map dimensions: {dims}")

    # check if min_bound and max_bound are inside dims
    for i in range(3):
        if min_bound_voxel[i] < 0 or max_bound_voxel[i] > dims[i]:
            raise ValueError(f"min_bound and max_bound must be within occupancy map dimensions {dims}, got min_bound_voxel={min_bound_voxel}, max_bound_voxel={max_bound_voxel}")
        
    # collect bounding boxes for each voxel type
    bboxs_occl = []
    bbox_unobserved = []
    bbox_hit = []
    
    for x in range(min_bound_voxel[0], max_bound_voxel[0]):
        for y in range(min_bound_voxel[1], max_bound_voxel[1]):
            for z in range(min_bound_voxel[2], max_bound_voxel[2]):
                min_bound_cube = min_bound + np.array([x, y, z]) * vox_dim
                max_bound_cube = min_bound_cube + vox_dim
                
                if occmap[x, y, z] == 3:  # Occluded
                    bboxs_occl.append(o3d.geometry.AxisAlignedBoundingBox(min_bound_cube, max_bound_cube))
                elif occmap[x, y, z] == 4:  # Unobserved
                    bbox_unobserved.append(o3d.geometry.AxisAlignedBoundingBox(min_bound_cube, max_bound_cube))
                elif occmap[x, y, z] == 1:  # Hit
                    bbox_hit.append(o3d.geometry.AxisAlignedBoundingBox(min_bound_cube, max_bound_cube))
    
    # create meshes
    plotter = pv.Plotter()
    
    if len(bboxs_occl) > 0:
        mesh_occl = batch_aabbs_to_mesh(bboxs_occl)
        pv_mesh_occl = o3d_mesh_to_pyvista(mesh_occl)
        plotter.add_mesh(pv_mesh_occl, opacity=opacity_occluded, color='red')
    
    if len(bbox_hit) > 0:
        mesh_hit = batch_aabbs_to_mesh(bbox_hit)
        pv_mesh_hit = o3d_mesh_to_pyvista(mesh_hit)
        plotter.add_mesh(pv_mesh_hit, opacity=opacity_hit, color='green')
    
    if len(bbox_unobserved) > 0:
        mesh_unobserved = batch_aabbs_to_mesh(bbox_unobserved)
        pv_mesh_unobserved = o3d_mesh_to_pyvista(mesh_unobserved)
        plotter.add_mesh(pv_mesh_unobserved, opacity=opacity_unobserved, color='blue')
    
    # add point cloud (cropped to visible region)
    if pointcloud_file is not None:
        las = laspy.read(pointcloud_file)
        points = np.vstack((las.x, las.y, las.z)).transpose()
        min_pc_crop = min_bound + min_bound_voxel * vox_dim
        max_pc_crop = min_bound + max_bound_voxel * vox_dim
        print(f"Crop bounds: {min_pc_crop} to {max_pc_crop}")
        
        mask = np.all((points >= min_pc_crop) & (points <= max_pc_crop), axis=1)
        points_in_crop = points[mask]
        if len(points_in_crop) == 0:
            print("Warning: No points in the point cloud are within the crop bounds.")
        else:
            point_cloud = pv.PolyData(points_in_crop)
            point_cloud["point_color"] = np.repeat(point_color, len(points_in_crop), axis=0).reshape(-1, 3)
            plotter.add_points(point_cloud, scalars='point_color', style="points", point_size=point_size, render_points_as_spheres=True)
            plotter.remove_scalar_bar()

    if return_plotter:
        return plotter
    else:
        plotter.show()

def vis_pv_rotating(occmap_file, 
                min_bound_voxel, 
                max_bound_voxel, 
                config_file,
                opath="occpy_pv_rotating.mp4",
                pointcloud_file=None,
                opacity_occluded=0.2,
                opacity_hit=0.1,
                opacity_unobserved=0.2,
                point_size=2,
                point_color=(0,0,0),
                framerate=30,
                n_frames=180,
                distance_factor=2,
                elevation=30):
    """
    Create a rotating visualization of an occupancy map with PyVista.

    Loads a voxel occupancy grid (.npy) and optionally point cloud (.las), crops a region based on max_bound and min_bound,
    and displays occluded (red), hit (green), and unobserved (blue) voxels as
    semi-transparent meshes, with optionally the point cloud overlaid.

    Ensure min_bound and max_bound are within the occupancy grid dimensions. Large regions may be quite slow to render.

    Parameters
    ----------
    occmap_file : str
        Path to the .npy file containing the occupancy map (3D array).
    min_bound_voxel : array-like of float, shape (3,)
        Minimum XYZ voxel coordinates to visualize
    max_bound_voxel : array-like of float, shape (3,)
        Maximum XYZ voxel coordinates to visualize
    config_file : str
        Path to config file containing occpy run parameters.
    opath: str, default "occpy_pv_rotating.mp4"
        Output path for the rotating video.
    pointcloud_file : str, default None
        If provided, visualize the point cloud from this .las file.
    opacity_occluded : float, default 0.2
        Opacity for occluded voxels (red).
    opacity_hit : float, default 0.1
        Opacity for hit voxels (green).
    opacity_unobserved : float, default 0.2
        Opacity for unobserved voxels (blue).
    point_size : float, default 2
        Point size for rendering the point cloud.
    point_color : tuple of float, default (0,0,0)
        RGB color for the point cloud points.
    framerate : int, default 30
        Frames per second for the output video.
    n_frames : int, default 180
        Number of frames in the rotation (e.g. 180 for a full 360 degree rotation at 2 degrees per frame).
    distance_factor: float, default 1.5
        Controls distance to scene for camera orbit. Higher values will show more of the scene but may reduce visibility of details.
    camera_elevation: int, default 30
        Elevation angle of the camera in degrees. Higher values will show more of the top-down view, lower values will be more level with the scene.
    """

    if not os.path.exists(occmap_file):
        raise FileNotFoundError(f"Occupancy map file not found: {occmap_file}")
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Config file not found: {config_file}")
    if pointcloud_file is not None and not os.path.exists(pointcloud_file):
        raise FileNotFoundError(f"Point cloud file not found: {pointcloud_file}")

    # read json config file
    with open(config_file) as file:
        settings = json.load(file)

    if "vox_dim" not in settings or "plot_dim" not in settings:
        raise ValueError(f"Config file must contain 'vox_dim' and 'plot_dim' keys. Found keys: {list(settings.keys())}")

    # check if min_bound_voxel and max_bound_voxel are valid
    if not all(isinstance(x, int) for x in min_bound_voxel) or not all(isinstance(x, int) for x in max_bound_voxel):
        raise ValueError(f"min_bound_voxel and max_bound_voxel must be array-like of int. Got min_bound_voxel={min_bound_voxel}, max_bound_voxel={max_bound_voxel}")
    min_bound_voxel = np.array(min_bound_voxel)
    max_bound_voxel = np.array(max_bound_voxel)

    vox_dim = settings["vox_dim"]
    plot_dim = settings["plot_dim"]
    # get min and max bounds
    min_bound = np.array(plot_dim[:3])
    max_bound = np.array(plot_dim[3:6])

    occmap = np.load(occmap_file)
    dims = occmap.shape
    
    # check if min_bound and max_bound are inside dims
    for i in range(3):
        if min_bound_voxel[i] < 0 or max_bound_voxel[i] > dims[i]:
            raise ValueError(f"min_bound and max_bound must be within occupancy map dimensions {dims}, got min_bound_voxel={min_bound_voxel}, max_bound_voxel={max_bound_voxel}")

    # collect bounding boxes for each voxel type
    bboxs_occl = []
    bbox_unobserved = []
    bbox_hit = []
    
    for x in range(min_bound_voxel[0], max_bound_voxel[0]):
        for y in range(min_bound_voxel[1], max_bound_voxel[1]):
            for z in range(min_bound_voxel[2], max_bound_voxel[2]):
                min_bound_cube = min_bound + np.array([x, y, z]) * vox_dim
                max_bound_cube = min_bound_cube + vox_dim
                
                if occmap[x, y, z] == 3:  # Occluded
                    bboxs_occl.append(o3d.geometry.AxisAlignedBoundingBox(min_bound_cube, max_bound_cube))
                elif occmap[x, y, z] == 4:  # Unobserved
                    bbox_unobserved.append(o3d.geometry.AxisAlignedBoundingBox(min_bound_cube, max_bound_cube))
                elif occmap[x, y, z] == 1:  # Hit
                    bbox_hit.append(o3d.geometry.AxisAlignedBoundingBox(min_bound_cube, max_bound_cube))
    
    # create meshes
    plotter = pv.Plotter(notebook=False, off_screen=True, window_size=[2048, 2048])
    
    if len(bboxs_occl) > 0:
        mesh_occl = batch_aabbs_to_mesh(bboxs_occl)
        pv_mesh_occl = o3d_mesh_to_pyvista(mesh_occl)
        plotter.add_mesh(pv_mesh_occl, opacity=opacity_occluded, color='red')
    
    if len(bbox_hit) > 0:
        mesh_hit = batch_aabbs_to_mesh(bbox_hit)
        pv_mesh_hit = o3d_mesh_to_pyvista(mesh_hit)
        plotter.add_mesh(pv_mesh_hit, opacity=opacity_hit, color='green')
    
    if len(bbox_unobserved) > 0:
        mesh_unobserved = batch_aabbs_to_mesh(bbox_unobserved)
        pv_mesh_unobserved = o3d_mesh_to_pyvista(mesh_unobserved)
        plotter.add_mesh(pv_mesh_unobserved, opacity=opacity_unobserved, color='blue')
    
    # read and add point cloud if given
    if pointcloud_file is not None:
        las = laspy.read(pointcloud_file)
        points = np.vstack((las.x, las.y, las.z)).transpose()
        min_pc_crop = min_bound + min_bound_voxel * vox_dim
        max_pc_crop = min_bound + max_bound_voxel * vox_dim
        print(f"Crop bounds: {min_pc_crop} to {max_pc_crop}")
        
        mask = np.all((points >= min_pc_crop) & (points <= max_pc_crop), axis=1)
        points_in_crop = points[mask]
        if len(points_in_crop) == 0:
            print("Warning: No points in the point cloud are within the crop bounds.")
        else:
            point_cloud = pv.PolyData(points_in_crop)
            point_cloud["point_color"] = np.repeat(point_color, len(points_in_crop), axis=0).reshape(-1, 3)
            plotter.add_points(point_cloud, scalars='point_color', style="points", point_size=point_size, render_points_as_spheres=True)
            plotter.remove_scalar_bar()

    # calculate center and radius for camera orbit
    min_bound_voxel = np.array(min_bound_voxel)
    max_bound_voxel = np.array(max_bound_voxel)
    center_voxel = (min_bound_voxel + max_bound_voxel) / 2
    bounds_size_voxel = max_bound_voxel - min_bound_voxel
    center_world = min_bound + center_voxel * vox_dim
    radius = np.linalg.norm(bounds_size_voxel * vox_dim) * distance_factor
    plotter.open_movie(opath, framerate=framerate)

    for i in range(n_frames):
        azimuth = i * (360 / n_frames)
        
        plotter.camera.position = (
            center_world[0] + radius * np.cos(np.radians(azimuth)) * np.cos(np.radians(elevation)),
            center_world[1] + radius * np.sin(np.radians(azimuth)) * np.cos(np.radians(elevation)),
            center_world[2] + radius * np.sin(np.radians(elevation))
        )
        
        plotter.camera.focal_point = center_world
        plotter.camera.up = (0, 0, 1)
        
        plotter.write_frame()

    plotter.close()

def vis_pv_interactive(occmap_file, 
                       min_bound_voxel, 
                       max_bound_voxel, 
                       config_file,
                       pointcloud_file=None,
                       opacity_occluded=0.2,
                       opacity_hit=0.1,
                       opacity_unobserved=0.2,
                       point_size=2,
                       point_color=(0,0,0),
                       return_plotter=True):
    """
    Visualize an occupancy map with PyVista in interactive mode.

    Loads a voxel occupancy grid (.npy) and optionally point cloud (.las), crops a region based on max_bound and min_bound,
    and displays occluded (red), hit (green), and unobserved (blue) voxels as
    semi-transparent meshes, with optionally the point cloud overlaid.

    Ensure min_bound and max_bound are within the occupancy grid dimensions. Large regions may be quite slow to render.

    Parameters
    ----------
    occmap_file : str
        Path to the .npy file containing the occupancy map (3D array).
    min_bound_voxel : array-like of int, shape (3,)
        Minimum XYZ voxel coordinates to visualize
    max_bound_voxel : array-like of int, shape (3,)
        Maximum XYZ voxel coordinates to visualize
    config_file : str
        Path to config file containing occpy run parameters.
    pointcloud_file : str, default None
        If provided, visualize the point cloud from this .las file.
    opacity_occluded : float, default 0.2
        Opacity for occluded voxels (red).
    opacity_hit : float, default 0.1
        Opacity for hit voxels (green).
    opacity_unobserved : float, default 0.2
        Opacity for unobserved voxels (blue).
    point_size : float, default 2
        Point size for rendering the point cloud.
    point_color : tuple of float, default (0,0,0)
        RGB color for the point cloud points.
    return_plotter : bool, default False
        If True, return the PyVista plotter object for further manipulation instead of showing the plot
    """

    occmap = np.load(occmap_file)
    print("occmap shape:", occmap.shape)

    # read json config file
    with open(config_file) as file:
        settings = json.load(file)

    vox_dim = settings["vox_dim"]
    plot_dim = settings["plot_dim"]
    # get min and max bounds
    min_bound = np.array(plot_dim[:3])
    max_bound = np.array(plot_dim[3:6])

    class_colors = {
        1: np.array([0, 255, 0], dtype=np.uint8),   # hit / green
        3: np.array([255, 0, 0], dtype=np.uint8),   # occluded / red
        4: np.array([0, 0, 255], dtype=np.uint8),   # unobserved / blue
    }
    class_alpha = {
        1: int(np.clip(round(opacity_hit * 255), 0, 255)),
        3: int(np.clip(round(opacity_occluded * 255), 0, 255)),
        4: int(np.clip(round(opacity_unobserved * 255), 0, 255)),
    }

    nx, ny, nz = occmap.shape

    # build full-scale meshes once and keep voxel indices for updates
    def build_mesh_for_class(target_class):
        verts, faces, cell_colors = [], [], []
        voxel_coords = []
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    if occmap[x, y, z] != target_class:
                        continue
                    voxel_coords.append((x, y, z))
                    base = np.array([x, y, z], dtype=float) * vox_dim + min_bound
                    cube_verts = base + vox_dim * np.array([
                        [0,0,0],[1,0,0],[1,1,0],[0,1,0],
                        [0,0,1],[1,0,1],[1,1,1],[0,1,1]
                    ])
                    idx0 = len(verts)
                    verts.extend(cube_verts)
                    cube_faces = [
                        [4, idx0, idx0+1, idx0+2, idx0+3],
                        [4, idx0+4, idx0+5, idx0+6, idx0+7],
                        [4, idx0, idx0+1, idx0+5, idx0+4],
                        [4, idx0+1, idx0+2, idx0+6, idx0+5],
                        [4, idx0+2, idx0+3, idx0+7, idx0+6],
                        [4, idx0+3, idx0, idx0+4, idx0+7],
                    ]
                    faces.extend(cube_faces)
                    cell_colors.extend([class_colors[target_class]]*6)
        verts = np.array(verts)
        faces = np.hstack(faces)
        cell_colors = np.array(cell_colors, dtype=np.uint8)
        # init alpha to 0, window updates will toggle visibility.
        rgba = np.hstack([cell_colors, np.zeros((cell_colors.shape[0],1), dtype=np.uint8)])
        mesh = pv.PolyData(verts, faces)
        mesh.cell_data["rgba"] = rgba
        return mesh, np.array(voxel_coords, dtype=np.int32)
    
    # build full meshes (slow)
    unique_values = np.unique(occmap)
    mesh_entries = []

    def _add_mesh_entry(mesh, voxels, alpha_on):
        if mesh is None:
            return
        n_cubes = voxels.shape[0]
        axis_bins = [[], [], []]
        for axis, size in enumerate((nx, ny, nz)):
            axis_vals = voxels[:, axis]
            axis_bins[axis] = [np.where(axis_vals == idx)[0] for idx in range(size)]
        mesh_entries.append({
            "mesh": mesh,
            "voxels": voxels,
            "n_cubes": n_cubes,
            "axis_bins": axis_bins,
            "alpha_on": alpha_on,
        })

    print("Constructing meshes, can be slow for large grids")

    if 1 in unique_values:
        mesh_hit, vox_hit = build_mesh_for_class(1)
        _add_mesh_entry(mesh_hit, vox_hit, class_alpha[1])
    else:
        mesh_hit = None
    if 3 in unique_values:
        mesh_occl, vox_occl = build_mesh_for_class(3)
        _add_mesh_entry(mesh_occl, vox_occl, class_alpha[3])
    else:
        mesh_occl = None
    if 4 in unique_values:  
        mesh_unobs, vox_unobs = build_mesh_for_class(4)
        _add_mesh_entry(mesh_unobs, vox_unobs, class_alpha[4])
    else:
        mesh_unobs = None

    print("Mesh construction done")
    
    # setup plotter
    plotter = pv.Plotter()
    if mesh_hit is not None:
        actor_hit = plotter.add_mesh(mesh_hit, scalars="rgba", rgb=True, lighting=False)
    if mesh_occl is not None:
        actor_occl = plotter.add_mesh(mesh_occl, scalars="rgba", rgb=True, lighting=False)
    if mesh_unobs is not None:
        actor_unobs = plotter.add_mesh(mesh_unobs, scalars="rgba", rgb=True, lighting=False)

    # add point cloud
    point_cloud = None
    point_actor = None
    point_rgba = None
    point_alpha_on = 255
    all_points = None
    points_vox_float = None
    point_axis_bins = None
    point_visible = None

    if pointcloud_file is not None:
        las = laspy.read(pointcloud_file)
        all_points = np.vstack((las.x, las.y, las.z)).transpose()
        points_vox_float = (all_points - min_bound) / vox_dim

        # build bins for fast updates
        points_vox_int = np.floor(points_vox_float).astype(np.int32)
        point_axis_bins = [[], [], []]
        for axis, size in enumerate((nx, ny, nz)):
            axis_vals = points_vox_int[:, axis]
            point_axis_bins[axis] = [np.where(axis_vals == idx)[0] for idx in range(size)]

        point_visible = np.zeros(all_points.shape[0], dtype=bool)

        min_pc_crop = min_bound + np.array([int(min_bound_voxel[0]), int(min_bound_voxel[1]), int(min_bound_voxel[2])]) * vox_dim
        max_pc_crop = min_bound + np.array([int(max_bound_voxel[0]), int(max_bound_voxel[1]), int(max_bound_voxel[2])]) * vox_dim
        init_mask = np.all((all_points >= min_pc_crop) & (all_points <= max_pc_crop), axis=1)
        point_visible[:] = init_mask

        # same logic as mesh: pre allocate points and set visibility with alpha channel
        point_cloud = pv.PolyData(all_points)
        point_rgba = np.zeros((all_points.shape[0], 4), dtype=np.uint8)
        point_rgba[:, :3] = point_color
        point_rgba[:, 3] = np.where(point_visible, point_alpha_on, 0).astype(np.uint8)
        point_cloud["rgba"] = point_rgba
        point_actor = plotter.add_points(
            point_cloud,
            scalars="rgba",
            rgb=True,
            style="points",
            point_size=point_size,
            render_points_as_spheres=True,
        )

    # setup window and updating
    window_x = [int(min_bound_voxel[0]), int(max_bound_voxel[0])]
    window_y = [int(min_bound_voxel[1]), int(max_bound_voxel[1])]
    window_z = [int(min_bound_voxel[2]), int(max_bound_voxel[2])]

    pref_size_x = max(1, window_x[1] - window_x[0])
    pref_size_y = max(1, window_y[1] - window_y[0])
    pref_size_z = max(1, window_z[1] - window_z[0])

    slider_widgets = {}

    bounds_text_actor = None

    def _update_window_text():
        nonlocal bounds_text_actor
        msg = (
            f"Window bounds: X [{window_x[0]}, {window_x[1]}) | "
            f"Y [{window_y[0]}, {window_y[1]}) | "
            f"Z [{window_z[0]}, {window_z[1]})"
        )
        if bounds_text_actor is None:
            bounds_text_actor = plotter.add_text(msg, position="upper_right", font_size=10, name="window_bounds")
        else:
            # update existing actor text
            try:
                bounds_text_actor.SetInput(msg)
            except AttributeError:
                plotter.add_text(msg, position="upper_right", font_size=10, name="window_bounds")

    def _window_mask_for_voxels(voxels):
        return (
            (window_x[0] <= voxels[:, 0]) & (voxels[:, 0] < window_x[1])
            & (window_y[0] <= voxels[:, 1]) & (voxels[:, 1] < window_y[1])
            & (window_z[0] <= voxels[:, 2]) & (voxels[:, 2] < window_z[1])
        )

    def update_window_full():
        nonlocal point_rgba
        for entry in mesh_entries:
            mesh = entry["mesh"]
            voxels = entry["voxels"]
            alpha_on = entry["alpha_on"]
            rgba = mesh.cell_data["rgba"].copy()
            alpha = rgba[:, 3].reshape(-1, 6)
            alpha[:] = 0
            mask = _window_mask_for_voxels(voxels)
            if np.any(mask):
                alpha[mask, :] = alpha_on
            mesh.cell_data["rgba"] = rgba

        if all_points is not None:
            mask = (
                (window_x[0] <= points_vox_float[:, 0]) & (points_vox_float[:, 0] < window_x[1])
                & (window_y[0] <= points_vox_float[:, 1]) & (points_vox_float[:, 1] < window_y[1])
                & (window_z[0] <= points_vox_float[:, 2]) & (points_vox_float[:, 2] < window_z[1])
            )
            point_visible[:] = mask
            point_rgba[:, 3] = np.where(point_visible, point_alpha_on, 0).astype(np.uint8)
            point_cloud["rgba"] = point_rgba

        _update_window_text()

        plotter.render()

    def update_window_incremental(axis, leaving_idx=None, entering_idx=None):
        nonlocal point_rgba
        for entry in mesh_entries:
            mesh = entry["mesh"]
            voxels = entry["voxels"]
            bins = entry["axis_bins"][axis]
            alpha_on = entry["alpha_on"]
            rgba = mesh.cell_data["rgba"].copy()
            alpha = rgba[:, 3].reshape(-1, 6)

            leaving_cubes = bins[leaving_idx] if leaving_idx is not None and 0 <= leaving_idx < len(bins) else np.empty(0, dtype=int)
            if leaving_cubes.size:
                alpha[leaving_cubes, :] = 0

            entering_cubes = bins[entering_idx] if entering_idx is not None and 0 <= entering_idx < len(bins) else np.empty(0, dtype=int)
            if entering_cubes.size:
                entering_voxels = voxels[entering_cubes]
                visible_mask = _window_mask_for_voxels(entering_voxels)
                if np.any(visible_mask):
                    alpha[entering_cubes[visible_mask], :] = alpha_on

            mesh.cell_data["rgba"] = rgba

        if all_points is not None:
            if leaving_idx is not None and 0 <= leaving_idx < len(point_axis_bins[axis]):
                leaving_points = point_axis_bins[axis][leaving_idx]
                if leaving_points.size:
                    point_visible[leaving_points] = False

            if entering_idx is not None and 0 <= entering_idx < len(point_axis_bins[axis]):
                entering_points = point_axis_bins[axis][entering_idx]
                if entering_points.size:
                    v = points_vox_float[entering_points]
                    keep = (
                        (window_x[0] <= v[:, 0]) & (v[:, 0] < window_x[1])
                        & (window_y[0] <= v[:, 1]) & (v[:, 1] < window_y[1])
                        & (window_z[0] <= v[:, 2]) & (v[:, 2] < window_z[1])
                    )
                    point_visible[entering_points] = keep

            point_rgba[:, 3] = np.where(point_visible, point_alpha_on, 0).astype(np.uint8)
            point_cloud["rgba"] = point_rgba

        _update_window_text()

        plotter.render()

    def _axis_window_and_dim(axis):
        if axis == 0:
            return window_x, nx
        if axis == 1:
            return window_y, ny
        return window_z, nz

    def _set_slider_representation(axis, value):
        widget = slider_widgets.get(axis)
        if widget is None:
            return
        rep = widget.GetSliderRepresentation()
        rep.SetValue(int(value))
        rep.SetLabelFormat('%.0f')

    def _set_window_start(axis, start_value):
        win, dim = _axis_window_and_dim(axis)
        if axis == 0:
            pref_size = pref_size_x
        elif axis == 1:
            pref_size = pref_size_y
        else:
            pref_size = pref_size_z

        # Sliders represent absolute indices [0, dim-1].
        # shrink at boundaries
        start = int(round(start_value))
        start = max(0, min(start, dim - 1))

        new_start = start
        new_end = min(new_start + pref_size, dim)

        old_start = win[0]
        old_end = win[1]

        if new_start == old_start and new_end == old_end:
            _set_slider_representation(axis, new_start)
            return
        
        win[0] = new_start
        win[1] = new_end

        # incremental shift using bins if possible
        if new_start == old_start + 1 and new_end == old_end + 1:
            update_window_incremental(axis=axis, leaving_idx=old_start, entering_idx=old_end)
        elif new_start == old_start + 1 and new_end == old_end:
            update_window_incremental(axis=axis, leaving_idx=old_start, entering_idx=None)
        elif new_start == old_start - 1 and new_end == old_end - 1:
            update_window_incremental(axis=axis, leaving_idx=old_end - 1, entering_idx=old_start - 1)
        elif new_start == old_start - 1 and new_end == old_end:
            update_window_incremental(axis=axis, leaving_idx=None, entering_idx=old_start - 1)
        else:
            update_window_full()

        _set_slider_representation(axis, new_start)

    def move_window(axis, step):
        win, _ = _axis_window_and_dim(axis)
        _set_window_start(axis, win[0] + int(step))

    # on-screen controls for notebook where key events don't work
    win_size_x = window_x[1] - window_x[0]
    win_size_y = window_y[1] - window_y[0]
    win_size_z = window_z[1] - window_z[0]

    plotter.add_text(
        "Window Controls: sliders (x/y/z) | Keys: x/b, y/n, z/m",
        position="upper_left",
        font_size=10,
    )

    def _slider_callback(value, widget, axis):
        discrete_value = int(round(value))
        _set_window_start(axis, discrete_value)
        rep = widget.GetSliderRepresentation()
        rep.SetValue(int(round(discrete_value)))
        rep.SetLabelFormat('%.0f')

    slider_widgets[0] = plotter.add_slider_widget(
        callback=partial(_slider_callback, axis=0),
        rng=[0, max(0, nx - 1)],
        value=window_x[0],
        title="X start",
        pointa=(0.03, 0.04),
        pointb=(0.33, 0.04),
        pass_widget=True,
        fmt="%.0f",
    )
    slider_widgets[1] = plotter.add_slider_widget(
        callback=partial(_slider_callback, axis=1),
        rng=[0, max(0, ny - 1)],
        value=window_y[0],
        title="Y start",
        pointa=(0.35, 0.04),
        pointb=(0.65, 0.04),
        pass_widget=True,
        fmt="%.0f",
    )
    slider_widgets[2] = plotter.add_slider_widget(
        callback=partial(_slider_callback, axis=2),
        rng=[0, max(0, nz - 1)],
        value=window_z[0],
        title="Z start",
        pointa=(0.67, 0.04),
        pointb=(0.97, 0.04),
        pass_widget=True,
        fmt="%.0f",
    )

    # key callbacks
    plotter.add_key_event("x", lambda: move_window(axis=0, step=1))
    plotter.add_key_event("b", lambda: move_window(axis=0, step=-1))
    plotter.add_key_event("y", lambda: move_window(axis=1, step=1))
    plotter.add_key_event("n", lambda: move_window(axis=1, step=-1))
    plotter.add_key_event("z", lambda: move_window(axis=2, step=1))
    plotter.add_key_event("m", lambda: move_window(axis=2, step=-1))
    
    # init window
    update_window_full()

    if return_plotter:
        return plotter
    else:
        plotter.show()


# util
def batch_aabbs_to_mesh(aabbs):
    """
    Convert a list of Open3D axis-aligned bounding boxes to a single triangle mesh.

    Parameters
    ----------
    aabbs : list of open3d.geometry.AxisAlignedBoundingBox
        List of bounding boxes to convert.

    Returns
    -------
    open3d.geometry.TriangleMesh
        Combined mesh of all bounding boxes.
    """
    
    all_vertices = []
    all_triangles = []
    offset = 0

    F = np.array([
        [0, 1, 2], [0, 2, 3],       # bottom (z=min)
        [4, 5, 6], [4, 6, 7],       # top (z=max)
        [0, 1, 5], [0, 5, 4],       # y=min face
        [3, 2, 6], [3, 6, 7],       # y=max face
        [0, 3, 7], [0, 7, 4],       # x=min face
        [1, 2, 6], [1, 6, 5],       # x=max face
    ], dtype=np.int32)

    for aabb in aabbs:
        min_x, min_y, min_z = aabb.get_min_bound()
        max_x, max_y, max_z = aabb.get_max_bound()

        V = np.array([
            [min_x, min_y, min_z],
            [max_x, min_y, min_z],
            [max_x, max_y, min_z],
            [min_x, max_y, min_z],
            [min_x, min_y, max_z],
            [max_x, min_y, max_z],
            [max_x, max_y, max_z],
            [min_x, max_y, max_z],
        ], dtype=np.float64)

        all_vertices.append(V)
        all_triangles.append(F + offset)
        offset += V.shape[0]

    all_vertices = np.vstack(all_vertices)
    all_triangles = np.vstack(all_triangles)

    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(all_vertices)
    mesh.triangles = o3d.utility.Vector3iVector(all_triangles)

    return mesh

def o3d_mesh_to_pyvista(o3d_mesh):
    """
    Convert an Open3D triangle mesh to a PyVista PolyData mesh.

    Parameters
    ----------
    o3d_mesh : open3d.geometry.TriangleMesh
        Open3D mesh to convert.

    Returns
    -------
    pyvista.PolyData
        PyVista mesh with vertices, faces, and optionally vertex colors.
    """
    
    vertices = np.asarray(o3d_mesh.vertices)
    triangles = np.asarray(o3d_mesh.triangles)
    faces = np.hstack([np.full((triangles.shape[0], 1), 3), triangles]).astype(np.int64).ravel()
    pv_mesh = pv.PolyData(vertices, faces)
    
    if o3d_mesh.has_vertex_colors():
        colors = np.asarray(o3d_mesh.vertex_colors)
        pv_mesh.point_data["Colors"] = (colors * 255).astype(np.uint8)
    
    return pv_mesh
