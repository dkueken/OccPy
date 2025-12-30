import os

import numpy as np
import pandas as pd

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


# TODO: TEMP TEST

# if __name__ == '__main__':
    
    # tree = "COL"
    # DATA_DIR = f"/Stor1/wout/data/occpy_barbara/{tree}"
    # AXIS = 1
    # odir = os.path.join(DATA_DIR, "output", "all")

    # interactive_figure(odir, axis=AXIS)