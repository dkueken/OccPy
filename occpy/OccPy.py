import time
import os
import glob

import numpy as np
import pandas as pd
import rasterio
from rasterio.fill import fillnodata
import laspy
import OSToolBox as ost
from tqdm import tqdm

# plotting functions
import matplotlib
try:
    matplotlib.use('TkAgg')
except ImportError:
    print("couldn't change matplotlib backend")
    
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import to_rgb
import seaborn as sns

from raytr import PyRaytracer
from occpy.TerrainModel import TerrainModel
from occpy.util import prepare_ply, read_trajectory_file, read_sensorpos_file, interpolate_traj, last_nonzero

# TODO: change print statements to logging like in occpyRIEGL

is_sorted = lambda a: np.all(a[:-1] <= a[1:])

def normalize_occlusion_output(input_folder, PlotDim, vox_dim, dtm_file, dsm_file=None, lower_threshold=0, output_voxels=False):
    """
    normalize_occlusion_output normalizes all occlusion output grids (Nhit, Nmiss, Nocc, Classification) with the specified DTM
    This function also calculates occlusion statistics for the total canopy volume (defined by the volume between DTM
    and DSM). Currently only binary occlusion is analysed at the moment (TODO: implement also fractional occlusion),
    i.e. only voxels that are completely occluded (Nhit==0 and Nmiss==0 and Nocc >0)
    :param input_folder: directory to the output of the raytracing algorithm
    :param PlotDim: Plot Dimensions defined as a list: (minX, minY, minZ, maxX, maxY, maxZ)
    :param dtm_file: DTM file (.tif) of the area of interest. Currently, both dimensions and pixel size should match the output grids
    :param dsm_file: DSM file (.tif) of the area of interest. Currently, both dimensions and pixel size should match the output grids
    :param lower_threshold: minimum Z coordinate to cut off lower part of the canopy, i.e. Voxels lieing at or below DTM. default=0
    :param output_voxels: if the voxel grids should be outputted as a ply file. default=False -> not yet working properly, recommend leaving this to False!
    :return:
    """

    Nhit = np.load(f"{input_folder}/Nhit.npy")
    Nmiss = np.load(f"{input_folder}/Nmiss.npy")
    Nocc = np.load(f"{input_folder}/Nocc.npy")
    Classification = np.load(f"{input_folder}/Classification.npy")

    # Get extent of voxel grid
    extent_voxgrid = (PlotDim[0], PlotDim[4], PlotDim[3], PlotDim[1])

    # This is a bit of a quick and dirty solution to check on the compatibility of voxel size and pixel size of terrain models. TODO: improve that!
    dtm = TerrainModel(dtm_file)
    gt = dtm.dtm.res
    pix_size = gt[0]

    dtm_fname = os.path.basename(dtm_file)

    if pix_size != vox_dim:
        dtm.crop2extent(
            extent=extent_voxgrid,
            out_file=f"{input_folder}\\{dtm_fname[:-4]}_resc_{vox_dim}.tif",
            res=vox_dim)

    ext = dtm.get_extent()
    extent_dtm = (ext.left, ext.top, ext.right, ext.bottom)

    if extent_dtm != extent_voxgrid:
        dtm.crop2extent(extent=extent_voxgrid,
                        out_file=f"{dtm.get_terrainmodel_path()[:-4]}_clipped.tif",
                        res=vox_dim)

    dtm_file = dtm.get_terrainmodel_path()

    with rasterio.open(dtm_file, 'r') as dtm_src:
        dtm = dtm_src.read(1)
        # TODO: check if this is still needed!
        dtm = np.flipud(dtm)  # we need to flip the terrain models in order to make them compatible with the Occlusion output
        # fill in data gaps in dtm
        dtm = fillnodata(dtm, mask=dtm != dtm_src.get_nodatavals()[0])

    if dsm_file is not None:

        dsm = TerrainModel(dsm_file)
        dsm_fname = os.path.basename(dsm_file)
        # check on pixel size
        gt = dsm.dtm.res
        pix_size = gt[0]

        if pix_size != vox_dim:
            dsm.crop2extent(
                extent=extent_voxgrid,
                out_file=f"{input_folder}\\{dsm_fname[:-4]}_resc_{vox_dim}.tif",
                res=vox_dim)

        ext = dsm.get_extent()
        extent_dsm = (ext.left, ext.top, ext.right, ext.bottom)
        if extent_dsm != extent_voxgrid:
            dsm.crop2extent(extent=extent_voxgrid,
                            out_file=f"{dsm.get_terrainmodel_path()[:-4]}_clipped.tif",
                            res=vox_dim)

        dsm_file = dsm.get_terrainmodel_path()

        with rasterio.open(dsm_file, 'r') as dsm_src:
            dsm = dsm_src.read(1)
            # TODO: check if this flip is still needed!
            dsm = np.flipud(dsm)  # we need to flip the terrain models in order to make them compatible with the Occlusion output
            dsm = fillnodata(dsm, mask=dsm != dsm_src.get_nodatavals()[0])

        chm = dsm - dtm

        Nhit_norm = np.zeros(
            (dtm.shape[1], dtm.shape[0], int(np.ceil(np.amax(chm) / vox_dim))), dtype=int)
        Nmiss_norm = np.zeros_like(Nhit_norm)
        Nocc_norm = np.zeros_like(Nhit_norm)
        Classification_norm = np.zeros_like(Nhit_norm)

        OcclFrac2D = np.zeros((dtm.shape[1], dtm.shape[0]))
        for y in range(0, dsm.shape[0], 1):
            for x in range(0, dsm.shape[1], 1):
                # get zind where DTM is located in grid at x,y
                zind_dtm = int(np.floor((dtm[y, x] - PlotDim[2]) / vox_dim))
                zind_dsm = int(np.floor((dsm[y, x] - PlotDim[2]) / vox_dim))
                # extract profile from grids
                prof_class = Classification[x, y, zind_dtm:zind_dsm]
                prof_class_buf = Classification[x, y,
                                 zind_dtm + int(np.ceil(lower_threshold / vox_dim)):zind_dsm]

                Classification_norm[x, y, 0:len(prof_class)] = prof_class
                # Calculate occlusion fraction for z profile
                num_occl = sum(prof_class_buf == 3)

                if len(prof_class_buf) == 0:
                    OcclFrac2D[x, y] = 0
                else:
                    OcclFrac2D[x, y] = num_occl / len(prof_class_buf)

                Nhit_norm[x, y, 0:len(prof_class)] = Nhit[x, y, zind_dtm:zind_dsm]
                Nmiss_norm[x, y, 0:len(prof_class)] = Nmiss[x, y, zind_dtm:zind_dsm]
                Nocc_norm[x, y, 0:len(prof_class)] = Nocc[x, y, zind_dtm:zind_dsm]


    else:
        # as we do not know the height of the scene a priori, we will initialize a 3 D grid with the same dimensions
        # as the unnormalized grids, introducing quite some overhead...
        Nhit_norm = np.zeros(Nhit.shape, dtype=int)
        Nmiss_norm = np.zeros(Nmiss.shape, dtype=int)
        Nocc_norm = np.zeros(Nocc.shape, dtype=int)
        Classification_norm = np.zeros(Classification.shape, dtype=int)

        OcclFrac2D = np.zeros((dtm.shape[1], dtm.shape[0]))
        chm = np.zeros(dtm.shape)

        max_len_prof = 0
        for y in range(0, dtm.shape[0], 1):
            for x in range(0, dtm.shape[1], 1):
                # get zind where DTM is located in grid at x,y
                zind_dtm = int(np.floor((dtm[y, x] - PlotDim[2]) / vox_dim))
                zind_dsm = last_nonzero(Nhit[x, y, :], axis=0)

                # If no dsm is provided, we take the DSM from the same
                # acquisition. This will introduce an under estimation of occlusion for ground based acquisitions
                # as occlusion on top of canopy is not counted.
                chm[y, x] = (zind_dsm - zind_dtm) * vox_dim

                # extract profile from grids
                prof_class = Classification[x, y, zind_dtm:zind_dsm]
                prof_class_buf = Classification[x, y,
                                 zind_dtm + int(np.ceil(lower_threshold / vox_dim)):zind_dsm]

                Classification_norm[x, y, 0:len(prof_class)] = prof_class
                # Calculate occlusion fraction for z profile
                num_occl = sum(prof_class_buf == 3)

                if len(prof_class_buf) == 0:
                    OcclFrac2D[x, y] = 0
                else:
                    OcclFrac2D[x, y] = num_occl / len(prof_class_buf)

                if len(prof_class) > max_len_prof:
                    max_len_prof = len(prof_class)

                Classification_norm[y, x, 0:len(prof_class)] = Classification[y, x, zind_dtm:zind_dsm]
                Nhit_norm[x, y, 0:len(prof_class)] = Nhit[x, y, zind_dtm:zind_dsm]
                Nmiss_norm[x, y, 0:len(prof_class)] = Nmiss[x, y, zind_dtm:zind_dsm]
                Nocc_norm[x, y, 0:len(prof_class)] = Nocc[x, y, zind_dtm:zind_dsm]


        # get rid of the excessive height of the grid
        Classification_norm = Classification_norm[:, :, 0:max_len_prof]
        Nhit_norm = Nhit_norm[:, :, 0:max_len_prof]
        Nmiss_norm = Nmiss_norm[:, :, 0:max_len_prof]
        Nocc_norm = Nocc_norm[:, :, 0:max_len_prof]

    print(f"Saving normalized output files into directory as .npy...")
    np.save(f"{input_folder}/Nhit_norm.npy", Nhit_norm)
    np.save(f"{input_folder}/Nmiss_norm.npy", Nmiss_norm)
    np.save(f"{input_folder}/Nocc_norm.npy", Nocc_norm)
    np.save(f"{input_folder}/Classification_norm.npy", Classification_norm)

    # write ply file TODO: This seems to not be working for me!
    if output_voxels:
        print(f"Saving normalized output files into directory as .ply...")
        tic = time.time()
        verts, faces = prepare_ply(vox_dim, PlotDim, Nhit_norm)
        ost.write_ply(f"{input_folder}/Nhit_norm.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
        verts, faces = prepare_ply(vox_dim, PlotDim, Nmiss_norm)
        ost.write_ply(f"{input_folder}/Nmiss_norm.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
        verts, faces = prepare_ply(vox_dim, PlotDim, Nocc_norm)
        ost.write_ply(f"{input_folder}/Nocc_norm.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
        verts, faces = prepare_ply(vox_dim, PlotDim, Classification_norm)
        ost.write_ply(f"{input_folder}/Classification_norm.ply", verts, ['X', 'Y', 'Z', 'data'],
                      triangular_faces=faces)
        toc = time.time()
        print("Elapsed Time: " + str(toc - tic) + " seconds")

    return Nhit_norm, Nmiss_norm, Nocc_norm, Classification_norm, chm

def slot_bbox(i, slot_width, gap, y_pos_axes, slot_height):
    start_x = i * (slot_width + gap)
    return (start_x, y_pos_axes, slot_width, slot_height)
def get_Occl_TransectFigure(Nhit, Classification, OcclFrac, plot_dim, vox_dim, out_dir, start_ind=None, end_ind=None, axis=0, chm=None, vertBuffer=0, fig_prop=None, show_plots=False):

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
            f"{out_dir}/Occlusion_Slice_YZ_{start_ind}_{end_ind}_voxels.{fig_prop['out_format']}",
            dpi=300, format=fig_prop['out_format'])
    elif axis == 1:
        plt.savefig(
            f"{out_dir}/Occlusion_Slice_XZ_{start_ind}_{end_ind}_voxels.{fig_prop['out_format']}",
            dpi=300, format=fig_prop['out_format'])
    else:
        plt.savefig(
            f"{out_dir}/Occlusion_Slice_XY_{start_ind}_{end_ind}_voxels.{fig_prop['out_format']}",
            dpi=300, format=fig_prop['out_format'])

    if show_plots:
        plt.show(block=True)
    else:
        plt.close()

def get_Occl_TransectFigure_BinaryOcclusion(Nhit, Classification, plot_dim, vox_dim, out_dir, start_ind=None, end_ind=None, axis=0, chm=None, vertBuffer=0, nhit_max=100000, nhit_min=1, fig_prop=None, show_plots=False):

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
            f"{out_dir}/Occlusion_Slice_YZ_{start_ind}_{end_ind}_voxels_binary.{fig_prop['out_format']}",
            dpi=300, format=fig_prop['out_format'])
    elif axis == 1:
        plt.savefig(
            f"{out_dir}/Occlusion_Slice_XZ_{start_ind}_{end_ind}_voxels_binary.{fig_prop['out_format']}",
            dpi=300, format=fig_prop['out_format'])
    else:
        plt.savefig(
            f"{out_dir}/Occlusion_Slice_XY_{start_ind}_{end_ind}_voxels_binary.{fig_prop['out_format']}",
            dpi=300, format=fig_prop['out_format'])

    if show_plots:
        plt.show(block=True)
    else:
        plt.close()

# Function to darken an RGB color
def darken_color(color, amount=0.6):
    r, g, b = to_rgb(color)
    return (r * amount, g * amount, b * amount)

def get_Occlusion_ProfileFigure(Classification, plot_dim, vox_dim, out_dir, low_thresh=0, vertBuffer=0, max_percentage=100, fig_prop=None, show_plots=False):

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

    plt.savefig(f"{out_dir}/OcclusionVertProf.{fig_prop['out_format']}", dpi=300, format=fig_prop['out_format'])
    if show_plots:
        plt.show(block=True)
    else:
        plt.close()




class OccPy:
    def __init__(self, laz_in, out_dir, vox_dim=0.1, lower_threshold=1, points_per_iter=10000000, plot_dim=None, output_voxels=False):
        self.laz_in_f = laz_in
        self.out_dir = out_dir
        os.makedirs(out_dir, exist_ok=True)
        self.vox_dim = vox_dim  # voxel dimension (cubic) in meters TODO: maybe implement non-cubic voxels?
        self.lower_threshold = lower_threshold  # lower threshold above ground to exclude TODO: check if necessary
        self.points_per_iter = points_per_iter  # number of points read in from laz file in each iteration
        self.output_voxels = output_voxels
        self.is_mobile = False
        self.traj_f = None
        self.traj = None
        self.senspos_f = None
        self.senspos = None
        self.single_return = None
        # TODO: find a better way to link scan position id between laz file and scan pos file!
        self.scan_pos_id_stridx = 0
        self.scan_pos_id_endstridx = 0

        # some parameters that will be filled during function calls
        self.TotalVolume = 0
        self.Volume0_3 = 0
        self.Volume3_10 = 0
        self.Volume10_max = 0
        self.TotalOcclusion = 0
        self.Occlusion0_3 = 0
        self.Occlusion3_10 = 0
        self.Occlusion10_max = 0

        self.OcclFrac2D = None


        # 3DGrids
        self.Nhit = None
        self.Nmiss = None
        self.Nocc = None
        self.Classification = None

        # 2D Grids
        self.dtm = None
        self.dsm = None
        self.chm = None

        if plot_dim is None:
            # TODO: Test if this works!
            # TODO: this assumes laz_in is single file?
            with laspy.open(laz_in) as file:
                hdr = file.header
                self.PlotDim = dict(minX=hdr.x_min,
                                    maxX=hdr.x_max,
                                    minY=hdr.y_min,
                                    maxY=hdr.y_max,
                                    minZ=hdr.z_min,
                                    maxZ=hdr.z_max)

        else:
            # we expect the format of plot_dim to be [minX, minY, minZ, maxX, maxY, maxZ]
            self.PlotDim = dict(minX=plot_dim[0],
                                maxX=plot_dim[3],
                                minY=plot_dim[1],
                                maxY=plot_dim[4],
                                minZ=plot_dim[2],
                                maxZ=plot_dim[5])

        self.grid_dim = dict(nx=int((self.PlotDim['maxX'] - self.PlotDim['minX']) / self.vox_dim),
                             ny=int((self.PlotDim['maxY'] - self.PlotDim['minY']) / self.vox_dim),
                             nz=int((self.PlotDim['maxZ'] - self.PlotDim['minZ']) / self.vox_dim))

        # initialize RayTr Object
        self.RayTr = PyRaytracer()

        # Define Grid
        minBound = np.array([self.PlotDim['minX'], self.PlotDim['minY'], self.PlotDim['minZ']])
        maxBound = np.array([self.PlotDim['maxX'], self.PlotDim['maxY'], self.PlotDim['maxZ']])
        self.RayTr.defineGrid(minBound, maxBound, self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz'],
                              self.vox_dim)

    def define_sensor_pos(self, path2file, is_mobile, single_return=None, delimiter=" ", hdr_time='%time', hdr_scanpos_id='', hdr_x='x', hdr_y='y', hdr_z='z', sens_pos_id_offset=0, str_idx_ScanPosID=0, str_end_idx_ScanPosID=0):
        """

        :param path2file: [mandatory] path to csv file with sensor position information
        :param is_mobile: [mandatory] True or False whether platform is mobile (MLS, ULS) or static (TLS)
        :param single_return: [adviced] True or False whether data is single return or multi return data
        :param delimiter: csv delimiter [default: " "]
        :param hdr_time: column header for time (only needed for mobile acquisition)
        :param hdr_scanpos_id: column header for scan pos id (only needed for static acquisitions -> equivalent to hdr_time in mobile acquisitions
        :param hdr_x: column header for x coordinates [default 'x']
        :param hdr_y: column header for y coordinates [default 'y']
        :param hdr_z: column header for z coordinates [default 'z']
        :param sens_pos_id_offset: Very specific use case where Scan Pos ID in position file does not correspond with Scan Pos ID in LAZ files and we need to add an offset
        :param str_idx_ScanPosID: string index of where the scan position identifier is written in the laz file name TODO: find a better way to handle this!
        :param str_end_idx_ScanPosID: string index of where the scan position identifier ends in the laz file name TODO: find a better way to handle this!
        :return:
        """
        self.is_mobile = is_mobile
        self.single_return = single_return

        if is_mobile: # case of mobile acquisitions (MLS, ULS)
            self.traj = read_trajectory_file(path2traj=path2file, delimiter=delimiter, hdr_time=hdr_time, hdr_x=hdr_x, hdr_y=hdr_y, hdr_z = hdr_z)
        else: # case of static acquisition (TLS)
            self.senspos = read_sensorpos_file(path2senspos=path2file, delimiter=delimiter, hdr_scanpos_id=hdr_scanpos_id, hdr_x=hdr_x, hdr_y=hdr_y, hdr_z=hdr_z, sens_pos_id_offset=sens_pos_id_offset)
            self.scan_pos_id_stridx = str_idx_ScanPosID
            self.scan_pos_id_endstridx = str_end_idx_ScanPosID

    # TODO: change to structure of occpyRIEGL: read input las files and link to scan positions/ trajectory first (how to do this for MLS/ULS?) -> allows us to error out/log early when position link is not properly found

    # TODO: ugly workaround for the case where a single laz file from a single TLS position should be run
    def define_sensor_pos_singlePos(self, scan_pos_id, x, y, z):
        d = {'ScanPos': scan_pos_id,
             'sensor_x': x, 'sensor_y': y, 'sensor_z': z}

        senspos = pd.DataFrame(data=d, index=[0])

        self.senspos = senspos


    def do_raytracing(self):
        """
        Perform ray tracing.

        This method processes either a directory of LAZ files (for TLS with known scan positions) or a single LAZ file 
        (Single TLS position or MLS/ULS with a trajectory, depends on self.is_mobile). 
        In the case of TLS, for each LAZ file, it extracts point positions and sensor positions, and 
        then performs ray tracing, accounting for single or multi-return pulse information. Multi-return handling 
        supports on-the-fly processing if the data is sorted by GPS time; otherwise, the full dataset must be loaded first.

        Raises
        ------
        FileNotFoundError
            If `self.laz_in_f` is not a valid file or directory.

        RuntimeWarning
            If multi-return data is detected but the LAZ file is not sorted by GPS time.

        """
        run_raytraycing_after_loading = False
        if os.path.isdir(self.laz_in_f):
            ## get list of laz files in input directory
            fCont = glob.glob(f"{self.laz_in_f}/*.laz")

            for f in fCont:
                # get scan position
                # lastdir_ind = f.rfind("/")
                scan_name = os.path.basename(f)
                # TODO: This is very specific to the test data and needs to be made generic!
                #end_of_scanID_idx = scan_name.rfind("_")
                scan_id = int(scan_name[self.scan_pos_id_stridx:self.scan_pos_id_endstridx])

                print(f"###############################")
                print(f"##### Processing {scan_name}...")
                print(f"###############################")

                scanpos_X = self.senspos.loc[self.senspos['ScanPos'] == scan_id, 'sensor_x'].values[0]
                scanpos_Y = self.senspos.loc[self.senspos['ScanPos'] == scan_id, 'sensor_y'].values[0]
                scanpos_Z = self.senspos.loc[self.senspos['ScanPos'] == scan_id, 'sensor_z'].values[0]



                # read in laz file
                tic = time.time()
                with laspy.open(f) as file:

                    count = 0
                    for points in file.chunk_iterator(self.points_per_iter):

                        # if self.single_return is not set, check data to find out if there are multiple returns stored per pulse
                        if self.single_return == None:
                            print(
                                f"TIPP: it was not specified if the point cloud stores single return data. to avoid any "
                                f"issues where there are only single return pulses in the first chunk, put there are "
                                f"multi return pulses in the whole dataset, please set the single_return variable within "
                                f"the define_sensor_pos function")
                            max_ret_num = np.max(points.return_number)
                            if max_ret_num > 1:
                                self.single_return = False
                            else:
                                self.single_return = True

                        if self.single_return:
                            sorted = True

                            x = points.x.copy()
                            y = points.y.copy()
                            z = points.z.copy()

                            sensor_x = np.ones(x.shape) * scanpos_X
                            sensor_y = np.ones(x.shape) * scanpos_Y
                            sensor_z = np.ones(x.shape) * scanpos_Z

                            gps_time = np.linspace(start=count + 1, stop=count + len(x), num=len(x), endpoint=True)

                            # run raytracing algorithme using singleReturnPulses version
                            #print("Do raytracing with all pulses in batch")
                            tic_r = time.time()
                            self.RayTr.doRaytracing_singleReturnPulses(x, y, z, sensor_x, sensor_y,
                                                                  sensor_z, gps_time)
                            toc_r = time.time()
                            #print("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))

                        else:
                            x = points.x.copy()
                            y = points.y.copy()
                            z = points.z.copy()
                            gps_time = points.gps_time.copy()
                            return_number = points.return_number.copy()
                            number_of_returns = points.number_of_returns.copy()

                            # check if gps_time is sorted
                            if count == 0:  # only check sort state in the first iteration of the for loop
                                if not is_sorted(gps_time):
                                    print(
                                        f"!!!!! input laz file is not sorted along gps_time. The algorithm will still run. However, the "
                                        f"performance will be greatly decreased as the entire content of the laz file has to be read into "
                                        f"the system memory. If you have multi return data, consider sorting your laz data first, e.g. using "
                                        f"LASTools lassort: lassort -i laz_in -gps_time -return_number -odix _sort -olaz -v !!!!")
                                    sorted = False
                                else:
                                    sorted = True

                            sensor_x = np.ones(gps_time.shape) * scanpos_X
                            sensor_y = np.ones(gps_time.shape) * scanpos_Y
                            sensor_z = np.ones(gps_time.shape) * scanpos_Z

                            self.RayTr.addPointData(x, y, z, sensor_x, sensor_y, sensor_z, gps_time, return_number,
                                               number_of_returns)

                            if sorted:  # only if pulses are sorted run raytracing now. Otherwise we have to read in the entire dataset first!
                                # Get report on pulse dataset - comment this out once everythin is working or TODO: add a verbose flag!
                                #self.RayTr.getPulseDatasetReport()

                                # run raytracing on added points
                                #print("Do raytracing with stored pulses")
                                tic_r = time.time()
                                self.RayTr.doRaytracing()
                                toc_r = time.time()
                                #print("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))

                                self.RayTr.clearPulseDataset()

                                # Check if traversed pulses have been deleted from map - comment this out once everything is working or TODO: add a verbose flag!
                                # RayTr.getPulseDatasetReport()

                        count = count + len(gps_time)

            toc = time.time()
            if sorted and not self.single_return:
                # optional: incomplete pulses can occur if the data has been filtered (either actively or during black box processing
                # of the processing software. We could actively turn the incomplete pulses into complete ones and do the raytracing
                # for them!
                print("convert incomplete pulses to complete ones - be cautious with that!")
                # RayTr.getPulseDatasetReport()
                self.RayTr.cleanUpPulseDataset()
                # RayTr.getPulseDatasetReport()
                print("Run raytracing for incomplete pulses")
                tic_r = time.time()
                self.RayTr.doRaytracing()
                toc_r = time.time()
                print("Time elapsed for raytracing incomplete pulses: {:.2f} seconds".format(toc_r - tic_r))
                print("Time elapsed for reading and raytracing entire data: {:.2f} seconds".format(toc_r - tic))
            elif not sorted and not self.single_return:
                print("Time elapsed for reading in data: {:.2f} seconds".format(toc - tic))

                # RayTr.getPulseDatasetReport()

                print("Clean up pulse dataset in order to handle incomplete pulses")
                self.RayTr.cleanUpPulseDataset()

                self.RayTr.getPulseDatasetReport()

                print("Do actual raytracing with all pulses")
                tic = time.time()
                self.RayTr.doRaytracing()
                toc = time.time()
                print("Time elapsed for raytracing: {:.2f} seconds".format(toc - tic))

        else: # if input is a single laz file
            with laspy.open(self.laz_in_f) as file:
                count = 0
                with tqdm(total=file.header.point_count, desc="Tracing Pulses...", unit="pulses") as pbar:
                    for points in file.chunk_iterator(points_per_iteration=self.points_per_iter):

                        # For performance we need to use copy
                        # so that the underlying arrays are contiguous
                        x = points.x.copy()
                        y = points.y.copy()
                        z = points.z.copy()
                        gps_time = points.gps_time.copy()
                        return_number = points.return_number.copy()
                        number_of_returns = points.number_of_returns.copy()

                        if np.max(
                                return_number) == 0:  # a not very nice hack for the special case where return_number and number_of_returns are all 0 for Horizon measurements - TODO: figure out why!
                            return_number[:] = 1
                            number_of_returns[:] = 1

                        # for the case of mobile acquisitions, inerpolate trajectory for gps_time
                        if self.is_mobile:
                            # call interpolate function for trajectory to extract sensor position for each gps_time
                            SensorPos = interpolate_traj(self.traj['time'], self.traj['sensor_x'], self.traj['sensor_y'],
                                                              self.traj['sensor_z'], gps_time)

                        else:
                            SensorPos = self.senspos

                            # TODO: figure out, why we have to repeat scan position to the shape of input point cloud and try to implement it, so that we only have to pass one position
                            SensorPos = pd.DataFrame(data={'ScanPos': np.ones(gps_time.shape) * self.senspos['ScanPos'].values[0],
                                                           'sensor_x': np.ones(gps_time.shape) * self.senspos['sensor_x'].values[0],
                                                           'sensor_y': np.ones(gps_time.shape) * self.senspos['sensor_y'].values[0],
                                                           'sensor_z': np.ones(gps_time.shape) * self.senspos['sensor_z'].values[0]})




                        if np.max(number_of_returns) == 1 or np.max(return_number) == 1:
                            run_raytraycing_after_loading = False
                            self.RayTr.doRaytracing_singleReturnPulses(x, y, z,  SensorPos['sensor_x'], SensorPos['sensor_y'],
                                                                       SensorPos['sensor_z'], gps_time)
                        else:
                            # check if gps_time is sorted
                            if count == 0:  # only check sort state in the first iteration of the for loop
                                if not is_sorted(gps_time):
                                    print(
                                        f"!!!!! input laz file is not sorted along gps_time. The algorithm will still run. However, the "
                                        f"performance will be greatly decreased as the entire content of the laz file has to be read into "
                                        f"the system memory. If you have multi return data, consider sorting your laz data first, e.g. using "
                                        f"LASTools lassort: lassort -i laz_in -gps_time -return_number -odix _sort -olaz -v !!!!")
                                    sorted = False
                                    run_raytraycing_after_loading=True


                                else:
                                    sorted = True
                                    run_raytraycing_after_loading=False


                            self.RayTr.addPointData(x, y, z, SensorPos['sensor_x'], SensorPos['sensor_y'],
                                                    SensorPos['sensor_z'],
                                                    gps_time, return_number, number_of_returns)

                            if sorted:  # only if pulses are sorted run raytracing now. Otherwise we have to read in the entire dataset first!
                                # Get report on pulse dataset - comment this out once everythin is working or TODO: add a verbose flag!
                                # self.RayTr.getPulseDatasetReport()

                                # run raytracing on added points
                                # print("Do raytracing with stored pulses")
                                tic_r = time.time()
                                self.RayTr.doRaytracing()
                                toc_r = time.time()
                                # print("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))

                                self.RayTr.clearPulseDataset() # clear out data that have been traced.

                                # Check if traversed pulses have been deleted from map - comment this out once everything is working or TODO: add a verbose flag!
                                # self.RayTr.getPulseDatasetReport()

                        count = count + len(gps_time)
                        pbar.update(len(points))

        if run_raytraycing_after_loading:
            self.RayTr.doRaytracing()

        self.get_raytracing_report()
        self.save_raytracing_output()

    def get_raytracing_report(self):
        """
        Print or log report on occlusion mapping statistics.
        """
        # Get report on traversal
        self.RayTr.reportOnTraversal()

    def save_raytracing_output(self):
        """
        Extract and save the outputs of the ray tracing process.

        This method performs the following steps:  
        1. Extracts the voxel-wise hit (`Nhit`), miss (`Nmiss`), and occlusion (`Nocc`) voxelgrids and saves as .npy in self.out_dir  
        2. Creates a voxel classification grid based on the `Nhit`, `Nmiss`, and `Nocc` values:  
        - 1 = observed (hit > 0)  
        - 2 = empty (miss > 0, hit == 0)  
        - 3 = occluded (occlusion > 0, hit == 0, miss == 0)  
        - 4 = unobserved (all three == 0)  
        3. Writes `.ply` files for all voxel outputs if `self.output_voxels` is True (takes long and creates large files)
        """
        print("Extracting Nhit")
        tic = time.time()
        self.Nhit = self.RayTr.getNhit()
        self.Nhit = np.array(self.Nhit, dtype=np.int32)

        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        print("Extracting Nocc")
        tic = time.time()
        self.Nocc = self.RayTr.getNocc()
        self.Nocc = np.array(self.Nocc, dtype=np.int32)

        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        print("Extracting Nmiss")
        tic = time.time()
        self.Nmiss = self.RayTr.getNmiss()
        self.Nmiss = np.array(self.Nmiss, dtype=np.int32)

        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        print("Saving Occlusion Outputs As .npy")
        tic = time.time()
        np.save(f"{self.out_dir}/Nhit.npy", self.Nhit)
        np.save(f"{self.out_dir}/Nmiss.npy", self.Nmiss)
        np.save(f"{self.out_dir}/Nocc.npy", self.Nocc)
        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        # Create Classification grid
        print("Classify Grid")
        tic = time.time()
        self.Classification = np.zeros((self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz']), dtype=int)

        self.Classification[np.logical_and.reduce((self.Nhit > 0, self.Nmiss >= 0, self.Nocc >= 0))] = 1  # voxels that were observed
        self.Classification[np.logical_and.reduce((self.Nhit == 0, self.Nmiss > 0, self.Nocc >= 0))] = 2  # voxels that are empty
        self.Classification[
            np.logical_and.reduce((self.Nhit == 0, self.Nmiss == 0, self.Nocc > 0))] = 3  # voxels that are hidden (occluded)
        self.Classification[np.logical_and.reduce((self.Nhit == 0, self.Nmiss == 0,
                                              self.Nocc == 0))] = 4  # voxels that were not observed # TODO: Figure out, why this overwrites voxels that are classified as occluded! -> this was because np.logical_and only takes in 2 arrays as input, not 3! use np.logical_and.reduce() for that!

        np.save(f"{self.out_dir}/Classification.npy", self.Classification)
        toc = time.time()
        print("Elapsed Time: " + str(toc - tic) + " seconds")

        # write ply file
        if self.output_voxels:
            print("Saving Occlusion Outputs As .ply")
            tic = time.time()
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nhit)
            ost.write_ply(f"{self.out_dir}/Nhit.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nmiss)
            ost.write_ply(f"{self.out_dir}/Nmiss.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nocc)
            ost.write_ply(f"{self.out_dir}/Nocc.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Classification)
            ost.write_ply(f"{self.out_dir}/Classification.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            self.occl = np.zeros(shape=self.Classification.shape)
            x4, y4, z4 = np.where(self.Classification == 4)
            self.occl[x4, y4, z4] = self.Classification[x4, y4, z4]
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.occl)
            ost.write_ply(f"{self.out_dir}/Occl.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            toc = time.time()
            print("Elapsed Time: " + str(toc - tic) + " seconds")


    def get_chm(self):
        if self.chm is None:
            print("No CHM was defined. To define CHM ")
        return self.chm

    def clean_up_RayTr(self):
        """
        Free up raytracer memory
        """
        del self.RayTr