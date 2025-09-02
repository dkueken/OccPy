import os

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# TODO: temp, give as parameter
VOX_DIM = 0.1

def lineplot_plusplus(orientation = "horizontal", **kwargs):
    """
    Create an enhanced seaborn line plot with rotated axes.

    The function applies an affine transformation that rotates the plot by 90 degrees
    and flips the y-axis. It swaps the x- and y-axis labels accordingly.

    Parameters
    ----------
    orientation : str, optional
        Orientation of the plot (default is "horizontal").
    **kwargs
        Additional keyword arguments passed to seaborn.lineplot.

    Returns
    -------
    matplotlib.axes.Axes
        The transformed seaborn line plot axes.
    """
    line = sns.lineplot(**kwargs)

    r = Affine2D().scale(sx=1, sy=-1).rotate_deg(90)
    for x in line.images + line.lines + line.collections:
        trans = x.get_transform()
        x.set_transform(r+trans)
        if isinstance(x, PathCollection):
            transoff = x.get_offset_transform()
            x._transOffset = r+transoff

    old = line.axis()
    line.axis(old[2:4] + old[0:2])
    xlabel = line.get_xlabel()
    line.set_xlabel(line.get_ylabel())
    line.set_ylabel(xlabel)

    return line

def interactive_figure(output_dir, axis=0):

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