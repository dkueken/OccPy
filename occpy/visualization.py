import os

import numpy as np
import pandas as pd
import open3d as o3d
import pyvista as pv
import laspy
import json

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# TODO: Test!
def mpl_occpy_interactive(output_dir, axis=0):
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

    # TODO: temp, give as parameter
    VOX_DIM = 0.1

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


def occmap_vis_pyvista(occmap_file, 
                       min_bound_voxel, 
                       max_bound_voxel, 
                       config_file,
                       pointcloud_file=None,
                       opacity_occluded=0.2,
                       opacity_hit=0.1,
                       opacity_unobserved=0.2,
                       point_size=2):
    
    # TODO : test!!
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
    min_bound_voxel : array-like of float, shape (3,)
        Minimum XYZ voxel coordinates to visualize
    max_bound_voxel : array-like of float, shape (3,)
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
    """
    
    # Read occupancy map
    occmap = np.load(occmap_file)

    # read json config file
    with open(config_file) as file:
        settings = json.load(file)

    vox_dim = settings["voxel_size"]
    plot_dim = settings["plot_dim"]

    # get min and max bounds
    min_bound = plot_dim[:3]
    max_bound = plot_dim[3:6]
    
    # Read point cloud
    if pointcloud_file is not None:
        las = laspy.read(pointcloud_file)
        points = np.vstack((las.x, las.y, las.z)).transpose()

    dims = occmap.shape
    
    # check if min_bound and max_bound are inside dims
    for i in range(3):
        if min_bound[i] < 0 or max_bound[i] > dims[i]:
            raise ValueError(f"min_bound and max_bound must be within occupancy map dimensions {dims}, got min_bound={min_bound}, max_bound={max_bound}")
    
    # Collect bounding boxes for each voxel type
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
    
    # Create meshes
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
    
    # Add point cloud (cropped to visible region)
    if pointcloud_file is not None:
        min_pc_crop = min_bound + min_bound_voxel * vox_dim
        max_pc_crop = min_bound + max_bound_voxel * vox_dim
        print(f"Crop bounds: {min_pc_crop} to {max_pc_crop}")
        
        mask = np.all((points >= min_pc_crop) & (points <= max_pc_crop), axis=1)
        points_in_crop = points[mask]
        point_cloud = pv.PolyData(points_in_crop)
        plotter.add_points(point_cloud, style="points", point_size=point_size)
    
    plotter.show()


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


def pyvista_static_figure():
    """
    Placeholder for future PyVista static figure implementation.
    """
    pass