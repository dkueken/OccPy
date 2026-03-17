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
    
    vox_dim = settings["vox_dim"]
    plot_dim = settings["plot_dim"]
    # get min and max bounds
    min_bound = plot_dim[:3]
    max_bound = plot_dim[3:6]

    occmap = np.load(occmap_file)
    dims = occmap.shape
    # check if min_bound and max_bound are inside dims
    for i in range(3):
        if min_bound[i] < 0 or max_bound[i] > dims[i]:
            raise ValueError(f"min_bound and max_bound must be within occupancy map dimensions {dims}, got min_bound={min_bound}, max_bound={max_bound}")
    
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
        point_cloud = pv.PolyData(points_in_crop)
        plotter.add_points(point_cloud, style="points", point_size=point_size)

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

    vox_dim = settings["vox_dim"]
    plot_dim = settings["plot_dim"]
    # get min and max bounds
    min_bound = np.array(plot_dim[:3])
    max_bound = np.array(plot_dim[3:6])

    occmap = np.load(occmap_file)
    dims = occmap.shape
    # check if min_bound and max_bound are inside dims
    for i in range(3):
        if min_bound[i] < 0 or max_bound[i] > dims[i]:
            raise ValueError(f"min_bound and max_bound must be within occupancy map dimensions {dims}, got min_bound={min_bound}, max_bound={max_bound}")

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
        point_cloud = pv.PolyData(points_in_crop)
        plotter.add_points(point_cloud, style="points", point_size=point_size)

    # calculate center and radius for camera orbit
    min_bound_voxel = np.array(min_bound_voxel)
    max_bound_voxel = np.array(max_bound_voxel)
    center = (min_bound_voxel + max_bound_voxel) / 2
    bounds_size = max_bound_voxel - min_bound_voxel
    radius = np.linalg.norm(bounds_size) * distance_factor
    plotter.open_movie(opath, framerate=framerate)

    for i in range(n_frames):
        azimuth = i * (360 / n_frames)
        
        plotter.camera.position = (
            center[0] + radius * np.cos(np.radians(azimuth)) * np.cos(np.radians(elevation)),
            center[1] + radius * np.sin(np.radians(azimuth)) * np.cos(np.radians(elevation)),
            center[2] + radius * np.sin(np.radians(elevation))
        )
        
        plotter.camera.focal_point = center
        plotter.camera.up = (0, 0, 1)
        
        plotter.write_frame()

    plotter.close()



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
