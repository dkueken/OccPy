
import numpy as np
import pandas as pd
import laspy
from scipy import interpolate
from tqdm import tqdm

from raytr import PyRaytracer

# raytracing helpers

def last_nonzero(arr, axis, invalid_val=-1):
    """
    Find the index of the last non-zero element along a specified axis.

    Parameters
    ----------
    arr : np.ndarray
        Input array to search.
    axis : int
        Axis along which to find the last non-zero element.
    invalid_val : int, optional
        Value to return if no non-zero elements are found (default is -1).

    Returns
    -------
    np.ndarray
        Indices of the last non-zero element along the specified axis.
        If none found, returns `invalid_val`.
    """

    mask = arr!=0
    val = arr.shape[axis] - np.flip(mask, axis=axis).argmax(axis=axis) - 1
    return np.where(mask.any(axis=axis), val, invalid_val)

def interpolate_traj(traj_time, traj_x, traj_y, traj_z, pts_gpstime):
    """
    Interpolate trajectory coordinates at specified GPS timestamps.

    Parameters
    ----------
    traj_time : array-like
        Known trajectory time stamps.
    traj_x : array-like
        Known x-coordinates of the trajectory.
    traj_y : array-like
        Known y-coordinates of the trajectory.
    traj_z : array-like
        Known z-coordinates of the trajectory.
    pts_gpstime : array-like
        GPS timestamps at which to interpolate the trajectory.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing interpolated 'time', 'sensor_x', 'sensor_y', and 'sensor_z' columns.
    """

    f_x = interpolate.interp1d(traj_time, traj_x, kind='linear', fill_value="extrapolate")
    sensor_x = f_x(pts_gpstime)
    f_y = interpolate.interp1d(traj_time, traj_y, kind='linear', fill_value="extrapolate")
    sensor_y = f_y(pts_gpstime)
    f_z = interpolate.interp1d(traj_time, traj_z, kind='linear', fill_value="extrapolate")
    sensor_z = f_z(pts_gpstime)

    d = {'time': pts_gpstime, 'sensor_x': sensor_x, 'sensor_y': sensor_y, 'sensor_z': sensor_z}

    df = pd.DataFrame(data=d)

    return df

def read_trajectory_file(path2traj, delimiter=" ", hdr_time='%time', hdr_x='x', hdr_y='y', hdr_z='z'):
    """
    Read a trajectory CSV file and extract trajectory data.

    Parameters
    ----------
    path2traj : str
        Path to the trajectory CSV file.
    delimiter : str, optional
        Delimiter used in the CSV file (default is space).
    hdr_time : str, optional
        Column name for time data (default is '%time').
    hdr_x : str, optional
        Column name for x-coordinate data (default is 'x').
    hdr_y : str, optional
        Column name for y-coordinate data (default is 'y').
    hdr_z : str, optional
        Column name for z-coordinate data (default is 'z').

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ['time', 'sensor_x', 'sensor_y', 'sensor_z'] representing the trajectory.
    
    Notes
    -----
    Defaults correspond to the GeoSLAM ZebHorizon scanner trajectory format.
    """

    traj_in = pd.read_csv(path2traj, sep=delimiter)

    d = {'time': traj_in[hdr_time], 'sensor_x': traj_in[hdr_x], 'sensor_y': traj_in[hdr_y],
         'sensor_z': traj_in[hdr_z]}

    traj = pd.DataFrame(data=d)

    return traj # retunr trajectory file if necessary

def read_sensorpos_file(path2senspos, delimiter=" ", hdr_scanpos_id='', hdr_x='', hdr_y='', hdr_z='', sens_pos_id_offset=0):
    """
    Read a sensor position CSV file and extract sensor position data.

    Parameters
    ----------
    path2senspos : str
        Path to the sensor position CSV file.
    delimiter : str, optional
        Delimiter used in the CSV file (default is space).
    hdr_scanpos_id : str, optional
        Column name for scan position ID (default is '').
    hdr_x : str, optional
        Column name for x-coordinate (default is '').
    hdr_y : str, optional
        Column name for y-coordinate (default is '').
    hdr_z : str, optional
        Column name for z-coordinate (default is '').
    sens_pos_id_offset : int, optional
        Offset to add to scan position IDs (default is 0).

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ['ScanPos', 'sensor_x', 'sensor_y', 'sensor_z'] representing sensor positions.
    """

    sens_pos_in = pd.read_csv(path2senspos, sep=delimiter)

    d = {'ScanPos': sens_pos_in[hdr_scanpos_id]+sens_pos_id_offset,
         'sensor_x': sens_pos_in[hdr_x], 'sensor_y': sens_pos_in[hdr_y], 'sensor_z': sens_pos_in[hdr_z]}

    senspos = pd.DataFrame(data=d)

    return senspos

def filterPointsIntersectingBox(laz_in, laz_out, min_bound, max_bound,sensor_pos=None, traj_in=None, points_per_iter=100000):
    """
    Filter points from a LAS/LAZ file whose pulses intersect a defined 3D bounding box.

    This filter considers the sensor or trajectory position to compute which pulses intersect
    the bounding box defined by `min_bound` and `max_bound`.

    Parameters
    ----------
    laz_in : str
        Path to the input LAS/LAZ file.
    laz_out : str
        Path to the output LAS/LAZ file where filtered points will be written.
    min_bound : tuple or list of float
        Minimum bounding box coordinates in the order (min_y, min_x, min_z).
    max_bound : tuple or list of float
        Maximum bounding box coordinates in the order (max_y, max_x, max_z).
    sensor_pos : pandas.DataFrame or None, optional
        Sensor position DataFrame with columns ['ScanPos', 'sensor_x', 'sensor_y', 'sensor_z'].
        Required if `traj_in` is not provided.
    traj_in : pandas.DataFrame or None, optional
        Trajectory DataFrame with columns ['time', 'sensor_x', 'sensor_y', 'sensor_z'].
        If provided, assumes a mobile platform.
    points_per_iter : int, optional
        Number of points to process per chunk iteration (default is 100000).

    Returns
    -------
    list
        List of GPS times for pulses intersecting the bounding box.

    Notes
    -----
    This filter does not work as intended for solid-state scanners (e.g., GeoSLAM ZebHorizon, FARO ORBIS)
    where GPSTime is not a unique pulse identifier. To properly handle these,
    pulse shooting directions should also be accounted for. This is partially implemented
    on the C++ side but not yet passed to Python.
    """

    pulses_intersecting = []
    #TODO: check if data has already been loaded from an initial do_raytracing run (and pulse dataset is complete).
    # If this is the case, the data does not need to be loaded and could be used directly
    if traj_in is not None:
        is_mobile = True
        traj = traj_in
    else:
        is_mobile = False
        sens_pos = sensor_pos

    with laspy.open(laz_in) as file:
        hdr = file.header

        writer = laspy.open(laz_out, mode="w", header=hdr) # create output laz file

        with tqdm(total=hdr.point_count, desc="Filtering points", unit="points") as pbar:
            for points in file.chunk_iterator(points_per_iteration=points_per_iter):

                # TODO: find another solution for this.
                if np.max(
                        points.return_number) == 0:  # a not very nice hack for the special case where return_number and number_of_returns are all 0 for Horizon measurements - TODO: figure out why!
                    points.return_number[:] = 1
                    points.number_of_returns[:] = 1

                # we only need first returns
                first_ret_ind = points.return_number == 1

                gps_time = points.gps_time.copy()[first_ret_ind]
                x = points.x.copy()[first_ret_ind]
                y = points.y.copy()[first_ret_ind]
                z = points.z.copy()[first_ret_ind]

                if is_mobile:
                    SensorPos = interpolate_traj(traj['time'], traj['sensor_x'],
                                                      traj['sensor_y'],
                                                      traj['sensor_z'], gps_time)

                    time = SensorPos['time'].to_numpy()
                    sensor_x = SensorPos['sensor_x'].to_numpy()
                    sensor_y = SensorPos['sensor_y'].to_numpy()
                    sensor_z = SensorPos['sensor_z'].to_numpy()

                else:
                    time = gps_time
                    sensor_x = np.ones(time.shape) * sens_pos['sensor_x']
                    sensor_y = np.ones(time.shape) * sens_pos['sensor_y']
                    sensor_z = np.ones(time.shape) * sens_pos['sensor_z']

                vmin = min_bound
                vmax = max_bound


                raytr = PyRaytracer()

                GPSTimes_intersect = raytr.getPulsesIntersectingBox(x, y, z, sensor_x, sensor_y, sensor_z, time, vmin, vmax)
                pulses_intersecting.extend(GPSTimes_intersect)
                # print(GPSTimes_intersect)
                del raytr # to free up some memory

                pulses2take = np.isin(time, GPSTimes_intersect)



                # write points with flag=True to new laz file
                writer.write_points(points[pulses2take])

                pbar.update(len(points))

        writer.close()

    return pulses_intersecting



# voxelgrid to ply helpers

def prepare_ply(vox_dim, PlotDim, data):
    """
    Convert a 3D voxel grid into mesh data for PLY export.

    Parameters
    ----------
    vox_dim : float
        Size of a single voxel in meters.
    PlotDim : dict
        Dictionary with plot bounding box coordinates. Must contain keys:
        'minX', 'maxX', 'minY', 'maxY', 'minZ', 'maxZ'.
    data : np.ndarray
        3D voxel grid (e.g., Nhit, Nmiss, Nocc, Classification).

    Returns
    -------
    verts : np.ndarray
        Array of mesh vertices with shape (N, 4), where the last column is the voxel value.
    faces : np.ndarray
        Array of triangular face indices with shape (M, 3).
    """

    data = array3Dto2D(data, vox_dim, PlotDim)
    # Generate mesh data
    verts, faces = generate_mesh_data(data, vox_dim, PlotDim)

    return verts, faces

def array3Dto2D(data, vox_dim, PlotDim):
    """
    Convert a 3D voxel grid into a 2D array of voxel center coordinates and values.

    Parameters
    ----------
    data : np.ndarray
        3D voxel grid (e.g., occupancy or classification values).
    vox_dim : float
        Size of each voxel in meters.
    PlotDim : dict
        Dictionary defining the spatial extent of the grid.
        Must include 'minX', 'maxX', 'minY', 'maxY', 'minZ', 'maxZ'.

    Returns
    -------
    np.ndarray
        2D array of shape (N, 4), where each row contains [x, y, z, value]
        for voxels with values greater than zero.
    """

    x = np.arange(PlotDim['minX'], PlotDim['maxX'], vox_dim)
    y = np.arange(PlotDim['minY'], PlotDim['maxY'], vox_dim)
    z = np.arange(PlotDim['minZ'], PlotDim['maxZ'], vox_dim)

    x, y, z = np.meshgrid(x, y, z)
    x, y, z, data = x.flatten(), y.flatten(), z.flatten(), data.flatten()
    data = np.c_[x,y,z,data]
    mask = np.ma.masked_where(data[:,-1]>0, data[:,-1]).mask
    data = data[mask]
    return data

def calculate_voxel_corners(data, vox_dim):
    """
    Compute the 3D corner coordinates of each voxel for mesh generation.

    Parameters
    ----------
    data : np.ndarray
        2D array of shape (N, 4), where each row contains [x, y, z, value]
        representing the center of a voxel and its scalar value.
    vox_dim : float
        Length of the voxel edge in meters (assumed cubic).
    Returns
    -------
    np.ndarray
        Array of shape (N * 8, 4), where each group of 8 rows corresponds to the
        corners of one voxel, and each row contains [x, y, z, value].
    """

    offset = vox_dim/2.
    xyz_offsets = np.array([
        [-offset, -offset, -offset],
        [offset, -offset, -offset],
        [offset, offset, -offset],
        [-offset, offset, -offset],
        [-offset, -offset, offset],
        [offset, -offset, offset],
        [offset, offset, offset],
        [-offset, offset, offset]
    ])

    # Fetch scalars array and repeat each line the number of lines in xyz_offsets
    scalars = data[:,3:]
    scalars = np.repeat(scalars, xyz_offsets.shape[0], axis=0)

    # Add offsets for each points
    xyz = data[:,:3][:, np.newaxis] + xyz_offsets
    xyz = np.vstack(xyz)

    # Concatenate the points with their respective scalars
    xyz_scalars = np.c_[xyz, scalars]

    return xyz_scalars

def generate_mesh_data(data, vox_dim, PlotDim):
    """
    Generate mesh vertex and face data from a 2D array of voxel center coordinates and values.

    Parameters
    ----------
    data : np.ndarray
        2D array of shape (N, 4), where each row contains [x, y, z, value]
        representing the center and value of a filled voxel.
    vox_dim : float
        Size of each voxel in meters.
    PlotDim : dict
        Dictionary defining the spatial extent of the voxel grid.
        Must include 'minX', 'maxX', 'minY', 'maxY', 'minZ', 'maxZ'.

    Returns
    -------
    verts : np.ndarray
        Array of vertex coordinates and associated voxel values.
    faces : np.ndarray
        Array of face indices defining triangles for visualization (PLY format).
    """
    
    # Calculate voxel corners
    corners = calculate_voxel_corners(data, vox_dim)

    # Define the face indices for a cube
    face_indices_0 = np.array([
        [0, 1, 2], [0, 2, 3],  # Front face
        [4, 5, 6], [4, 6, 7],  # Back face
        [0, 1, 5], [0, 5, 4],  # Bottom face
        [2, 3, 7], [2, 7, 6],  # Top face
        [1, 2, 6], [1, 6, 5],  # Right face
        [3, 0, 4], [3, 4, 7]   # Left face
    ])

    # Total number of voxels
    num_voxels = int(corners.shape[0] / 8)

    # Repeat face indices for each voxel
    faces = np.tile(face_indices_0, (num_voxels, 1))

    # Offset face indices for each voxel
    offset = np.arange(0, num_voxels * 8, 8).repeat(12)
    faces += offset[:, np.newaxis]

    return corners, faces
