#!/usr/bin/env python3
"""
Occpy

Drivers for handling RIEGL rdbx and rxp files

Written for pylidar_tls_canopy by John Armston, University of Maryland
Adapted for OccPy by Wout Cherlet, Ghent University
"""

try:
    import riegl_rdb
except ImportError:
    print('RIEGL RDBlib is not available')

try:
    import riegl_rxp
except ImportError:
    print('RIEGL RiVlib is not available')

import re
import sys
import json
import numpy as np


class RDBFile:
    """
    RIEGL RDBX point reader.

    This class wraps cpp riegl_rdb reader to provide acces .rdbx files. Can apply transforms and 
    filter points.
    """

    def __init__(self, filename, transform_file=None, pose_file=None, query_str=None):
        """
        Initialize RDBFile.

        Parameters
        ----------
        filename : str
            Path to .rdbx file.
        transform_file : str, optional
            Path to text file containing a 4x4 transform matrix.
            If provided, applied to coordinates.
        pose_file : str, optional
            Path to a JSON file with keys 'pitch', 'roll', 'yaw' (degrees) to
            construct a rotation matrix. Ignored if `transform_file` is provided.
        query_str : str or list of str, optional
            Filter expression(s) applied to point attributes prior to aligning with
            pulses. See `run_query` for syntax.
        """

        self.filename = filename
        if transform_file is not None:
            self.transform = read_transform_file(transform_file)
        elif pose_file is not None:
            with open(pose_file,'r') as f:
                pose = json.load(f)
            self.transform = calc_transform_matrix(pose['pitch'], pose['roll'], pose['yaw'])
        else:
            self.transform = None
        self.query_str = query_str
        self.read_file()

    def __enter__(self):
        return self
        
    def __exit__(self, type, value, traceback):
        pass
            
    def run_query(self, points):
        """
        Build a boolean mask by evaluating the configured query over a structured array.
        - Supported operators: <, >, ==, >=, <=, !=.
        - Query can be a string or a list of strings; all are ANDed.

        Parameters
        ----------
        points : numpy.ndarray
            Structured array of point attributes as returned by riegl_rdb.

        Returns
        -------
        numpy.ndarray of bool
            Boolean mask of shape (N,)
        """
        
        if not isinstance(self.query_str, list):
            self.query_str = [self.query_str]

        r = re.compile(r'((?:\(|))\s*([a-z_]+)\s*(<|>|==|>=|<=|!=)\s*([-+]?\d+(?:\.\d+)?)\s*((?:\)|))')
        valid = np.ones(points.shape[0], dtype=bool)
        for q in self.query_str:
            m = r.match(q)
            if m is not None:
                if m.group(2) in points.dtype.names:
                    q_str = f'(points[\'{m.group(2)}\'] {m.group(3)} {m.group(4)})'
                else:
                    msg = f'{m.group(2)} is not a valid point attribute name'
                    print(msg)
                    return
                try:
                    valid &= eval(q_str)
                except SyntaxError:
                    msg = f'{q} is not a valid query string'
                    print(msg)

        return valid

    def read_file(self):
        """
        Read file and get global stats
        """

        self.meta, points = riegl_rdb.readFile(self.filename)
        
        self.points = {}
        if self.query_str is not None:
            idx = np.lexsort((points['target_index'],points['scanline_idx'],points['scanline']))
            points = points[idx]
            valid = self.run_query(points)
            points = points[valid]
            self.points['target_index'],self.points['target_count'] = reindex_targets(points['target_index'],
                points['target_count'], points['scanline'], points['scanline_idx'])

        if self.transform is not None:
            x_t,y_t,z_t = apply_transformation(points['x'], points['y'], points['z'],
                points['x'].shape[0], self.transform)
            self.points['x'] = x_t + self.transform[3,0] 
            self.points['y'] = y_t + self.transform[3,1]
            self.points['z'] = z_t + self.transform[3,2]
            _, self.points['zenith'], self.points['azimuth'] = xyz2rza(x_t, y_t, z_t)
        else:
            _, self.points['zenith'], self.points['azimuth'] = xyz2rza(points['x'], points['y'], points['z']) 
        
        self.points['valid'] = np.ones(points.shape[0], dtype=bool)
        for name in points.dtype.names:
            if name not in self.points:
                self.points[name] = points[name]

        self.minc = 0
        self.maxc = np.max(self.points['scanline'])
        self.minr = 0
        self.maxr = np.max(self.points['scanline_idx'])
        self.max_range = np.max(self.points['range'])
        self.max_target_count = np.max(self.points['target_count'])

    def get_data(self, name):
        """
        Return a point attribute filtered by the current validity mask.
        Exits the process if the attribute is not available.

        Parameters
        ----------
        name : str
            Attribute name present in `self.points`.

        Returns
        -------
        numpy.ndarray
            The requested attribute array, filtered by `self.points['valid']`.

        """

        if name in self.points:
            data = self.points[name]
            valid = self.points['valid']
        else:
            print(f'{name:} is not a point attribute')
            sys.exit()

        return data[valid]

    def get_meta(self, key):
        """
        Return individual metadata item by key.

        Possible keys: 'riegl.pose_estimation', 'riegl.geo_tag', 'riegl.notch_filter', 'riegl.window_echo_correction',
        'riegl.detection_probability', 'riegl.angular_notch_filter', 'riegl.pulse_position_modulation',
        'riegl.noise_estimates', 'riegl.device', 'riegl.atmosphere', 'riegl.near_range_correction', 
        'riegl.time_base', 'riegl.scan_pattern', 'riegl.device_geometry', 'riegl.range_statistics', 
        'riegl.beam_geometry', 'riegl.reflectance_calculation', 'riegl.window_analysis',
        'riegl.mta_settings', 'riegl.point_attribute_groups'  

        Parameters
        ----------
        key : str
            Metadata key present in `self.meta`.

        Returns
        -------
        Any
            The JSON-parsed value for the given key.
        """

        return json.loads(self.meta[key])

class RXPFile:
    """
    RIEGL RXP point reader.

    This class wraps cpp riegl_rxp reader to provide access to .rxp files. Can apply transforms and 
    filter points and pulses.
    """

    def __init__(self, filename, transform_file=None, pose_file=None, query_str=None):
        """
        Initialize an RXPFile reader.

        Parameters
        ----------
        filename : str
            Path to .rxp file.
        transform_file : str, optional
            Path to text file containing a 4x4 transform matrix.
            If provided, applied to beam origins, directions, and points.
        pose_file : str, optional
            Path to a JSON file with keys 'pitch', 'roll', 'yaw' (degrees) to
            construct a rotation matrix. Ignored if `transform_file` is provided.
        query_str : str or list of str, optional
            Filter expression(s) applied to point attributes prior to aligning with
            pulses. See `run_query` for syntax.
        """

        self.filename = filename
        if transform_file is not None:
            self.transform = read_transform_file(transform_file)
        elif pose_file is not None:
            with open(pose_file,'r') as f:
                pose = json.load(f)
            self.transform = calc_transform_matrix(pose['pitch'], pose['roll'], pose['yaw'])
        else:
            self.transform = None
        self.query_str = query_str
        self.read_file()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def run_query(self, points):
        """
        Build a boolean mask by evaluating the configured query over a structured array.
        - Supported operators: <, >, ==, >=, <=, !=.
        - Query can be a string or a list of strings; all are ANDed.

        Parameters
        ----------
        points : numpy.ndarray
            Structured array of point attributes as returned by riegl_rdb.

        Returns
        -------
        numpy.ndarray of bool
            Boolean mask of shape (N,)
        """

        if not isinstance(self.query_str, list):
            self.query_str = [self.query_str]

        r = re.compile(r'((?:\(|))\s*([a-z_]+)\s*(<|>|==|>=|<=|!=)\s*([-+]?\d+(?:\.\d+)?)\s*((?:\)|))')
        valid = np.ones(points.shape[0], dtype=bool)
        for q in self.query_str:
            m = r.match(q)
            if m is not None:
                if m.group(2) in points.dtype.names:
                    q_str = f'(points[\'{m.group(2)}\'] {m.group(3)} {m.group(4)})' 
                else:
                    msg = f'{m.group(2)} is not a valid point attribute name'
                    print(msg)
                    return
                try:
                    valid &= eval(q_str)
                except SyntaxError:
                    msg = f'{q} is not a valid query string'
                    print(msg)
        
        return valid

    def read_file(self):
        """
        Read file and get global stats
        """

        self.meta, points, pulses = riegl_rxp.readFile(self.filename)

        self.points = {}
        self.pulses = {}
        if self.query_str is not None:
    
            valid = self.run_query(points)
            points = points[valid]
           
            target_count = np.repeat(pulses['target_count'], pulses['target_count'])
            scanline = np.repeat(pulses['scanline'], pulses['target_count'])
            scanline_idx = np.repeat(pulses['scanline_idx'], pulses['target_count']) 
           
            self.points['target_index'],new_target_count = reindex_targets(points['target_index'], 
                target_count[valid], scanline[valid], scanline_idx[valid])

            max_scanline_idx = np.max(pulses['scanline_idx'])
            pulse_id = pulses['scanline'] * max_scanline_idx + pulses['scanline_idx']
            point_id = scanline[valid] * max_scanline_idx + scanline_idx[valid]
            
            pulse_sort_idx = np.argsort(pulse_id)
            point_sort_idx = np.argsort(point_id)
            first = (self.points['target_index'][point_sort_idx] == 1)
            idx = np.searchsorted(pulse_id, point_id[point_sort_idx][first], sorter=pulse_sort_idx)
            self.pulses['target_count'] = np.zeros(pulses.shape[0], dtype=np.uint8)         
            self.pulses['target_count'][idx] = new_target_count[point_sort_idx][first]
             
            self.points['valid'] = np.ones(np.count_nonzero(valid), dtype=bool)
        else:
            self.pulses['target_count'] = pulses['target_count']
            self.points['valid'] = np.ones(points.shape[0], dtype=bool)

        if self.transform is None:
            if 'PITCH' in self.meta:
                self.transform = calc_transform_matrix(self.meta['PITCH'], self.meta['ROLL'], self.meta['YAW'])
        
        if self.transform is not None:
            x_t,y_t,z_t = apply_transformation(pulses['beam_direction_x'], pulses['beam_direction_y'], 
                pulses['beam_direction_z'], pulses['beam_direction_x'].shape[0], self.transform)
            self.pulses['beam_direction_x'] = x_t
            self.pulses['beam_direction_y'] = y_t
            self.pulses['beam_direction_z'] = z_t
            
            x_t_o,y_t_o,z_t_o = apply_transformation(pulses['beam_origin_x'], pulses['beam_origin_y'], 
                pulses['beam_origin_z'], pulses['beam_origin_x'].shape[0], self.transform)
            self.pulses['beam_origin_x'] = x_t_o + self.transform[3,0] 
            self.pulses['beam_origin_y'] = y_t_o + self.transform[3,1]
            self.pulses['beam_origin_z'] = z_t_o + self.transform[3,2]

            _, self.pulses['zenith'], self.pulses['azimuth'] = xyz2rza(x_t, y_t, z_t)
        else:
            _, self.pulses['zenith'], self.pulses['azimuth'] = xyz2rza(pulses['beam_direction_x'],
                pulses['beam_direction_y'], pulses['beam_direction_z'])
        self.pulses['valid'] = np.ones(pulses.shape[0], dtype=bool)

        if self.transform is not None:
            x_t,y_t,z_t = apply_transformation(points['x'], points['y'], points['z'], 
                points['x'].shape[0], self.transform)
            self.points['x'] = x_t + self.transform[3,0]
            self.points['y'] = y_t + self.transform[3,1]
            self.points['z'] = z_t + self.transform[3,2]

        for name in pulses.dtype.names:
            if name not in self.pulses:
                self.pulses[name] = pulses[name]
            
        for name in points.dtype.names:
            if name not in self.points:
                self.points[name] = points[name]

        self.minc = 0
        self.maxc = np.max(self.pulses['scanline'])
        self.minr = 0
        self.maxr = np.max(self.pulses['scanline_idx'])
        self.max_range = np.max(self.points['range'])
        self.max_target_count = np.max(self.pulses['target_count'])

    def get_points_by_pulse(self, names):
        """
        Reshape data as a number of pulses by max_target_count array
        Multiple point attributes are handled using a structured array

        Parameters
        ----------
        names : list of str
            Point attribute names to include.
        """

        dtype_list = []
        for name in names:
            t = self.points[name].dtype.str
            dtype_list.append((str(name), t, self.max_target_count))
        
        pulse_id = np.repeat(self.pulses['pulse_id'] - 1, self.pulses['target_count'])
        
        npulses = self.pulses['pulse_id'].shape[0]
        data = np.empty(npulses, dtype=dtype_list)
        for i in range(self.max_target_count):
            point_idx = self.points['target_index'] == i + 1
            for name,t,s in dtype_list:
                idx = pulse_id[point_idx]
                data[name][idx,i] = self.points[name][point_idx]

        return data[self.pulses['valid']]

    def get_data(self, name, return_as_point_attribute=False):
        """
        Return a pulse or point attribute, optionally broadcast to points.

        Parameters
        ----------
        name : str
            Attribute name present in `self.pulses` or `self.points`.
        return_as_point_attribute : bool, default False
            When `True` and `name` is a pulse attribute, repeat values per target
            count and return a point-aligned array. Ignored for point attributes.

        Returns
        -------
        numpy.ndarray
            The attribute array
        """

        if name in self.pulses:
            if return_as_point_attribute:
                data = np.repeat(self.pulses[name], self.pulses['target_count'])
                valid = self.points['valid']
            else:
                data = self.pulses[name]
                valid = self.pulses['valid']
        elif name in self.points:
            data = self.points[name]
            valid = self.points['valid']
        else:
            print(f'{name:} is not a pulse or point attribute')
            sys.exit()
        
        return data[valid]


def calc_transform_matrix(pitch, roll, yaw):
    """
    Construct a 4x4 rotation matrix from pitch, roll, yaw.

    Parameters
    ----------
    pitch : float
        Pitch angle in degrees.
    roll : float
        Roll angle in degrees.
    yaw : float
        Yaw angle in degrees. If NaN, treated as 0.0.

    Returns
    -------
    numpy.ndarray
        4x4 rotation matrix.
    """
    pitch = np.radians(pitch)
    pitch_mat = np.identity(4)
    pitch_mat[0,0] = np.cos(pitch)
    pitch_mat[0,2] = np.sin(pitch)
    pitch_mat[2,0] = -np.sin(pitch)
    pitch_mat[2,2] = np.cos(pitch)

    roll = np.radians(roll)
    roll_mat = np.identity(4)
    roll_mat[1,1] = np.cos(roll)
    roll_mat[1,2] = -np.sin(roll)
    roll_mat[2,1] = np.sin(roll)
    roll_mat[2,2] = np.cos(roll)

    yaw = np.radians(yaw)
    if np.isnan(yaw):
        yaw = 0.0
    yaw_mat = np.identity(4)
    yaw_mat[0,0] = np.cos(yaw)
    yaw_mat[0,1] = -np.sin(yaw)
    yaw_mat[1,0] = np.sin(yaw)
    yaw_mat[1,1] = np.cos(yaw)

    tmp_mat = yaw_mat.dot(pitch_mat)
    transform = tmp_mat.dot(roll_mat)

    return transform

def apply_transformation(x, y, z, size, transform_matrix, translate=False):
    """
    Apply a 4x4 transformation to 3D coordinates.

    Parameters
    ----------
    x, y, z : array-like
        Coordinate components of equal length N.
    size : int
        Number of points
    transform_matrix : numpy.ndarray
        4x4 transform matrix.
    translate : bool, default False
        If True, apply translation on top of rotation.

    Returns
    -------
    tuple of numpy.ndarray
        Transformed (x, y, z) arrays.
    """

    xyz = np.vstack((x, y, z)).T
    if translate:
        t = np.ones((size,1))
    else:
        t = np.zeros((size,1))

    xyz = np.concatenate((xyz, t), 1)
    xyz_t = np.dot(xyz, transform_matrix)

    return xyz_t[:,0],xyz_t[:,1],xyz_t[:,2]

def xyz2rza(x, y, z):
    """
    Convert Cartesian coordinates to spherical (range, zenith, azimuth).

    Parameters
    ----------
    x, y, z : array-like
        Cartesian components.

    Returns
    -------
    tuple of numpy.ndarray
        r : range (radius)
        theta : zenith angle in radians (0 at +Z)
        phi : azimuth angle in radians, wrapped to [0, 2π).
    """

    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    phi = np.arctan2(x, y)
    np.add(phi, 2*np.pi, out=phi, where=x < 0)

    return r, theta, phi

def read_transform_file(fn):
    """
    Read a 4x4 transform matrix from a text file.

    Parameters
    ----------
    fn : str
        Path to the transform file (whitespace-separated floats, row-major).

    Returns
    -------
    numpy.ndarray
        4x4 transform matrix (transposed to match expected orientation).
    """
    
    with open(fn, 'rb') as f:
        transform = np.loadtxt(f, delimiter=' ', dtype=np.float32)
    return transform.T

def reindex_targets(target_index, target_count, scanline, scanline_idx):
    """
    Reindex the target index and count
    Assumes the input data are time-sequential

    Parameters
    ----------
    target_index : numpy.ndarray
        Original per-point target index array of shape (N,).
    target_count : numpy.ndarray
        Per-point target count array of shape (N,).
    scanline : numpy.ndarray
        Per-point scanline identifier of shape (N,).
    scanline_idx : numpy.ndarray
        Per-point scanline index within the scanline of shape (N,).

    Returns
    -------
    tuple of numpy.ndarray
        new_target_index : numpy.ndarray
            Reindexed target indices, same shape and dtype as `target_index`.
        new_target_count : numpy.ndarray
            Recomputed per-point target counts, same shape and dtype as `target_count`.
    """
    
    new_target_index = np.ones_like(target_index)
    new_target_count = np.ones_like(target_count)

    n = 1
    for i in range(1, target_index.shape[0], 1):
        same_pulse = (scanline[i] == scanline[i-1]) & (scanline_idx[i] == scanline_idx[i-1])
        if same_pulse:
            new_target_index[i] = new_target_index[i-1] + 1 
            new_target_count[i-n:i+1] = n + 1
            n += 1
        else:
            n = 1
        
    return new_target_index,new_target_count