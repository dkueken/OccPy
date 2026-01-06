"""
GPU-Accelerated 3D Voxel Tracer
===============================

A GPU-based implementation of the Amanatides & Woo ray traversal algorithm,
optimized for massive parallel ray tracing in 3D voxel grids.

This module is designed to be used as a standalone single-file library.
Simply copy this script into your project and import it.

Information
-----------
Author   : Junpei Kariyazono, Hideki Kobayashi, and Tomo'omi Kumagai
License  : MIT License
Version  : 0.1.1
DOI      : 10.5281/zenodo.18112352
github  : https://github.com/junpei-kariyazono/gpu-voxel-tracer

Citation
--------
If you use this software in your research, please cite it using the DOI above.
"""

import numpy as np
import math
from numba import cuda


# ===============================
# CUDA Kernel (Device Function)
# ===============================
@cuda.jit
def _voxel_traversal_kernel(starts, ends, matrix_shape, bounds, voxel_size,
                            point_counts, beam_counts, path_lengths):
    """
    Core Kernel: Executes the Amanatides & Woo algorithm on GPU.
    """
    i = cuda.grid(1)
    if i >= ends.shape[0]:
        return

    # Unpack bounds
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    sx, sy, sz = voxel_size

    # Ray Parameters
    x_position = starts[i, 0]
    y_position = starts[i, 1]
    z_position = starts[i, 2]

    dx = ends[i, 0] - starts[i, 0]
    dy = ends[i, 1] - starts[i, 1]
    dz = ends[i, 2] - starts[i, 2]

    norm = math.sqrt(dx * dx + dy * dy + dz * dz)
    if norm == 0:
        return

    # Direction steps
    step_x = 1 if dx > 0 else -1
    step_y = 1 if dy > 0 else -1
    step_z = 1 if dz > 0 else -1

    # Initial Voxel Indices
    x_index = int(math.floor((x_position - x_min) / sx))
    y_index = int(math.floor((y_position - y_min) / sy))
    z_index = int(math.floor((z_position - z_min) / sz))

    # Ray Parameter (t) Initialization
    if dx != 0:
        x_init_id = math.floor((x_position - x_min) / sx) if step_x > 0 else math.ceil((x_position - x_min) / sx)
        t_max_x = ((x_init_id + step_x) * sx + x_min - x_position) / dx
        t_delta_x = abs(sx / dx)
    else:
        t_max_x = float('inf')
        t_delta_x = float('inf')

    if dy != 0:
        y_init_id = math.floor((y_position - y_min) / sy) if step_y > 0 else math.ceil((y_position - y_min) / sy)
        t_max_y = ((y_init_id + step_y) * sy + y_min - y_position) / dy
        t_delta_y = abs(sy / dy)
    else:
        t_max_y = float('inf')
        t_delta_y = float('inf')

    if dz != 0:
        z_init_id = math.floor((z_position - z_min) / sz) if step_z > 0 else math.ceil((z_position - z_min) / sz)
        t_max_z = ((z_init_id + step_z) * sz + z_min - z_position) / dz
        t_delta_z = abs(sz / dz)
    else:
        t_max_z = float('inf')
        t_delta_z = float('inf')

    t_current = 0.0

    # DDA Loop
    while t_current < 1.0:
        # Determine next voxel boundary
        if t_max_x < t_max_y:
            if t_max_x < t_max_z:
                target_axis = 0
                dt = (1 - t_current) if t_max_x > 1 else (t_max_x - t_current)
            else:
                target_axis = 2
                dt = (1 - t_current) if t_max_z > 1 else (t_max_z - t_current)
        else:
            if t_max_y < t_max_z:
                target_axis = 1
                dt = (1 - t_current) if t_max_y > 1 else (t_max_y - t_current)
            else:
                target_axis = 2
                dt = (1 - t_current) if t_max_z > 1 else (t_max_z - t_current)

        # Update Grids if within bounds
        if ((0 <= x_index < matrix_shape[0])
                and (0 <= y_index < matrix_shape[1])
                and (0 <= z_index < matrix_shape[2])):

            cuda.atomic.add(beam_counts, (x_index, y_index, z_index), 1)
            cuda.atomic.add(path_lengths, (x_index, y_index, z_index), dt * norm)

            # Check if ray ends in this voxel (approximate)
            if min(t_max_x, t_max_y, t_max_z) >= 1.0:
                cuda.atomic.add(point_counts, (x_index, y_index, z_index), 1)

        # Advance
        if target_axis == 0:
            t_current = t_max_x
            x_index += step_x
            t_max_x += t_delta_x
        elif target_axis == 1:
            t_current = t_max_y
            y_index += step_y
            t_max_y += t_delta_y
        else:
            t_current = t_max_z
            z_index += step_z
            t_max_z += t_delta_z


# ===============================
# Host Class (Interface)
# ===============================
class Tracer:
    def __init__(self, bounds, voxel_size):
        """
        Initialize the Voxel Tracer.
        :param bounds: Tuple (xmin, xmax, ymin, ymax, zmin, zmax)
        :param voxel_size: Float or Tuple (sx, sy, sz)
        """
        self.bounds = bounds
        self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax = bounds

        if isinstance(voxel_size, (int, float)):
            self.voxel_size = (float(voxel_size),) * 3
        else:
            self.voxel_size = tuple(voxel_size)

        self.sx, self.sy, self.sz = self.voxel_size

        # Calculate Grid Dimensions
        self.dims = (
            int((self.xmax - self.xmin) / self.sx),
            int((self.ymax - self.ymin) / self.sy),
            int((self.zmax - self.zmin) / self.sz)
        )
        print(f"Initialized Tracer Grid: {self.dims}")

    def run(self, starts, ends):
        """
        Execute ray traversal on GPU.
        :param starts: numpy array [N, 3]
        :param ends: numpy array [N, 3]
        :return: dict containing 'point_counts', 'beam_counts', 'path_lengths' (numpy arrays)
        """
        n_rays = len(starts)

        # 1. Host -> Device
        d_starts = cuda.to_device(starts.astype(np.float64))
        d_ends = cuda.to_device(ends.astype(np.float64))

        # Output grids (Device)
        d_point_counts = cuda.to_device(np.zeros(self.dims, dtype=np.int32))
        d_beam_counts = cuda.to_device(np.zeros(self.dims, dtype=np.int32))
        d_path_lengths = cuda.to_device(np.zeros(self.dims, dtype=np.float64))

        # 2. Kernel Configuration
        threads_per_block = 128
        blocks_per_grid = (n_rays + threads_per_block - 1) // threads_per_block

        # 3. Launch Kernel
        _voxel_traversal_kernel[blocks_per_grid, threads_per_block](
            d_starts, d_ends, self.dims,
            self.bounds, self.voxel_size,
            d_point_counts, d_beam_counts, d_path_lengths
        )
        cuda.synchronize()

        # 4. Device -> Host
        return {
            "point_counts": d_point_counts.copy_to_host(),
            "beam_counts": d_beam_counts.copy_to_host(),
            "path_lengths": d_path_lengths.copy_to_host()
        }
