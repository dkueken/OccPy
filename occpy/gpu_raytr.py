import numpy as np
import math
from occpy.gpu_voxel_tracer import Tracer


def _distance(p1, p2):
    return math.sqrt(sum((a - b)**2 for a, b in zip(p1, p2)))


class GPURaytracerWrapper:
    def __init__(self):
        self.tracer = None

        # Buffers for Ray data
        self.buffer_points = []  # (x, y, z) - used for Hit calculation
        self.buffer_sensor = []  # (sx, sy, sz)
        self.buffer_meta = []  # (gps_time, return_num, num_returns)

        # Grid definition
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.min_bound = np.zeros(3)
        self.voxel_size = 1.0
        self.accum_initialized = False

        # Accumulation grids (stored on CPU)
        self.Nhit_accum = None
        self.Nmiss_accum = None
        self.Nocc_accum = None

        # extension constant from ends to far-points for occlusion calculation
        self.max_distance = 100


    def defineGrid(self, minBound, maxBound, nx, ny, nz, vox_dim):
        print(f"[GPU] Initializing Grid: {nx}x{ny}x{nz}, Voxel: {vox_dim}m")
        self.nx, self.ny, self.nz = nx, ny, nz
        self.voxel_size = vox_dim
        self.min_bound = np.array(minBound, dtype=np.float64)

        # Initialize Tracer (used for Miss/Occ traversal)
        bounds = (minBound[0], maxBound[0], minBound[1], maxBound[1], minBound[2], maxBound[2])

        self.max_distance = _distance((minBound[0],minBound[1], minBound[2]), (maxBound[0],maxBound[1], maxBound[2]))

        self.tracer = Tracer(bounds, vox_dim)

        self._init_accumulators()


    def _init_accumulators(self):
        self.Nhit_accum = np.zeros((self.nx, self.ny, self.nz), dtype=np.int32)
        self.Nmiss_accum = np.zeros((self.nx, self.ny, self.nz), dtype=np.int32)
        self.Nocc_accum = np.zeros((self.nx, self.ny, self.nz), dtype=np.int32)
        self.accum_initialized = True


    def addPointData(self, x, y, z, sensor_x, sensor_y, sensor_z, gps_time, return_number, number_of_returns):
        """
        Adds data to the buffer. Kept as lists and concatenated during execution.
        """
        # Point coordinates
        pts = np.stack([x, y, z], axis=1).astype(np.float64)
        self.buffer_points.append(pts)

        # Sensor coordinates
        sens = np.stack([sensor_x, sensor_y, sensor_z], axis=1).astype(np.float64)
        self.buffer_sensor.append(sens)

        # Metadata (required for pulse identification)
        meta = np.stack([gps_time, return_number, number_of_returns], axis=1).astype(np.float64)
        self.buffer_meta.append(meta)


    def doRaytracing(self):

        if not self.buffer_points:
            return

        # 1. Concatenate and prepare data
        points_all = np.concatenate(self.buffer_points, axis=0)
        sensors_all = np.concatenate(self.buffer_sensor, axis=0)
        meta_all = np.concatenate(self.buffer_meta, axis=0)

        # Clear buffers
        self.clearPulseDataset()

        print(f"[GPU] Processing batch: {len(points_all)} points...")


        # -------------------------------------------------
        # Step 1: Calculate Nhit (Point Rasterization)
        # -------------------------------------------------

        # 1. Convert coords to index
        #    idx = floor((pos - min) / voxel_size)
        norm_coords = (points_all - self.min_bound) / self.voxel_size
        voxel_indices = np.floor(norm_coords).astype(np.int32)

        # 2. Filter out-of-bounds
        mask = (
                (voxel_indices[:, 0] >= 0) & (voxel_indices[:, 0] < self.nx) &
                (voxel_indices[:, 1] >= 0) & (voxel_indices[:, 1] < self.ny) &
                (voxel_indices[:, 2] >= 0) & (voxel_indices[:, 2] < self.nz)
        )
        valid_indices = voxel_indices[mask]

        # 3. Accumulate
        current_nhit = np.zeros((self.nx, self.ny, self.nz), dtype=np.int32)
        np.add.at(current_nhit, (valid_indices[:, 0], valid_indices[:, 1], valid_indices[:, 2]), 1)

        self.Nhit_accum += current_nhit


        # -------------------------------------------------
        # Step 2: Define Rays per Pulse (Sensor -> LastReturn)
        # -------------------------------------------------
        # Filter for "Last Return" only.
        # Assumption: Using (return_number == number_of_returns) to identify the last return.

        # gps_time = meta_all[:, 0] # Not used for filtering in this version
        ret_num = meta_all[:, 1]
        num_rets = meta_all[:, 2]

        # Index mask for Last Returns
        is_last = (ret_num == num_rets)

        # Ray start/end points for Miss/Occ calculation
        ray_starts = sensors_all[is_last]
        ray_ends = points_all[is_last]

        if len(ray_starts) == 0:
            return

        print(f"[GPU] Tracing {len(ray_starts)} rays (unique pulses)...")

        # -------------------------------------------------
        # Step 3: Nmiss (Sensor -> LastReturn)
        # -------------------------------------------------
        # Trace from Sensor to LastReturn
        # Note: Tracer counts all traversed voxels including the end voxel.
        results_miss = self.tracer.run(ray_starts, ray_ends)
        grid_traversed = results_miss['beam_counts']

        # We only correct based on the hits from the current batch to approximate the ray path logic.
        current_nmiss = grid_traversed - current_nhit
        current_nmiss[current_nmiss < 0] = 0

        self.Nmiss_accum += current_nmiss.astype(np.int32)

        # -------------------------------------------------
        # Step 4: Nocc (LastReturn -> Boundary) //TODO: this is not efficient, we should limit the raytrace within bounds.
        # -------------------------------------------------
        # Extend rays from LastReturn further out
        directions = ray_ends - ray_starts
        # Extend significantly to reach grid boundary (e.g., *1000)
        ray_occ_ends = ray_ends + directions * self.max_distance

        # Trace LastReturn -> Infinity
        results_occ = self.tracer.run(ray_ends, ray_occ_ends)
        grid_occ = results_occ['beam_counts']


        current_nocc = grid_occ.astype(np.int32)
        # Mask out voxels that have hits (approximating "start after last return")
        current_nocc[current_nhit > 0] = 0

        self.Nocc_accum += current_nocc


    def doRaytracing_singleReturnPulses(self, x, y, z, sensor_x, sensor_y, sensor_z, gps_time):
        # Prepare data
        starts = np.stack([sensor_x, sensor_y, sensor_z], axis=1).astype(np.float64)
        ends = np.stack([x, y, z], axis=1).astype(np.float64)

        # execute trace
        results = self.tracer.run(starts, ends)
        current_nhit = results['point_counts']
        self.Nhit_accum += current_nhit

        current_nmiss = results['beam_counts'] - current_nhit
        current_nmiss[current_nmiss < 0] = 0
        self.Nmiss_accum += current_nmiss.astype(np.int32)

        # 3. Occ calculation (Hit -> Far) //TODO: this is not efficient, we should limit the raytrace within bounds.
        directions = ends - starts
        far_ends = ends + directions * self.max_distance
        results_occ = self.tracer.run(ends, far_ends)
        current_nocc = results_occ['beam_counts'].astype(np.int32)
        current_nocc[current_nhit > 0] = 0  # Do not mark Hit voxels as Occ
        self.Nocc_accum += current_nocc


    # ---------------------------------------------------------
    # Data Management & Getters (Compatible with C++ bindings)
    # ---------------------------------------------------------
    def clearPulseDataset(self):
        self.buffer_points = []
        self.buffer_sensor = []
        self.buffer_meta = []

    def cleanUpPulseDataset(self):
        self.clearPulseDataset()

    def getPulseDatasetReport(self):
        pass

    def reportOnTraversal(self):
        print(f"[GPU Report] Hits: {np.sum(self.Nhit_accum)}, Misses: {np.sum(self.Nmiss_accum)}, Occ: {np.sum(self.Nocc_accum)}")

    def getNhit(self): return self.Nhit_accum
    def getNmiss(self): return self.Nmiss_accum
    def getNocc(self): return self.Nocc_accum

    # Dummy methods for compatibility
    def getGridDimensions(self): return [0, self.nx, 0, self.ny]
    def getGridOrigin(self): return [0, 0, 0]
    def get_num_traversed_pulses(self): return 0
    def get_total_pulses_in_dataset(self): return 0
    def get_num_registered_hits(self): return np.sum(self.Nhit_accum)
    def get_num_echoes_outside(self): return 0
    def get_num_missing_returns(self): return 0
    def get_num_pulses_no_intersection(self): return 0