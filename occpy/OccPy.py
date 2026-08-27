import logging
import time
import os
import glob
import json
from contextlib import contextmanager

import numpy as np
import pandas as pd
import laspy
import OSToolBox as ost
from tqdm import tqdm

from raytr import PyRaytracer
from occpy.util import prepare_ply, read_trajectory_file, read_sensorpos_file, interpolate_traj
from occpy.pulse_util import (
    ReturnMode,
    ScanJob,
    TraceMode,
    chunk_sortedness,
    constant_sensor_positions,
    detect_return_mode,
    normalise_chunk
)


class TqdmLoggingHandler(logging.StreamHandler):
    """
    Log handler that plays nicely with a live tqdm progress bar.

    A plain StreamHandler writes straight to the stream while tqdm is redrawing
    on the same one, so log lines land in the middle of the bar and the bar is
    reprinted on the next update. ``tqdm.write`` clears the bar, emits the line,
    and redraws it, which keeps both readable.

    Falls back to normal StreamHandler behaviour when no bar is active, so this
    is safe for non-chunked runs and for logging outside the read loop.
    """
    def emit(self, record):
        try:
            tqdm.write(self.format(record), file=self.stream)
            self.flush()
        except RecursionError:
            raise
        except Exception:
            self.handleError(record)


class OccPy:

    def __init__(self, config=None, config_file=None):
        """
        initialize OccPy object

        Parameters in config file:  
        Must include:  
            - 'laz_in' : path to single .las/.laz file or directory with multiple .las/.laz files
            - 'vox_dim' : voxel size in meters  
            - 'plot_dim': grid for occlusion mapping: [minX, minY, minZ, maxX, maxY, maxZ]  
        Optional parameters:  
            - 'out_dir' : output directory (default: ./output)
            - 'output_voxels' : whether to export `.ply` voxel grids (default: False)
            - 'verbose': set logging level  (default: False)
            - 'debug': set logging level (default: False)
            - 'lower_threshold': lower threshold above ground to exclude from occlusion mapping in voxels (default: 0)
            - 'points_per_iter': number of points read in from las/laz file in each iteration (default: 10000000)
            - 'delimiter': csv delimiter for scan position file (default: ",")
            - 'root_folder': if given, will assume other paths are relative to this root folder and will prepend it to the paths (default: None)
            - 'is_mobile': whether the acquisition is mobile (MLS/ULS) or static (TLS) (default: False)
            - 'single_return': whether the data is single return or multi return data. If omitted (default), it is detected from the data before ray tracing starts: from the LAS header when trustworthy, otherwise by scanning the return-number fields.
            - 'check_returns_all_files': when laz_in is a directory, probe every file for its return mode instead of only the first (default: False)
            - 'str_idxs_ScanPosID': string indices of where the scan position identifier is written in the las/laz file name. If not given, will use file name as ID (without extension) (default: None)
            - 'cleanup_incomplete_pulses': whether incomplete pulses should be cleaned so there are no missing returns (e.g. a pulse with number_of_return==3 only has 2 returns with the same GPSTime. If set to True, return number and number of return values for this pulse is manually changed to make it appear complete. (default: False)) CAUTION: This can result in unexpected behavior. If set to True we recommend to set 'verbose' and 'debug' to True
            - 'move_senspos_to_collinearity': There are occasions where the sensor position is not on a line built up by all returns (if multiple returns). In this case one could force collinearity with this flag. Be aware of potential caveats when using this flag. default=False
        Parameters
        ----------
        config : dict, optional
            Configuration dictionary containing processing parameters.
        config_file : str, optional
            Path to a JSON configuration file containing processing parameters.
            If both config and config_file are provided, config takes precedence.
        """

        if config is None and config_file is None:
            raise ValueError("Either 'config' or 'config_file' must be provided.")

        if config is None:
            with open(config_file) as f:
                config = json.load(f)
        elif not isinstance(config, dict):
            raise TypeError("'config' must be a dict when provided.")

        # Keep an internal copy of the input config for record keeping.
        config = dict(config)

        necessary_args = ["laz_in", "vox_dim", "plot_dim"]

        missing = []
        for key in necessary_args:
            if key not in config:
                missing.append(key)

        if len(missing) > 0:
            raise ValueError(f"Missing necessary arguments in config file: {missing}")
        
        self.laz_in = config["laz_in"]
        self.vox_dim = config["vox_dim"]
        self.plot_dim = config["plot_dim"]
        
        optional_args = ["out_dir", "is_mobile", "output_voxels", "verbose", "debug", "lower_threshold", "points_per_iter", "delimiter", "root_folder", "single_return", "str_idxs_ScanPosID", "cleanup_incomplete_pulses", "move_senspos_to_collinearity", "check_returns_all_files"]
        
        print(f"INFO: optional arguments: {optional_args}")

        # set optional arguments with defaults if not provided
        self.out_dir = config.get("out_dir", os.path.join(os.getcwd(), "output"))
        self.output_voxels = config.get("output_voxels", False)
        self.verbose = config.get("verbose", False)
        self.debug = config.get("debug", False)
        self.lower_threshold = config.get("lower_threshold", 0)
        self.points_per_iter = config.get("points_per_iter", 10000000)
        self.root_folder = config.get("root_folder", None)
        self.is_mobile = config.get("is_mobile", False)
        self.single_return = config.get("single_return", None) # None -> detect from data
        self.str_idxs_ScanPosID = config.get("str_idxs_ScanPosID", None)
        self.cleanup_incomplete_pulses = config.get("cleanup_incomplete_pulses", False)
        self.move_senspos_to_collinearity = config.get("move_senspos_to_collinearity", False)
        self.check_returns_all_files = config.get("check_returns_all_files", False)
        # Resolved at run time, not here; recorded in run_summary.json.
        self.returns_all_zero = None
        self.return_mode_result = None
        self.return_mode_probed_files = None
        self.single_return_source = None
        self.trace_mode = None
        self.traversal_stats = {}

        # config logging 
        if self.debug:
            logging_level = logging.DEBUG
        elif self.verbose:
            logging_level = logging.INFO
        else:
            logging_level = logging.WARNING
        self.logger = logging.getLogger('occpy_logger')
        self.logger.setLevel(logging_level)
        self.logger.propagate = False
        if self.logger.handlers:
            self.logger.handlers.clear()
        console_handler = TqdmLoggingHandler()
        console_handler.setLevel(logging_level)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

        # config paths
        if self.root_folder is not None and self.root_folder != "":
            self.root_folder = os.path.abspath(self.root_folder)
            if not os.path.exists(self.root_folder):
                raise ValueError(f"Root folder {self.root_folder} does not exist.")
            self.logger.info(f"Prepending root folder {self.root_folder} to input and output paths.")
            self.out_dir = os.path.join(self.root_folder, self.out_dir)
            self.laz_in = os.path.join(self.root_folder, self.laz_in)

        if not os.path.exists(self.laz_in):
            raise ValueError(f"Input path {self.laz_in} does not exist.")

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        # configure occpy input
        self.TotalVolume = 0
        self.TotalOcclusion = 0
        self.OcclFrac2D = None
        self.Nhit = None
        self.Nmiss = None
        self.Nocc = None
        self.Classification = None
        self.dtm = None
        self.dsm = None
        self.chm = None
        self.sens_pos_initialized = False
        self.scans_linked = None

        # ensure extents are exactly divisible by vox_dim by extending max bounds if needed.
        self.plot_dim, warnings = self.align_plot_dim_to_voxel_size(self.plot_dim, self.vox_dim)
        for msg in warnings:
            self.logger.warning(msg)
        if len(warnings) > 0:
            # write the adapted max bound to the config for record keeping
            config["adapted_max_bound"] = self.plot_dim[3:6]

        # initialize RayTr Object and define grid
        self.RayTr = PyRaytracer()        
        self.PlotDim = dict(minX=self.plot_dim[0],
                                maxX=self.plot_dim[3],
                                minY=self.plot_dim[1],
                                maxY=self.plot_dim[4],
                                minZ=self.plot_dim[2],
                                maxZ=self.plot_dim[5])
        
        self.grid_dim = dict(nx=int((self.PlotDim['maxX'] - self.PlotDim['minX']) / self.vox_dim),
                     ny=int((self.PlotDim['maxY'] - self.PlotDim['minY']) / self.vox_dim),
                     nz=int((self.PlotDim['maxZ'] - self.PlotDim['minZ']) / self.vox_dim))
                     
        minBound = np.asarray([self.PlotDim['minX'], self.PlotDim['minY'], self.PlotDim['minZ']], dtype=np.float64)
        maxBound = np.asarray([self.PlotDim['maxX'], self.PlotDim['maxY'], self.PlotDim['maxZ']], dtype=np.float64)
        self.RayTr.defineGrid(minBound, maxBound, self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz'], self.vox_dim)
        
        # Write config to output directory for record keeping.
        with open(os.path.join(self.out_dir, "config.json"), "w") as to:
            json.dump(config, to)

        return

    def define_sensor_pos(self, path2file, delimiter=" ", hdr_time='%time', hdr_scanpos_id='', hdr_x='x', hdr_y='y', hdr_z='z', sens_pos_id_offset=0):
        """
        defines sensor position based on the provided csv file. CSV file needs to include
        Parameters
        ----------
        path2file: str [mandatory]
            path to csv file with sensor position information
        delimiter: str [default: " "]
            csv delimiter
        hdr_time: str
            column header for time (only needed for mobile acquisitions)
        hdr_scanpos_id str
            column header for scan pos id (only needed for static acquisitions -> equivalent to hdr_time in mobile acquisitions
        hdr_x: str [default 'x']
            column header for x coordinates
        hdr_y: str [default 'y']
            column header for y coordinates
        hdr_z: str [default 'z']
            column header for z coordinates
        sens_pos_id_offset: int
            Very specific use case where Scan Pos ID in position file does not correspond with Scan Pos ID in LAZ files and we need to add an offset

        """
        # resolve path
        if not os.path.isabs(path2file):
            self.logger.info(f"Provided path to sensor position file {path2file} is not absolute. Resolving relative to root folder {self.root_folder}.")
            if self.root_folder is not None:
                path2file = os.path.join(self.root_folder, path2file)
            else:
                path2file = os.path.abspath(path2file)
        if not os.path.exists(path2file):
            raise ValueError(f"Provided path to sensor position file {path2file} does not exist. Please check the path and try again.")

        if self.is_mobile: # case of mobile acquisitions (MLS, ULS)
            self.logger.info("Defining sensor positions for mobile acquisition. Interpolating trajectory for each pulse based on provided trajectory file.")
            self.traj = read_trajectory_file(path2traj=path2file, delimiter=delimiter, hdr_time=hdr_time, hdr_x=hdr_x, hdr_y=hdr_y, hdr_z = hdr_z)
        else: # case of static acquisition (TLS)
            self.logger.info("Defining sensor positions for static acquisition.")
            self.senspos = read_sensorpos_file(path2senspos=path2file, delimiter=delimiter, hdr_scanpos_id=hdr_scanpos_id, hdr_x=hdr_x, hdr_y=hdr_y, hdr_z=hdr_z, sens_pos_id_offset=sens_pos_id_offset)
        self.sens_pos_initialized = True

    def define_sensor_pos_singlePos(self, scan_pos_id, x, y, z):    # TODO: ugly workaround for the case where a single laz file from a single TLS position should be run
        """
        defines the scanner position of a single TLS scan position. This is currently just a work-around where we have
        the case of a single laz file and a single position without a text file defining e.g. multiple scan positions.
        Writes scanner position into self.senspos

        Parameters
        ----------
        scan_pos_id: int
            Scan Position Identificaiton number
        x: float
            X-Coordinates of scanner position
        y: float
            Y-Coordinates of scanner position
        z: float
            Z-Coordinates of scanner position

        """
        d = {'ScanPos': scan_pos_id,
             'sensor_x': x, 'sensor_y': y, 'sensor_z': z}

        senspos = pd.DataFrame(data=d, index=[0])

        self.senspos = senspos

        self.sens_pos_initialized = True

    def do_raytracing(self):
        """
        Perform ray tracing over all configured input.

        Handles three input shapes uniformly:

        * a directory of TLS LAS/LAZ files, each linked to its own scan position,
        * a single TLS LAS/LAZ file with one fixed scan position
        * a single mobile (MLS/ULS) LAS/LAZ file with a trajectory

        The return mode (single vs. multi) is resolved once, up front, for the whole dataset.
        Multi-return data is traced on the fly when it is sorted by GPS time;
        if unsorted data is detected -- at any point, not only in the first chunk -- processing
        switches to deferred mode and the remaining dataset is traced in single pass at the end.
        """
        if not self.sens_pos_initialized:
            raise ValueError(
                "Sensor positions not defined. Please call define_sensor_pos or "
                "define_sensor_pos_singlePos before running ray tracing."
            )

        started = time.perf_counter()

        jobs = self._build_scan_jobs()
        self.check_multi_return_handling(jobs)

        self.logger.info(f"Ray tracing {len(jobs)} scan(s).")
        mode = TraceMode.SINGLE_RETURN if self.single_return else TraceMode.STREAMING

        with self._timed("Reading and ray tracing all input"):
            for i, job in enumerate(jobs, start=1):
                mode = self._process_scan(job, mode, position=i, total_scans=len(jobs))

        self._finalize(mode)
        self.trace_mode = mode

        elapsed = time.perf_counter() - started
        self.get_raytracing_report(elapsed_seconds=elapsed)
        self.save_raytracing_output()
        self._write_run_summary(jobs, elapsed_seconds=elapsed)

    # ------------------------------------------------------------------
    # input assembly
    # ------------------------------------------------------------------

    def _build_scan_jobs(self):
        """
        Turn the three supported input shapes into one list of scan jobs.
        """
        if os.path.isdir(self.laz_in):
            if self.is_mobile:
                raise NotImplementedError(
                    "Mobile acquisition with multiple LAZ files is not yet implemented. "
                    "Provide a single LAZ file for mobile acquisitions or set "
                    "is_mobile to False for TLS. If you have multiple input LAZ files for "
                    "a mobile acquisition (e.g. multiple flight campaigns, multiple scan pattern), "
                    "please run OccPy once per acquisition."
                )
            if self.scans_linked is None:
                self.link_positions_to_laz_files()

            return [
                ScanJob(
                    laz_file=scan["laz_file"],
                    name=scan["scan_name"],
                    sensor_xyz=(scan["sensor_x"], scan["sensor_y"], scan["sensor_z"]),
                )
                for scan in self.scans_linked
            ]

        name = os.path.basename(self.laz_in)

        if self.is_mobile:
            return [ScanJob(laz_file=self.laz_in, name=name, sensor_xyz=None)]

        if len(self.senspos) > 1:
            self.logger.warning(
                f"{len(self.senspos)} sensor positions are defined but laz_in is a single file. "
                f"Using the first position (ScanPos={self.senspos['ScanPos'].values[0]}) for all pulses."
            )
        row = self.senspos.iloc[0]
        return [
            ScanJob(
                laz_file=self.laz_in,
                name=name,
                sensor_xyz=(
                    float(row["sensor_x"]),
                    float(row["sensor_y"]),
                    float(row["sensor_z"]),
                ),
            )
        ]

    # ------------------------------------------------------------------
    # per-scan processing
    # ------------------------------------------------------------------
    def _process_scan(self, job, mode, position=None, total_scans=None):
        """
        Process one LAS/LAZ file. Returns the (possibly downgraded) trace mode, so
        that a switch to DEFERRED persists for the rest of the run.

        Per-chunk ray-tracing times are accumulated and reported once, after the
        progress bar has closed, rather than logged per chunk. The per-chunk
        numbers are still available at DEBUG level.

        Parameters
        ----------
        job
        mode
        position
        total_scans
        """
        self.logger.info(f"===== Processing {job.name} =====")

        prev_gps_max = None
        traced_seconds = 0.0
        traced_batches = 0

        with self._timed(f"Reading and tracing {job.name}"):
            for chunk in self._iter_chunks(job, position=position, total_scans=total_scans):

                if mode is TraceMode.STREAMING:
                    is_sorted, reason = chunk_sortedness(chunk, prev_gps_max)
                    if not is_sorted:
                        self.logger.warning(
                            f"{job.name} is not sorted by gps_time ({reason}). Switching to "
                            f"deferred mode: the remaining data has to be held in memory "
                            f"until the end of the run, which is considerably slower. "
                            f"Consider sorting first, e.g. LAStools: "
                            f"lassort -i laz_in -gps_time -return_number -odix _sort -olaz -v"
                            f" or PDAL: "
                            f'pdal translate in.laz out.laz sort --filters.sort.dimensions=""GpsTime"" --filters.sort.order=""ASC"" --filters.sort.algorithm=""STABLE"" --writers.las.forward=all --writers.las.extra_dims=all'
                        )
                        mode = TraceMode.DEFERRED
                    prev_gps_max = float(chunk.gps_time[-1])

                sensor_x, sensor_y, sensor_z = self._sensor_positions(chunk, job)

                if mode is TraceMode.SINGLE_RETURN:
                    tic = time.perf_counter()
                    self.RayTr.doRaytracing_singleReturnPulses(
                        chunk.x, chunk.y, chunk.z,
                        sensor_x, sensor_y, sensor_z,
                        chunk.gps_time,
                    )
                    elapsed = time.perf_counter() - tic
                    self.logger.debug(f"Ray tracing batch (single-return): {elapsed:.2f} seconds")
                    traced_seconds += elapsed
                    traced_batches += 1
                    continue

                self.RayTr.addPointData(
                    chunk.x, chunk.y, chunk.z,
                    sensor_x, sensor_y, sensor_z,
                    chunk.gps_time, chunk.return_number, chunk.number_of_returns,
                )

                if mode is TraceMode.STREAMING:
                    # Sorted data: trace what we have, then release the traced
                    # pulses. Incomplete pulses straddling the chunk boundary are
                    # retained by clearPulseDataset and complete by the next chunk.
                    traced_seconds += self._trace_pulse_dataset(
                        "stored pulses", clear_after=True, level=logging.DEBUG
                    )
                    traced_batches += 1

        if traced_batches:
            self.logger.info(
                f"{job.name}: ray traced {traced_batches} chunk(s) in "
                f"{traced_seconds:.2f} seconds."
            )

        return mode

    def _iter_chunks(self, job, position=None, total_scans=None):
        """
        Yield normalised chunks from a LAS/LAZ file, with a progress bar.
        Parameters
        ----------
        job
        """
        require_gps_time = not self.single_return

        desc = job.name
        if position is not None and total_scans is not None and total_scans > 1:
            desc = f"[{position}/{total_scans}] {job.name}"

        with laspy.open(job.laz_file) as file:
            with tqdm(
                    total=file.header.point_count,
                    desc=desc,
                    unit="pts",
                    unit_scale=True,
                    disable=False,
            ) as pbar:
                for points in file.chunk_iterator(self.points_per_iter):
                    n = len(points)
                    if n == 0:
                        continue
                    yield normalise_chunk(points, job.name, require_gps_time)
                    pbar.update(n)

    def _sensor_positions(self, chunk, job):
        """
        Per-echo sensor positions as contiguous float64 arrays.

        Parameters
        ----------
        chunk
        job
        """
        if job.is_mobile:
            pos = interpolate_traj(
                self.traj["time"],
                self.traj["sensor_x"],
                self.traj["sensor_y"],
                self.traj["sensor_z"],
                chunk.gps_time
            )
            return (
                np.ascontiguousarray(pos["sensor_x"], dtype=np.float64),
                np.ascontiguousarray(pos["sensor_y"], dtype=np.float64),
                np.ascontiguousarray(pos["sensor_z"], dtype=np.float64),
            )

        return constant_sensor_positions(len(chunk), job.sensor_xyz)

    # ------------------------------------------------------------------
    # tracing and finalisation
    # ------------------------------------------------------------------

    def _trace_pulse_dataset(self, label, clear_after, level=logging.INFO):
        """
        The single place where the accumulated pulse dataset is traced.

        Centralising this guarantees that move_senspos_to_collinearity, the bug reports and the timing
        are applied consistently.

        ``level`` controls how loudly the timing and the collinearity notice are
        reported. The streaming path calls this once per chunk and passes DEBUG,
        so a live progress bar is not interrupted; the caller aggregates the
        returned times into a single INFO line per file.

        Returns the ray-tracing time in seconds.

        Parameters
        ----------
        label
        clear_after
        """
        self._pulse_report(f"before ray tracing {label}")
        if self.move_senspos_to_collinearity:
            self.logger.log(level,f"Moving sensor positions to force collinearity ({label}).")
            self.RayTr.moveSensorPos2Collinearity()

        tic = time.perf_counter()
        self.RayTr.doRaytracing()
        elapsed = time.perf_counter() - tic
        self.logger.log(level, f"Ray tracing {label}: {elapsed:.2f} seconds")

        if clear_after:
            self.RayTr.clearPulseDataset()

        self._pulse_report(f"after ray tracing {label}")

        return elapsed

    def _finalize(self, mode):
        """
        Run what is left in dataset.
        If deffered mode, this is full pulse dataset.
        Otherwise, if set, incomplete pulses are cleaned up and traced.

        Parameters
        ----------
        mode
        """
        if mode is TraceMode.SINGLE_RETURN:
            # Nothing was ever added to the pulse dataset, so there is nothing to flush and nothing to clean up.
            if self.cleanup_incomplete_pulses:
                self.logger.info(
                    "cleanup_incomplete_pulses has no effect for single-return data; skipping."
                )
            return

        if mode is TraceMode.DEFERRED:
            self.logger.info("Tracing the full accumulated pulse dataset.")
            self._trace_pulse_dataset("deferred pulses", clear_after=True)

        if self.cleanup_incomplete_pulses:
            # Incomplete pulses arise when data has been filtered, actively or
            # inside the vendor's processing. Completing them artificially is deliberately opt-in.
            self.logger.info(
                "Converting incomplete pulses into complete ones -- be cautious with this."
            )
            self._pulse_report("before cleaning up incomplete pulses")
            self.RayTr.cleanUpPulseDataset()
            self._pulse_report("after cleaning up incomplete pulses")

            self._trace_pulse_dataset("cleaned-up incomplete pulses", clear_after=True)

    # ------------------------------------------------------------------
    # small helpers
    # ------------------------------------------------------------------

    def _write_run_summary(self, jobs, elapsed_seconds):
        """
        Record the facts that were resolved at run time, not at construction.

        ``config.json`` says what was requested. This says what actually
        happened: which return mode was used and where that answer came from,
        whether the return fields were degenerate, and whether the run had to
        fall back from streaming to deferred tracing. None of that is
        recoverable from the config or the output grids, and all of it changes
        how the grids should be interpreted.

        Parameters
        ----------
        jobs:
        elapsed_seconds: int
        """
        summary = {
            "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "elapsed_seconds": round(elapsed_seconds, 2),
            "scans_processed": len(jobs),
            "single_return": self.single_return,
            "single_return_source": self.single_return_source,
            "returns_all_zero": self.returns_all_zero,
            "trace_mode": self.trace_mode.value if self.trace_mode else None,
            "cleanup_incomplete_pulses": self.cleanup_incomplete_pulses,
            "move_senspos_to_collinearity": self.move_senspos_to_collinearity,
        }

        if self.return_mode_result is not None:
            summary["return_mode"] = self.return_mode_result.mode.value
            summary["return_mode_source"] = self.return_mode_result.source
            summary["return_mode_detail"] = self.return_mode_result.detail
            summary["return_mode_probed_files"] = self.return_mode_probed_files

        if self.trace_mode is TraceMode.DEFERRED and not self.single_return:
            summary["note"] = (
                "Unsorted gps_time was detected, so tracing fell back to deferred mode. "
                "Results are unaffected, but the run held the pulse dataset in memory; "
                "sorting the input by gps_time would be considerably faster."
            )

        # Traversal counters, read once by get_raytracing_report so that the
        # console report and this file cannot disagree.
        stats = self.traversal_stats if self.traversal_stats else self._traversal_stats()
        if stats:
            summary["traversal"] = stats

        path = os.path.join(self.out_dir, "run_summary.json")
        try:
            with open(path, "w") as to:
                json.dump(summary, to, indent=2)
            self.logger.info(f"Wrote run summary to {path}")
        except OSError as exc:
            self.logger.warning(f"Could not write the run summary to {path}: {exc}")

        return summary

    @contextmanager
    def _timed(self, label, level=logging.INFO):
        tic = time.perf_counter()
        try:
            yield
        finally:
            self.logger.log(level, f"{label}: {time.perf_counter() - tic:.2f} seconds")

    def _pulse_report(self, when):
        if self.debug:
            self.logger.debug(f"Pulse dataset report {when}")
            self.RayTr.getPulseDatasetReport()

    def _traversal_stats(self):
        """
        Read the traversal counters from the C++ side.

        Each getter is wrapped individually so that a signature change in raytr
        costs one field rather than the whole report or the run summary. The
        result is cached on ``self.traversal_stats`` so the report and the run
        summary cannot disagree.
        """
        getters = {
            "total_pulses_in_dataset": self.getTotalNumPulses,
            "traversed_pulses": self.getNumTraversedPulses,
            "registered_hits": self.getNumRegisteredHits,
            "echoes_outside_grid": self.getNumEchoesOutside,
            "missing_returns": self.getNumMissingReturns,
            "pulses_not_intersecting_grid": self.getNumNonIntersectPulses,
        }
        stats = {}
        for name, getter in getters.items():
            try:
                stats[name] = int(getter())
            except Exception as exc:  # noqa: BLE001 - a report must not break a finished run
                self.logger.debug(f"Could not read {name}: {exc}")

        self.traversal_stats = stats
        return stats

    def get_raytracing_report(self, elapsed_seconds=None):
        """
        Print a report on the voxel traversal.

        The numbers are read from the C++ counters on the Python side and
        written out here, rather than relying on ``reportOnTraversal`` printing
        from C++. That call writes to the C++ standard output stream, which is
        buffered separately from Python's: in a terminal it usually still
        appears, but it goes missing in notebooks and under output capture, and
        it ignores the logger configuration entirely.

        The report is written with ``tqdm.write`` rather than through the
        logger, because it is a result rather than a log message: it should
        appear even when the run is not verbose, and it must not collide with a
        progress bar.

        Parameters
        ----------
        elapsed_seconds : float, optional
            Wall time of the ray-tracing run, included in the report if given.

        Returns
        -------
        dict
            The traversal counters, also stored on ``self.traversal_stats``.
        """
        stats = self._traversal_stats()

        if self.debug:
            # The C++ report, for cross-checking against the Python one.
            self.RayTr.reportOnTraversal()

        if not stats:
            self.logger.warning(
                "Could not read any traversal counters, so no report can be produced. "
                "The occlusion grids themselves are unaffected."
            )
            return stats

        def fmt(key):
            return f"{stats[key]:,}" if key in stats else "n/a"

        def pct(key, of_key):
            if key not in stats or of_key not in stats or not stats[of_key]:
                return ""
            return f"  ({100.0 * stats[key] / stats[of_key]:.1f}%)"

        total = "total_pulses_in_dataset"
        lines = [
            "",
            "==================== Occlusion mapping report ====================",
            f"  Pulses traversed          {fmt('traversed_pulses'):>15}"
            f" of {fmt(total)}{pct('traversed_pulses', total)}",
            f"  Pulses missing the grid   {fmt('pulses_not_intersecting_grid'):>15}"
            f"{'':<{len(fmt(total)) + 4}}{pct('pulses_not_intersecting_grid', total)}",
            f"  Returns registered        {fmt('registered_hits'):>15}",
            f"  Returns outside the grid  {fmt('echoes_outside_grid'):>15}",
            f"  Returns missed            {fmt('missing_returns'):>15}",
        ]
        if elapsed_seconds is not None:
            lines.append(f"  Ray tracing time          {elapsed_seconds:>15.2f} s")
        lines.append("=" * 66)

        tqdm.write("\n".join(lines))

        if stats.get("missing_returns", 0) > 0:
            self.logger.warning(
                f"{stats['missing_returns']:,} returns were missed during traversal: they did "
                f"not fall in a voxel crossed by their own pulse, which means the returns are "
                f"not collinear with the sensor position. Setting "
                f"move_senspos_to_collinearity=True can recover them, with the caveats noted "
                f"in the config documentation."
            )

        return stats

    def save_raytracing_output(self):
        """
        Extract and save the outputs of the ray-tracing process.
        This method performs the following steps:
        1. Extracts the voxel-wise hit (`Nhit`), miss (`Nmiss`), and occlusion (`Nocc`) voxelgrids and saves as .npy in self.out_dir
        2. Creates a voxel classification grid based on the `Nhit`, `Nmiss`, and `Nocc` values:
        - 1 = observed (hit > 0)
        - 2 = empty (miss > 0, hit == 0)
        - 3 = occluded (occlusion > 0, hit == 0, miss == 0)
        - 4 = unobserved (all three == 0)
        3. Writes `.ply` files for all voxel outputs if `self.output_voxels` is True (takes long and creates large files)

        """
        self.logger.info("Extracting Nhit")
        tic = time.time()
        self.Nhit = self.RayTr.getNhit()
        self.Nhit = np.array(self.Nhit, dtype=np.int32)

        toc = time.time()
        self.logger.info("Elapsed Time: {:.2f} seconds".format(toc - tic))

        self.logger.info("Extracting Nocc")
        tic = time.time()
        self.Nocc = self.RayTr.getNocc()
        self.Nocc = np.array(self.Nocc, dtype=np.int32)

        toc = time.time()
        self.logger.info("Elapsed Time: {:.2f} seconds".format(toc - tic))

        self.logger.info("Extracting Nmiss")
        tic = time.time()
        self.Nmiss = self.RayTr.getNmiss()
        self.Nmiss = np.array(self.Nmiss, dtype=np.int32)

        toc = time.time()
        self.logger.info("Elapsed Time: {:.2f} seconds".format(toc - tic))

        self.logger.info("Saving Occlusion Outputs As .npy")
        tic = time.time()
        np.save(os.path.join(self.out_dir, "Nhit.npy"), self.Nhit)
        np.save(os.path.join(self.out_dir, "Nmiss.npy"), self.Nmiss)
        np.save(os.path.join(self.out_dir, "Nocc.npy"), self.Nocc)
        toc = time.time()
        self.logger.info("Elapsed Time: {:.2f} seconds".format(toc - tic))

        # Create Classification grid
        self.logger.info("Classify Grid")
        tic = time.time()
        self.Classification = np.zeros((self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz']), dtype=np.uint8)

        self.Classification[np.logical_and.reduce((self.Nhit > 0, self.Nmiss >= 0, self.Nocc >= 0))] = 1  # voxels that were observed
        self.Classification[np.logical_and.reduce((self.Nhit == 0, self.Nmiss > 0, self.Nocc >= 0))] = 2  # voxels that are empty
        self.Classification[
            np.logical_and.reduce((self.Nhit == 0, self.Nmiss == 0, self.Nocc > 0))] = 3  # voxels that are hidden (occluded)
        self.Classification[np.logical_and.reduce((self.Nhit == 0, self.Nmiss == 0,
                                              self.Nocc == 0))] = 4  # voxels that were not observed # TODO: Figure out, why this overwrites voxels that are classified as occluded! -> this was because np.logical_and only takes in 2 arrays as input, not 3! use np.logical_and.reduce() for that!

        np.save(os.path.join(self.out_dir, "Classification.npy"), self.Classification)
        toc = time.time()
        self.logger.info("Elapsed Time: " + str(toc - tic) + " seconds")

        # write ply file
        if self.output_voxels:
            self.logger.info("Saving Occlusion Outputs As .ply")
            tic = time.time()
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nhit)
            ost.write_ply(os.path.join(self.out_dir, "Nhit.ply"), verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nmiss)
            ost.write_ply(os.path.join(self.out_dir, "Nmiss.ply"), verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nocc)
            ost.write_ply(os.path.join(self.out_dir, "Nocc.ply"), verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Classification)
            ost.write_ply(os.path.join(self.out_dir, "Classification.ply"), verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            self.occl = np.zeros(shape=self.Classification.shape)
            x4, y4, z4 = np.where(self.Classification == 4)
            self.occl[x4, y4, z4] = self.Classification[x4, y4, z4]
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.occl)
            ost.write_ply(os.path.join(self.out_dir, "Occl.ply"), verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            toc = time.time()
            self.logger.info("Elapsed Time: " + str(toc - tic) + " seconds")


    def link_positions_to_laz_files(self):
        """
        Link TLS LAS/LAZ files from a directory input to scanner positions before processing.

        This method is only applicable for TLS runs where ``self.laz_in`` is a directory.
        It validates that each LAS/LAZ file can be linked to exactly one sensor position and
        stores the links for re-use in ``do_raytracing``.

        Returns
        -------
        pandas.DataFrame
            Table with columns ["laz_file", "scan_name", "scan_id", "sensor_x", "sensor_y", "sensor_z"].
        """

        if self.is_mobile:
            raise ValueError("link_positions_to_laz_files is only valid for TLS (is_mobile=False).")

        if not os.path.isdir(self.laz_in):
            raise ValueError("link_positions_to_laz_files requires laz_in to be a directory containing TLS LAS/LAZ files.")

        if not self.sens_pos_initialized:
            raise ValueError("Sensor positions not defined. Please call define_sensor_pos first.")

        laz_files = sorted(glob.glob(os.path.join(self.laz_in, "*.laz")))
        las_files = sorted(glob.glob(os.path.join(self.laz_in, "*.las")))
        laz_files.extend(las_files)
        if len(laz_files) == 0:
            raise ValueError(f"No .las or .laz files found in input directory: {self.laz_in}")

        links = []
        self.logger.debug(f"Sensor position file: {self.senspos}")
        for laz_file in laz_files:
            scan_name = os.path.basename(laz_file)
            scan_id = os.path.splitext(scan_name)[0]

            if self.str_idxs_ScanPosID is not None:
                scan_id = scan_name[self.str_idxs_ScanPosID[0]:self.str_idxs_ScanPosID[1]]
                self.logger.debug(f"Using str_idxs_ScanPosID {self.str_idxs_ScanPosID} to extract scan ID from file name, resulting in scan ID {scan_id}.")

            self.logger.debug(f"Using scan ID {scan_id} for file {scan_name}.")
            matches = self.senspos.loc[self.senspos['ScanPos'] == scan_id]

            if matches.empty:
                self.logger.warning(f"No sensor position found for scan ID '{scan_id}' (file '{scan_name}'). This scan will be skipped.")
                continue

            if len(matches) > 1:
                raise ValueError(
                    f"Multiple sensor positions found for scan ID '{scan_id}' (file '{scan_name}'). "
                    f"ScanPos IDs must be unique."
                )
            
            self.logger.debug(f"Found matching sensor position for scan ID {scan_id}: (x={matches['sensor_x'].values[0]:.3f}, y={matches['sensor_y'].values[0]:.3f}, z={matches['sensor_z'].values[0]:.3f}).")

            links.append({
                'laz_file': laz_file,
                'scan_name': scan_name,
                'scan_id': scan_id,
                'sensor_x': matches['sensor_x'].values[0],
                'sensor_y': matches['sensor_y'].values[0],
                'sensor_z': matches['sensor_z'].values[0],
            })

        if len(links) == 0:
            raise ValueError(
                f"No sensor position found for any of the {len(laz_files)} LAS/LAZ files in "
                f"{self.laz_in}. Please check str_idxs_ScanPosID and the scan position file."
            )

        self.scans_linked = links

        self.logger.info(f"Linked {len(links)} TLS LAS/LAZ files to scan positions.")
        return pd.DataFrame(links)

    def check_multi_return_handling(self, jobs=None):
        """
       Resolve the return mode for the whole dataset, once, before any tracing.

        The header field ``number_of_points_by_return`` answers this for free
        when it is trustworthy; otherwise the points are scanned, keeping only
        the return-number fields, with an early exit as soon as a second return
        appears. See ``occpy.pulse_util.detect_return_mode``.

        For a directory of LAS/LAZ files only the first file is probed. Multi-station
        TLS acquisitions write one file per scan position with identical sensor
        settings, so the return mode is a property of the acquisition rather than
        of the individual file. Set ``check_returns_all_files`` in the config to
        probe every file instead.

        If ``single_return`` was set explicitly in the config, detection
        validates it rather than overriding it: a conflict that would corrupt the
        results raises, one that only costs performance warns.
        Parameters
        ----------
        jobs: list of ScanJob, optional
            Scan jobs to probe. Built from the config when omitted.

        Returns
        ----------
        bool
            the resolved value of ``self.single_return``.
        """

        if jobs is None:
            jobs = self._build_scan_jobs()
        if len(jobs) == 0:
            raise ValueError("No input files to check for return mode.")

        probes = jobs if self.check_returns_all_files else jobs[:1]
        if len(jobs) > 1 and not self.check_returns_all_files:
            self.logger.info(
                f"Probing the return mode on {jobs[0].name} only and applying it to all "
                f"{len(jobs)} files. Set check_returns_all_files=True to probe every file."
            )

        results = []
        for job in probes:
            result = detect_return_mode(job.laz_file, self.points_per_iter, progress=self.verbose)
            self.logger.info(
                f"{job.name}: return mode {result.mode.value} (from {result.source}) -- {result.detail}"
            )
            results.append((job, result))

        # Kept for the run summary: which file was probed and whether the answer
        # came from the header or from scanning the points
        self.return_mode_result = results[0][1]
        self.return_mode_probed_files = [job.name for job, _ in results]

        conflicting = {r.mode for _, r in results}
        if len(conflicting) > 1:
            summary = ", ".join(f"{j.name}={r.mode.value}" for j, r in results)
            raise ValueError(
                f"Input files disagree on the return mode ({summary}). Process them in "
                f"separate runs, or fix the acquisition metadata."
            )

        mode = results[0][1].mode
        detected_single = mode.implies_single_return
        self.returns_all_zero = mode is ReturnMode.DEGENERATE_ZERO

        if self.returns_all_zero:
            self.logger.warning(
                "return_number and number_of_returns are 0 for every echo. The return "
                "fields carry no information, so the data is processed as single-return."
            )
            if self.single_return is False:
                raise ValueError(
                    "single_return is set to False, but the return-number fields are 0 "
                    "everywhere, so echoes cannot be grouped into pulses. Set "
                    "single_return to True, or repair the return fields."
                )

        if self.single_return is None:
            self.single_return = detected_single
            self.single_return_source = "detected"
            self.logger.info(f"Resolved single_return={self.single_return} from the data.")
            return self.single_return

        self.single_return_source = "config"

        if self.single_return and mode is ReturnMode.MULTI:
            raise ValueError(
                f"single_return is set to True, but the data contains multiple returns per "
                f"pulse ({results[0][1].detail}). This would silently corrupt the occlusion "
                f"result: every echo would be traced as an independent pulse."
            )

        if not self.single_return and detected_single:
            self.logger.warning(
                "single_return is set to False, but the data contains only single returns. "
                "The result is unaffected, but setting single_return=True avoids converting "
                "the point cloud into a pulse dataset and is considerably faster."
            )

        return self.single_return

    def get_chm(self):
        """
        returns chm.
        TODO: check if this is still needed!
        Returns
        -------
        reference to canopy height model

        """
        if self.chm is None:
            self.logger.warning("No CHM was defined. To define CHM ")
        return self.chm

    def clean_up_RayTr(self):
        """
        Free up memory after raytracing.
        """
        del self.RayTr

    def getGridDimensions(self):
        """
        Get the grid dimensions
        Returns
        -------
        gridDim: [int, int, int, int]
            vector with grid diemsnions [minx, maxx, miny, maxy]

        """
        gridDim = self.RayTr.getGridDimensions()
        gridDim = np.asarray(gridDim, dtype=np.int32)

        return gridDim

    def getGridOrigin(self):
        """
        Get the grid origin
        Returns
        -------
        origin: [float, float, float]
            with grid origins [minx, miny, minz]

        """
        origin = np.asarray(self.RayTr.getGridOrigin(), dtype=np.float64)
        return origin

    def getNumTraversedPulses(self):
        """
        Get the number of traversed pulses
        Returns
        -------
        num_trav_pulses: int
            total number of pulses traversed by the algorithm

        """
        return np.int32(self.RayTr.get_num_traversed_pulses())

    def getTotalNumPulses(self):
        """
        Get the total number of pulses loaded into the pulse dataset on the c++ side
        Returns
        -------
        total_num_pulses: int
            total number of pulses stored in the pulse dataset

        """
        return np.int32(self.RayTr.get_total_pulses_in_dataset())

    def getNumRegisteredHits(self):
        """
        Get the number of registered hits
        Returns
        -------
        num_registered_hits: int
            Number of registered hits during voxel traversal

        """
        return np.int32(self.RayTr.get_num_registered_hits())

    def getNumEchoesOutside(self):
        """
        Get the number of echoes registered outside the voxel grid
        Returns
        -------
        num_echoes_outside: int
            number of echoes outside the voxel grid

        """
        return np.int32(self.RayTr.get_num_echoes_outside())

    def getNumMissingReturns(self):
        """
        Get the number of missing returns. Missing returns can occur, if the laser returns do not fall into an exact line
        that is defined by the pulse origin and the last return. If the laser return falls outside of the voxel which is
        traversed by the pulse, this return is counted to num_missing_returns.
        Returns
        -------
        num_missing_returns: int
            number of returns that were not registered during voxel traversal.

        """
        return np.int32(self.RayTr.get_num_missing_returns())

    def getNumNonIntersectPulses(self):
        """
        Get the number of pulses that do not intersect the defined voxel grid
        Returns
        -------
        num_pulses_not_intersect: int
            number of pulses that do not intersect the defined voxel grid.

        """
        return np.int32(self.RayTr.get_num_pulses_no_intersection())

    @staticmethod
    def align_plot_dim_to_voxel_size(plot_dim, vox_dim):
        """extend max bounds so each axis extent is divisible by vox_dim."""

        adjusted = [float(v) for v in plot_dim]
        messages = []
        tol = 1e-9
        axes = (("X", 0, 3), ("Y", 1, 4), ("Z", 2, 5))

        for axis_name, min_idx, max_idx in axes:
            min_bound = adjusted[min_idx]
            max_bound = adjusted[max_idx]
            extent = max_bound - min_bound

            if extent <= 0:
                raise ValueError(
                    f"Invalid plot_dim on axis {axis_name}: max ({max_bound}) must be greater than min ({min_bound})."
                )

            n_voxels = int(np.ceil((extent / vox_dim) - tol))
            adjusted_extent = n_voxels * vox_dim

            if not np.isclose(extent, adjusted_extent, rtol=0.0, atol=tol):
                new_max = min_bound + adjusted_extent
                messages.append(
                    f"Axis {axis_name}: extent {extent:.12g} is not divisible by vox_dim {vox_dim:.12g}. "
                    f"Extending max bound from {max_bound:.12g} to {new_max:.12g}."
                )
                adjusted[max_idx] = new_max

        return adjusted, messages

