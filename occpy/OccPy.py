import logging
import time
import os
import glob
import json

import numpy as np
import pandas as pd
import laspy
import OSToolBox as ost
from tqdm import tqdm

from raytr import PyRaytracer
from occpy.util import prepare_ply, read_trajectory_file, read_sensorpos_file, interpolate_traj, last_nonzero


class OccPy:

    def __init__(self, config=None, config_file=None):
        """
        initialize OccPy object

        Parameters in config file:  
        Must include:  
            - 'laz_in' : path to single .laz file or directory with multiple .laz files
            - 'vox_dim' : voxel size in meters  
            - 'plot_dim': grid for occlusion mapping: [minX, minY, minZ, maxX, maxY, maxZ]  
        Optional parameters:  
            - 'out_dir' : output directory (default: ./output)
            - 'output_voxels' : whether to export `.ply` voxel grids (default: False)
            - 'verbose': set logging level  (default: False)
            - 'debug': set logging level (default: False)
            - 'lower_threshold': lower threshold above ground to exclude from occlusion mapping in voxels (default: 0)
            - 'points_per_iter': number of points read in from laz file in each iteration (default: 10000000)
            - 'delimiter': csv delimiter for scan position file (default: ",")
            - 'root_folder': if given, will assume other paths are relative to this root folder and will prepend it to the paths (default: None)
            - 'is_mobile': whether the acquisition is mobile (MLS/ULS) or static (TLS) (default: False)
            - 'single_return': whether the data is single return or multi return data (default: False)
            - 'str_idxs_ScanPosID': string indices of where the scan position identifier is written in the laz file name. If not given, will use file name as ID (without extension) (default: None)
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
        
        optional_args = ["out_dir", "output_voxels", "verbose", "debug", "lower_threshold", "points_per_iter", "delimiter", "root_folder", "single_return", "str_idxs_ScanPosID"]
        
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
        self.single_return = config.get("single_return", False)
        self.str_idxs_ScanPosID = config.get("str_idxs_ScanPosID", None)

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
        console_handler = logging.StreamHandler()
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

        Returns
        -------

        """
        d = {'ScanPos': scan_pos_id,
             'sensor_x': x, 'sensor_y': y, 'sensor_z': z}

        senspos = pd.DataFrame(data=d, index=[0])

        self.senspos = senspos

        self.sens_pos_initialized = True


    def do_raytracing(self):
        """
        Perform ray tracing.

        This method processes either a directory of LAZ files (for TLS with known scan positions) or a single LAZ file
        (Single TLS position or MLS/ULS with a trajectory, depends on self.is_mobile).
        In the case of TLS, for each LAZ file, it extracts point positions and sensor positions, and
        then performs ray tracing, accounting for single or multi-return pulse information. Multi-return handling
        supports on-the-fly processing if the data is sorted by GPS time; otherwise, the full dataset must be loaded first.
        """

        if not self.sens_pos_initialized:
            raise ValueError("Sensor positions not defined. Please call define_sensor_pos or define_sensor_pos_singlePos before running ray tracing.")
        
        is_sorted = lambda a: np.all(a[:-1] <= a[1:])

        # TODO: link positions to laz files first (in case of TLS at least), and warn/error if one/more positions cant be linked before processing

        run_raytracing_after_loading = False
        if os.path.isdir(self.laz_in): # multi-pos TLS 

            if self.is_mobile: # TODO: could it occur that we have mobile acquisition with multiple laz files corresponding to one traj file? if so, should implement
                raise NotImplementedError("The case of mobile acquisition with multiple laz files is currently not implemented. Please provide a single laz file for mobile acquisitions or set is_mobile to False for TLS.")
            
            if self.scans_linked is None:
                self.link_positions_to_laz_files()

            for scan in self.scans_linked:
                f = scan['laz_file']
                scan_name = scan['scan_name']

                print(f"###############################")
                print(f"##### Processing {scan_name}...")
                print(f"###############################")

                scanpos_X = scan['sensor_x']
                scanpos_Y = scan['sensor_y']
                scanpos_Z = scan['sensor_z']

                # read in laz file
                tic = time.time()
                with laspy.open(f) as file:

                    count = 0
                    for points in file.chunk_iterator(self.points_per_iter):

                        self.check_multi_return_handling(points, scan_name)

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
                            self.logger.info("Do raytracing with all pulses in batch")
                            tic_r = time.time()
                            self.RayTr.doRaytracing_singleReturnPulses(x, y, z, sensor_x, sensor_y,
                                                                  sensor_z, gps_time)
                            toc_r = time.time()
                            self.logger.info("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))

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
                                    self.logger.warning(
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
                                # Get report on pulse datase
                                if self.debug:
                                    self.RayTr.getPulseDatasetReport()

                                # run raytracing on added points
                                self.logger.info("Do raytracing with stored pulses")
                                tic_r = time.time()
                                self.RayTr.doRaytracing()
                                toc_r = time.time()
                                self.logger.info("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))

                                self.RayTr.clearPulseDataset()

                                # Check if traversed pulses have been deleted from map
                                if self.debug:
                                    self.RayTr.getPulseDatasetReport()

                        count = count + len(gps_time)

            toc = time.time()
            if sorted and not self.single_return:
                # optional: incomplete pulses can occur if the data has been filtered (either actively or during black box processing
                # of the processing software. We could actively turn the incomplete pulses into complete ones and do the raytracing
                # for them!
                self.logger.info("convert incomplete pulses to complete ones - be cautious with that!")
                if self.debug:
                    self.RayTr.getPulseDatasetReport()
                self.RayTr.cleanUpPulseDataset()
                if self.debug:
                    self.RayTr.getPulseDatasetReport()
                self.logger.info("Run raytracing for incomplete pulses")
                tic_r = time.time()
                self.RayTr.doRaytracing()
                toc_r = time.time()
                self.logger.info("Time elapsed for raytracing incomplete pulses: {:.2f} seconds".format(toc_r - tic_r))
                self.logger.info("Time elapsed for reading and raytracing entire data: {:.2f} seconds".format(toc_r - tic))
            elif not sorted and not self.single_return:
                self.logger.info("Time elapsed for reading in data: {:.2f} seconds".format(toc - tic))

                if self.debug:
                    self.RayTr.getPulseDatasetReport()

                self.logger.info("Clean up pulse dataset in order to handle incomplete pulses")
                self.RayTr.cleanUpPulseDataset()

                if self.debug:
                    self.RayTr.getPulseDatasetReport()

                self.logger.info("Do actual raytracing with all pulses")
                tic = time.time()
                self.RayTr.doRaytracing()
                toc = time.time()
                self.logger.info("Time elapsed for raytracing: {:.2f} seconds".format(toc - tic))

        else: # if input is a single laz file, TLS with single scan position or MLS/ULS with trajectory
            with laspy.open(self.laz_in) as file:
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

                        self.check_multi_return_handling(points, self.laz_in)

                        # for the case of mobile acquisitions, inerpolate trajectory for gps_time
                        if self.is_mobile:
                            # call interpolate function for trajectory to extract sensor position for each gps_time
                            SensorPos = interpolate_traj(self.traj['time'], self.traj['sensor_x'], self.traj['sensor_y'],
                                                              self.traj['sensor_z'], gps_time)

                        else:
                            SensorPos = self.senspos
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
                                    self.logger.warning(
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
                                # Get report on pulse dataset
                                if self.debug:
                                    self.RayTr.getPulseDatasetReport()

                                # run raytracing on added points
                                self.logger.info("Do raytracing with stored pulses")
                                tic_r = time.time()
                                self.RayTr.doRaytracing()
                                toc_r = time.time()
                                self.logger.info("Time elapsed for raytracing batch: {:.2f} seconds".format(toc_r - tic_r))

                                self.RayTr.clearPulseDataset() # clear out data that have been traced.

                                # Check if traversed pulses have been deleted from map
                                
                                if self.debug:
                                    self.RayTr.getPulseDatasetReport()

                        count = count + len(gps_time)
                        pbar.update(len(points))

        if run_raytracing_after_loading:
            self.RayTr.doRaytracing()

        self.get_raytracing_report()
        self.save_raytracing_output()

    def get_raytracing_report(self):
        """
        prints a report on the voxel traversal to the console

        Returns
        -------

        """
        # Get report on traversal
        self.RayTr.reportOnTraversal()

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

        Returns
        -------

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
        self.Classification = np.zeros((self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz']), dtype=int)

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
        Link TLS LAZ files from a directory input to scanner positions before processing.

        This method is only applicable for TLS runs where ``self.laz_in`` is a directory.
        It validates that each LAZ file can be linked to exactly one sensor position and
        stores the links for re-use in ``do_raytracing``.

        Returns
        -------
        pandas.DataFrame
            Table with columns ["laz_file", "scan_name", "scan_id", "sensor_x", "sensor_y", "sensor_z"].
        """

        if self.is_mobile:
            raise ValueError("link_positions_to_laz_files is only valid for TLS (is_mobile=False).")

        if not os.path.isdir(self.laz_in):
            raise ValueError("link_positions_to_laz_files requires laz_in to be a directory containing TLS LAZ files.")

        if not self.sens_pos_initialized:
            raise ValueError("Sensor positions not defined. Please call define_sensor_pos first.")

        laz_files = sorted(glob.glob(os.path.join(self.laz_in, "*.laz")))
        if len(laz_files) == 0:
            raise ValueError(f"No LAZ files found in input directory: {self.laz_in}")

        links = []
        for laz_file in laz_files:
            scan_name = os.path.basename(laz_file)
            scan_id = os.path.splitext(scan_name)[0]

            if self.str_idxs_ScanPosID is not None:
                # TODO: we are forcing integer ids, this might not always be te case, but for now ok.
                scan_id = int(scan_name[self.str_idxs_ScanPosID[0]:self.str_idxs_ScanPosID[1]])

            self.logger.debug(f"Using scan ID {scan_id} for file {scan_name}.")
            matches = self.senspos.loc[self.senspos['ScanPos'] == scan_id]

            if matches.empty:
                raise ValueError(
                    f"No sensor position found for scan ID '{scan_id}' (file '{scan_name}'). "
                    f"Please check str_idxs_ScanPosID and the scan position file."
                )
            if len(matches) > 1:
                raise ValueError(
                    f"Multiple sensor positions found for scan ID '{scan_id}' (file '{scan_name}'). "
                    f"ScanPos IDs must be unique."
                )

            links.append({
                'laz_file': laz_file,
                'scan_name': scan_name,
                'scan_id': scan_id,
                'sensor_x': matches['sensor_x'].values[0],
                'sensor_y': matches['sensor_y'].values[0],
                'sensor_z': matches['sensor_z'].values[0],
            })

        self.scans_linked = links

        self.logger.info(f"Linked {len(links)} TLS LAZ files to scan positions.")
        return pd.DataFrame(links)

    def check_multi_return_handling(self, points, scan_name):
        """
        Check if the input laz data contains multiple returns per pulse and set self.single_return accordingly if not already set by the user.
        If self.single_return is already set by the user, this function will check if the data is consistent with the provided setting and will raise a warning if not.

        Parameters
        ----------
        points: laspy points object
            The points object read from a laz file chunk, containing point attributes such as x, y, z, gps_time, return_number, number_of_returns, etc.

        """

        if self.single_return:
            if np.any(points.return_number > 1):
                raise ValueError(f"Data appears to contain multiple returns per pulse (detected in {scan_name}), but single_return is set to True. This leads to unexpected behavior.")
        else:
            if not np.any(points.return_number > 1) and not np.any(points.number_of_returns > 1):
                self.logger.warning(f"Data appears to contain only single returns per pulse (detected in {scan_name}), but single_return is set to False. Consider setting single_return to True for more efficient processing of single return data.")

        return


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

        Returns
        -------

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

