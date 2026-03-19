import os
import glob
import logging
import cv2
import time
import json
from random import randrange

import pandas as pd
import numpy as np
import OSToolBox as ost

from occpy import riegl_io
from occpy.util import prepare_ply
from occpy.visualization import plot_riegl_grid
from raytr import PyRaytracer


class OccPyRIEGL:
    """
    Voxel-based occlusion mapping from RIEGL terrestrial laser scanning data.

    This class handles RIEGL-specific data workflows, including reading `.rdbx`, `.rxp`, and transformation files,
    optionally modeling empty pulses using preview images, performing voxel-based ray tracing, and saving the output
    grids. It provides preprocessing, collinearity checks, and output saving in `.npy` and `.ply` formats.
    """

    def __init__(self, config_file):
        """
        Initialize an OccPyRIEGL instance for RIEGL TLS occlusion mapping.

        Parameters in config file:  
        Must include:  
            - 'proj_folder' : path to RIEGL project directory with `.rdbx` and `.rxp` files  
            - 'riscan_folder' : path to RiSCAN PRO scan directory with transformation files  
            - 'vox_dim' : voxel size in meters  
            - 'plot_dim': grid for occlusion mapping: [minX, minY, minZ, maxX, maxY, maxZ]  
        Optional parameters:  
            - 'out_dir' : output directory (default: ./output)
            - 'buffer' : spatial buffer around point cloud   
            - 'output_voxels' : whether to export `.ply` voxel grids  
            - 'model_empty_pulses' : whether to model empty pulses  
            - 'verbose': set logging level  
            - 'debug': set logging level, run extra checks and output debug files
            - 'auto_dim': not implemented yet
            - 'buffer': not implemented yet
            - 'exclude_scan_pattern': string pattern, to exclude scans containing this pattern (in name of rdbx and transform files)

        Parameters
        ----------
        config_file : str
            Path to a JSON configuration file containing processing parameters.
        """

        with open(config_file) as f:
            config = json.load(f)

        necessary_args = ["proj_folder", "riscan_folder", "vox_dim", "plot_dim"]
        missing = []
        for key in necessary_args:
            if key not in config:
                missing.append(key)

        if len(missing) > 0:
            raise ValueError(f"Missing necessary arguments in config file: {missing}")

        optional_args = ["model_empty_pulses", "verbose", "debug", "output_voxels", "lower_threshold", "out_dir", "auto_dim", "buffer", "exclude_scan_pattern"]

        print(f"INFO: optional arguments: {optional_args}")

        self.riscan_folder = config["riscan_folder"]
        self.proj_folder = config["proj_folder"]
        self.vox_dim = config["vox_dim"]

        self.model_empty_pulses = config.get("model_empty_pulses", False)
        self.verbose = config.get("verbose", False)
        self.debug = config.get("debug", False)
        self.output_voxels = config.get("output_voxels", False)
        self.lower_threshold = config.get("lower_threshold", 0)
        self.exclude_scan_pattern = config.get("exclude_scan_pattern", None)
        self.out_dir = config.get("out_dir", os.path.join(os.getcwd(), "output"))

        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir, exist_ok=True)
        # copy config file for future reference
        with open(os.path.join(self.out_dir, "config.json"), "w") as to:
            json.dump(config, to)

        # -- config logging 
        if self.debug:
            logging_level = logging.DEBUG
        elif self.verbose:
            logging_level = logging.INFO
        else:
            logging_level = logging.WARNING
        self.logger = logging.getLogger('occpy_logger')
        self.logger.setLevel(logging_level)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging_level)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

        self.logger.info(f"config:")
        self.logger.info(f"{config}")

        # --- prepare input

        self.prepare_input()

        # --- init dimensions and raytracer

        if "auto_dim" in config and config["auto_dim"] is True:
            if "buffer" in config:
                buffer = config["buffer"]
            else:
                buffer = [0,0,0,0]
            # TODO: implement
            plot_dim = self.determine_grid(buffer)
        else:
            if "plot_dim" not in config:
                raise ValueError("plot_dim must be provided if auto_dim is not set to True")
            
            plot_dim = config["plot_dim"]
            self.plot_dim = dict(minX=plot_dim["minX"],
                                maxX=plot_dim["maxX"],
                                minY=plot_dim["minY"],
                                maxY=plot_dim["maxY"],
                                minZ=plot_dim["minZ"],
                                maxZ=plot_dim["maxZ"])
        
        self.RayTr = PyRaytracer()
        # Define Grid
        self.grid_dim = dict(nx=int((self.plot_dim['maxX'] - self.plot_dim['minX']) / self.vox_dim),
                             ny=int((self.plot_dim['maxY'] - self.plot_dim['minY']) / self.vox_dim),
                             nz=int((self.plot_dim['maxZ'] - self.plot_dim['minZ']) / self.vox_dim))
        min_bound = np.array([self.plot_dim['minX'], self.plot_dim['minY'], self.plot_dim['minZ']])
        max_bound = np.array([self.plot_dim['maxX'], self.plot_dim['maxY'], self.plot_dim['maxZ']])
        self.RayTr.defineGrid(min_bound, max_bound, self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz'],
                              self.vox_dim)

    def prepare_input(self):
        """
        Prepare input file mappings from project folder and scan directory.

        Scans the RIEGL project folder for `.rdbx`, `.rxp`, `.DAT`, and `.png` files,
        and builds dictionaries that map scan IDs to each file type to prepare raytracing.
        """

        # get all rdbx
        self.rdbx_scans = {}
        self.scan_pos_to_name = {}

        rdbx_scan_list = glob.glob(os.path.join(self.riscan_folder, "project.rdb", "SCANS", "**"))
        if len(rdbx_scan_list) == 0:
            raise ValueError(f'No rdbx files found in riscan folder {self.riscan_folder}. Please check the path and ensure it contains RIEGL scan data with .rdbx files. Path checked: {os.path.join(self.riscan_folder, "project.rdb", "SCANS", "**")}')

        for folder in rdbx_scan_list:
            if "@" in folder:
                # indicates deleted folder
                self.logger.info(f"Skipping deleted folder {folder}")
                continue

            folder_name = os.path.basename(folder)
            scan_pos_base = folder_name[:10]

            if self.exclude_scan_pattern is not None and self.exclude_scan_pattern in folder_name:
                self.logger.info(f"Excluding scan {folder_name} based on exclude_scan_pattern in config")
                continue

            rdbx_folders = glob.glob(os.path.join(folder, "SINGLESCANS", "**"))
            rdbx_folders = [folder for folder in rdbx_folders if "residual" not in folder]
            # skip deleted scans as well
            rdbx_folders = [folder for folder in rdbx_folders if "@" not in folder]

            # if multiple scans were left, we take latest one as this is usually the best one
            # not sure how to make this customizable, should warn in documentation maybe
            if len(rdbx_folders) > 1:
                self.logger.debug(f"multiple rdbx folders found for position {scan_pos_base}, taking latest one ( {rdbx_folders} )")
                rdbx_folders = sorted(rdbx_folders)
            if len(rdbx_folders) == 0:
                self.logger.warning(f"no rdbx folder found for position {scan_pos_base}, skipping.")
                self.logger.debug(f"Path checked: {os.path.join(folder, 'SINGLESCANS', '**')}")
                continue

            rdbx_folder_final = rdbx_folders[-1]
            self.scan_pos_to_name[scan_pos_base] = os.path.basename(rdbx_folder_final)

            rdbx_files = glob.glob(os.path.join(rdbx_folder_final, "*.rdbx"))
            rdbx_files = [file for file in rdbx_files if "residual" not in file]

            if len(rdbx_files) > 1:
                self.logger.warning(f"multiple rdbx files for single scan found ({scan_pos_base}), skipping.")
                self.logger.debug(f"Path checked: {os.path.join(rdbx_folder_final, '*.rdbx')}")
                continue
            if len(rdbx_files) == 0:
                self.logger.warning(f"no rdbx files for scan found ({scan_pos_base}), skipping.")
                self.logger.debug(f"Path checked: {os.path.join(rdbx_folder_final, '*.rdbx')}")
                continue

            rdbx_final = rdbx_files[0]
            self.rdbx_scans[scan_pos_base] = rdbx_final
        
        # get transform files
        # TODO: option for csv? not priority
        self.transform_files = {}
        for scan_pos_name in self.rdbx_scans:
            transform_file = glob.glob(os.path.join(self.riscan_folder, f'{scan_pos_name}.DAT'))
            if len(transform_file) > 1:
                # should never happen
                self.logger.warning(f"Multiple DAT files found for {scan_pos_name}, taking random one.") 
                self.logger.debug(f"Path checked: {os.path.join(self.riscan_folder, f'{scan_pos_name}.DAT')}")
            if len(transform_file) == 0:
                self.logger.warning(f"No DAT files found for {scan_pos_name}, checking if any dat files contain scan position name.") 
                self.logger.debug(f"Path checked: {os.path.join(self.riscan_folder, f'{scan_pos_name}.DAT')}")
                # also check contains
                transform_file = glob.glob(os.path.join(self.riscan_folder, f'*{scan_pos_name}*.DAT'))
                if len(transform_file) == 0:
                    self.logger.warning(f"Cant find DAT file for scan {scan_pos_name}, will not be processed")
                    continue
                self.logger.warning(f"Alternative path used for DAT file: {transform_file[0]}")

            self.transform_files[scan_pos_name] = transform_file[0]

        # get rxp's and optionally previews
        self.rxp_scans = {}
        if self.model_empty_pulses:
            self.png_scans = {}
        
        # look for rxp files matching rdbx files
        for pos in self.rdbx_scans:
            rxp_folder = os.path.join(self.proj_folder, f"{pos}.SCNPOS")

            if not os.path.exists(rxp_folder):
                self.logger.warning(f"rxp folder not found for position {pos}, skipping")
                self.logger.debug(f"Path checked: {rxp_folder}")
                continue
            
            # search for file with exact name of rdbx scan
            scan_name = self.scan_pos_to_name[pos]
            rxp_file = os.path.join(rxp_folder, "scans", f"{scan_name}.rxp")

            if not os.path.exists(rxp_file):
                self.logger.warning(f"rxp file not found for position {pos}, skipping")
                self.logger.debug(f"Path checked: {rxp_file}")
                continue
            
            self.rxp_scans[pos] = rxp_file

            if self.model_empty_pulses:
                png_file = rxp_file[:-4] + ".png"

                if not os.path.exists(png_file):
                    self.logger.warning(f"preview not found for position {scan_pos_base}, skipping")
                    self.logger.debug(f"Path checked: {png_file}")
                    continue

                self.png_scans[pos] = png_file

    def rdbx_rxp_to_df(self, rdbx, rxp):
        """
        Convert RDBX and RXP binary scan files into DataFrames with relevant fields.

        Parameters
        ----------
        rdbx : str
            Path to the `.rdbx` file containing return information.
        rxp : str
            Path to the `.rxp` file containing pulse information.

        Returns
        -------
        tuple of pandas.DataFrame  
            A tuple of two DataFrames:  
            - df_rdbx : Returns with coordinates and beam info.  
            - df_rxp : Pulses with origin and beam direction.  
        """
        columns_rxp = ["beam_origin_x", "beam_origin_y", "beam_origin_z", "beam_direction_x", "beam_direction_y", "beam_direction_z", "scanline", "scanline_idx", "timestamp"]
        subset_rxp = {k: rxp.pulses[k] for k in columns_rxp}
        df_rxp = pd.DataFrame.from_dict(subset_rxp)

        columns_rdbx = ["x", "y", "z", "scanline", "scanline_idx", "reflectance", "target_index", "target_count"]
        subset_rdbx = {k: rdbx.points[k] for k in columns_rdbx}
        df_rdbx = pd.DataFrame.from_dict(subset_rdbx)

        min_scanline = df_rdbx[["scanline"]].to_numpy().min()
        if min_scanline < -1:
            # scanline in rdbx is in reverse (not sure if this is due to rotation of scanner or just bug)
            # shift (for vis)
            df_rdbx[["scanline"]] = df_rdbx[["scanline"]] + abs(df_rdbx[["scanline"]].min())
            max_scanline = df_rdbx[["scanline"]].to_numpy().max()
            # then invert rxp scanline + shift
            df_rxp[["scanline"]] = df_rxp[["scanline"]]*(-1)
            df_rxp[["scanline"]] = df_rxp[["scanline"]] + max_scanline

            if df_rxp[["scanline"]].to_numpy().max() > df_rdbx[["scanline"]].to_numpy().max():
                # drop last column
                df_rxp.drop(df_rxp.loc[df_rxp['scanline'] > df_rdbx["scanline"].to_numpy().max()].index, inplace=True)

        return df_rdbx, df_rxp

    def merge_df_rdbx_rxp(self, df_rdbx, df_rxp):
        """
        Merge return and pulse DataFrames to separate hits and empty pulses.

        Parameters
        ----------
        df_rdbx : pandas.DataFrame
            Return data from RDBX file.
        df_rxp : pandas.DataFrame
            Pulse data from RXP file.

        Returns
        -------
        tuple of pandas.DataFrame
            A tuple of two DataFrames:  
            - df_hit : Pulses with associated returns.  
            - df_empty : Pulses without any return.  
        """
        merged_df = df_rxp.merge(df_rdbx, how="left", on=["scanline", "scanline_idx"])

        na_mask = merged_df["reflectance"].isna()
        point_df = merged_df[~na_mask]
        empty_pulse_df = merged_df[na_mask]

        # for points, keep x,y,z, beam_origin and reflectance, timestamp
        if not self.debug:
            point_df = point_df.drop(["scanline", "scanline_idx", "beam_direction_x", "beam_direction_y", "beam_direction_z"], axis=1)
        else:
            # if debug: also keep beam_direction to check colinearity
            point_df = point_df.drop(["scanline", "scanline_idx"], axis=1)

        # for empty pulses, just keep beam_origin and beam_direction and timestamp
        empty_pulse_df = empty_pulse_df.drop(["x", "y", "z", "reflectance"], axis=1)

        return point_df, empty_pulse_df

    def mask_empty_pulses_preview(self, df_empty, preview_png, max_scanline_idx, max_scanline):
        """
        Mask out empty pulses that fall within occluded regions using preview image.

        Parameters
        ----------
        df_empty : pandas.DataFrame
            DataFrame of empty pulses.
        preview_png : str
            Path to preview PNG image showing occlusion pattern.
        max_scanline_idx : int
            Maximum scanline index in current scan.
        max_scanline : int
            Maximum scanline in current scan.

        Returns
        -------
        pandas.DataFrame
            Filtered DataFrame with only empty pulses in non-occluded regions.
        """
        image = cv2.imread(preview_png)
        image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)  # Convert to RGB

        blue_lower = np.array([0, 0, 255])
        blue_upper = np.array([0, 0, 255])

        blue_mask = cv2.inRange(image, blue_lower, blue_upper)
        blue_mask = blue_mask > 0

        h,w = blue_mask.shape

        scanline_idx_mask_idxs, scanline_mask_idxs = blue_mask.nonzero()

        h_scale = (max_scanline_idx+1)/h
        w_scale= (max_scanline+1)/w

        scanline_idx_lower_h = np.ceil(scanline_idx_mask_idxs*h_scale).astype(int)
        scanline_idx_upper_h = np.floor((scanline_idx_mask_idxs+1)*h_scale).astype(int)
        scanline_lower_h = np.ceil(scanline_mask_idxs*w_scale).astype(int)
        scanline_upper_h = np.floor((scanline_mask_idxs+1)*w_scale).astype(int)

        img_upscaled_manual = np.zeros(shape=(max_scanline_idx+1, max_scanline+1), dtype=bool)

        for i in range(len(scanline_lower_h)):
            img_upscaled_manual[scanline_idx_lower_h[i]:scanline_idx_upper_h[i], scanline_lower_h[i]:scanline_upper_h[i]] = 1

        # filter empty pulses df
        mask = np.any(
            (df_empty["scanline"].values[:, None] >= scanline_lower_h) & 
            (df_empty["scanline"].values[:, None] <= scanline_upper_h) & 
            (df_empty["scanline_idx"].values[:, None] >= scanline_idx_lower_h) & 
            (df_empty["scanline_idx"].values[:, None] <= scanline_idx_upper_h), 
            axis=1
        )

        # Apply the mask to keep only rows that match
        df_filtered = df_empty[mask]

        # also mask out the lower few lines
        # on upright scans, this consistently seems to be the first 5 scanlines
        # on tilts, sometimes up to 50-60
        # -> filter out lowest 80 to be sure
        df_filtered_lower = df_filtered[df_filtered["scanline_idx"]<max_scanline_idx-80]

        if self.debug:
            # write image to output folder
            ofolder = os.path.join(self.odir, "debug")
            if not os.path.exists(ofolder):
                os.makedirs(ofolder)
            opath = os.path.join(ofolder, f"preview_mask_{os.path.basename(preview_png)}")
            self.logger.debug(f"Filtering done using preview {preview_png}, saved at {opath}")

            plot_riegl_grid(df_empty, max_scanline, max_scanline_idx, blue_mask, out_path=opath)

        return df_filtered_lower

    def test_colinearity(self, point_df, n_points=None):
        """
        Check geometric collinearity between pulse origin, beam direction, and return point.

        Parameters
        ----------
        point_df : pandas.DataFrame
            DataFrame containing pulse origin, direction, and return coordinates.
        n_points : int or None, optional
            If specified, tests only the first `n_points` entries.

        Returns
        -------
        int
            number of points failing collinearity check
        """
        def check_parallel(beam_origin, beam_direction, point, epsilon=1e-6):
            vector_point_origin = point - beam_origin
            return (np.dot(beam_direction, vector_point_origin))/(np.linalg.norm(vector_point_origin)*np.linalg.norm(beam_direction)) > 1 - epsilon

        beam_origin = point_df[["beam_origin_x", "beam_origin_y", "beam_origin_z"]].to_numpy()
        beam_direction = point_df[["beam_direction_x", "beam_direction_y", "beam_direction_z"]].to_numpy()
        point = point_df[["x", "y", "z"]].to_numpy()

        count = 0
        if n_points is None:
            n_points = len(point_df)
        for i in range(n_points):
            if n_points is None:
                # check all points
                idx = i
            else:
                # generate random index to check
                idx = randrange(len(point_df))
            if not check_parallel(beam_origin[idx,:], beam_direction[idx,:], point[idx,:]):
                count += 1
        return count
        

    def do_raytracing(self):
        """
        Perform voxel-based ray tracing for all scan positions.

        For each scan position:  
        - Reads RDBX and RXP data.  
        - Optionally models empty pulses using preview image.  
        - Adds pulses to ray tracer and performs traversal.  
        - Clears memory after each scan.  
        """

        for i, scan in enumerate(self.rdbx_scans):
            # read rdbx file for point data

            self.logger.info(f"Processing {scan} ({i+1}/{len(self.rdbx_scans)})")

            if scan not in self.transform_files:
                self.logger.info(f"Transform file not found for pos {scan}, skipping.")
                continue

            self.logger.info(f"Reading RDBX and RXP")


            # read rdbx and optionally rxp
            rdbx = riegl_io.RDBFile(self.rdbx_scans[scan], transform_file=self.transform_files[scan])
            if scan in self.rxp_scans:
                rxp = riegl_io.RXPFile(self.rxp_scans[scan], transform_file=self.transform_files[scan])
            else:
                self.logger.warning(f"RXP not found for pos {scan}, skipping.")
                continue

            if self.debug:
                matrix = rdbx.transform
                self.logger.debug(f"Transformation file: {self.transform_files[scan]}")
                self.logger.debug(f"Transformation matrix: {matrix}")
            
            self.logger.info(f"Merging RDBX and RXP")

            # merge rxp and rdbx
            df_rdbx, df_rxp = self.rdbx_rxp_to_df(rdbx, rxp)
            point_df, empty_pulse_df = self.merge_df_rdbx_rxp(df_rdbx, df_rxp)
            # get max scanline and idx
            max_scanline = max(df_rdbx[["scanline"]].to_numpy().max(), df_rxp[["scanline"]].to_numpy().max())
            max_scanline_idx= max(df_rdbx[["scanline_idx"]].to_numpy().max(), df_rxp[["scanline_idx"]].to_numpy().max())

            # test colinearity to see if rdbx and rxp merge succesfull
            if self.debug:
                # n_tested = len(point_df) # TODO: adapt
                n_tested = 100000
                self.logger.debug(f"Testing collinearity for {n_tested} points")
                n = self.test_colinearity(point_df, n_points=n_tested)
                if n > 0:
                    self.logger.warning(f"Collinearity test for {scan} returned {n} non-colinear points out of {n_tested} tested")

            if self.model_empty_pulses:
                self.logger.info("Reading preview and masking empty pulses")
                if scan not in self.png_scans:
                    raise ValueError(f"Model empty pulses on but scan preview not found for {scan}, exiting.")
                
                empty_pulse_df = self.mask_empty_pulses_preview(empty_pulse_df, self.png_scans[scan], max_scanline_idx, max_scanline)

            self.logger.info("Adding point data")

            x = point_df["x"].to_numpy()
            y = point_df["y"].to_numpy()
            z = point_df["z"].to_numpy()
            gps_time = point_df["timestamp"].to_numpy()
            self.logger.debug(f"Number of points in rdbx: {len(df_rdbx)}")
            unique_gps = np.unique(gps_time)
            self.logger.debug(f"Number of unique gps entries in point_df: {len(unique_gps)}")
            self.logger.debug(f"Number of pulses in rxp: {len(df_rxp)}, number of empty pulses: {len(empty_pulse_df)}, diff: {len(df_rxp)-len(empty_pulse_df)}")

            return_number = point_df["target_index"].to_numpy()
            number_of_returns = point_df["target_count"].to_numpy()

            sensor_x = point_df["beam_origin_x"].to_numpy()
            sensor_y = point_df["beam_origin_y"].to_numpy()
            sensor_z = point_df["beam_origin_z"].to_numpy()

            self.RayTr.addPointData(x, y, z, sensor_x, sensor_y, sensor_z, gps_time, return_number, number_of_returns)

            self.logger.info("Fixing incomplete pulses (number of returns not correct, likely due to filtered points)")

            self.RayTr.cleanUpPulseDataset()

            self.RayTr.getPulseDatasetReport()

            self.logger.info("Perform raytracing")
            tic = time.time()
            self.RayTr.doRaytracing()
            toc = time.time()
            self.logger.info("Time elapsed for raytracing: {:.2f} seconds".format(toc - tic))

            if self.model_empty_pulses:

                self.logger.info("Running raytracing for empty pulses")
                tic = time.time()

                sensor_x = empty_pulse_df["beam_origin_x"].to_numpy()
                sensor_y = empty_pulse_df["beam_origin_y"].to_numpy()
                sensor_z = empty_pulse_df["beam_origin_z"].to_numpy()
                direction_x = empty_pulse_df["beam_direction_x"].to_numpy()
                direction_y = empty_pulse_df["beam_direction_y"].to_numpy()
                direction_z = empty_pulse_df["beam_direction_z"].to_numpy()
                gps_time = empty_pulse_df["timestamp"].to_numpy()

                self.RayTr.addEmptyPulseData(sensor_x, sensor_y, sensor_z, direction_x, direction_y, direction_z, gps_time)
                self.RayTr.doRaytracingEmptyPulses()
                toc = time.time()
                self.logger.info("Time elapsed for raytracing empty pulses: {:.2f} seconds".format(toc - tic))

            self.RayTr.clearPulseDataset()

        self.logger.info("Report on traversal:")
        self.RayTr.reportOnTraversal()
        
        return
    
    def determine_grid(self, buffer):
        # TODO: automatically derive grid from scan_positions
        # Take predefined buffer in x,y,z direction (configurable)
        # then find x,y,z extremes of scan positions and define grid based on this
        # can have auto_dim flag in config
        # if auto_dim, xyz_buf must be defined
        # if not, plot_dim must be defined
        


        # should probably be util/helper function in seperate module
        raise NotImplementedError

    def save_raytracing_output(self):
        """
        Save ray tracing results (`Nhit`, `Nmiss`, `Nocc`, `Classification`) to disk.

        Saves the voxel outputs as `.npy` arrays, and optionally writes `.ply` files
        for visualization if `self.output_voxels` is True.
        """
        self.logger.info("Saving output")
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
        np.save(f"{self.odir}/Nhit.npy", self.Nhit)
        np.save(f"{self.odir}/Nmiss.npy", self.Nmiss)
        np.save(f"{self.odir}/Nocc.npy", self.Nocc)
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

        np.save(f"{self.odir}/Classification.npy", self.Classification)
        toc = time.time()
        print("Elapsed Time: " + str(toc - tic) + " seconds")

        # write ply file
        if self.output_voxels:
            print("Saving Occlusion Outputs As .ply")
            tic = time.time()
            verts, faces = prepare_ply(self.vox_dim, self.plot_dim, self.Nhit)
            ost.write_ply(f"{self.odir}/Nhit.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.plot_dim, self.Nmiss)
            ost.write_ply(f"{self.odir}/Nmiss.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.plot_dim, self.Nocc)
            ost.write_ply(f"{self.odir}/Nocc.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            # TODO: TEMP disable: these take up a lot of space, so disable for now
            # verts, faces = prepare_ply(self.vox_dim, self.plot_dim, self.Classification)
            # ost.write_ply(f"{self.odir}/Classification.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            # self.occl = np.zeros(shape=self.Classification.shape)
            # x4, y4, z4 = np.where(self.Classification == 4)
            # self.occl[x4, y4, z4] = self.Classification[x4, y4, z4]
            # verts, faces = prepare_ply(self.vox_dim, self.plot_dim, self.occl)
            # ost.write_ply(f"{self.odir}/Occl.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            toc = time.time()
            print("Elapsed Time: " + str(toc - tic) + " seconds")
