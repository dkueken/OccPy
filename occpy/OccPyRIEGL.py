import os
import glob
import logging
import cv2
import time

import pandas as pd
import numpy as np
import OSToolBox as ost

from occpy import riegl_io
from occpy.PreparePly import prepare_ply

from raytr import PyRaytracer

# TODO: inherit from OccPy class? only save and viz functions are the same
# alternatively, viz functions could/should be in seperate module imo

class OccPyRIEGL:
    def __init__(self, riscan_folder, proj_folder, model_empty_pulses=False, verbose=False, debug=False, odir=None, output_voxels=False):

        # TODO: other inputs as config file or parameters
        self.vox_dim = 0.1
        self.lower_threshold = 1
        self.buffer = [10,10,10,10,5,40] # [x_min, x_max, y_min, y_max, z_min, z_max]
        self.plot_dim = dict(minX=-10,
                                maxX=10,
                                minY=-10,
                                maxY=10,
                                minZ=-5,
                                maxZ=30)

        self.riscan_folder = riscan_folder
        self.proj_folder = proj_folder
        self.model_empty_pulses = model_empty_pulses

        self.verbose = verbose
        self.debug = debug

        self.output_voxels = output_voxels

        if odir is None:
            odir = os.path.join(os.getcwd(), "out")
            if not os.path.exists(odir):
                os.makedirs(odir, exist_ok=True)
        self.out_dir = odir 

        # -- config logging 
        if debug:
            logging_level = logging.DEBUG
        elif verbose:
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

        # --- prepare input

        self.prepare_input()

        # --- init dimensions and raytracer

        # TODO: TEST
        # if self.auto_dim:
        #     if csv_positions is None or self.buffer is None:
        #         self.logger.error("Must provide csv file and buffer in config if auto_dim is false")
        #         os._exit(1)
        #     plot_dim = self.determine_grid(csv_positions)
        # else:
        #     if plot_dim is None:
        #         self.logger.error("must provide plot dimensions in config if auto_dim is false")
        #         os._exit(1)
        #     # format of plot_dim in config: [minX, minY, minZ, maxX, maxY, maxZ]
        #     self.plot_dim = dict(minX=plot_dim[0],
        #                         maxX=plot_dim[3],
        #                         minY=plot_dim[1],
        #                         maxY=plot_dim[4],
        #                         minZ=plot_dim[2],
        #                         maxZ=plot_dim[5])
        
        self.grid_dim = dict(nx=int((self.plot_dim['maxX'] - self.plot_dim['minX']) / self.vox_dim),
                             ny=int((self.plot_dim['maxY'] - self.plot_dim['minY']) / self.vox_dim),
                             nz=int((self.plot_dim['maxZ'] - self.plot_dim['minZ']) / self.vox_dim))
        
        self.RayTr = PyRaytracer()

        # Define Grid
        min_bound = np.array([self.plot_dim['minY'], self.plot_dim['minX'], self.plot_dim['minZ']])
        max_bound = np.array([self.plot_dim['maxY'], self.plot_dim['maxX'], self.plot_dim['maxZ']])
        self.RayTr.defineGrid(min_bound, max_bound, self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz'],
                              self.vox_dim)

    def prepare_input(self):

        # get all rdbx
        self.rdbx_scans = {}
        self.scan_pos_mapping = {}

        rdbx_scan_list = glob.glob(os.path.join(self.riscan_folder, "project.rdb", "SCANS", "**"))
        if len(rdbx_scan_list) == 0:
            self.logger.error(f"No rdbx files found in riscan folder {self.riscan_folder}")
            self.logger.debug(f'path checked: {os.path.join(self.riscan_folder, "project.rdb", "SCANS", "**")}')
            os._exit(1)

        for folder in rdbx_scan_list:
            folder_name = os.path.basename(folder)
            scan_pos_base = folder_name[:10]
            if scan_pos_base != folder_name:
                self.scan_pos_mapping[scan_pos_base] = folder_name

            rdbx_folders = glob.glob(os.path.join(folder, "SINGLESCANS", "**"))
            rdbx_folders = [folder for folder in rdbx_folders if "residual" not in folder]

            if len(rdbx_folders) > 1:
                rdbx_folders = sorted(rdbx_folders)
            if len(rdbx_folders) == 0:
                self.logger.warning(f"no rdbx folder found for position {scan_pos_base}, skipping.")
                self.logger.debug(f"Path checked: {os.path.join(folder, 'SINGLESCANS', '**')}")
                continue

            rdbx_folder_final = rdbx_folders[-1]

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
        
        # get transform files for rxp's
        # we assume rdbx already transformed (SOP backup applied in riscan pro) TODO: could make this customizable as well, not priority
        self.transform_files = {}
        for scan_pos_name in self.rdbx_scans:
            transform_file = glob.glob(os.path.join(self.riscan_folder, f'{scan_pos_base}.DAT'))
            if len(transform_file) > 1:
                # should never happen
                self.logger.warning(f"Multiple DAT files found for {scan_pos_name}, skipping.") 
                self.logger.debug(f"Path checked: {os.path.join(self.riscan_folder, f'{scan_pos_base}.DAT')}")
            if len(transform_file) == 0:
                # TODO: what is desired behaviour here? report, but then use untransformed in raytracing? or exit out
                self.logger.warning(f"No DAT files found for {scan_pos_name}, skipping.") 
                self.logger.debug(f"Path checked: {os.path.join(self.riscan_folder, f'{scan_pos_base}.DAT')}")
            self.transform_files[scan_pos_name] = transform_file[0]

        # get rxp's and optionally previews
        self.rxp_scans = {}
        if self.model_empty_pulses:
            self.png_scans = {}
        for folder in glob.glob(os.path.join(self.proj_folder, "**.SCNPOS")):
            folder_name = os.path.basename(folder)
            scan_pos_base = folder_name[:10]
            rxp_files = glob.glob(os.path.join(folder, "scans", "*.rxp"))
            # filter out .residuals and .mon
            rxp_files = [file for file in rxp_files if not (file.endswith(".residual.rxp") or file.endswith(".mon.rxp")) ]

            if len(rxp_files) > 1:
                # keep only latest scan (TODO: make customizable? in our data latest is always best scan)
                rxp_files = sorted(rxp_files)

            if len(rxp_files) == 0:
                self.logger.warning(f"no rxp file found for position {scan_pos_base}, skipping")
                self.logger.debug(f"Path checked: {os.path.join(folder, 'scans', '*.rxp')}")
                continue

            scan_final = rxp_files[-1]
            self.rxp_scans[scan_pos_base] = scan_final

            if self.model_empty_pulses:
                png_file = scan_final[:-4] + ".png"

                if not os.path.exists(png_file):
                    self.logger.warning(f"preview not found for position {scan_pos_base}, skipping")
                    self.logger.debug(f"Path checked: {png_file}")
                    continue

                self.png_scans[scan_pos_base] = png_file

    def rdbx_rxp_to_df(self, rdbx, rxp):
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

        # TODO: debug: save fig of mask to debug folder

        return df_filtered

    def test_colinearity(self, point_df, n_points=None):
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
            if not check_parallel(beam_origin[i,:], beam_direction[i,:], point[i,:]):
                count += 1
        return count
        

    def do_raytracing(self):

        # TODO:
        # tests to do:
            # 5. single rdbx + rxp with empty pulse modelling (requires editing raytracer further)
            # 6. multiple rdbx's + rxp's with empty pulse modelling

        for scan in self.rdbx_scans:
            # read rdbx file for point data
            
            # TODO: test
            scan_test = ["ScanPos001", "ScanPos002"]
            if scan not in scan_test:
                continue

            self.logger.info(f"Processing {scan}")

            self.logger.info(f"Reading RDBX and RXP")
            
            # TODO: test: do we actually need to give the transform file here? as the rdb points should already be transformed!!
            # can test with colinearity test, do for all files (will take long but worth it)
            rdbx = riegl_io.RDBFile(self.rdbx_scans[scan], self.transform_files[scan])
            max_scanline = rdbx.maxc
            max_scanline_idx = rdbx.maxr

            if scan in self.rxp_scans:
                rxp = riegl_io.RXPFile(self.rxp_scans[scan], self.transform_files[scan])
            else:
                self.logger.warning(f"RXP not found for pos {scan}, skipping.")
                continue
            
            self.logger.info(f"Merging RDBX and RXP")

            df_rdbx, df_rxp = self.rdbx_rxp_to_df(rdbx, rxp)
            point_df, empty_pulse_df = self.merge_df_rdbx_rxp(df_rdbx, df_rxp)

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
                    self.logger.error(f"Model empty pulses on but scan preview not found for {scan}, exiting.")
                    os._exit(1)
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

            self.RayTr.clearPulseDataset()

        self.logger.info("Report on traversal:")
        self.RayTr.reportOnTraversal()



            # TODO: if model_empty_pulses, do raytracing with empty pulses (likely have to modify Raytracer.cpp)
        
        return
    
    def determine_grid(self, csv_positions, buffer):
        # TODO: automatically derive grid from scan_positions
        # Take predefined buffer in x,y,z direction (configurable)
        # then find x,y,z extremes of scan positions and define grid based on this
        # can have auto_dim flag in config
        # if auto_dim, xyz_buf must be defined
        # if not, plot_dim must be defined
        


        # should probably be util/helper function in seperate module
        raise NotImplementedError

    def save_raytracing_output(self):
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
        np.save(f"{self.out_dir}/Nhit.npy", self.Nhit)
        np.save(f"{self.out_dir}/Nmiss.npy", self.Nmiss)
        np.save(f"{self.out_dir}/Nocc.npy", self.Nocc)
        toc = time.time()
        print("Elapsed Time: {:.2f} seconds".format(toc - tic))

        # Create Classification grid
        print("Classify Grid")
        tic = time.time()
        self.Classification = np.zeros((self.grid_dim['ny'], self.grid_dim['nx'], self.grid_dim['nz']), dtype=int)

        self.Classification[np.logical_and.reduce((self.Nhit > 0, self.Nmiss >= 0, self.Nocc >= 0))] = 1  # voxels that were observed
        self.Classification[np.logical_and.reduce((self.Nhit == 0, self.Nmiss > 0, self.Nocc >= 0))] = 2  # voxels that are empty
        self.Classification[
            np.logical_and.reduce((self.Nhit == 0, self.Nmiss == 0, self.Nocc > 0))] = 3  # voxels that are hidden (occluded)
        self.Classification[np.logical_and.reduce((self.Nhit == 0, self.Nmiss == 0,
                                            self.Nocc == 0))] = 4  # voxels that were not observed # TODO: Figure out, why this overwrites voxels that are classified as occluded! -> this was because np.logical_and only takes in 2 arrays as input, not 3! use np.logical_and.reduce() for that!

        np.save(f"{self.out_dir}/Classification.npy", self.Classification)
        toc = time.time()
        print("Elapsed Time: " + str(toc - tic) + " seconds")

        # write ply file
        if self.output_voxels:
            print("Saving Occlusion Outputs As .ply")
            tic = time.time()
            verts, faces = prepare_ply(self.vox_dim, self.plot_dim, self.Nhit)
            ost.write_ply(f"{self.out_dir}/Nhit.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.plot_dim, self.Nmiss)
            ost.write_ply(f"{self.out_dir}/Nmiss.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.plot_dim, self.Nocc)
            ost.write_ply(f"{self.out_dir}/Nocc.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            # TODO: TEMP disable: these take up a lot of space, so disable for now
            # verts, faces = prepare_ply(self.vox_dim, self.plot_dim, self.Classification)
            # ost.write_ply(f"{self.out_dir}/Classification.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            # self.occl = np.zeros(shape=self.Classification.shape)
            # x4, y4, z4 = np.where(self.Classification == 4)
            # self.occl[x4, y4, z4] = self.Classification[x4, y4, z4]
            # verts, faces = prepare_ply(self.vox_dim, self.plot_dim, self.occl)
            # ost.write_ply(f"{self.out_dir}/Occl.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            toc = time.time()
            print("Elapsed Time: " + str(toc - tic) + " seconds")


# TODO: TEMP TEST

if __name__ == "__main__":
    proj_folder = "/Stor1/wout/occlusion/cls_raw/Ficus/2023-04-30_LOT_peru2.PROJ"
    riscan_folder  = "/Stor1/wout/occlusion/oxa_occpy_test.RiSCAN"

    # odir = "/Stor1/wout/occlusion/output_test/OXA"
    odir = "./test_out/OXA/2_pos_origin"
    if not os.path.exists(odir):
        os.makedirs(odir, exist_ok=True)

    OUTPUT_VOXELS = True

    occpy_riegl = OccPyRIEGL(riscan_folder, proj_folder, model_empty_pulses=False, debug=True, odir=odir, output_voxels=OUTPUT_VOXELS)

    occpy_riegl.do_raytracing()

    occpy_riegl.save_raytracing_output()




