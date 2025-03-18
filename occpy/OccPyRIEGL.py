import os
import glob
import logging

import pandas as pd
import numpy as np

from occpy import riegl_io

from raytr import PyRaytracer


# TODO:
# 4. custom configexception? for when config is not correct, like with autodim and plot_dim

class OccPyRIEGL:
    def __init__(self, riscan_folder, proj_folder, model_empty_pulses=False, csv_positions=None, verbose=False, debug=False):

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
        # ---

        # TODO: TEMP: read csv with scan positions
        if csv_positions is not None:
            self.csv_positions = csv_positions

        self.prepare_input()

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
        # TODO: test rdbx reading and SOP reading
        # get all rdbx
        self.rdbx_scans = {}
        self.scan_pos_mapping = {}
        for folder in glob.glob(os.path.join(self.riscan_folder, "project.rdb", "SCANS", "**")):
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

    def do_raytracing(self):

        # TODO: implement raytracing similar to occpy class

        # tests to do:
            # 1. 1 rdbx with single scan position and grid around scan
            # 2. multiple rdbx's with multiple scan position
            # 3. 1 rdbx+rxp with pulse from beam_origin to point (requires editing raytracer)
            # 4. multiple rdbx's + rxp's
            # 5. single rdbx + rxp with empty pulse modelling (requires editing raytracer further)
            # 6. multiple rdbx's + rxp's with empty pulse modelling

        for scan in self.rdbx_scans:
            # read rdbx file for point data

            scan_test = ["ScanPos001"]
            if scan not in scan_test:
                continue

            self.logger.info(f"Processing {scan}")



            # RDBX_file = riegl_io.RDBFile()
            

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




# TODO: TEMP TEST

if __name__ == "__main__":
    proj_folder = "/Stor1/wout/occlusion/cls_raw/Ficus/2023-04-30_LOT_peru2.PROJ"

    occpy_riegl = OccPyRIEGL("", proj_folder, model_empty_pulses=True, debug=True)

    # occpy_riegl.read_input()




