import os
import glob

import pandas as pd
import numpy as np

from occpy import riegl_io



class OccPyRIEGL:
    def __init__(self, riscan_folder, proj_folder, model_empty_pulses=False):

        # TODO: other inputs as config file or parameters
        vox_dim = 0.1
        lower_threshold = 1

        self.riscan_folder = riscan_folder
        self.proj_folder = proj_folder
        self.model_empty_pulses = model_empty_pulses


    def prepare_input(self):
        # TODO: test rdbx reading and SOP reading
        # get all rdbx
        rdbx_scans = {}
        scan_pos_mapping = {}
        for folder in glob.glob(os.path.join(self.riscan_folder, "project.rdb", "SCANS", "**")):
            folder_name = os.path.basename(folder)
            scan_pos_base = folder_name[:10]
            if scan_pos_base != folder_name:
                scan_pos_mapping[scan_pos_base] = folder_name

            rdbx_folders = glob.glob(os.path.join(folder, "SINGLESCANS", "**"))
            rdbx_folders = [folder for folder in rdbx_folders if "residual" not in folder]

            if len(rdbx_folders) > 1:
                rdbx_folders = sorted(rdbx_folders)
            if len(rdbx_folders) == 0:
                print(f"WARNING: no rdbx folder found for position {scan_pos_base}, skipping.")
                print(f"Path checked: {os.path.join(folder, 'SINGLESCANS', '**')}")
                continue

            rdbx_folder_final = rdbx_folders[-1]

            rdbx_files = glob.glob(os.path.join(rdbx_folder_final, "*.rdbx"))
            rdbx_files = [file for file in rdbx_files if "residual" not in file]

            if len(rdbx_files) > 1:
                print(f"Warning: multiple rdbx files for single scan found ({scan_pos_base}), skipping.")
                print(f"Path checked: {os.path.join(rdbx_folder_final, "*.rdbx")}")
                continue
            if len(rdbx_files) == 0:
                print(f"Warning: no rdbx files for scan found ({scan_pos_base}), skipping.")
                print(f"Path checked: {os.path.join(rdbx_folder_final, "*.rdbx")}")
                continue

            rdbx_final = rdbx_files[0]
            rdbx_scans[scan_pos_base] = rdbx_final
        
        # get transform files for rxp's
        # we assume rdbx already transformed (SOP backup applied in riscan pro) TODO: could make this customizable as well, not priority
        transform_files = {}
        for scan_pos_name in rdbx_scans:
            transform_file = glob.glob(os.path.join(self.riscan_folder, f'{scan_pos_base}.DAT'))
            if len(transform_file) > 1:
                # should never happen
                print(f"WARNING: Multiple DAT files found for {scan_pos_name}, skipping.") 
                print(f"Path checked: {os.path.join(self.riscan_folder, f'{scan_pos_base}.DAT')}")
            if len(transform_file) == 0:
                # TODO: what is desired behaviour here? report, but then use untransformed in raytracing? or exit out
                print(f"WARNING: No DAT files found for {scan_pos_name}, skipping.") 
                print(f"Path checked: {os.path.join(self.riscan_folder, f'{scan_pos_base}.DAT')}")
            transform_files[scan_pos_name] = transform_file[0]


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
                print(f"WARNING: no rxp file found for position {scan_pos_base}, skipping")
                print(f"Path checked: {os.path.join(folder, 'scans', '*.rxp')}")
                continue

            scan_final = rxp_files[-1]
            self.rxp_scans[scan_pos_base] = scan_final

            if self.model_empty_pulses:
                png_file = scan_final[:-4] + ".png"

                if not os.path.exists(png_file):
                    print(f"WARNING: preview not found for position {scan_pos_base}, skipping")
                    print(f"Path checked: {png_file}")
                    continue

                self.png_scans[scan_pos_base] = png_file
        

        




# TODO: TEMP TEST

if __name__ == "__main__":
    proj_folder = "/Stor1/wout/occlusion/cls_raw/Ficus/2023-04-30_LOT_peru2.PROJ"

    occpy_riegl = OccPyRIEGL("", proj_folder, model_empty_pulses=True)

    occpy_riegl.read_input()




