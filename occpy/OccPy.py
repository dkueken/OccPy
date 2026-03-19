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

# plotting functions
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import to_rgb
import seaborn as sns

from raytr import PyRaytracer
from occpy.TerrainModel import TerrainModel
from occpy.util import prepare_ply, read_trajectory_file, read_sensorpos_file, interpolate_traj, last_nonzero

is_sorted = lambda a: np.all(a[:-1] <= a[1:])

# Function to darken an RGB color
def darken_color(color, amount=0.6):
    r, g, b = to_rgb(color)
    return (r * amount, g * amount, b * amount)

""" Moved to visualization module. delete once everything is working.
def get_Occlusion_ProfileFigure(Classification, plot_dim, vox_dim, out_dir, low_thresh=0, vertBuffer=0, max_percentage=100, fig_prop=None, show_plots=False):

    grid_dim = (int((plot_dim[3] - plot_dim[0]) / vox_dim), int((plot_dim[4] - plot_dim[1]) / vox_dim),
                int((plot_dim[5] - plot_dim[2]) / vox_dim))

    vert_vect = np.arange(start=low_thresh, stop=Classification.shape[2] * vox_dim, step=vox_dim)
    Classification = Classification[:,:,int(low_thresh / vox_dim):]
    # a hack to make sure that vert_vect is of the same length as OcclVertProf TODO: this has to be checked if it is generic!


    OcclVertProf = np.sum(Classification == 3, axis=0)
    OcclVertProf = np.sum(OcclVertProf, axis=0)
    OcclVertProf_Rel = OcclVertProf / ((grid_dim[0]) * (grid_dim[1]))

    FilledVertProf = np.sum(Classification == 1, axis=0)
    FilledVertProf = np.sum(FilledVertProf, axis=0)
    FilledVertProf_Rel = FilledVertProf / (grid_dim[0] * grid_dim[1])

    EmptyVertProf = np.sum(np.logical_or(Classification == 2, Classification==0), axis=0)
    EmptyVertProf = np.sum(EmptyVertProf, axis=0)
    EmptyVertProf_Rel = EmptyVertProf / (grid_dim[0] * grid_dim[1])

    heights = vert_vect[0:len(OcclVertProf)]

    percentages = np.column_stack([FilledVertProf_Rel*100, OcclVertProf_Rel*100, EmptyVertProf_Rel*100])
    categories = ['Filled', 'Occluded', 'Empty']
    colors = ['skyblue', 'salmon', 'lightgreen']

    # Compute cumulative percentages for stacking
    cumulative = np.cumsum(percentages, axis=1)

    palette = sns.color_palette('colorblind', n_colors=len(categories))
    palette[2] = (1.0, 1.0, 1.0) # white for empty

    fig, ax = plt.subplots(figsize=fig_prop['fig_size'])

    for i, cat in enumerate(categories):
        left = cumulative[:, i - 1] if i > 0 else np.zeros_like(heights)
        face_color = palette[i]
        edge_color = darken_color(face_color, 0.8)  # slightly darker for lines

        # Fill area
        ax.fill_betweenx(
            heights, left, cumulative[:, i],
            color=face_color, alpha=0.6
        )
        # Outline
        ax.plot(cumulative[:, i], heights, color=edge_color, linewidth=1.5, label="_nolegend_")

    ax.set_xlabel('Percentage of voxels [%]', fontsize=fig_prop['label_size'])
    ax.set_ylabel('Height above ground [m]', fontsize=fig_prop['label_size'])
    ax.set_xlim(0.1,max_percentage)
    ax.set_ylim(0,np.max(heights) + vertBuffer)
    plt.xticks(fontsize=fig_prop['label_size_ticks'])
    plt.yticks(fontsize=fig_prop['label_size_ticks'])
    ax.legend(categories[0:2], fontsize=fig_prop['label_size_ticks'])
    plt.tight_layout()

    plt.savefig(f"{out_dir}/OcclusionVertProf.{fig_prop['out_format']}", dpi=300, format=fig_prop['out_format'])
    if show_plots:
        plt.show(block=True)
    else:
        plt.close()
"""




class OccPy:

    def __init__(self, config_file):
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
        config_file : str
            Path to a JSON configuration file containing processing parameters.
        """

        with open(config_file) as f:
            config = json.load(f)

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
        self.out_dir = config.get("out_dir", "./output")
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

        # copy config file to output directory for record keeping
        with open(os.path.join(self.out_dir, "config.json"), "w") as to:
            json.dump(config, to)


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
        
        
        # initialize RayTr Object and define grid
        self.RayTr = PyRaytracer()
        self.PlotDim = dict(minX=self.plot_dim[0],
                                maxX=self.plot_dim[3],
                                minY=self.plot_dim[1],
                                maxY=self.plot_dim[4],
                                minZ=self.plot_dim[2],
                                maxZ=self.plot_dim[5])
        # TODO: here, we should watch out with non-integer bounds. nx is rounded to nearest multiple of vox_dim, which can lead to edge effects when using non-integer bounds. 
        # desired behaviour would be 
        #       a. not require integer bounds -> then voxel indexing logic needs to be rewritten on c++ side to not use max_bound-min_bound but min_bound + index
        #           in this case, the voxelgrid is started at min_bound always, and the amount of voxels is determined by how many voxels fit into the plot_dim. Max_bound in practice is then not equal to max_bound if the extent is not a multiple of vox_dim, so user needs to be warned about this. But this shouldnt be a huge problem.
        #           some untested code for this case is already ready
        #       b. require integer bounds, round here -> then we need to clearly warn the user, and save the actual bounds used to the output directory in some way
        self.grid_dim = dict(nx=int((self.PlotDim['maxX'] - self.PlotDim['minX']) / self.vox_dim),
                             ny=int((self.PlotDim['maxY'] - self.PlotDim['minY']) / self.vox_dim),
                             nz=int((self.PlotDim['maxZ'] - self.PlotDim['minZ']) / self.vox_dim))
        minBound = np.array([self.PlotDim['minX'], self.PlotDim['minY'], self.PlotDim['minZ']])
        maxBound = np.array([self.PlotDim['maxX'], self.PlotDim['maxY'], self.PlotDim['maxZ']])
        self.RayTr.defineGrid(minBound, maxBound, self.grid_dim['nx'], self.grid_dim['ny'], self.grid_dim['nz'],
                              self.vox_dim)
        
        return


    def define_sensor_pos(self, path2file, delimiter=" ", hdr_time='%time', hdr_scanpos_id='', hdr_x='x', hdr_y='y', hdr_z='z', sens_pos_id_offset=0):
        """
        defines sensor position based on the provided csv file. CSV file needs to include
        Parameters
        ----------
        path2file: str [mandatory]
            path to csv file with sensor position information
        is_mobile: bool [mandatory]
            True or False whether platform is mobile (MLS, ULS) or static (TLS)
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
            self.logger.warning(f"Provided path to sensor position file {path2file} is not absolute. Resolving relative to root folder {self.root_folder}.")
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

    # TODO: ugly workaround for the case where a single laz file from a single TLS position should be run
    def define_sensor_pos_singlePos(self, scan_pos_id, x, y, z):
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

    def prepare_input(self):
        """
        Prepare input data for ray tracing. Checks for each scan file if a link can be made to the scan positions file.
        Can be called manually to check if linking works, or is called by do_raytracing if not manually called before.
        """
        raise NotImplementedError("TODO!")

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

        # TODO: link positions to laz files first (in case of TLS at least), and warn/error if one/more positions cant be linked before processing

        run_raytraycing_after_loading = False
        if os.path.isdir(self.laz_in): # TLS
            ## get list of laz files in input directory
            fCont = glob.glob(f"{self.laz_in}/*.laz")

            for f in fCont:
                # TODO: this shouldnt be done here but before, see TODO above
                # get scan position ID
                scan_name = os.path.basename(f)
                scan_id = scan_name.split(".")[0] # default: use file name without extension as scan ID
                if self.str_idxs_ScanPosID is not None: # if given, use sting indices
                    scan_id = int(scan_name[self.str_idxs_ScanPosID[0]:self.str_idxs_ScanPosID[1]])

                self.logger.debug(f"Using scan ID {scan_id} for file {scan_name}.")

                print(f"###############################")
                print(f"##### Processing {scan_name}...")
                print(f"###############################")

                scanpos_X = self.senspos.loc[self.senspos['ScanPos'] == scan_id, 'sensor_x'].values[0]
                scanpos_Y = self.senspos.loc[self.senspos['ScanPos'] == scan_id, 'sensor_y'].values[0]
                scanpos_Z = self.senspos.loc[self.senspos['ScanPos'] == scan_id, 'sensor_z'].values[0]

                # read in laz file
                tic = time.time()
                with laspy.open(f) as file:

                    count = 0
                    for points in file.chunk_iterator(self.points_per_iter):

                        # if self.single_return is not set, check data to find out if there are multiple returns stored per pulse
                        if self.single_return == None:
                            self.logger.warning(
                                f"TIPP: it was not specified if the point cloud stores single return data. to avoid any "
                                f"issues where there are only single return pulses in the first chunk, put there are "
                                f"multi return pulses in the whole dataset, please set the single_return variable within "
                                f"the define_sensor_pos function")
                            max_ret_num = np.max(points.return_number)
                            if max_ret_num > 1:
                                self.single_return = False
                            else:
                                self.single_return = True

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

                        # for the case of mobile acquisitions, inerpolate trajectory for gps_time
                        if self.is_mobile:
                            # call interpolate function for trajectory to extract sensor position for each gps_time
                            SensorPos = interpolate_traj(self.traj['time'], self.traj['sensor_x'], self.traj['sensor_y'],
                                                              self.traj['sensor_z'], gps_time)

                        else:
                            SensorPos = self.senspos

                            # TODO: figure out, why we have to repeat scan position to the shape of input point cloud and try to implement it, so that we only have to pass one position
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

        if run_raytraycing_after_loading:
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
        np.save(f"{self.out_dir}/Nhit.npy", self.Nhit)
        np.save(f"{self.out_dir}/Nmiss.npy", self.Nmiss)
        np.save(f"{self.out_dir}/Nocc.npy", self.Nocc)
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

        np.save(f"{self.out_dir}/Classification.npy", self.Classification)
        toc = time.time()
        self.logger.info("Elapsed Time: " + str(toc - tic) + " seconds")

        # write ply file
        if self.output_voxels:
            self.logger.info("Saving Occlusion Outputs As .ply")
            tic = time.time()
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nhit)
            ost.write_ply(f"{self.out_dir}/Nhit.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nmiss)
            ost.write_ply(f"{self.out_dir}/Nmiss.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Nocc)
            ost.write_ply(f"{self.out_dir}/Nocc.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.Classification)
            ost.write_ply(f"{self.out_dir}/Classification.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            self.occl = np.zeros(shape=self.Classification.shape)
            x4, y4, z4 = np.where(self.Classification == 4)
            self.occl[x4, y4, z4] = self.Classification[x4, y4, z4]
            verts, faces = prepare_ply(self.vox_dim, self.PlotDim, self.occl)
            ost.write_ply(f"{self.out_dir}/Occl.ply", verts, ['X', 'Y', 'Z', 'data'], triangular_faces=faces)
            toc = time.time()
            self.logger.info("Elapsed Time: " + str(toc - tic) + " seconds")


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
        origin: [int, int, int]
            with grid origins [minx, miny, minz]

        """
        origin = np.asarray(self.RayTr.getGridOrigin(), dtype=np.int32)
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

