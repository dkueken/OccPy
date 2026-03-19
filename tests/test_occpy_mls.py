import numpy as np
import json
import os
import pooch
from occpy.OccPy import OccPy
from test_occpy_tls import find_project_root

def test_occpy_run_mls():

    # Download test data
    p = pooch.create(
        path=pooch.os_cache("occpy_test_data_mls"),
        base_url="https://zenodo.org/records/18095821/files/",
        registry={"MLS_demo.zip": "md5:f217fe4753e97ab9ccf40a2a2a351834"}
    )
    p.fetch("MLS_demo.zip", processor=pooch.Unzip(members=["MLS_demo"]), progressbar=True)
    data_path = os.path.join(p.path, "MLS_demo.zip.unzip", "MLS_demo")

    occpy_root = find_project_root()

    # load the json config file
    with open(os.path.join(occpy_root, 'config', 'settings_MLS_tutorial.JSON'), "r") as file:
        settings_orig = json.load(file)

    # copy original settings file
    settings = settings_orig.copy()
    # alter root folder
    settings['root_folder'] = os.path.abspath(occpy_root)
    # update input paths to be relative to the root folder
    settings['laz_in'] = os.path.join(data_path, 'LAZ', 'MLS_TestData_20perc_FP10_2025.laz')
    settings['tif_in']['DTM'] = os.path.join(data_path, 'Grids', 'Ramerenwald_DTM_20250305.tif')
    settings['tif_in']['DSM'] = os.path.join(data_path, 'Grids', 'Ramerenwald_DSM_20250305.tif')
    settings['ScanPos'] = os.path.join(data_path, 'ScanPos', 'MLS_TestData_traj_FP10_2025.txt')
    settings['out_dir'] = os.path.join(settings['root_folder'], 'output', 'MLS')

    # initiate OccPy object from config dict
    test = OccPy(config=settings)

    # Define sensor position
    test.define_sensor_pos(path2file=settings['ScanPos'],  # Path to the trajectory file
                           delimiter=" ",  # delimiter used in the trajectory file
                           hdr_time='//world_time',  # column header for the time information in the trajectory file
                           hdr_x='x',  # column header for the x coordinate in the trajectory file
                           hdr_y='y',  # column header for the y coordinate in the trajectory file
                           hdr_z='z', )  # column header for the z coordinate in the trajectory file


    # run occpy and measure the time
    import time
    tic = time.time()
    test.do_raytracing()
    toc = time.time()
    print(f"Raytracing took {toc - tic:.2f} seconds.")

    # load Classification npy file
    Classification = np.load(os.path.join(settings['out_dir'], 'Classification.npy'))

    # Load the expected Classification file from the repository
    Classification_expected = np.load(os.path.join(data_path, 'Output', 'Classification.npy'))

    np.testing.assert_array_equal(Classification, Classification_expected)

    # Normalize occlusion mapping output
    from occpy.util import normalize_occlusion_output
    Nhit_norm, Nmiss_norm, Nocc_norm, Classification_norm, chm = normalize_occlusion_output(
        input_folder=settings['out_dir'],
        PlotDim=settings['plot_dim'],
        vox_dim=settings['vox_dim'],
        dtm_file=settings['tif_in']['DTM'],
        dsm_file=settings['tif_in']['DSM'],
        lower_threshold=settings['lower_threshold'],
        output_voxels=settings['output_voxels'])

    # Load expected normalized Classification output
    Classification_norm_expected = np.load(os.path.join(data_path, 'Output', 'Classification_norm.npy'))

    np.testing.assert_array_equal(Classification_norm, Classification_norm_expected)
