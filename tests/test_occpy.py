import numpy as np
import json
import os
from pathlib import Path
import pooch

def test_grid_definition():
    """

    Returns
    -------

    """
    from occpy.OccPy import OccPy

    minX = 10
    minY = 20
    minZ = 30
    maxX = 100
    maxY = 200
    maxZ = 300

    vox_size = 1
    expected_extent = [int((maxX - minX)/vox_size),
                       int((maxY - minY)/vox_size),
                       int((maxZ - minZ)/vox_size)]

    expected_origin = [minX, minY, minZ]

    test = OccPy(laz_in=f"dummy_path/to_dummy_file.laz",
                 out_dir=f"dummy_path/to_dummy_output_dir",
                 vox_dim=1,
                 lower_threshold=1,
                 points_per_iter=1,
                 plot_dim=[minX, minY, minZ, maxX, maxY, maxZ])

    extent = test.getGridDimensions()
    np.testing.assert_array_equal(extent, expected_extent)

    origin = test.getGridOrigin()
    np.testing.assert_array_equal(origin, expected_origin)


def find_project_root(markers=('README.md', '.gitignore', 'environment.yml', 'setup.py')):
    try:
        import ipynbname
        current = ipynbname.path().parent
    except Exception:
        current = Path(os.getcwd())

    for parent in [current, *current.parents]:
        if any((parent / m).exists() for m in markers):
            return parent

    return current  # fallback if nothing found

def test_occpy_run():
    from occpy.OccPy import OccPy

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

    # initiate Occpy object
    test = OccPy(laz_in=f"{settings['laz_in']}",
             out_dir=settings['out_dir'],
             vox_dim=settings['vox_dim'],
             lower_threshold=settings['lower_threshold'],
             points_per_iter=settings['points_per_iter'],
             plot_dim=settings['plot_dim'])

    # Define sensor position
    test.define_sensor_pos(path2file=settings['ScanPos'],  # Path to the trajectory file
                           is_mobile=True,  # whether acquisition is mobile. Always true for MLS or ULS
                           single_return=True,  # wheter the data is single or multi return
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




