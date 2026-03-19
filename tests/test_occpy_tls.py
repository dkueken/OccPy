import numpy as np
import json
import os
from pathlib import Path
import pooch

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

def test_occpy_run_tls():
    from occpy.OccPy import OccPy

    # Download test data
    p = pooch.create(
        path=pooch.os_cache("occpy_test_data"),
        base_url="https://zenodo.org/records/18095821/files/",
        registry={"TLS_demo.zip": "md5:ae9a8de9ff7110595fffe613b2dc83f2"},
    )
    p.fetch("TLS_demo.zip", processor=pooch.Unzip(members=["TLS_demo"]), progressbar=True)
    data_path = os.path.join(p.path, "TLS_demo.zip.unzip", "TLS_demo")

    occpy_root = find_project_root()

    # Load the json settings file
    with open(os.path.join(occpy_root, 'config', 'settings_TLS_tutorial.JSON'), "r") as file:
        settings_orig = json.load(file)

    # copy original settings into new variable
    settings = settings_orig.copy()
    settings['root_folder'] = os.path.abspath(occpy_root)
    settings['laz_in'] = os.path.join(data_path, 'LAZ')
    settings['tif_in']['DTM'] = os.path.join(data_path, 'Grids', 'Ramerenwald_DTM_20250305.tif')
    settings['tif_in']['DSM'] = os.path.join(data_path, 'Grids', 'Ramerenwald_DSM_20250305.tif')
    settings['ScanPos'] = os.path.join(data_path, 'ScanPos', 'ScanPositions.txt')
    settings['out_dir'] = os.path.join(settings['root_folder'], 'output', 'TLS')

    # Initialize OccPy object from config dict
    test = OccPy(config=settings)

    # Define sensor positions
    test.define_sensor_pos(path2file=settings['ScanPos'],
                           delimiter=",",
                           hdr_scanpos_id='ID',
                           hdr_x='X',
                           hdr_y='Y',
                           hdr_z='Z',
                           sens_pos_id_offset=0)

    # Run occpy
    import time
    tic = time.time()
    test.do_raytracing()
    toc = time.time()

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