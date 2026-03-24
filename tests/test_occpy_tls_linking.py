import glob
import json
import os
from pathlib import Path

import pooch
import pytest


def find_project_root(markers=("README.md", ".gitignore", "environment.yml", "setup.py")):
    try:
        import ipynbname

        current = ipynbname.path().parent
    except Exception:
        current = Path(os.getcwd())

    for parent in [current, *current.parents]:
        if any((parent / m).exists() for m in markers):
            return parent

    return current


def get_tls_demo_data_path():
    p = pooch.create(
        path=pooch.os_cache("occpy_test_data"),
        base_url="https://zenodo.org/records/18095821/files/",
        registry={"TLS_demo.zip": "md5:ae9a8de9ff7110595fffe613b2dc83f2"},
    )
    p.fetch("TLS_demo.zip", processor=pooch.Unzip(members=["TLS_demo"]), progressbar=True)
    return os.path.join(p.path, "TLS_demo.zip.unzip", "TLS_demo")


def load_tls_settings(data_path, out_dir, str_idxs=None):
    occpy_root = find_project_root()
    with open(os.path.join(occpy_root, "config", "settings_TLS_tutorial.JSON"), "r") as file:
        settings = json.load(file)

    settings["root_folder"] = os.path.abspath(occpy_root)
    settings["laz_in"] = os.path.join(data_path, "LAZ")
    settings["tif_in"]["DTM"] = os.path.join(data_path, "Grids", "Ramerenwald_DTM_20250305.tif")
    settings["tif_in"]["DSM"] = os.path.join(data_path, "Grids", "Ramerenwald_DSM_20250305.tif")
    settings["ScanPos"] = os.path.join(data_path, "ScanPos", "ScanPositions.txt")
    settings["out_dir"] = str(out_dir)
    if str_idxs is not None:
        settings["str_idxs_ScanPosID"] = str_idxs

    return settings


def _init_tls_occpy(settings):
    from occpy.OccPy import OccPy

    test = OccPy(config=settings)
    test.define_sensor_pos(
        path2file=settings["ScanPos"],
        delimiter=",",
        hdr_scanpos_id="ID",
        hdr_x="X",
        hdr_y="Y",
        hdr_z="Z",
        sens_pos_id_offset=0,
    )
    return test


def test_link_positions_to_laz_files_tls_wrong_indices(tmp_path):
    data_path = get_tls_demo_data_path()
    settings = load_tls_settings(data_path=data_path, out_dir=tmp_path / "TLS_wrong", str_idxs=[7, 9])

    test = _init_tls_occpy(settings)

    with pytest.raises(ValueError, match="No sensor position found"):
        test.link_positions_to_laz_files()


def test_link_positions_to_laz_files_tls_correct_indices(tmp_path):
    data_path = get_tls_demo_data_path()
    settings = load_tls_settings(data_path=data_path, out_dir=tmp_path / "TLS_ok")

    test = _init_tls_occpy(settings)
    links = test.link_positions_to_laz_files()

    assert not links.empty
    assert {"laz_file", "scan_name", "scan_id", "sensor_x", "sensor_y", "sensor_z"}.issubset(links.columns)
    assert len(links) == len(glob.glob(os.path.join(settings["laz_in"], "*.laz")))
