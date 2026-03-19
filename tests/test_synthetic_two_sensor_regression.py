from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pytest

from occpy.OccPy import OccPy


FIXTURE_ROOT = Path(__file__).resolve().parent / "fixtures" / "two_sensor_case"


def _required_expected_files(expected_dir: Path) -> list[Path]:
    return [
        expected_dir / "Nhit.npy",
        expected_dir / "Nmiss.npy",
        expected_dir / "Nocc.npy",
        expected_dir / "Classification.npy",
    ]

def test_occpy_two_sensor_regression_against_fixture_expected(tmp_path: Path) -> None:
    """Run deterministic 2-sensor synthetic input and compare to expected grids."""
    cfg_path = FIXTURE_ROOT / "case_config.json"
    if not cfg_path.exists():
        pytest.skip("Fixture config missing. Run generator script first.")

    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))

    input_laz_dir = FIXTURE_ROOT / cfg["laz_dir"]
    senspos_file = FIXTURE_ROOT / cfg["sensor_positions"]
    expected_dir = FIXTURE_ROOT / cfg["expected_dir"]

    missing_input = [p for p in [input_laz_dir, senspos_file] if not p.exists()]
    if missing_input:
        pytest.skip(
            "Fixture input missing. Run tests/fixtures/two_sensor_case/generate_two_sensor_case.py first."
        )

    required_expected = _required_expected_files(expected_dir)
    missing_expected = [p for p in required_expected if not p.exists()]
    if missing_expected:
        pytest.skip(
            "Expected baseline arrays are missing in tests/fixtures/two_sensor_case/expected/. "
            "Generate and commit Nhit/Nmiss/Nocc/Classification first."
        )

    out_dir = tmp_path / "occpy_out"

    occpy_cfg = {
        "laz_in": str(input_laz_dir),
        "out_dir": str(out_dir),
        "vox_dim": float(cfg["vox_dim"]),
        "lower_threshold": 0,
        "points_per_iter": 100000,
        "plot_dim": cfg["plot_dim"],
        "output_voxels": False,
        "is_mobile": False,
        "single_return": True,
        "str_idxs_ScanPosID": [
            int(cfg["scan_pos_id_stridx"]),
            int(cfg["scan_pos_id_endstridx"]),
        ],
    }
    occpy = OccPy(config=occpy_cfg)

    occpy.define_sensor_pos(
        path2file=str(senspos_file),
        delimiter=" ",
        hdr_scanpos_id="ScanPos",
        hdr_x="x",
        hdr_y="y",
        hdr_z="z",
    )

    occpy.do_raytracing()

    actual_nhit = np.load(out_dir / "Nhit.npy")
    actual_nmiss = np.load(out_dir / "Nmiss.npy")
    actual_nocc = np.load(out_dir / "Nocc.npy")
    actual_classification = np.load(out_dir / "Classification.npy")

    expected_nhit = np.load(expected_dir / "Nhit.npy")
    expected_nmiss = np.load(expected_dir / "Nmiss.npy")
    expected_nocc = np.load(expected_dir / "Nocc.npy")
    expected_classification = np.load(expected_dir / "Classification.npy")

    np.testing.assert_array_equal(actual_nhit, expected_nhit)
    np.testing.assert_array_equal(actual_nmiss, expected_nmiss)
    np.testing.assert_array_equal(actual_nocc, expected_nocc)
    np.testing.assert_array_equal(actual_classification, expected_classification)
