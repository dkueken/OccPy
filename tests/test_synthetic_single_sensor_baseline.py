from __future__ import annotations
from pathlib import Path

import laspy
import numpy as np
from occpy.OccPy import OccPy


def _create_synthetic_single_sensor_laz(laz_path: Path) -> dict[str, np.ndarray]:
    """Create a small deterministic LAZ with 100 single-return pulses.

    All returns are aligned above a single sensor position, which makes the
    expected traversal easy to derive analytically.
    """

    n_pulses = 100

    # Single sensor / single column setup.
    sensor_x = 5.5
    sensor_y = 5.5
    sensor_z = 0.5

    x = np.full(n_pulses, sensor_x, dtype=np.float64)
    y = np.full(n_pulses, sensor_y, dtype=np.float64)

    # Keep returns safely inside the integer bbox [0, 10] to avoid boundary ties.
    # This creates z voxel indices 1..8 (at vox_dim=1), repeated deterministically.
    z = 1.5 + (np.arange(n_pulses, dtype=np.float64) % 8)

    gps_time = np.arange(1, n_pulses + 1, dtype=np.float64)
    return_number = np.ones(n_pulses, dtype=np.uint8)
    number_of_returns = np.ones(n_pulses, dtype=np.uint8)

    header = laspy.LasHeader(point_format=3, version="1.2")
    header.offsets = np.array([0.0, 0.0, 0.0])
    header.scales = np.array([0.001, 0.001, 0.001])

    las = laspy.LasData(header)
    las.x = x
    las.y = y
    las.z = z
    las.gps_time = gps_time
    las.return_number = return_number
    las.number_of_returns = number_of_returns
    las.write(laz_path)

    return {
        "sensor": np.array([sensor_x, sensor_y, sensor_z], dtype=np.float64),
        "x": x,
        "y": y,
        "z": z,
    }


def _expected_grids_for_vertical_single_return_case(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    grid_shape: tuple[int, int, int],
    min_bound: tuple[float, float, float],
    vox_dim: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute expected Nhit/Nmiss/Nocc/Classification for this synthetic setup."""
    
    nx, ny, nz = grid_shape
    min_x, min_y, min_z = min_bound

    nhit = np.zeros(grid_shape, dtype=np.int32)
    nmiss = np.zeros(grid_shape, dtype=np.int32)
    nocc = np.zeros(grid_shape, dtype=np.int32)

    for px, py, pz in zip(x, y, z):
        ix = int(np.floor((px - min_x) / vox_dim))
        iy = int(np.floor((py - min_y) / vox_dim))
        iz = int(np.floor((pz - min_z) / vox_dim))

        # Traversal for this synthetic case is purely vertical in one column:
        # misses until return voxel, hit at return voxel, occluded above.
        nmiss[ix, iy, :iz] += 1
        nhit[ix, iy, iz] += 1
        nocc[ix, iy, iz + 1 : nz] += 1

    classification = np.zeros(grid_shape, dtype=np.int32)
    classification[np.logical_and.reduce((nhit > 0, nmiss >= 0, nocc >= 0))] = 1
    classification[np.logical_and.reduce((nhit == 0, nmiss > 0, nocc >= 0))] = 2
    classification[np.logical_and.reduce((nhit == 0, nmiss == 0, nocc > 0))] = 3
    classification[np.logical_and.reduce((nhit == 0, nmiss == 0, nocc == 0))] = 4

    return nhit, nmiss, nocc, classification


def test_occpy_single_sensor_synthetic_baseline(tmp_path: Path) -> None:
    """Regression baseline for a tiny deterministic OccPy run.

    This test intentionally uses integer plot bounds to avoid bbox float->int
    edge effects. It can be reused later to validate bbox rounding fixes.
    """
    laz_path = tmp_path / "synthetic_single_sensor.laz"
    out_dir = tmp_path / "occpy_out"
    expected_dir = tmp_path / "expected_out"
    expected_dir.mkdir(parents=True, exist_ok=True)

    synthetic = _create_synthetic_single_sensor_laz(laz_path)

    plot_dim = [0, 0, 0, 10, 10, 10]  # integer bbox by design
    vox_dim = 1.0

    cfg = {
        "laz_in": str(laz_path),
        "out_dir": str(out_dir),
        "vox_dim": vox_dim,
        "lower_threshold": 0,
        "points_per_iter": 1000,
        "plot_dim": plot_dim,
        "output_voxels": False,
        "is_mobile": False,
        "single_return": True,
    }
    occpy = OccPy(config=cfg)
    occpy.define_sensor_pos_singlePos(
        scan_pos_id=1,
        x=float(synthetic["sensor"][0]),
        y=float(synthetic["sensor"][1]),
        z=float(synthetic["sensor"][2]),
    )
    occpy.do_raytracing()

    actual_nhit = np.load(out_dir / "Nhit.npy")
    actual_nmiss = np.load(out_dir / "Nmiss.npy")
    actual_nocc = np.load(out_dir / "Nocc.npy")
    actual_classification = np.load(out_dir / "Classification.npy")

    expected_nhit, expected_nmiss, expected_nocc, expected_classification = (
        _expected_grids_for_vertical_single_return_case(
            x=synthetic["x"],
            y=synthetic["y"],
            z=synthetic["z"],
            grid_shape=(10, 10, 10),
            min_bound=(0.0, 0.0, 0.0),
            vox_dim=vox_dim,
        )
    )

    # Save expected arrays to make the baseline directly inspectable.
    np.save(expected_dir / "expected_Nhit.npy", expected_nhit)
    np.save(expected_dir / "expected_Nmiss.npy", expected_nmiss)
    np.save(expected_dir / "expected_Nocc.npy", expected_nocc)
    np.save(expected_dir / "expected_Classification.npy", expected_classification)

    np.testing.assert_array_equal(actual_nhit, expected_nhit)
    np.testing.assert_array_equal(actual_nmiss, expected_nmiss)
    np.testing.assert_array_equal(actual_nocc, expected_nocc)
    np.testing.assert_array_equal(actual_classification, expected_classification)
