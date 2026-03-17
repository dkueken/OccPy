from __future__ import annotations

from pathlib import Path
import json

import laspy
import numpy as np


def _build_pulse_endpoints(sensor_xyz: np.ndarray) -> np.ndarray:
    """Deterministic set of pulse endpoints covering full azimuth directions."""
    azimuth_deg = np.arange(0, 360, 15, dtype=np.float64)  # 24 directions around
    elevation_deg = np.array([-25.0, -10.0, 10.0, 25.0], dtype=np.float64)

    azimuth = np.deg2rad(azimuth_deg)
    elevation = np.deg2rad(elevation_deg)

    points = []
    k = 0
    for elev in elevation:
        cos_e = np.cos(elev)
        sin_e = np.sin(elev)
        for az in azimuth:
            # Vary distance deterministically to avoid repeated points.
            distance = 4.5 + 0.25 * (k % 6)
            direction = np.array([
                cos_e * np.cos(az),
                cos_e * np.sin(az),
                sin_e,
            ])
            points.append(sensor_xyz + distance * direction)
            k += 1

    return np.array(points, dtype=np.float64)


def _write_laz(points_xyz: np.ndarray, gps_time: np.ndarray, out_file: Path) -> None:
    n = points_xyz.shape[0]
    header = laspy.LasHeader(point_format=3, version="1.2")
    header.offsets = np.array([0.0, 0.0, 0.0])
    header.scales = np.array([0.001, 0.001, 0.001])

    las = laspy.LasData(header)
    las.x = points_xyz[:, 0]
    las.y = points_xyz[:, 1]
    las.z = points_xyz[:, 2]
    las.gps_time = gps_time
    las.return_number = np.ones(n, dtype=np.uint8)
    las.number_of_returns = np.ones(n, dtype=np.uint8)

    las.write(out_file)


def generate_case(base_dir: Path) -> None:
    """Generate deterministic 2-sensor synthetic test input files."""
    input_dir = base_dir / "input"
    laz_dir = input_dir / "laz"
    laz_dir.mkdir(parents=True, exist_ok=True)

    # Integer-bounded grid setup (intentionally avoids bbox float edge behavior).
    plot_dim = [0, 0, 0, 20, 20, 20]
    vox_dim = 1.0

    sensors = {
        1: np.array([6.0, 6.0, 6.0], dtype=np.float64),
        2: np.array([14.0, 14.0, 6.0], dtype=np.float64),
    }

    # Build and write two LAZ files, one per static sensor position.
    for scan_id, sensor_xyz in sensors.items():
        points = _build_pulse_endpoints(sensor_xyz)

        # Keep points away from exact grid faces for deterministic voxel indexing.
        points = np.clip(points, 0.5, 19.5)

        if scan_id == 1:
            gps = np.arange(1, points.shape[0] + 1, dtype=np.float64)
        else:
            gps = np.arange(10001, 10001 + points.shape[0], dtype=np.float64)

        out_laz = laz_dir / f"{scan_id:03d}_synthetic.laz"
        _write_laz(points, gps, out_laz)

    # Sensor positions table used by OccPy.define_sensor_pos(..., is_mobile=False).
    sensor_table = np.array([
        [1, sensors[1][0], sensors[1][1], sensors[1][2]],
        [2, sensors[2][0], sensors[2][1], sensors[2][2]],
    ])
    np.savetxt(
        input_dir / "sensor_positions.txt",
        sensor_table,
        fmt=["%03d", "%.3f", "%.3f", "%.3f"],
        delimiter=" ",
        header="ScanPos x y z",
        comments="",
    )

    config = {
        "plot_dim": plot_dim,
        "vox_dim": vox_dim,
        "scan_pos_id_stridx": 0,
        "scan_pos_id_endstridx": 3,
        "laz_dir": "input/laz",
        "sensor_positions": "input/sensor_positions.txt",
        "expected_dir": "expected",
    }
    (base_dir / "case_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")


if __name__ == "__main__":
    fixture_root = Path(__file__).resolve().parent
    generate_case(fixture_root)
    print(f"Generated deterministic two-sensor case under: {fixture_root}")
