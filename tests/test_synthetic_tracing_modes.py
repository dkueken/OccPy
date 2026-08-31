from __future__ import annotations

from pathlib import Path

import laspy
import numpy as np
from occpy.occpy import OccPy
from occpy.pulse_util import TraceMode


# ---------------------------------------------------------------------------
# Shared synthetic-data helpers
#
# All tests below use purely vertical pulses (sensor and every echo of a pulse
# share x/y) above a 10x10x10, vox_dim=1 grid over [0,10]x[0,10]x[0,10]. This
# keeps the expected traversal analytically obvious and independent of the
# ray-tracer implementation: for one column (ix, iy) and one pulse with
# echoes at increasing voxel-z indices, everything before the first echo is a
# miss, every echo voxel is a hit, gaps between consecutive echoes are misses
# (the pulse passed through empty space to reach the next echo), and
# everything after the last echo is occluded.
# ---------------------------------------------------------------------------


def _expected_grids_for_pulses(pulses, grid_shape):
    """Analytically derive Nhit/Nmiss/Nocc for a set of purely vertical pulses.

    Parameters
    ----------
    pulses : list of (ix, iy, iz_list)
        One entry per pulse: the voxel column it is fired into, and the voxel-z
        indices of its echoes in travel order (nearest echo first).
    grid_shape : tuple of int
    """
    nx, ny, nz = grid_shape
    nhit = np.zeros(grid_shape, dtype=np.int32)
    nmiss = np.zeros(grid_shape, dtype=np.int32)
    nocc = np.zeros(grid_shape, dtype=np.int32)

    for ix, iy, iz_list in pulses:
        iz_sorted = sorted(iz_list)
        nmiss[ix, iy, : iz_sorted[0]] += 1
        for k, iz in enumerate(iz_sorted):
            nhit[ix, iy, iz] += 1
            if k < len(iz_sorted) - 1:
                nxt = iz_sorted[k + 1]
                if nxt > iz + 1:
                    nmiss[ix, iy, iz + 1 : nxt] += 1
        nocc[ix, iy, iz_sorted[-1] + 1 : nz] += 1

    return nhit, nmiss, nocc


def _write_laz(path: Path, rows: list[dict]) -> None:
    """Write a minimal LAZ file from a list of per-echo row dicts."""
    header = laspy.LasHeader(point_format=3, version="1.2")
    header.offsets = np.array([0.0, 0.0, 0.0])
    header.scales = np.array([0.001, 0.001, 0.001])

    las = laspy.LasData(header)
    las.x = np.array([r["x"] for r in rows], dtype=np.float64)
    las.y = np.array([r["y"] for r in rows], dtype=np.float64)
    las.z = np.array([r["z"] for r in rows], dtype=np.float64)
    las.gps_time = np.array([r["gps_time"] for r in rows], dtype=np.float64)
    las.return_number = np.array([r["return_number"] for r in rows], dtype=np.uint8)
    las.number_of_returns = np.array([r["number_of_returns"] for r in rows], dtype=np.uint8)
    las.write(str(path))


def _rows_for_pulse(ix, iy, iz_list, gps, declared_n=None, return_numbers=None):
    """Rows for one vertical pulse in voxel column (ix, iy), echoes at iz_list."""
    x = ix + 0.5
    y = iy + 0.5
    iz_sorted = sorted(iz_list)
    n = declared_n if declared_n is not None else len(iz_sorted)
    rns = return_numbers if return_numbers is not None else list(range(1, len(iz_sorted) + 1))
    return [
        dict(x=x, y=y, z=iz + 0.5, gps_time=gps, return_number=rn, number_of_returns=n)
        for iz, rn in zip(iz_sorted, rns)
    ]


def _load_grids(out_dir):
    out_dir = Path(out_dir)
    return (
        np.load(out_dir / "Nhit.npy"),
        np.load(out_dir / "Nmiss.npy"),
        np.load(out_dir / "Nocc.npy"),
    )


_GRID_SHAPE = (10, 10, 10)
_PLOT_DIM = [0, 0, 0, 10, 10, 10]


def _run_occpy(laz_path, out_dir, extra_cfg=None):
    cfg = dict(
        laz_in=str(laz_path),
        out_dir=str(out_dir),
        vox_dim=1.0,
        lower_threshold=0,
        points_per_iter=1000,
        plot_dim=_PLOT_DIM,
        output_voxels=False,
        is_mobile=False,
        single_return=False,
    )
    if extra_cfg:
        cfg.update(extra_cfg)
    occpy = OccPy(config=cfg)
    occpy.define_sensor_pos_single(scan_pos_id=1, x=5.5, y=5.5, z=0.5)
    occpy.do_raytracing()
    return occpy


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_streaming_and_deferred_modes_agree_on_synthetic_multi_return_data(tmp_path: Path) -> None:
    """STREAMING (gps-sorted) and DEFERRED (unsorted) tracing of the same
    multi-return pulses must produce identical, analytically-correct grids.

    This targets the mode-selection/switching logic introduced by the
    tracing-mode refactor: whichever path is taken, the physical result must
    not depend on input ordering.
    """
    pulses = [
        (5, 5, [1, 3]),  # gap between the two echoes
        (5, 5, [4, 5]),  # adjacent echoes, no gap
        (5, 5, [6, 8]),  # gap between the two echoes
    ]
    expected_nhit, expected_nmiss, expected_nocc = _expected_grids_for_pulses(pulses, _GRID_SHAPE)

    gps = {0: 10.0, 1: 20.0, 2: 30.0}
    pulse_rows = [_rows_for_pulse(ix, iy, iz_list, gps[i]) for i, (ix, iy, iz_list) in enumerate(pulses)]

    sorted_rows = [row for rows in pulse_rows for row in rows]  # already gps-ascending
    shuffled_rows = [row for rows in (pulse_rows[2], pulse_rows[0], pulse_rows[1]) for row in rows]

    sorted_laz = tmp_path / "sorted.laz"
    shuffled_laz = tmp_path / "shuffled.laz"
    _write_laz(sorted_laz, sorted_rows)
    _write_laz(shuffled_laz, shuffled_rows)

    sorted_run = _run_occpy(sorted_laz, tmp_path / "sorted_out")
    shuffled_run = _run_occpy(shuffled_laz, tmp_path / "shuffled_out")

    assert sorted_run.trace_mode is TraceMode.STREAMING
    assert shuffled_run.trace_mode is TraceMode.DEFERRED

    for run in (sorted_run, shuffled_run):
        actual_nhit, actual_nmiss, actual_nocc = _load_grids(run.out_dir)
        np.testing.assert_array_equal(actual_nhit, expected_nhit)
        np.testing.assert_array_equal(actual_nmiss, expected_nmiss)
        np.testing.assert_array_equal(actual_nocc, expected_nocc)


def test_cleanup_incomplete_pulses_traces_each_pulse_exactly_once(tmp_path: Path) -> None:
    """Regression test for the bug fixed in commit 92036ff: the deferred trace
    and the cleaned-up-incomplete-pulses trace must each clear the pulse
    dataset afterwards, or pulses traced in the first pass get traced again
    (and double-counted) in the second.

    The synthetic file has three complete multi-return pulses and one
    incomplete pulse (declared 3 returns, only 1st and 3rd present), written
    out of gps_time order to force deferred tracing.
    """
    complete_pulses = [
        (5, 5, [1, 3]),
        (5, 5, [4, 5]),
        (5, 5, [6, 8]),
    ]
    incomplete_pulse = (5, 5, [2, 7])  # after cleanup: echoes renumbered to 1 and 2

    gps = {0: 10.0, 1: 20.0, 2: 30.0}
    complete_rows = [
        _rows_for_pulse(ix, iy, iz_list, gps[i]) for i, (ix, iy, iz_list) in enumerate(complete_pulses)
    ]
    # number_of_returns says 3, but only returns 1 and 3 are present (return 2 missing).
    incomplete_rows = _rows_for_pulse(
        *incomplete_pulse, gps=40.0, declared_n=3, return_numbers=[1, 3]
    )

    # Non-monotonic gps_time forces deferred tracing for the whole file.
    rows = complete_rows[2] + complete_rows[0] + incomplete_rows + complete_rows[1]

    laz_path = tmp_path / "unsorted_with_incomplete.laz"
    _write_laz(laz_path, rows)

    run_without_cleanup = _run_occpy(
        laz_path, tmp_path / "no_cleanup_out", extra_cfg={"cleanup_incomplete_pulses": False}
    )
    run_with_cleanup = _run_occpy(
        laz_path, tmp_path / "with_cleanup_out", extra_cfg={"cleanup_incomplete_pulses": True}
    )

    assert run_without_cleanup.trace_mode is TraceMode.DEFERRED
    assert run_with_cleanup.trace_mode is TraceMode.DEFERRED

    # Without cleanup, the incomplete pulse is silently dropped.
    expected_no_cleanup = _expected_grids_for_pulses(complete_pulses, _GRID_SHAPE)
    # With cleanup, the incomplete pulse is repaired and traced too -- but each
    # pulse (including the three already traced in the deferred pass) must be
    # counted exactly once.
    expected_with_cleanup = _expected_grids_for_pulses(
        complete_pulses + [incomplete_pulse], _GRID_SHAPE
    )

    actual_no_cleanup = _load_grids(run_without_cleanup.out_dir)
    actual_with_cleanup = _load_grids(run_with_cleanup.out_dir)

    for actual, expected in zip(actual_no_cleanup, expected_no_cleanup):
        np.testing.assert_array_equal(actual, expected)
    for actual, expected in zip(actual_with_cleanup, expected_with_cleanup):
        np.testing.assert_array_equal(actual, expected)


def test_mobile_trajectory_interpolation_matches_analytic_expectation(tmp_path: Path) -> None:
    """A moving sensor with a known linear trajectory should place each
    single-return pulse's origin exactly where linear interpolation predicts,
    which fully determines the expected occlusion grids.
    """
    traj_csv = tmp_path / "trajectory.csv"
    traj_csv.write_text("time,x,y,z\n0.0,0.5,5.5,0.5\n7.0,7.5,5.5,0.5\n")

    n_pulses = 8
    times = np.arange(n_pulses, dtype=np.float64)  # 0..7
    sensor_x = 0.5 + times  # linear interpolation between the two waypoints above

    pulses = []
    rows = []
    for i, t in enumerate(times):
        ix = int(np.floor(sensor_x[i]))
        iy = 5
        iz = 1 + i
        pulses.append((ix, iy, [iz]))
        rows.append(
            dict(
                x=sensor_x[i], y=5.5, z=iz + 0.5, gps_time=t,
                return_number=1, number_of_returns=1,
            )
        )

    expected_nhit, expected_nmiss, expected_nocc = _expected_grids_for_pulses(pulses, _GRID_SHAPE)

    laz_path = tmp_path / "mobile_synthetic.laz"
    _write_laz(laz_path, rows)

    cfg = dict(
        laz_in=str(laz_path),
        out_dir=str(tmp_path / "mobile_out"),
        vox_dim=1.0,
        lower_threshold=0,
        points_per_iter=1000,
        plot_dim=_PLOT_DIM,
        output_voxels=False,
        is_mobile=True,
        single_return=True,
    )
    occpy = OccPy(config=cfg)
    occpy.define_sensor_pos(
        path2file=str(traj_csv), delimiter=",",
        hdr_time="time", hdr_x="x", hdr_y="y", hdr_z="z",
    )
    occpy.do_raytracing()

    actual_nhit, actual_nmiss, actual_nocc = _load_grids(occpy.out_dir)
    np.testing.assert_array_equal(actual_nhit, expected_nhit)
    np.testing.assert_array_equal(actual_nmiss, expected_nmiss)
    np.testing.assert_array_equal(actual_nocc, expected_nocc)
