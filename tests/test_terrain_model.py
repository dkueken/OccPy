from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pytest
import rasterio

from occpy.TerrainModel import (
    NODATA_VALUE,
    get_dsm_from_nhit_grid,
    get_dtm_from_nhit_grid,
)


def _make_synthetic_nhit_grid() -> tuple[np.ndarray, dict]:
    """3x3x6 Nhit grid: a flat ground return at z_idx=1 in every column, plus a
    single higher canopy return at z_idx=4 in the center column (x=1, y=1)."""

    nx, ny, nz = 3, 3, 6
    vox_dim = 1.0
    plot_dim = [0, 0, 0, nx * vox_dim, ny * vox_dim, nz * vox_dim]

    nhit = np.zeros((nx, ny, nz), dtype=np.int32)
    nhit[:, :, 1] = 1
    nhit[1, 1, 4] = 1

    config = {"plot_dim": plot_dim, "vox_dim": vox_dim}
    return nhit, config

def test_get_dtm_and_chm_from_nhit_grid_raster_properties(tmp_path: Path) -> None:
    nhit, config = _make_synthetic_nhit_grid()
    plot_dim = config["plot_dim"]
    vox_dim = config["vox_dim"]

    dtm_file = tmp_path / "dtm.tif"
    dsm_file = tmp_path / "dsm.tif"

    dtm = get_dtm_from_nhit_grid(nhit, config, str(dtm_file))
    dsm = get_dsm_from_nhit_grid(nhit, config, str(dsm_file))

    # Ground is flat at z_idx=1 -> voxel-center elevation 1.5 everywhere.
    np.testing.assert_allclose(dtm, 1.5)

    # Canopy top is at z_idx=1 everywhere except the center column (z_idx=4 -> 4.5).
    expected_dsm = np.full((3, 3), 1.5)
    expected_dsm[1, 1] = 4.5
    np.testing.assert_allclose(dsm, expected_dsm)

    for path in (dtm_file, dsm_file):
        with rasterio.open(path) as src:
            assert src.width == 3
            assert src.height == 3
            assert src.res == (vox_dim, vox_dim)
            assert src.nodata == NODATA_VALUE
            # north-up: origin at (minX, maxY), row 0 = maxY, col 0 = minX
            assert src.bounds == (plot_dim[0], plot_dim[1], plot_dim[3], plot_dim[4])
            assert src.transform.c == plot_dim[0]
            assert src.transform.f == plot_dim[4]

        # The center Nhit column is (x_idx=1, y_idx=1); with row 0 = maxY, that
        # lands on raster row (ny - 1 - y_idx) = 1, column x_idx = 1.
        with rasterio.open(dsm_file) as src:
            band = src.read(1)
            assert band[1, 1] == 4.5


def test_get_dtm_from_nhit_grid_fills_nodata(tmp_path: Path) -> None:
    nx, ny, nz = 2, 2, 4
    vox_dim = 1.0
    config = {"plot_dim": [0, 0, 0, nx * vox_dim, ny * vox_dim, nz * vox_dim], "vox_dim": vox_dim}

    nhit = np.zeros((nx, ny, nz), dtype=np.int32)
    nhit[0, 0, 1] = 1  # only one column has any hit at all

    dtm_file = tmp_path / "dtm_sparse.tif"
    dtm = get_dtm_from_nhit_grid(nhit, config, str(dtm_file))

    assert (dtm == NODATA_VALUE).sum() == 3
    assert (dtm == 1.5).sum() == 1

    with rasterio.open(dtm_file) as src:
        assert src.nodata == NODATA_VALUE
        band = src.read(1)
        assert (band == NODATA_VALUE).sum() == 3


def test_get_dtm_from_nhit_grid_smoothing_sigma_reduces_spike_noise(tmp_path: Path) -> None:
    """A flat ground plane with a single spurious low outlier in the center
    column should, once smoothed, land close to its (flat) neighbors instead
    of the raw, deep pit."""

    nx, ny, nz = 9, 9, 8
    vox_dim = 1.0
    config = {"plot_dim": [0, 0, 0, nx * vox_dim, ny * vox_dim, nz * vox_dim], "vox_dim": vox_dim}

    nhit = np.zeros((nx, ny, nz), dtype=np.int32)
    nhit[:, :, 2] = 1  # flat ground at z_idx=2 (elevation 2.5) everywhere
    nhit[4, 4, 2] = 0
    nhit[4, 4, 0] = 1  # spurious outlier hit near the very bottom in one column

    dtm_raw = get_dtm_from_nhit_grid(nhit, config, str(tmp_path / "dtm_raw.tif"))
    assert dtm_raw[4, 4] == 0.5  # unsmoothed: full depth of the spurious pit

    dtm_smooth = get_dtm_from_nhit_grid(
        nhit, config, str(tmp_path / "dtm_smooth.tif"), smoothing_sigma=1.0
    )

    # Smoothing should pull the outlier much closer to its flat (2.5m) neighbors...
    assert abs(dtm_smooth[4, 4] - 2.5) < abs(dtm_raw[4, 4] - 2.5)
    # ...without disturbing columns far away from the outlier.
    assert dtm_smooth[0, 0] == pytest.approx(2.5, abs=1e-6)


def test_dtm_and_dsm_from_nhit_grid_normalize_occlusion_output(tmp_path: Path) -> None:
    from occpy.util import normalize_occlusion_output

    nhit, config = _make_synthetic_nhit_grid()
    plot_dim = config["plot_dim"]
    vox_dim = config["vox_dim"]
    nx, ny, _ = nhit.shape

    nmiss = np.zeros_like(nhit)
    nmiss[1, 1, 2] = 1
    nmiss[1, 1, 3] = 1

    nocc = np.zeros_like(nhit)

    classification = np.zeros_like(nhit)
    classification[nhit > 0] = 1
    classification[(nhit == 0) & (nmiss > 0)] = 2
    classification[(nhit == 0) & (nmiss == 0) & (nocc == 0)] = 4

    input_folder = tmp_path / "raytracing_out"
    input_folder.mkdir()
    np.save(input_folder / "Nhit.npy", nhit)
    np.save(input_folder / "Nmiss.npy", nmiss)
    np.save(input_folder / "Nocc.npy", nocc)
    np.save(input_folder / "Classification.npy", classification)

    dtm_file = tmp_path / "dtm.tif"
    dsm_file = tmp_path / "dsm.tif"
    get_dtm_from_nhit_grid(nhit, config, str(dtm_file))
    get_dsm_from_nhit_grid(nhit, config, str(dsm_file))

    Nhit_norm, Nmiss_norm, Nocc_norm, Classification_norm, chm = normalize_occlusion_output(
        input_folder=str(input_folder),
        PlotDim=plot_dim,
        vox_dim=vox_dim,
        dtm_file=str(dtm_file),
        dsm_file=str(dsm_file),
    )

    # chm = dsm - dtm, indexed [y, x]: 3m of canopy above ground at the center
    # column (x=1, y=1), flat ground (0m) everywhere else.
    expected_chm = np.zeros((ny, nx))
    expected_chm[1, 1] = 3.0
    np.testing.assert_allclose(chm, expected_chm)

    # Center column's normalized profile starts at the ground voxel (hit) and
    # covers the two empty voxels below the canopy return.
    np.testing.assert_array_equal(Nhit_norm[1, 1, 0:3], [1, 0, 0])
    np.testing.assert_array_equal(Nmiss_norm[1, 1, 0:3], [0, 1, 1])
    np.testing.assert_array_equal(Classification_norm[1, 1, 0:3], [1, 2, 2])

    # Columns without canopy have dtm == dsm, so nothing to normalize.
    np.testing.assert_array_equal(Nhit_norm[0, 0, :], np.zeros(Nhit_norm.shape[2]))


def test_plot_terrain_models_runs_and_saves_figure(tmp_path: Path) -> None:
    from occpy.visualization import plot_terrain_models

    nhit, config = _make_synthetic_nhit_grid()
    plot_dim = config["plot_dim"]

    dtm = get_dtm_from_nhit_grid(nhit, config, str(tmp_path / "dtm.tif"))
    dsm = get_dsm_from_nhit_grid(nhit, config, str(tmp_path / "dsm.tif"))

    fig = plot_terrain_models(dtm, dsm, plot_dim, out_dir=str(tmp_path))

    # 3 map panels + 3 colorbars
    assert len(fig.axes) == 6
    titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
    assert titles == ["DTM", "DSM", "CHM (DSM - DTM)"]
    out_file = tmp_path / "TerrainModels_debug.png"
    assert out_file.exists()
