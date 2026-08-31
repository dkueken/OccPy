# Changelog

All notable changes to this project are documented in this file.

Releases before 0.3.0 were pre-release candidates (`0.2.0rc1`-`0.2.0rc5`) and are not
documented here; see the git history for those.

## [0.3.0] - 2026-08-28

This release includes a major refactor of the main `OccPy`class, plus a
package-wide alignment of names with PEP 8. Several fixes and enhancements to usability are included.

### Added

- **Automatic return-mode detection.** `single_return` is now resolved from the data
  when it is not set in the config: the LAS header's `number_of_points_by_return` is
  used when it can be trusted, otherwise the return-number fields are scanned.
- **`check_returns_all_files` config option** (default `False`) to probe every file in a
  directory for its return mode rather than only the first.
- **`run_summary.json`**, written to the output directory alongside `config.json`,
  recording resolved return mode, trace mode, elapsed time, and the traversal counters.
- **Python-side traversal report.** `OccPy.get_raytracing_report()` now formats and
  returns the traversal counters itself, and accepts an optional `elapsed_seconds`. The
  previous report was printed from C++ and was lost in notebooks and under output
  capture.
- **`occpy.pulse_util` module** with following classes and helpers:
  `ReturnMode`, `TraceMode`, `ReturnModeResult`, `ScanJob`, `Chunk`,
  `detect_return_mode`, `return_mode_from_header`, `scan_return_mode`,
  `normalise_chunk`, `chunk_sortedness` and `constant_sensor_positions`.
- **`occpy.terrainmodel.nhit_to_dtm` and `nhit_to_dsm`**: derive a DTM or
  canopy-surface GeoTIFF directly from an `Nhit` grid, with optional Gaussian smoothing.
- **`occpy.util.normalize_occlusion_grids`**: height-normalize grids already    held in memory. `normalize_occlusion_output` now delegates to it.
- **`occpy.visualization.plot_terrain_models`** -- side-by-side DTM, DSM and derived CHM
  debug figure.
- **Figure-property defaults and validation**: `check_fig_prop` plus the
  `TRANSECT_FIG_PROP_DEFAULTS`, `PROFILE_FIG_PROP_DEFAULTS` and
  `TERRAIN_FIG_PROP_DEFAULTS` dictionaries.
- **Progress bar for TLS directory runs**, which previously had none, and a
  `TqdmLoggingHandler` so log records no longer break a live progress bar.
- **Tests**: `tests/test_synthetic_tracing_modes.py` (streaming/deferred equivalence,
  the double-count regression, and mobile trajectory interpolation, all against
  analytically derived grids) and `tests/test_terrain_model.py`.

### Changed

- **`OccPy` main class refactored as a pipeline.** All supported input shapes are now normalized into one list of scan jobs and processed by a single
  pipeline , replacing a lot of near-duplicate branches and overlapping code blocks.
- **`single_return` defaults to auto-detect** instead of `False`.
- **Peak memory for unsorted multi-scan TLS is bounded to a single scan.** In deferred
  mode the accumulated pulse dataset is now flushed at each scan boundary instead of
  once at the very end.
- **GPS-time sortedness is checked on every chunk** and across chunk boundaries, not only on the first chunk of each file.
- **`.las` files are accepted everywhere `.laz` files are**, including directory input, and extension matching is case-insensitive.
- **Unmatched scan IDs no longer abort the run.** `link_scan_positions` warns and skips a LAZ file with no matching sensor position.

### Changed -- naming (PEP 8)

Module, function and method names were brought in line with
[PEP 8](https://peps.python.org/pep-0008/#naming-conventions).

Modules:

| Old | New |
| --- | --- |
| `occpy/OccPy.py` | `occpy/occpy.py` |
| `occpy/TerrainModel.py` | `occpy/terrainmodel.py` |
| `occpy/OccPyRIEGL.py` | `occpy/occpyRIEGL.py` |


`occpy.visualization`:

| Old | New |
| --- | --- |
| `get_Occl_TransectFigure`, `get_Occl_TransectFigure_BinaryOcclusion` | `plot_transect` (merged, see below) |
| `get_Occlusion_ProfileFigure` | `plot_profile` |
| `vis_pv_static_bounds` | `plot_pyvista_static` |
| `vis_pv_rotating` | `plot_pyvista_rotating` |
| `vis_pv_interactive` | `plot_pyvista_interactive` |

`OccPy` methods:

| Old | New |
| --- | --- |
| `getGridDimensions` | `get_grid_dimensions` |
| `getGridOrigin` | `get_grid_origin` |
| `getNumTraversedPulses` | `get_num_traversed_pulses` |
| `getTotalNumPulses` | `get_total_num_pulses` |
| `getNumRegisteredHits` | `get_num_registered_hits` |
| `getNumEchoesOutside` | `get_num_echoes_outside` |
| `getNumMissingReturns` | `get_num_missing_returns` |
| `getNumNonIntersectPulses` | `get_num_non_intersecting_pulses` |
| `clean_up_RayTr` | `clean_up_raytracer` |
| `define_sensor_pos_singlePos` | `define_sensor_pos_single` |
| `link_positions_to_laz_files` | `link_scan_positions` |

Other modules:

| Old | New |
| --- | --- |
| `occpy.util.filterPointsIntersectingBox` | `filter_points_intersecting_box` |
| `occpy.util.array3Dto2D` | `array_3d_to_2d` |
| `occpy.terrainmodel.get_dtm_from_nhit_grid` | `nhit_to_dtm` |
| `occpy.terrainmodel.get_dsm_from_nhit_grid` | `nhit_to_dsm` |
| `OccPyRIEGL.test_colinearity` | `check_collinearity` |

- `get_Occl_TransectFigure` and `get_Occl_TransectFigure_BinaryOcclusion` were merged
  into `plot_transect(..., mode=...)`.
  `mode="binary"` (the default) reproduces the former
  `get_Occl_TransectFigure_BinaryOcclusion`; `mode="fraction"` reproduces
  `get_Occl_TransectFigure`, which additionally takes the `OcclFrac` grid.

[0.3.0]: https://github.com/dkueken/OccPy/releases/tag/v0.3.0
