import numpy as np
from occpy.OccPy import OccPy

def test_grid_definition(tmp_path):
    """

    Returns
    -------

    """

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

    laz_file = tmp_path / "dummy_input.laz"
    laz_file.touch()

    cfg = {
        "laz_in": str(laz_file),
        "out_dir": str(tmp_path / "dummy_output"),
        "vox_dim": 1,
        "lower_threshold": 1,
        "points_per_iter": 1,
        "plot_dim": [minX, minY, minZ, maxX, maxY, maxZ],
    }
    test = OccPy(config=cfg)

    extent = test.getGridDimensions()
    np.testing.assert_array_equal(extent, expected_extent)

    origin = test.getGridOrigin()
    np.testing.assert_array_equal(origin, expected_origin)


def test_non_integer_bounds_preserve_origin_and_extend_max(tmp_path):
    laz_file = tmp_path / "dummy_input.laz"
    laz_file.touch()

    cfg = {
        "laz_in": str(laz_file),
        "out_dir": str(tmp_path / "dummy_output"),
        "vox_dim": 0.3,
        "lower_threshold": 1,
        "points_per_iter": 1,
        "plot_dim": [0.25, 1.75, -0.1, 1.06, 2.76, 0.91],
    }

    test = OccPy(config=cfg)

    # max bounds are extended to the next voxel edge from min bound
    np.testing.assert_allclose(test.plot_dim, [0.25, 1.75, -0.1, 1.15, 2.95, 1.1])

    # origin stays exactly at provided min bound
    origin = test.getGridOrigin()
    np.testing.assert_allclose(origin, [0.25, 1.75, -0.1])

    # extents now map to integer voxel counts
    extent = test.getGridDimensions()
    np.testing.assert_array_equal(extent, [3, 4, 4])


def test_non_integer_bounds_already_divisible_unchanged(tmp_path):
    laz_file = tmp_path / "dummy_input.laz"
    laz_file.touch()

    cfg = {
        "laz_in": str(laz_file),
        "out_dir": str(tmp_path / "dummy_output"),
        "vox_dim": 0.25,
        "lower_threshold": 1,
        "points_per_iter": 1,
        "plot_dim": [0.125, -0.75, 0.5, 1.125, 0.25, 1.5],
    }

    test = OccPy(config=cfg)

    # no extension needed, dimensions remain unchanged
    np.testing.assert_allclose(test.plot_dim, cfg["plot_dim"])

    origin = test.getGridOrigin()
    np.testing.assert_allclose(origin, [0.125, -0.75, 0.5])

    extent = test.getGridDimensions()
    np.testing.assert_array_equal(extent, [4, 4, 4])

