import numpy as np

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



