import numpy as np

def test_grid_definition():
    from occpy.OccPy import OccPy

    minX = 10
    minY = 10
    minZ = 10
    maxX = 100
    maxY = 100
    maxZ = 100

    vox_size = 1
    expected_extent = [int((maxX - minX)/vox_size),
                       int((maxY - minY)/vox_size),
                       int((maxZ - minZ)/vox_size)]

    test = OccPy(laz_in=f"dummy_path/to_dummy_file.laz",
                 out_dir=f"dummy_path/to_dummy_output_dir",
                 vox_dim=1,
                 lower_threshold=1,
                 points_per_iter=1,
                 plot_dim=[minX, minY, minZ, maxX, maxY, maxZ])

    extent = test.getGridDimensions()
    np.testing.assert_array_equal(extent, expected_extent)



