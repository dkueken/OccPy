def test_import_occpy():
    import occpy
    assert occpy is not None

def test_imports_other():
    import raytr
    from occpy import OccPy
    from occpy import OccPyRIEGL
    from occpy import TerrainModel
    import occpy.visualization
    import occpy.util
    from raytr import PyRaytracer

    assert raytr is not None
    assert OccPy is not None
    assert OccPyRIEGL is not None
    assert TerrainModel is not None
    assert occpy.visualization is not None
    assert occpy.util is not None
    assert PyRaytracer is not None
