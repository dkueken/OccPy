import lexcube
import numpy as np


# define input
data = np.load(r"D:\_tmp_wdir\TestData_occpy\MLS\OcclusionMapping_OccPy_omp\Classification.npy")
w = lexcube.Cube3DWidget(data, cmap="prism", vmin=0, vmax=768)
w
