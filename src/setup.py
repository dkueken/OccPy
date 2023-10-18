# Cython compile instructions
# Use python setup.py build_ext --inplace
# to compile

import sys
import os
import shutil

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

# clean previous build
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if (name.startswith("raytr") and not(name.endswith(".pyx") or name.endswith(".pxd"))):
            os.remove(os.path.join(root, name))
    #for name in dirs:
    #    if (name == "build"):
    #        shutil.rmtree(name)

# build "raytr.so" python extension to be added to "PYTHONPATH" afterwards...
extensions = [
    Extension("raytr",
                  sources=["raytr.pyx",
                           "Raytracer.cpp",
                           "Pulse.cpp",
                           "Echo.cpp"
                       ],
                  libraries=[],  # refers to "liblas.2.3.0.dylib"
                  language="c++",  # remove this if C and not C++
                  include_dirs=["C:/Users/kueken/Miniconda3/envs/occPy/Library/include/"], #["C:/Users/kueken/Miniconda3/envs/OcclusionMapping_PDAL/Library/include/"],
                  #include_dirs=["/usr/local/Cellar/boost/1.75.0_2/include", "/usr/local/Cellar/pdal/2.2.0_3/include", "/usr/local/Cellar/laszip/3.4.3/include"],
                  library_dirs=["C:/Users/kueken/Miniconda3/envs/occPy/Library/lib"], #["C:/Users/kueken/Miniconda3/envs/OcclusionMapping_PDAL/Library/lib"],
                  #library_dirs=["/usr/local/Cellar/pdal/2.2.0_3/lib"],
                  requires=['Cython'],
                  #extra_compile_args=["-ferror-limit=0"]
                  #extra_compile_args=["-std=c++11"],
                  #requires=['Cython'],
                  #extra_compile_args=["-fopenmp", "-O3"],
                  #extra_link_args=["L.C:/Users/kueken/Miniconda3/envs/OcclusionMapping_PDAL/Library/lib/pdalcpp"]
                  #, "-DSOME_DEFINE_OPT","-L./some/extra/dependency/dir/"]
                  )
]
setup(
    name = 'raytr',
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize(extensions),
    install_requires=['Cython', 'numpy']
)
