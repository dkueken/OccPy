# Cython compile instructions
# Use python setup.py build_ext --inplace
# to compile

import sys
import os
import numpy
import platform
import shutil
import ctypes
import occpy

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

# clean previous build
for root, dirs, files in os.walk("", topdown=False):
    for name in files:
        if (name.startswith("raytr") and not(name.endswith(".pyx") or name.endswith(".pxd"))):
            os.remove(os.path.join(root, name))
    #for name in dirs:
    #    if (name == "build"):
    #        shutil.rmtree(name)


# Removes hardcoding by using Environment path to form the relevant directory paths. 
# exact paths platform-dependent, see below
env_path = sys.prefix

src_path = os.path.abspath("src")
sources=[os.path.join(src_path, "raytr.pyx"),
                        os.path.join(src_path, "Raytracer.cpp"),
                        os.path.join(src_path, "Pulse.cpp"),
                        os.path.join(src_path, "Echo.cpp")
                        ]

# build "raytr.so" python extension to be added to "PYTHONPATH" afterwards...
if platform.system() == 'Linux':
    include_path = os.path.join(env_path, "include")
    library_path = os.path.join(env_path, "lib")

    extensions = [
        Extension("raytr",
                sources=sources,
                libraries=[],  # refers to "liblas.2.3.0.dylib"
                language="c++",  # remove this if C and not C++
                include_dirs=[include_path],
                library_dirs=[library_path],
                requires=['Cython']
                )
    ]
else:
    include_path = os.path.join(env_path, "Library/include")
    library_path = os.path.join(env_path, "Library/lib")
    extensions = [
        Extension("raytr",
                sources=sources,
                libraries=[],  # refers to "liblas.2.3.0.dylib"
                language="c++",  # remove this if C and not C++
                include_dirs=[include_path],
                # ["C:/Users/kueken/Miniconda3/envs/OcclusionMapping_PDAL/Library/include/"],
                # include_dirs=["/usr/local/Cellar/boost/1.75.0_2/include", "/usr/local/Cellar/pdal/2.2.0_3/include", "/usr/local/Cellar/laszip/3.4.3/include"],
                library_dirs=[library_path],
                # ["C:/Users/kueken/Miniconda3/envs/OcclusionMapping_PDAL/Library/lib"],
                #library_dirs=["/usr/local/Cellar/pdal/2.2.0_3/lib"],
                requires=['Cython'],
                #extra_compile_args=["-ferror-limit=0"]
                #extra_compile_args=["-std=c++11"],
                #extra_compile_args=["-fopenmp", "-O3"],
                extra_compile_args=['-openmp'],  #Windows only
                extra_link_args=['-openmp']  #Windows only
                #extra_link_args=["L.C:/Users/kueken/Miniconda3/envs/OcclusionMapping_PDAL/Library/lib/pdalcpp"]
                #, "-DSOME_DEFINE_OPT","-L./some/extra/dependency/dir/"]
                )
    ]


# ------------ RIEGL IO ----------------#

NUMPY_MACROS = ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')

def getExtraCXXFlags():
    """
    Looks at the $PYLIDAR_CXX_FLAGS environment variable.
    If it exists, this function returns a list of flags
    to be passed as the extra_compile_args argument to
    the Extension constructor.
    Otherwise None.
    """
    if 'PYLIDAR_CXX_FLAGS' in os.environ:
        return os.environ['PYLIDAR_CXX_FLAGS'].split()
    else:
        return None

def addRieglRXPDriver(extModules, cxxFlags):
    """
    RIEGL RXP driver
    """
    if 'RIVLIB_ROOT' in os.environ:
        print('Building RIEGL RXP Extension...')
        rivlibRoot = os.environ['RIVLIB_ROOT']
        rivlibs = ['scanlib-mt', 'riboost_chrono-mt', 
                   'riboost_date_time-mt', 'riboost_filesystem-mt', 
                   'riboost_regex-mt', 'riboost_system-mt', 
                   'riboost_thread-mt']
        
        # on Windows the libs do not follow the normal naming convention
        # and start with 'lib'. On Linux the compiler prepends this automatically
        # but on Windows we need to do it manually
        if sys.platform == 'win32':
            rivlibs = ['lib' + name for name in rivlibs]
        elif sys.platform == 'darwin':
            print('RiVlib library is not available for MacOS')

        rieglModule = Extension(name='riegl_rxp', 
                define_macros=[NUMPY_MACROS],
                sources=['src/riegl_io/riegl_rxp.cpp', 'src/riegl_io/pylidar.c'],
                include_dirs=[os.path.join(rivlibRoot, 'include'), numpy.get_include()],
                extra_compile_args=cxxFlags,
                libraries=rivlibs,
                library_dirs=[os.path.join(rivlibRoot, 'lib')],
                runtime_library_dirs=[os.path.join(rivlibRoot, 'lib')])
                 
        extModules.append(rieglModule)
    else:
        print('RIEGL RXP libraries not found.')
        print('If installed set $RIVLIB_ROOT to the install location of RiVLib')

def addRieglRDBDriver(extModules, cxxFlags):
    """
    Decides if the Riegl RDB driver is to be built. If so 
    adds the Extension class to extModules.
    """
    if 'RDBLIB_ROOT' in os.environ:
        print('Building RIEGL RDB Extension...')

        rdblibRoot = os.environ['RDBLIB_ROOT']
        if sys.platform == 'win32':
            rdbLibName = 'rdblib'
        else:
            rdbLibName = 'rdb'

        defines = getRieglRDBLibVersion(rdblibRoot, rdbLibName)
        defines.extend([NUMPY_MACROS])

        rieglRDBModule = Extension(name='riegl_rdb',
                define_macros=defines,
                sources=['src/riegl_io/riegl_rdb.cpp', 'src/riegl_io/pylidar.c'],
                include_dirs=[os.path.join(rdblibRoot, 'interface', 'c'), numpy.get_include()],
                extra_compile_args=cxxFlags,
                libraries=[rdbLibName],
                library_dirs=[os.path.join(rdblibRoot, 'library')],
                runtime_library_dirs=[os.path.join(rdblibRoot, 'library')])

        extModules.append(rieglRDBModule)
    else:
        print('RIEGL RDB libraries not found.')
        print('If installed set $RDBLIB_ROOT to the install location of RDBLib')

def getRieglRDBLibVersion(rdbRoot, libname):
    """
    Because we cannot distribute the rdblib library, we need
    to check that the major version at compile time matches the 
    version the user has at runtime. We do this by getting the
    version now and setting it as a #define. The library can then
    use the #define to check at runtime.
    
    Unfortunately the headers don't give us this information.
    
    """
    if sys.platform == 'win32':
        libname = os.path.join(rdbRoot, 'library', libname + '.dll')
    elif sys.platform == 'darwin':
        libname = os.path.join(rdbRoot, 'library', 'lib' + libname + '.dylib')
    else:
        libname = os.path.join(rdbRoot, 'library', 'lib' + libname + '.so')
    rdb = ctypes.cdll.LoadLibrary(libname)

    context = ctypes.c_void_p()
    logLevel = ctypes.c_char_p(b"NONE")
    logPath = ctypes.c_char_p(b"")
    rdb.rdb_context_new(ctypes.byref(context), logLevel, logPath)

    version = ctypes.c_char_p()
    rdb.rdb_library_version(context, ctypes.byref(version))
    versionString = version.value
    if sys.version_info[0] >= 3:
        versionString = versionString.decode()

    rdb.rdb_context_delete(ctypes.byref(context))

    # versionString is quite specific  -something like:
    # 2.2.1-2094 (x86_64-linux, Jul 11 2019, 13:10:32)
    # we probably don't have to have the exact same string
    # so just extract major and minor version numbers
    arr = versionString.split('.')
    major = arr[0]
    minor = arr[1]

    return [("RIEGL_RDB_MAJOR", major),
            ("RIEGL_RDB_MINOR", minor)]

# get any C++ flags
cxxFlags = getExtraCXXFlags()

# External modules
pre_riegl_ext_len = len(extensions)
ext_cy = cythonize(extensions)
# addRieglRXPDriver(ext_cy, cxxFlags)
# addRieglRDBDriver(ext_cy, cxxFlags)

# if len(extensions) == pre_riegl_ext_len:
#     print('RiVLib and/or RDBLib not found. RIEGL I/O will not be supported.')

setup(
    name = 'occpy',
    version= occpy.__version__,
    packages = ["occpy"],
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_cy,
    install_requires=['Cython', 'numpy']
)
