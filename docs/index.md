# Welcome to OccPy


<figure markdown>
  ![OccPy logo](assets/occpy_logo_v3_white.png#only-light){ width="500" }
  ![OccPy logo](assets/occpy_logo_v3_trans.png#only-dark){ width="500" } 
</figure>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📑 What is OccPy

OccPy is a python tool to map occluded area from LiDAR data in 3D using a voxel traversal algorithm implemented in C++.
OccPy traces each laser pulse from its origin along the pulse vector through the defined voxel grid. Only the point cloud
(in las or laz format) and pulse origin information (e.g. scanner trajectories for mobile acquisitions or scanner position
for stationary TLS acquisitions) as well as the grid dimensions and voxel resolution are needed as input to OccPy (see 
[Requirements for successful occlusion mapping](#requirements-for-successful-occlusion-mapping) for more information on how to use OccPy depending on your input data)

Implementation of the voxel traversal in C++ allows for a fast voxel traversal, while the python side allows for a user-friendly
and easy implementation of the occlusion mapping approach in a python environment and subsequent analysis and visualization 
of the occlusion maps. The interface between C++ and python is established via cython, which enables the calling of the C++ 
functions directly from python via the provided OccPy class. While transferring data between python and C++ can result in 
some overhead, running of the voxel traversal algorithm on the C++ side still results in a significant boost in performance
compared to a pure python implementation.

Alternatively to using las or laz files, it is also possible to run occpy using Riegl RXP and RDBX files. This data format
has the advantage that information on empty pulses is available. Therefore, also empty pulses can be traced for occlusion 
mapping, which can be essential when using ground-base (i.e. MLS or TLS) acquisitions, where pulses can travel into the 
atmosphere without any interaction. These pulses are not occluded and can have implication on occlusion. At the moment, 
these empty pulses are not tracable when using las or laz input files.

## 🛠 Installation
OccPy can be easily installed via pip:
(TODO: replace with pypi version)

```commandline
 pip install occpy-ls
```

Pre-built wheels are available for Python versions 3.10, 3.11, 3.12, 3.13, on:
- Linux (x86_64)
- Windows
- Mac OS

A source distribution is also available, which will require a C++ environment with boost libraries installed for a working installation.

Note: Riegl libraries are not packaged, to use RXP and RDBX files, you must build from source.

### Build from source
Clone the repository

```commandline
git clone https://github.com/dkueken/OccPy.git
cd OccPy
```

Set up the environment

```commandline
conda env create -f environment.yml
conda activate occPy
```

Build extensions and install occpy using:
```commandline
pip install -v .
```

NOTE: if you want to use RIEGL .rdbx and .rxp files as input, make sure to set the RIVLIB_ROOT and RDBLIB_ROOT environment variables to the root path of the corresponding libraries before installing.
The paths to RIVLIB_ROOT and RDBLIB_ROOT can be defined within the environment YAML file [environment.yml](https://github.com/dkueken/OccPy/blob/master/environment.yml), where the necessary Variables are commented out per default.

## Usage
The behavior of OccPy can be configured using JSON setting files (see e.g. [settings_MLS_tutorial.JSON](https://github.com/dkueken/OccPy/blob/master/config/settings_MLS_tutorial.JSON) as an example). 
Alternatively, a dictionary with all necessary input variables can be passed. We recommend setting up setting files for later reference on the settings.

First we import OccPy into your python (or jupyter) script:
```python
from occpy.OccPy import OccPy  # this loads the OccPy class with the core functionality
from occpy.util import normalize_occlusion_output   # within occpy.util multiple additional utility functions can be loaded, e.g. to normalize occlusion output
```
Afterwards we initiate an OccPy object for voxel traversal using the configured settings file:
```python
test = OccPy(config_file="path/to/config_file.json")
```
In the next step we define the sensor position, either by defining the scanner position or by providing trajectory information for mobile acquisitions. In this example we show the sensor position definition based on a handheld MLS acquistion using a GeoSLAM ZebHorizon processed using FARO Connect processing facilities.

```python
test.define_sensor_pos(path2file='path/to/trajectory_file.txt', 
                       delimiter=' ',               # delimiter used in the trajectory file
                       hdr_time='//world_time',     # column header for the time information in the trajectory file
                       hdr_x='x',                   # column header for the x coordiante in the trajectory file
                       hdr_y='y',                   # column header for the y coordinate in the trajectory file
                       hdr_z='z'                    # column header for the z coordinate in the trajectory file
                       )
```
Note, sensor position definition differs slightly when using TLS scans. See [TLS jupyter notbooks](notebooks/TLS_notebook.ipynb) for more information on this.

Once OccPy is fully parameterized, we can run it by simply calling

```python
test.do_raytracing()
```

This will store four output files as .npy files into the defined output directory: Nhit.npy, Nmiss.npy, Nocc.npy and Classification.npy
These can be loaded into your python script like:
```python
import numpy as np
Nhit = np.load("output_dir/Nhit.npy")
Nmiss = np.load("output_dir/Nmiss.npy")
Nocc = np.load("output_dir/Nocc.npy")
Classification = np.load("output_dir/Classification.npy")
```
- Nhit.npy: 3D numpy array with the number of laser hits per voxel
- Nmiss.npy: 3D numpy array with the number of misses (i.e. pulses that have no laser return in the specific voxel, but last return has not yet been reached)
- Nocc.npy: 3D numpy array witht the number of occluded pulses (i.e. the number of pulses that have already reached the last return before traversing the specific voxel)
- Classification.npy: 3D numpy array stating the classifciation for each voxel into Hit, Empty, Occluded and Unobserved, with the class definition like:
  - 1: Observed voxel with at least on registered return
  - 2: Empty voxel which was observed by at least one pulse that has not yet reached its last return
  - 3: occluded voxel, where all pulses traversing it were occluded
  - 4: unobserved voxel, where no pulse was traversed through it.

The classification into these classes is performed on the python side by using the following code:
```python
Classification[np.logical_and.reduce((Nhit > 0, Nmiss >= 0, Nocc >= 0))] = 1  # voxels that were observed
Classification[np.logical_and.reduce((Nhit == 0, Nmiss > 0, Nocc >= 0))] = 2  # voxels that are empty
Classification[np.logical_and.reduce((Nhit == 0, Nmiss == 0, Nocc > 0))] = 3  # voxels that are hidden (occluded)
Classification[np.logical_and.reduce((Nhit == 0, Nmiss == 0, Nocc == 0))] = 4  # voxels that were not observed # 

```
Note, this Classification grid follows a binary definition of occlusion, i.e. a voxel is labelled as occluded, only if no pulse was labelled as miss or return in the specific voxel.
If you prefer to define your own threshold for occlusion, or you would like to assess fractional occlusion, you could calculate the occlusion fraction per voxel like this:

```python
occl_frac = Nocc.astype(float) / (Nhit.astype(float) + Nmiss.astype(float) + Nocc.astype(float))
```

If you would like to have the output grids height normalized, this can be performed using the normalize_occlusion_output function

```python
from occpy.util import normalize_occlusion_output

normalize_occlusion_output(input_folder='path/to/occpy_output_dir',
                           PlotDim=plot_dim,           # Plot dimensions, i.e. corner coordinates like [min_x, min_y, min_z, max_x, max_y, max_z]
                           vox_dim=vox_dim,            # Voxel dimesnions in meters
                           dtm_file='path/to/DTM.tif', # path to DTM tif 
                           dsm_file='path/to/DSM.tif', # optional path to DSM
                           lower_threshold=lower_threshold,    # if voxels close to DTM should be ignored
                           output_voxels=False
                           )
```

For more information and examples, please see the following example jupyter scripts.


## 🎮 Examples

<table>
  <tr>
    <td style="background-color:#ffffff; color:black; padding:10px;">
      Example 1: TLS <br><br>
      <a href="notebooks/TLS_notebook.html">
        <img src="img/TLS_notebook.png" width="400">
      </a>
    </td>
  </tr>
  <tr>
    <td style="background-color:#ffffff; color:black; padding:10px;">
      Example 2: MLS <br>
      <a href="notebooks/MLS_notebook.html">
        <img src="img/MLS_notebook.png" width="400">
      </a>
    </td>
  </tr>
  <tr>
    <td style="background-color:#ffffff; color:black; padding:10px;">
      Example 3: UAVLS <br>
      <a href="notebooks/ULS_notebook.html">
        <img src="img/UAVLS_notebook.png" width="400">
      </a>
    </td>
  </tr>
  <tr>
</table>

## Config file examples
### TLS
This is an example for a configuration file to define the behavior of OccPy. This configuration file is defined for a 
TLS acquisition campaign and is used to run the [TLS jupyter notebook](notebooks/TLS_notebook.ipynb).
```json
{
    "root_folder": "../..",
    "laz_in": "data_notebooks/TLS_demo/LAZ/",
    "out_dir": "output/TLS",
    "vox_dim": 0.1,
    "lower_threshold": 1,
    "points_per_iter": 1000000,
    "plot_dim": [
        2676515,
        1246063,
        545,
        2676525,
        1246113,
        590
    ],
    "output_voxels": false,
    "single_return": false,
    "is_mobile": false,
    "str_idxs_ScanPosID": [7,10],
    "debug": true,
    "verbose": true,
    "tif_in": {
        "DTM": "data_notebooks/TLS_demo/Grids/Ramerenwald_DTM_20250305.tif",
        "DSM": "data_notebooks/TLS_demo/Grids/Ramerenwald_DSM_20250305.tif"
    },
    "ScanPos": "data_notebooks/TLS_demo/ScanPos/ScanPositions.txt"
}
```
### MLS / ULS
Tis is an example for a configuration file for the use with mobile acquisition (MLS/ULS). This specific configuration file
is used in the [MLS jupyter notebook](notebooks/MLS_notebook.ipynb). Most importantly, *is_mobile* should be set to *true*
and the trajectory file should be defined in *ScanPos*.
```json
{
    "root_folder": "../..",
    "laz_in": "data_notebooks/MLS_demo/LAZ/MLS_TestData_20perc_FP10_2025.laz",
    "tif_in": {
        "DTM": "data_notebooks/MLS_demo/Grids/Ramerenwald_DTM_20250305.tif",
        "DSM": "data_notebooks/MLS_demo/Grids/Ramerenwald_DSM_20250305.tif"
    },
    "out_dir": "output/MLS",
    "vox_dim": 0.1,
    "lower_threshold": 1,
    "points_per_iter": 1000000,
    "plot_dim": [
        2676515,
        1246063,
        545,
        2676525,
        1246113,
        590
    ],
    "ScanPos": "data_notebooks/MLS_demo/ScanPos/MLS_TestData_traj_FP10_2025.txt",
    "output_voxels": false,
    "single_return": true,
    "is_mobile": true
}
```
Please see documentation of the OccPy class for further explanation on the various parameters.



## Requirements for successful occlusion mapping
In order for the occlusion mapping to work, several requirements on the input data have to be met. 
These are listed below specifically for the different flavors of LiDAR platforms.

### TLS
**Scan Positions**

Scan Position file (as txt), where position should be referring to the laser source position.

**LAZ File**

- 1 LAZ or LAS file per scan position, preferably not filtered.

Currently, OccPy only works for TLS acquisitions if each scan station is stored in a separate LAZ file. Most Processing 
software solutions are capable of exporting the individual scan stations separately. If not, point identification to the 
corresponding scan station is usually performed via the point_source_id field in the LAZ file. Use this field to split a
combined point cloud. OccPy currently expects the linkage between scan position information (see [Scan Position paragraph above](#tls))
via the laz-file name. E.g. if your scan position information file has the following information:
```csv
ScanID, X,  Y,  Z
ScanPos001,      1,  2,  3
ScanPos002,      5,  2,  4
```
We would expect the ScanID number to be present in the laz-file name like this:
```commandline
ScanPos_001.laz
```
The exact location of the scan position identification number can be specified when initializing the OccPy object using 
the *str_idxs_ScanPosID* parameter either defined in *config_file* or the *config* dictionary. A list with the start and 
end index of the scan position identifier is expected, i.e. in the case of 'ScanPos010', ```python str_idxs_ScanPosID=[7,10]```.
Please see the [TLS config file example](#config-file-examples) for an example on setting the scan position linkage.

Note that according to python convention, positional indices in a vector is zero-based. Therefore, in the example given 
above, the scan posdtion identification *001* starts with the first *0* at index 7. Also following python convention, subsetting
a string will exclude the ending index. Therefore we have to add *+1* to the actual ending index of the scan position identification
(i.e. in the example given above, the *1* is located at index 9, therefore we have to specify 10).

If a multi return TLS is used, you can improve performance by sorting the LAZ file according to GPS Time and return number, 
e.g. by using LASTools's lassort function:

```commandline
lassort -i in_laz -gps_time -return_number -odix _sort -olaz -cpu64 -v
```

### MLS
**Trajectory file**

A trajectory file is strongly needed for the algorithm to work with MLS data. 
The following data should be present in trajectory file:

- Time (usually GPS time in seconds) - be sure that the GPS time format corresponds to the one stored in the gps_time field of the laz file
- Position of the sensor in X, Y, Z coordinates

The pose of the sensor is currently not regarded (e.g. quaternions). We are expecting that the coordinates in the trajectory 
corresponds to the position of the laser source.

The gps time tags do not need to be exactly the same as found in the gps_time field of the laz file, as the exact position will be interpolated based on the gps_time. 
However, a higher frequency in positional readings of the trajectory file will result in more accurate interpolation of the scanner position and hence a more accurate occlusion map.

**LAZ file**

As stated before, the biggest requirement for the LAZ file is that gps_time field is corresponding to the gps time readings in the trajectory file.

### ULS
As in the case for MLS data, trajectories are a hard requirement for UAVLS data. Please refer to 
[MLS](#mls) section for the requirements on the trajectory file. Also check out [MLS_notebook.ipynb](notebooks/MLS_notebook.ipynb)  or 
[ULS_notebook.ipynb](notebooks/ULS_notebook.ipynb) jupyter notebooks for how to use this tool for occlusion mapping.

## 📚 Related publications
<details>
<summary>Journal</summary>

```
@article{Amanatides1987,
    author = {Amanatides, John and Woo, Andrew},
    year = {1987},
    month = {08},
    pages = {},
    title = {A Fast Voxel Traversal Algorithm for Ray Tracing},
    volume = {87},
    journal = {Proceedings of EuroGraphics}
}
```

```
@article{Bienert2010,
    author = {Bienert, Anne and Queck, Ronald and A, A. and Maas, Hans-Gerd},
    year = {2010},
    month = {01},
    pages = {92-97},
    title = {Voxel space analysis of terrestrial laser scans in forests for wind field modelling},
    volume = {XXXVIII, Part 5},
    journal = {International Archives of Photogrammetry, Remote Sensing and Spatial Information Sciences}
}
```

```
@article{KUKENBRINK2017424,
    title = {Quantification of hidden canopy volume of airborne laser scanning data using a voxel traversal algorithm},
    journal = {Remote Sensing of Environment},
    volume = {194},
    pages = {424-436},
    year = {2017},
    issn = {0034-4257},
    doi = {https://doi.org/10.1016/j.rse.2016.10.023},
    url = {https://www.sciencedirect.com/science/article/pii/S0034425716303959},
    author = {Daniel Kükenbrink and Fabian D. Schneider and Reik Leiterer and Michael E. Schaepman and Felix Morsdorf}}
```

```
@article{SCHNEIDER2019249,
    title = {Quantifying 3D structure and occlusion in dense tropical and temperate forests using close-range LiDAR},
    journal = {Agricultural and Forest Meteorology},
    volume = {268},
    pages = {249-257},
    year = {2019},
    issn = {0168-1923},
    doi = {https://doi.org/10.1016/j.agrformet.2019.01.033},
    url = {https://www.sciencedirect.com/science/article/pii/S0168192319300267},
    author = {Fabian D. Schneider and Daniel Kükenbrink and Michael E. Schaepman and David S. Schimel and Felix Morsdorf}}
```

</details>

## 📂 Credits

**How to cite**

Please cite the following studies when using OccPy.

```
@article{KUKENBRINK2017424,
    title = {Quantification of hidden canopy volume of airborne laser scanning data using a voxel traversal algorithm},
    journal = {Remote Sensing of Environment},
    volume = {194},
    pages = {424-436},
    year = {2017},
    issn = {0034-4257},
    doi = {https://doi.org/10.1016/j.rse.2016.10.023},
    url = {https://www.sciencedirect.com/science/article/pii/S0034425716303959},
    author = {Daniel Kükenbrink and Fabian D. Schneider and Reik Leiterer and Michael E. Schaepman and Felix Morsdorf}}
```

```
@article{SCHNEIDER2019249,
    title = {Quantifying 3D structure and occlusion in dense tropical and temperate forests using close-range LiDAR},
    journal = {Agricultural and Forest Meteorology},
    volume = {268},
    pages = {249-257},
    year = {2019},
    issn = {0168-1923},
    doi = {https://doi.org/10.1016/j.agrformet.2019.01.033},
    url = {https://www.sciencedirect.com/science/article/pii/S0168192319300267},
    author = {Fabian D. Schneider and Daniel Kükenbrink and Michael E. Schaepman and David S. Schimel and Felix Morsdorf}}
```


## Authors, Funding and  acknowledgements

The algorithm is strongly based on the initial publication of a voxel traversal algorithm as seen in Amanatides & Woo (1987). 
This algorithm has been used in the publication by Kükenbrink et al. (2017) to map occlusion in ALS data and is openly 
available as a Matlab code here: https://www.eufar.net/documents/6028 (user account needed). Big motivation for the 
development of this study came from the interesting paper by Bienert et al. (2010). This implementation is a substantial 
evolution to the Matlab implementation and should now be able to run for any lidar platform available, when requirements 
as stated in [Requirements section](#requirements-for-successful-occlusion-mapping) are met. Also, performance of this Cython implementation should be largely increased 
compared to the Matlab implementation.

Development of the initial Matlab implementation was performed during the PhD studies of Daniel Kükenbrink at the University of 
Zurich within the EUFAR JRA - HYLIGHT project (EUFAR2 contract no. 312609). The initial development of the Cython version 
has started during the same PhD and was used in the study by Schneider et al. (2019) to map occlusion from TLS and UAVLS 
acquisitions in a temperate and tropical forest. Substantial improvements and further development has been done at the 
Swiss Federal Institute WSL since then. The development is still ongoing also in the framework of the 3DForEcoTech COST 
action (working group 1).

Big thank you go out to all contributing to this code base since the beginning of my PhD,  Felix Morsdorf, Fabian Schneider, 
Meinrad Abegg, Ruedi Bösch, Christian Ginzler as well as to those pushing the code base towards the publication of OccPy as
a python package: William Albert, Wout Cherlet, Bernhard Höfle, and Jonas Wenk. 


## Contact / bugs / feature requests
Have you found a bug or have specific request for a new feature? Please open a new issue in the online code repository on <a href="https://github.com/dkueken/OccPy">Github</a>.

Scientific requests or questions can be directed to
<a href="https://www.wsl.ch/en/staff/kueken/">Daniel Kükenbrink</a>.

## Roadmap
Several open issues and improvements are currently worked on or planned for the future:

- Add support for reading in a DTM file into the voxel traversal, so the algorithm could stop, once the pulse reached the terrain.
- Substantial performance improvement by using multi core processing
- Add functionality for PAI/PAD calculation of each voxel (i.e. calculation of path length within voxel for each pulse)


## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License
This is licensed under the [MIT license](https://opensource.org/licenses/MIT).

## Project status
This tool is still under development and substantial testing with different datasets should be performed.
