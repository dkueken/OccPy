from sample_scripts.Test_UAVLS import laz_in

# Welcome to OccPy


<figure markdown>
  ![OccPy logo](assets/occpy_logo_v3_white.png#only-light){ width="500" }
  ![OccPy logo](assets/occpy_logo_v3_trans.png#only-dark){ width="500" } 
</figure>

<span style="color:red"> 
TOCHECK: Change license if not the good one
</span>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📑 What is OccPy

OccPy is a python tool to map occluded area from LiDAR data in 3D using a voxel traversal algorithm implemented in C++.
OccPy traces each laser pulse from its origin along the pulse vector through the defined voxel grid. Only the point cloud
(in las or laz format) and pulse origin information (e.g. scanner trajectories for mobile acquisitions or scanner position
for stationary TLS acquisitions) as well as the grid dimensions and voxel resolution are needed as input to OccPy (see 
[Requirements for successful occlusion mapping](#requirements-for-successful-occlusion-mapping) for more information on how to use OccPy depending on your input data)

Alternatively to using las or laz files, it is also possible to run occpy using Riegl RXP and RDBX files. This data format
has the advantage that information on empty pulses is available. Therefore, also empty pulses can be traced for occlusion 
mapping, which can be essential when using ground-base (i.e. MLS or TLS) acquisitions, where pulses can travel into the 
atmosphere without any interaction. These pulses are not occluded and can have implication on occlusion. At the moment, 
these empty pulses are not tracable when using las or laz input files.

## Installation
OccPy can be easily installed via pip:
(TODO: replace with pypi version)

```commandline
 pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple occpy_ls==0.1
```

Pre-built wheels are available for Python versions 3.10, 3.11, 3.12, 3.13, on:
- Linux (x86_64)
- Windows (TODO)

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
The paths to RIVLIB_ROOT and RDBLIB_ROOT can be defined within the environment YAML file [environment.yml](../environment.yml), where the necessary Variables are commented out per default.

## Usage
The behavior of OccPy can be configured using JSON setting files (see e.g. [settings_MLS_tutorial.JSON](../config/settings_MLS_tutorial.JSON) as an example). 
While this is not strictly necessary, we still recommend setting up such setting files for running OccPy for later reference on the settings.

First we import OccPy into your python (or jupyter) script:
```python
from occpy.OccPy import OccPy  # this loads the OccPy class with the core functionality
from occpy.util import normalize_occlusion_output   # within occpy.util multiple additional utility functions can be loaded, e.g. to normalize occlusion output
```
Afterwards we initiate an OccPy object for voxel traversal like the following:
```python
test = OccPy(laz_in="path/to/laz_file.laz",
             out_dir="path/to/output_dir",
             vox_dim=0.1,                   # voxel dimesnions in m (for the moment only cubic voxels are allowed)
             lower_threshold=1,             # lower threshold in meters, to reduce effects caused by terrain (only actually functional if Terrain Model is provided)
             points_per_iter=10000000,      # number of points to be loaded at once. Note, this is only active if point cloud is sorted along gps_time or it is single return. Otherwise entire dataset needs to be loaded at once
             plot_dim=[min_x, min_y, min_z, max_x, max_y, max_z] # Corner coordinates of the voxel grid. Note only integer values are currently supported here (do not define sub-meter corner coordiantes for the moment!)
             )
```
In the next step we define the sensor position, either by defining the scanner position or by providing trajectory information for mobile acquisitions. In this example we show the sensor position definition based on a handheld MLS acquistion using a GeoSLAM ZebHorizon processed using FARO Connect processing facilities.

```python
test.define_sensor_pos(path2file="path/to/trajectory_file.txt", 
                       is_mobile=True,              # whether acquisition is mobile. Always true for MLS or ULS
                       single_return=True,          # whether the data is single or multi return
                       delimiter=" ",               # delimiter used in the trajectory file
                       hdr_time='//world_time',     # column header for the time information in the trajectory file
                       hdr_x='x',                   # column header for the x coordiante in the trajectory file
                       hdr_y='y',                   # column header for the y coordinate in the trajectory file
                       hdr_z='z'                    # column header for the z coordinate in the trajectory file
                       )
```
Note, sensor position definition differs slightly when using TLS scans. See [TLS jupyter notbooks](notebooks/TLS_notebook.ipynb) for more information on this.

Once OccPy is fully parameterized, we can run it bz simply calling

```python
test.do_raytracing()
```

This will store four output files as .npy files into the defined output directory: Nhit.npy, Nmiss.npy, Nocc.npy and Classification.npy
These can be loaded into you python script like:
```python
import numpy as np
Nhit = np.load("output_dir/Nhit.npy")
Nmiss = np.load("output_dir/Nmiss.npy")
Nocc = np.load("output_dir/Nocc.npy")
Classification + np.load("output_dir/Classification.npy")
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

If you would like to have the output grids height normalized, this can be perfored using the normalize_occlusion_output function

```python
from occpy.util import normalize_occlusion_output

Nhit_norm, Nmiss_norm, Nocc_norm, Classification_norm, chm = normalize_occlusion_output(input_folder='out_dir',
                                                                                        PlotDim=plot_dim,           # Plot dimensions, i.e. corner coordinates like [min_x, min_y, min_z, max_x, max_y, max_z]
                                                                                        vox_dim=vox_dim,            # Voxel dimesnions in meters
                                                                                        dtm_file='path/to/DTM.tif', # path to DTM tif 
                                                                                        dsm_file='path/to/DSM.tif', # optional path to DSM
                                                                                        lower_threshold=lower_threshold,    # if voxels close to DTM should be ignored
                                                                                        output_voxels=False)
```

For more information and examples, plese see the following example jupytor scripts.


## 🎮 Examples

<span style="color:red"> 
TOCHECK: Update figures with more representative ones
</span>

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


## 🛠️ How to install

To install OccPy, several steps are required which may or may not go through easily. The tool has been tested on Windows 10. 
If you encounter any issues installing the OccPy, please open an issue on the GitHub repository.


### Clone repository

Clone the repository using git with the following command (or download the zip from gitlab):
```commandline
git clone https://github.com/dkueken/OccPy.git
```

cd into the cloned repository
```commandline
cd OccPy
```

### Setting up the environment
We expect you to have a working conda installation (either through Anaconda or miniconda)
Either setup a new environment using the following command
```commandline
conda env create -f environment.yml
```
activate the environment
```commandline
conda activate occPy
```
or within an existing environment install all the necessary packages:
```commandline
pip install -r requirements.txt
```

<span style="color:red">
<br>
TOCHECK: If you would prefer to have the further information ("Potential issues", "List of needed packages", and "Compile the c++ side of the OccPy tool") in a drop down menu like the "Related publications", it is possible with the 'details' and 'summary' tags. In that case, the formating (underline, titles, etc) would have to be switched to HTML instead of markdown.
</span>

## Requirements for successful occlusion mapping
In order for the occlusion mapping to work, several requirements on the input data have to be met. 
These are listed below specifically for the different flavors of LiDAR platforms.

### TLS
**Scan Positions**

Scan Position file (as txt), where position should be referring to the laser source position.

**LAZ File**

- 1 LAZ or LAS file per scan position, preferably not filtered.

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


**Funding / acknowledgements**

If any, add funding or acknowledgements here.


**Contact / bugs / feature requests**

Have you found a bug or have specific request for a new feature? Please open a new issue in the online code repository on <a href="https://github.com/dkueken/OccPy">Github</a>.

Scientific requests can be directed to
<a href="https://www.wsl.ch/en/staff/kueken/">Daniel Kükenbrink</a>.



**License**
<span style="color:red">  
TOCHECK: Change license if not the good one
</span>

This is licensed under the [MIT license](https://opensource.org/licenses/MIT).