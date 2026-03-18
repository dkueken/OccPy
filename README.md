![Occpy](assets/occpy_logo_v3_trans.png "Occpy logo")

---

OccPy is a python tool to map occluded area from LiDAR data in 3D using a voxel traversal algorithm implemented in C++. 

## Installation

Via pip:
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

## Usage
There are example notebooks provided for different flavors of LiDAR platforms (TLS, MLS, UAVLS), under docs/notebooks.

TODO: maybe add some quick start guide to the webpage?

## Requirements for a successful occlusion mapping
In order for the occlusion mapping to work, several requirements on the input data have to be met. These are listed below specifically for the different flavors of LiDAR platforms.

### TLS
**Scan Positions**
- Scan Position file (as txt), where position should be referring to the laser source position.

**LAZ File**
- 1 LAZ or LAS file per scan position, preferably not filtered. 

If a multi return TLS is used, you can improve
performance by sorting the LAZ file according to GPS Time and return number, e.g. by using LASTools's lassort function:

```commandline
lassort -i in_laz -gps_time -return_number -odix _sort -olaz -cpu64 -v
```
### MLS
**Trajectory file**

A trajectory file is strongly needed for the algorithm to work with MLS data. The following data should be present in 
trajectory file:
- Time (usually GPS time in seconds) - be sure that the GPS time format corresponds to the one stored in the gps_time field of the laz file
- Position of the sensor in X, Y, Z coordinates

The pose of the sensor is currently not regarded (e.g. quaternions) We are expecting that the coordinates in the trajectory
corresponds to the position of the laser source.

The gps time tags do not need to be exactly the same as found in the gps_time field of the laz file, as the exact position
will be interpolated based on the gps_time. However, a higher frequency in positional readings of the trajectory file will 
result in more accurate interpolation of the scanner position and hence a more accurate occlusion map.

**LAZ file**

As stated before, the biggest requirement for the LAZ file is that gps_time field is corresponding to the gps time readings in the trajectory file.

### UAVLS
**Trajectory file**

As in the case for MLS data, trajectories are a hard requirement for UAVLS data. Please refer to 
[MLS](MLS) section for the requirements on the trajectory file. Also check out [Test_MLS.py](Test_MLS.py) or 
[Test_UAVLS](Test_UAVLS.py) python scripts for how to use this tool for occlusion mapping.

**LAZ file**

As UAVLS data often come as multi-return data, it is again recommended to sort the LAZ file based on gps_time and 
return_number like:
```commandline
lassort -i laz_in -gps_time -return_number -odix _sort -olaz -cpu64 -v
```

Unsorted LAZ files will also work, however, there will be a substantial computational overhead, as the entire dataset 
needs to be read and stored at once.

## Visualization

Check the pyvista and pyvista_interactive notebooks to visualize occlusion map outputs in 3D.

![Pyvista demo](assets/pv_readme.webp)

## Support
For questions and support, please contact Daniel Kükenbrink via daniel.kuekenbrink@wsl.ch

## Roadmap
Several open issues and improvements are currently worked on or planned for the future:

- Improve (add) documentation of the different functions and example scripts
- Add example data which should be used in the example scripts
- There is currently still an issue with UAVLS data, where some (very few) LiDAR returns are not registered by the algorithm. The implications for that should be analysed and the problem mitigated. This could cause an underestimation of occlusion, as the e.g. the last return is never reached and the pulse will traverse further without declaring an voxels as occluded for that pulse. There is the possibility to overcome this issue by using the function ```RayTr.doRaytracing_singleReturnPulses(x, y, z, sensor_x, sensor_y, sensor_z, gps_time, return_number, number_of_returns)``` as used in the script _Test_MLS.py_, where the input data is not initially converted to a pulse dataset, but each return is basically treated as a single pulse. We would only recommend to use this approach, if you are confident about your trajectory information.
- Add support for reading in a DTM file, so the algorithm could stop, once the pulse reached the terrain.
- Add functionality for height normalisation of outputs
- Substantial performance improvement by using multi core processing
- Add functionality for PAI/PAD calculation of each voxel (i.e. calculation of path length within voxel for each pulse) 
- Add visualization solution like (potential idea: https://github.com/msoechting/lexcube)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## Authors and acknowledgment
The algorithm is strongly based on the initial publication of a voxel traversal algorithm as seen in Amanatides & Woo (1987). 
This algorithm has been used in the publication by Kükenbrink et al. (2017) to map occlusion in ALS data and 
is openly available as a Matlab code here: https://www.eufar.net/documents/6028 (user account needed). Big motivation for
the development of this study came from the interesting paper by Bienert et al. (2010).
This implementation is a substantial evolution to the Matlab implementation and should now be able to run 
for any lidar platform available, when requirements as stated in [Requirements section](#Requirements for a successful occlusion mapping) 
are met. Also performance of this Cython implementation should be largely increased compared to the Matlab implementation.

Development of the initial Matlab implementation was performed during the PhD studies of Daniel Kükenbrink at the University of Zurich 
within the EUFAR JRA - HYLIGHT project (EUFAR2 contract no. 312609). 
The initial development of the Cython version has started during the same PhD and was used in the study by 
Schneider et al. (2019) to map occlusion from TLS and UAVLS acquisitions in a temperate and tropical forest. 
Substantial improvements and further development has been done at the Swiss Federal Institute WSL since then. The development
is still ongoing also in the framework of the [3DForEcoTech COST action](https://3dforecotech.eu/) (working group 1). 

Big thank you go out to all contributing to this code base since the beginning of my PhD:
Felix Morsdorf, Fabian Schneider, Meinrad Abegg, Ruedi Bösch, Christian Ginzler

## Literature
Amanatides, J., Woo, A., 1987. A fast voxel traversal algorithm for ray tracing. Proc. EUROGRAPHICS 87, 3–10.

Bienert, A., Queck, R., Schmidt, A., 2010. Voxel Space Analysis of Terrestrial Laser Scans in Forests for Wind Field Modeling. Int. Arch. Photogramm. Remote Sens. Spat. Inf. Sci. - ISPRS Arch. 38, 92–97.

Kükenbrink, D., Schneider, F.D., Leiterer, R., Schaepman, M.E., Morsdorf, F., 2017. Quantification of hidden canopy volume of airborne laser scanning data using a voxel traversal algorithm. Remote Sens. Environ. 194, 424–436. https://doi.org/10.1016/j.rse.2016.10.023

Schneider, F.D., Kükenbrink, D., Schaepman, M.E., Schimel, D.S., Morsdorf, F., 2019. Quantifying 3D structure and occlusion in dense tropical and temperate forests using close-range LiDAR. Agric. For. Meteorol. 268. https://doi.org/10.1016/j.agrformet.2019.01.033

## How to cite
For now, please cite the following studies

Kükenbrink, D., Schneider, F.D., Leiterer, R., Schaepman, M.E., Morsdorf, F., 2017. Quantification of hidden canopy volume of airborne laser scanning data using a voxel traversal algorithm. Remote Sens. Environ. 194, 424–436. https://doi.org/10.1016/j.rse.2016.10.023

and

Schneider, F.D., Kükenbrink, D., Schaepman, M.E., Schimel, D.S., Morsdorf, F., 2019. Quantifying 3D structure and occlusion in dense tropical and temperate forests using close-range LiDAR. Agric. For. Meteorol. 268. https://doi.org/10.1016/j.agrformet.2019.01.033


## License
See [LICENSE](LICENSE).

## Project status
This tool is still under development and substantial testing with different datasets should be performed.
