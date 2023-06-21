# OccPy

<!---
## Name
Choose a self-explaining name for your project.
--->

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

<!---
## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.
--->

## Installation
To install OccPy, several steps are required which may or may not go through easily. The tool has been tested on Windows 10 and Mac OS X. Please let me know if you encounter any issues installing the tool.

### Seting up environment
We expect you to have a working conda installation (either through Anaconda or miniconda)
Either setup a new environment using the following command
```commandline
conda env create -f environment.yml
```
or within an existing environment install all the necessary packages:
```commandline
pip install -r requirements.txt
```

### Compile the c++ side of the OccPy tool
first we need to edit the _setup.py_ file in order for the compiler to find the necessary libraries.
Only two lines have to be changed as shown below (example for Windows machines). Adapt based on your system (paths will be different based on your conda installation):
```python
include_dirs=["C:/Users/_your_user_name_/Miniconda3/envs/_name_of_your_environment_/Library/include/"]
library_dirs=["C:/Users/_your_user_name_/Miniconda3/envs/_name_of_your_environment_/Library/lib"]
```
Once you have adapted _setup.py_ you can compile the C++ code using the following command:
```commandline
python setup.py build_ext --inplace
```
There will probably appear several warning messages. These can mostly be ignored (TODO: @kueken: check on these warnings!)
If compilation was successful, the tool should be ready to use. 

## Usage
There are three example scripts provided that should show how the tool can work for different flavors of LiDAR platforms (TLS, MLS, UAVLS)

- For TLS see Test_TLS.py
- For MLS see Test_MLS.py
- For UAVLS see Test_UAVLS.py

## Requirements for a successful occlusion mapping
In order for the occlusion mapping to work, several requirements on the input data have to be met. These are listed below specifically for the diffreent flavors of LiDAR platforms

### TLS

### MLS

### UAVLS

## Support
For questions and support, please contact Daniel Kükenbrink via daniel.kuekenbrink@wsl.ch

## Roadmap
Several open issues and improvements are currently worked on or planned for the future:

- [ ] Improve (add) documentation of the different functions and example scripts
- [ ] Add example data which should be used in the example scripts
- [ ] There is currently still an issue with UAVLS data, where some (very few) LiDAR returns are not registered by the algorithm. The implications for that should be analysed and the problem mitigated. This could cause an underestimation of occlusion, as the e.g. the last return is never reached and the pulse will traverse further without declaring an voxels as occluded for that pulse. There is the possibility to overcome this issue by using the function ```RayTr.doRaytracing_singleReturnPulses(x, y, z, sensor_x, sensor_y, sensor_z, gps_time, return_number, number_of_returns)``` as used in the script _Test_MLS.py_, where the input data is not initially converted to a pulse dataset, but each return is basically treated as a single pulse. We would only recommend to use this approach, if you are confident about your trajectory information.
- [ ] Add support for reading in a DTM file, so the algorithm could stop, once the pulse reached the terrain.
- [ ] substantial performance improvement by using multi core processing

<!---
## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.
--->

## Authors and acknowledgment
The algorithm is strongly based on the initial publication of the voxel traversal algorithm as seen in Amanatides & Woo (1987). 
This algorithm has been initially used in the publication by Kükenbrink et al. (2017) to map occlusion in ALS data and 
is openly available as a Matlab code here: https://www.eufar.net/documents/6028 (user account needed). Big motivation for
the development of this study came from the interesting paper by Bienert et al. (2010).
This implementation is a substantial evolution to the Matlab implementation and should now be able to run 
for any lidar platform available, when requirements as stated in [Requirements section](#Requirements for a successful occlusion mapping).
Also performance of this Cython implementation should be largely increased compared to the Matlab implementation.

Development of the initial Matlab implementation was performed during the PhD studies of Daniel Kükenbrink within the
EUFAR JRA - HYLIGHT project (EUFAR2 contract no. 312609). 
The initial development of the Cython version has started during the same PhD and was used in the study by 
Schneider et al. (2019) to map occlusion from TLS and UAVLS acquisitions in a temperate and tropical forest. 
Substantial improvements and further development has been done at the Swiss Federal Institute WSL since then. The development
is still ongoing also in the framework of the 3DForEcoTech COST action (working group 1). 

Big thank you go out to all contributing to this code base since the beginning of my PhD:
Felix Morsdorf, Fabian Schneider, Meinrad Abegg, Ruedi Bösch, Christian Ginzler

## Literature
Amanatides, J., Woo, A., 1987. A fast voxel traversal algorithm for ray tracing. Proc. EUROGRAPHICS 87, 3–10.

Bienert, A., Queck, R., Schmidt, A., 2010. Voxel Space Analysis of Terrestrial Laser Scans in Forests for Wind Field Modeling. Int. Arch. Photogramm. Remote Sens. Spat. Inf. Sci. - ISPRS Arch. 38, 92–97.

Kükenbrink, D., Schneider, F.D., Leiterer, R., Schaepman, M.E., Morsdorf, F., 2017. Quantification of hidden canopy volume of airborne laser scanning data using a voxel traversal algorithm. Remote Sens. Environ. 194, 424–436. https://doi.org/10.1016/j.rse.2016.10.023

Schneider, F.D., Kükenbrink, D., Schaepman, M.E., Schimel, D.S., Morsdorf, F., 2019. Quantifying 3D structure and occlusion in dense tropical and temperate forests using close-range LiDAR. Agric. For. Meteorol. 268. https://doi.org/10.1016/j.agrformet.2019.01.033

<!---
## License
For open source projects, say how it is licensed.
--->

## Project status
This tool is still under development and substantial testing with different datasets should be performed.