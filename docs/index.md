# Welcome to OccPy


<figure markdown>
  ![OccPy logo](img/occpy_logo_light.png#only-light){ width="500" }
  ![OccPy logo](img/occpy_logo_dark.png#only-dark){ width="500" } 
</figure>

<span style="color:red"> 
TOCHECK: Change license if not the good one
</span>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📜 History of OccPy

Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.

Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.



### Overall objective

Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.

Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.


### Applications

Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.

Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.


## 🎮 Examples

<span style="color:red"> 
TOCHECK: Update figures with more representative ones
</span>

<table>
  <tr>
    <td style="background-color:#ffffff; color:black; padding:10px;">
      Example 1: TLS <br><br>
      <a href="jupyter/TLS_notebook.html">
        <img src="img/TLS_notebook.png" width="400">
      </a>
    </td>
  </tr>
  <tr>
    <td style="background-color:#ffffff; color:black; padding:10px;">
      Example 2: MLS <br>
      <a href="jupyter/MLS_notebook.html">
        <img src="img/MLS_notebook.png" width="400">
      </a>
    </td>
  </tr>
  <tr>
    <td style="background-color:#ffffff; color:black; padding:10px;">
      Example 3: UAVLS <br>
      <a href="jupyter/UAVLS_notebook.html">
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

### Potential Issues
Currently, packages installed via pip are not listed in the environment.yml file (e.g. laspy). It could be that you need
to install these packages using pip:
e.g.
```commandline
pip install laspy[laszip]
```
the _[laszip]_ option enables the reading of .laz files 

or you could also install all packages with the following command (potential package conflicts between conda and pip not tested)

```commandline
pip install -r requirements.txt
```


### List of needed packages
This is a list of needed packages for the tool to run, if the instal via environment.yml or requirements.txt fails:

- boost
- cython
- numpy
- laspy*
- pandas
- scipy

'*' should be installed through pip (see below)

installing most current version should in theory work.
In order for laz compatibility, install _laspy_ with _laszip_ compatibility through
```commandline
pip install laspy[laszip]
```

### Compile the c++ side of the OccPy tool

Compile c++ side and install occpy package using:
pip install -v .

There will probably appear several warning messages. These can mostly be ignored (TODO: @kueken: check on these warnings!)
If compilation was successful, the tool can be imported using 

Navigate to the src folder where setup.py is located
```commandline
cd occpy/src
```
Then compile the code using the following command:
```commandline
python setup.py build_ext --inplace
```
 
However, if the compiled raytr file is not in the same directory as your python code calling the raytracer, you should add the following lines to the beginning of your python script
```python
import sys
sys.path.append(r".\src) # assuming you have your python file in the root directory of the tool. change adequately if this is not the case.
```


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