![Occpy](assets/occpy_logo_v3_trans.png "Occpy logo")

---

OccPy is a python tool to map occluded area from LiDAR data in 3D using a voxel traversal algorithm implemented in C++. 

## Installation

Via pip:

```commandline
 pip install occpy-ls
```

Pre-built wheels are available for Python versions 3.10, 3.11, 3.12, 3.13, on:
- Linux (x86_64)
- Windows
- MacOS
A source distribution is also available, which will require a C++ environment with boost libraries installed for a working installation.

Note: Proprietary RIEGL libraries are not packaged, to use RXP and RDBX files, you must build from source.

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
Head over to the Documentation page [TODO: add link to published github page] and have a look at the various jupyter notebooks
to find out how to configure and run OccPy.

### TLS examples
There are two ways to run OccPy on a multi-station TLS dataset, where the point cloud of each scan position is stored in a separate LAZ file:
- passing the input folder and let OccPy handle all files automatically: [TODO: Link to notebook]
- handling the separate laz files in the input folder individually yourself: [TODO: Link to notebook]. This approach can bring some performance benefits and provides more flexibility on how to treat the individual scans.
### Mobile examples
Both MLS and ULS acquisitions can be treated similarly, as shown in these two notebooks:
- MLS: [TODO: Link to notebook]
- ULS: [TODO: Link to notebook]

### Visualization
The above stated notebooks provide some inputs on creating 2D visualizations of the occlusion mapping outputs.
If you would like to visualize occlusion outputs in 3D, please check out the pyvista and pyvista_interactive notebooks:
- PyVista demo: [TODO: Link to notebook]
- PyVista interactive demo: [TODO: Link to notebook]

Here is an example 3D visualization of an occlusion map for a multi-station TLS campaign as provided by Wout Cherlet and 
shown at SilviLaser 2025 in Québec City.

![Pyvista demo](assets/pv_readme.webp)

## Support
For questions and support, please contact Daniel Kükenbrink via daniel.kuekenbrink@wsl.ch

## Roadmap
Several open issues and improvements are currently worked on or planned for the future:

- Add support for reading in a DTM file into the voxel traversal, so the algorithm could stop, once the pulse reached the terrain.
- Substantial performance improvement by using multi core processing
- Add functionality for PAI/PAD calculation of each voxel (i.e. calculation of path length within voxel for each pulse) 

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## Authors and acknowledgment
The algorithm is strongly based on the initial publication of a voxel traversal algorithm as seen in Amanatides & Woo (1987). 
This algorithm has been used in the publication by Kükenbrink et al. (2017) to map occlusion in ALS data and 
is openly available as a Matlab code here: https://www.eufar.net/documents/6028 (user account needed). Big motivation for
the development of this study came from the interesting paper by Bienert et al. (2010).
This implementation is a substantial evolution to the Matlab implementation and should now be able to run 
for any lidar platform available, when requirements as stated in [Requirements section](#requirements-for-a-successful-occlusion-mapping) 
are met. Also performance of this Cython implementation should be largely increased compared to the Matlab implementation.

Development of the initial Matlab implementation was performed during the PhD studies of Daniel Kükenbrink at the University of Zurich 
within the EUFAR JRA - HYLIGHT project (EUFAR2 contract no. 312609). 
The initial development of the Cython version has started during the same PhD and was used in the study by 
Schneider et al. (2019) to map occlusion from TLS and UAVLS acquisitions in a temperate and tropical forest. 
Substantial improvements and further development has been done at the Swiss Federal Institute WSL since then. The development
is still ongoing also in the framework of the [3DForEcoTech COST action](https://3dforecotech.eu/) (working group 1). 

Big thank you go out to all contributing to this code base since the beginning of my PhD,  Felix Morsdorf, Fabian Schneider, 
Meinrad Abegg, Ruedi Bösch, Christian Ginzler as well as to those pushing the code base towards the publication of OccPy as
a python package: William Albert, Wout Cherlet, Bernhard Höfle, and Jonas Wenk. 


## Literature

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
    author = {Daniel Kükenbrink and Fabian D. Schneider and Reik Leiterer and Michael E. Schaepman and Felix Morsdorf}
    }
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
    author = {Fabian D. Schneider and Daniel Kükenbrink and Michael E. Schaepman and David S. Schimel and Felix Morsdorf}
    }
```
## How to cite
For now, please cite the following studies

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
    author = {Daniel Kükenbrink and Fabian D. Schneider and Reik Leiterer and Michael E. Schaepman and Felix Morsdorf}
    }
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
    author = {Fabian D. Schneider and Daniel Kükenbrink and Michael E. Schaepman and David S. Schimel and Felix Morsdorf}
    }
```

## License
See [LICENSE](LICENSE).

## Project status
This tool is still under development and substantial testing with different datasets should be performed.
