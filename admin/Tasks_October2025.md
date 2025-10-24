# TODOs as of October 24, 2025

## Code cleanup, refactoring
- get rid of unnecessary lines, comments etc. [Daniel]
- OccPy.py still has a lot of functions on top that could be moved to uitl.py or visualization.py [Daniel]


## Documentation
- make function documentation template consistent, so that on the webpage all looks similar. See docstrings by Wout in OccPyRIEGL [Daniel]
- add function documentation for those functions that do not yet have one [Daniel]
- *Question:* What should we do with the c++ side? Should this be documented as well?
- Updated README [Daniel]

## Notebooks
- add tutorial data to zenodo or other repository and get data through URL, so data does not have to be stored alongside code [Daniel]
- Cleanup notebooks and make more concise [Daniel]

## Test scripts
- Cleanup and adapt to code changes [Daniel]
- Define which are necessary [Daniel]

## Test OccPyRIEGL
- Update to code base changes [Wout]
- Test everything [Wout]
- Add Section about OccPyRIEGL to documentation

# Installation
- cleanup environment.yml
  - get rid of unnecessary entries. 
  - test if environment.yml is working for any OS flavors
- test pip install . on various OS flavors
  - There is currently a deprecation warning which should be resolved 

# Testing Testing Testing
- Test on various OS
- Test various data sources
- Implement Testing routines (Unit tests etc.) **Ask Hannah!**

# GitHub Page
- Cleanup
- Design OccPy Logo [Wout?]

# Packaging
- How to implement Packages and Releases feature on github?

# Visualization
- add all visualization functions to visualization.py and test them [Daniel]
- test pyvista? [Wout]
- 

# Known bugs
- Issues with non-integer box coordinates, see Teams message by Wout. This has not yet been merged into the integration/masters branch

## Nice-to-have

- [Performance] Multi-core processing
- [Extension] Pulse footprint consideration
- [Performance] DTM on C++ side for early voxel traversal termination
- [Extension] PAD estimation
- [Extension] Dynamic visualization option?