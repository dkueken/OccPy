# Ideas for Tasks to be done
## Familiarize with tool
- create test datasets (smaller ones) for testing
- Test all kinds of lidar flavors and report issues (also different level of details to see if we run into problems with larger point clouds)

## Performance improvements
- currently only a very rudimentary implementation of multicore-processing is implemented (only for single return pulses so far)
	- This should be enhanced and tested
	- Also for multi return pulses! There is a different function for that, which first sets up a pulse dataset and then does the raytracing -> test if this can be done using the openmp library
	- check on transferability to other systems (openmp is only avaialable for windows?)
		- It would be fine to develop the tool only for windows, but if there are easy ways to make it compatible also with other systems that would be nice
- Progress report improvements when using multi core processing


## Adding Functionality
- improve visualization capabilities. Initial steps have been done but could be improved. Are there potential good packages around?
- Currently, only infinitesimally small pulses are assumed, which is a big simplification. See if we can include pulse footprint size into the raytracing -> might be a bit complex but I think Uni Heidelberg has already some knowledge on how to implement this
- potentially add DTM as an input to the C++ side, so that voxel traversal can stop once the pulse hit the terrain (not sure how big the performance gain can be from that, but in topographically challenging terrain this could be substantial)
- Add ways to estimate PAD from the voxel traversal output (most of the necessary information should be there)

# Packaging
- Find out ways how we could package everything and distribute via pip and/or conda
	- I have no experience in doing this, but you should have plenty and I hope you could help me here to figure out, what is needed and should be done to reach this goal
- See how compiling would work so that the user does not have to specify any paths in the setup.py file.
- I guess you will have some suggestions on how to improve the structure of the package
- setup a light weight testing dataset that could be used in testing the tool using notebooks (see below) and find ways how to distribute this (we could upload some dataset to a repository or use this dataset [https://www.doi.org/10.16904/envidat.383](https://www.doi.org/10.16904/envidat.383))

# Documentation
- improve documentation 
- setup notebooks showcasing the tool for different lidar flavors

# Code clean-up
- there is a lot of potential for cleaning up code 