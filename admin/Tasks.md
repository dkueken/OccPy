# TODOs
## Testing
- [ ] Setup light weight testing data for each LiDAR flavor (MLS, TLS, UAVLS, ALS etc. multi-return and single-return sensors). -> @Jonas
- [ ] Extensively test tool and report on issues

## Performance improvements
- [ ] improve implementation of multi core processing (e.g. using open-mp library)
  - [ ] see how open mp works if we iterate over a map (e.g. in function **doRaytracing**), where we do not iterate over every return but iterate over pulses that store all returns of each pulse.
  - [ ] check if output of multi-core version is the same as single-core version
- [ ] improve progress report on run-time when using open mp (i.e. a single progressbar for entire process)
- [ ] see if open mp would be platform transferable (linux, mac-os) or if there would be better alternatives out
  - [ ] It is not of high priority to provide the tool as a multi-platform tool. But if the efforts to achieve this are within reasonable limits, we should investigate this

## Adding Functionalities
- [ ] Enhance and build upon visualization capabilities of the occlusion outputs. Are there potentially good packages out there?
  - [x] add funcitonality to output figures directly   -> @Jonas
  - [x] Check on output format of figures (png vs pdf) -> @Jonas
- [ ] Implement footprint size for raytracing. Bernhard suggested to look into this notebook for generating subrays within the footprint: https://github.com/3dgeo-heidelberg/helios/blob/implementation_docu/doc/power_sep_interactive.ipynb
- [ ] Possibility to pass a DTM to the voxel traversal in order to stop the traversing of pulses once they reached the terrain
- [ ] *lower priority:* Add possibility to derive plant area density from the traversal output

## Documentation and Notebooks
- [ ] improve and enhance documentation (on code side as well as on readme)
- [ ] setup notebooks showcasing the tool for different lidar sensors (MLS, TLS, UAVLS; multi-return pulses, single return pulses) -> @Jonas
- [ ] add test scripts for single and multi tls scenario
- [ ] integrate lexcube into notebooks
- [ ] 

## Code clean-up
- [ ] there is a lot of potential for cleaning up code   -> @Jonas
- [ ] remove hard-coding of paths etc. (also on setup.py) -> @Jonas

## Packaging
- [ ] How to package and distribute tool via pip and/or conda?
- [ ] Improve structure of the package before distribution
- [ ] how could we add sample data to the distribution package?
- [x] update requirements.txt    -> @Jonas

## Ubuntu
- [ ] make code run on ubuntu -> @William
