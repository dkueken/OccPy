
# Files to check:

## cpp side: 

[x] Raytracer.cpp+hpp

output variables order changed:
    - Nhit
    - Nmiss
    - Nocc

[x] Pulse.cpp+hpp  
[x] Echo.cpp+hpp  
[x] vector3.h  
[x] riegl_io (don't think any changes needed)  

[x] raytr.pyx

## python side

[] Occpy.py  
[] occpyriegl.py  
[] prepare_trajectory.py  
[] prepareply.py  
[] riegl_io.py  
[] terrainmodel.py  
[] visualization.py  

## notebooks

[] MLS  
[] TLS  
[] UAVLS  


TODO checks:

 - rayboxintersection: check (numNoGridIntersection) before and after flip
 -> use getPulsesIntersectingBox() function



