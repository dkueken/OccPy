# distutils: language = c++
# distutils: sources = Raytracer.cpp

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

# c++ interface to cython
cdef extern from "Raytracer.hpp":
  cdef cppclass Raytracer:
        Raytracer() except +
        vector[vector[vector[int]]] Nhit
        vector[vector[vector[int]]] Nmiss
        vector[vector[vector[int]]] Nocc

        void defineGrid(vector[int], vector[int], int, int, int, float)

        void doRaytracing(vector[double] x, vector[double] y, vector[double] z, vector[double] sensor_x, vector[double] sensor_y, vector[double] sensor_z, vector[double] gps_time)

        void rayBoxIntersection(vector[double], vector[double], vector[int], vector[int], int &, double &);

        #get functions
        vector[vector[vector[int]]]& getNhit()
        vector[vector[vector[int]]]& getNmiss()
        vector[vector[vector[int]]]& getNocc()
        vector[int]& getGridDimensions()


# creating a cython wrapper class
cdef class PyRaytracer:
    cdef Raytracer *thisptr     # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = new Raytracer()
    def __dealloc__(self):
        del self.thisptr
    def rayBoxIntersection(self, origin,  direction, vmin, vmax, flag, tmin):
        self.thisptr.rayBoxIntersection(origin, direction, vmin, vmax, flag, tmin)
    def defineGrid(self, minBound, maxBound, nx, ny, nz, voxSize):
        self.thisptr.defineGrid(minBound, maxBound, nx, ny, nz, voxSize)
    def doRaytracing(self, X, Y, Z, sensor_x, sensor_y, sensor_z, gps_time):
        self.thisptr.doRaytracing(X, Y, Z, sensor_x, sensor_y, sensor_z, gps_time)
    def getNhit(self):
        return self.thisptr.getNhit()
    def getNmiss(self):
        return self.thisptr.getNmiss()
    def getNocc(self):
        return self.thisptr.getNocc()
    def getGridDimensions(self):
        return self.thisptr.getGridDimensions()
