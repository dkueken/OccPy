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

        void addPointData(vector[double], vector[double], vector[double], vector[double], vector[double], vector[double], vector[double], vector[int], vector[int])

        void doRaytracing()

        void rayBoxIntersection(vector[double], vector[double], vector[int], vector[int], int &, double &);

        void cleanUpPulseDataset();

        #get functions
        vector[vector[vector[int]]]& getNhit()
        vector[vector[vector[int]]]& getNmiss()
        vector[vector[vector[int]]]& getNocc()
        vector[int]& getGridDimensions()
        void getPulseDatasetReport()


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
    def addPointData(self, X, Y, Z, sensor_x, sensor_y, sensor_z, gps_time, return_number, number_of_returns):
        self.thisptr.addPointData(X, Y, Z, sensor_x, sensor_y, sensor_z, gps_time, return_number, number_of_returns)
    def doRaytracing(self):
        self.thisptr.doRaytracing()
    def cleanUpPulseDataset(self):
        self.thisptr.cleanUpPulseDataset()
    def getNhit(self):
        return self.thisptr.getNhit()
    def getNmiss(self):
        return self.thisptr.getNmiss()
    def getNocc(self):
        return self.thisptr.getNocc()
    def getGridDimensions(self):
        return self.thisptr.getGridDimensions()
    def getPulseDatasetReport(self):
        self.thisptr.getPulseDatasetReport()
