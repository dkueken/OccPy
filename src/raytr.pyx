# distutils: language = c++

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

        void defineGrid(vector[double], vector[double], int, int, int, float)

        void addPointData(vector[double], vector[double], vector[double], vector[double], vector[double], vector[double], vector[double], vector[int], vector[int])

        void doRaytracing()
        void doRaytracing_singleReturnPulses(vector[double] X, vector[double] Y, vector[double] Z, vector[double] sensor_x, vector[double] sensor_y, vector[double] sensor_z, vector[double] gps_time)

        void rayBoxIntersection(vector[double], vector[double], vector[double], vector[double], int &, double &)

        void cleanUpPulseDataset()

        void clearPulseDataset()

        void moveSensorPos2Collinearity()

        string reportOnTraversal()

        #get functions
        vector[vector[vector[int]]]& getNhit()
        vector[vector[vector[int]]]& getNmiss()
        vector[vector[vector[int]]]& getNocc()
        vector[int] getGridDimensions()
        vector[double] getGridOrigin()
        void getPulseDatasetReport()
        vector[vector[double]]& reportSensorShifts()
        vector[double] getPulsesIntersectingBox(vector[double], vector[double], vector[double], vector[double], vector[double], vector[double], vector[double], vector[double], vector[double])

        void addEmptyPulseData(vector[double] sensor_x, vector[double] sensor_y, vector[double] sensor_z, vector[double] direction_x, vector[double] direction_y, vector[double] direction_z, vector[double] gps_time)
        void doRaytracingEmptyPulses()

        int get_num_traversed_pulses() const
        int get_total_pulses_in_dataset() const
        int get_num_registered_hits() const
        int get_num_echoes_outside() const
        int get_num_missing_returns() const
        int get_num_pulses_no_intersection() const



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
    def doRaytracing_singleReturnPulses(self, X, Y, Z, sensor_x, sensor_y, sensor_z, gps_time):
        self.thisptr.doRaytracing_singleReturnPulses(X, Y, Z, sensor_x, sensor_y, sensor_z, gps_time)
    def cleanUpPulseDataset(self):
        self.thisptr.cleanUpPulseDataset()
    def clearPulseDataset(self):
        self.thisptr.clearPulseDataset()
    def moveSensorPos2Collinearity(self):
        self.thisptr.moveSensorPos2Collinearity()
    def reportSensorShifts(self):
        return self.thisptr.reportSensorShifts()
    def getNhit(self):
        return self.thisptr.getNhit()
    def getNmiss(self):
        return self.thisptr.getNmiss()
    def getNocc(self):
        return self.thisptr.getNocc()
    def getGridDimensions(self):
        return self.thisptr.getGridDimensions()
    def getGridOrigin(self):
        return self.thisptr.getGridOrigin()
    def getPulseDatasetReport(self):
        self.thisptr.getPulseDatasetReport()
    def reportOnTraversal(self):
        return self.thisptr.reportOnTraversal()
    def getPulsesIntersectingBox(self, x, y, z, sensor_x, sensor_y, sensor_z, gps_time, vmin, vmax):
        return self.thisptr.getPulsesIntersectingBox(x, y, z, sensor_x, sensor_y, sensor_z, gps_time, vmin, vmax)
    def addEmptyPulseData(self, sensor_x, sensor_y, sensor_z, direction_x, direction_y, direction_z, gps_time):
        return self.thisptr.addEmptyPulseData(sensor_x, sensor_y, sensor_z, direction_x, direction_y, direction_z, gps_time)
    def doRaytracingEmptyPulses(self):
        return self.thisptr.doRaytracingEmptyPulses()
    def get_num_traversed_pulses(self) -> int:
        return self.thisptr.get_num_traversed_pulses()
    def get_total_pulses_in_dataset(self) -> int:
        return self.thisptr.get_total_pulses_in_dataset()
    def get_num_registered_hits(self) -> int:
        return self.thisptr.get_num_registered_hits()
    def get_num_echoes_outside(self) -> int:
        return self.thisptr.get_num_echoes_outside()
    def get_num_missing_returns(self) -> int:
        return self.thisptr.get_num_missing_returns()
    def get_num_pulses_no_intersection(self) -> int:
        return self.thisptr.get_num_pulses_no_intersection()


