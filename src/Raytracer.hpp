//
//  Raytracer.hpp
//  AmanatidesWooAlgorithm_Test
//
//  Created by dkuekenb on 27/04/23.
//  Copyright © 2023 Daniel Kükenbrink. All rights reserved.
//

#define NOMINMAX
#define _USE_MATH_DEFINES

#include <cstring>
#include <iostream>
#include <memory>
#include <iomanip>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include "vector3.h"

#include "Pulse.hpp"

using namespace std;

struct dimensions {
    int nx;
    int ny;
    int nz;
    vector<int> minBound;
    vector<int> maxBound;
    vector<int> origin;
    float voxSize;
};

class Raytracer {

    map<double,boost::shared_ptr<Pulse> > pulsedataset;
    map<double,boost::shared_ptr<Pulse> > incompletePulses;
    map<double,boost::shared_ptr<Pulse> > emptypulsedataset;
    //TODO: Check if these flags are necessary
    int invalidPoints = 0;
    int storedPoints = 0;
    int duplicateReturns = 0;
    int numUnusedReturns = 0;
    /*
    vector<double> X;
    vector<double> Y;
    vector<double> Z;
    vector<double> SensorX;
    vector<double> SensorY;
    vector<double> SensorZ;
    vector<double> GPSTime;
    */

    dimensions gridDim;
    bool hasDTM;

    vector<vector<vector<int > > > Nhit;
    vector<vector<vector<int > > > Nmiss;
    vector<vector<vector<int > > > Nocc;
    vector<vector< double > > DTMind;

    vector<vector< double > > SensorShifts;

    // variables to report back about voxel traversal
    int traversedPulses = 0;            // number of pulses successfully traversed
    int totalPulsesInDataset = 0;       // number of pulses stored and potentially traversed (if voxel grid is intersected)
    int regHit = 0;                     // number of returns registered by the voxel traversal
    int echoesOutside = 0;               // number of returns outside the voxel grid
    int numMissingReturns = 0;           // number of missing returns
    int numNoGridIntersection = 0;      // number of pulses that did not intersect the grid.

    public:
        Raytracer();
        ~Raytracer();

        void rayBoxIntersection (vector<double> origin, vector<double> direction, vector<int> vmin, vector<int> vmax, int & flag, double & tmin);

        void addPointData(vector<double> X, vector<double> Y, vector<double> Z,
                          vector<double> sensor_x, vector<double> sensor_y, vector<double> sensor_z,
                          vector<double> gps_time, vector<int> return_number, vector<int> number_of_returns);


        void doRaytracing();
        void doRaytracing_singleReturnPulses(vector<double> X, vector<double> Y, vector<double> Z, vector<double> sensor_x, vector<double> sensor_y, vector<double> sensor_z, vector<double> gps_time);

        void defineGrid(vector<int> minBound, vector<int> maxBound, int nx, int ny, int nz, float voxSize);

        vector<vector<vector<int > > >& getNhit();
        vector<vector<vector<int > > >& getNmiss();
        vector<vector<vector<int > > >& getNocc();
        vector<vector<double > >& reportSensorShifts();
        vector<int >& getGridDimensions();
        vector<double> getPulsesIntersectingBox(vector<double> x, vector<double> y, vector<double> z, vector<double> sensor_x, vector<double> sensor_y, vector<double> sensor_z, vector<double> gps_time, vector<int> vmin, vector<int> vmax);

        void moveSensorPos2Collinearity();

        bool linePlaneIntersection(vector3& contact, vector3 ray, vector3 rayOrigin, vector3 normal, vector3 coord);

        void getPulseDatasetReport();

        void cleanUpPulseDataset();

        void clearPulseDataset();

        string reportOnTraversal();

        void addEmptyPulseData(vector<double> sensor_x, vector<double> sensor_y, vector<double> sensor_z, vector<double> direction_x, vector<double> direction_y, vector<double> direction_z, vector<double> gps_time);
        void doRaytracingEmptyPulses();
};