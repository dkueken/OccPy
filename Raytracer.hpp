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
#include <map>
#include "vector3.h"

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

    vector<double> X;
    vector<double> Y;
    vector<double> Z;
    vector<double> SensorX;
    vector<double> SensorY;
    vector<double> SensorZ;
    vector<double> GPSTime;

    dimensions gridDim;
    bool hasDTM;

    vector<vector<vector<int > > > Nhit;
    vector<vector<vector<int > > > Nmiss;
    vector<vector<vector<int > > > Nocc;
    vector<vector< double > > DTMind;

    public:
        Raytracer();
        ~Raytracer();

        void rayBoxIntersection (vector<double> origin, vector<double> direction, vector<int> vmin, vector<int> vmax, int & flag, double & tmin);

        //void readMLSData(vector<double> x, vector<double> y, vector<double> z, vector<double> sensor_x, vector<double> sensor_y, vector<double> sensor_z, vector<double> gps_time);

        //void readDTM(vector<vectory double > > DTM);

        void doRaytracing(vector<double> X, vector<double> Y, vector<double> Z, vector<double> sensor_x, vector<double> sensor_y, vector<double> sensor_z, vector<double> gps_time);

        void defineGrid(vector<int> minBound, vector<int> maxBound, int nx, int ny, int nz, float voxSize);

        vector<vector<vector<int > > >& getNhit();
        vector<vector<vector<int > > >& getNmiss();
        vector<vector<vector<int > > >& getNocc();
        vector<int >& getGridDimensions();

        bool linePlaneIntersection(vector3& contact, vector3 ray, vector3 rayOrigin, vector3 normal, vector3 coord);
};