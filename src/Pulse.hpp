//
//  Pulse.hpp
//  AmanatidesWooAlgorithm_OO
//
//  Created by Daniel Kükenbrink on 30.10.15.
//  Copyright © 2015 Daniel Kükenbrink. All rights reserved.
//

#ifndef Pulse_hpp
#define Pulse_hpp


#include <cstring>
#include <iostream>
#include <vector>
#include "Echo.hpp"
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>



using namespace std;

class Pulse {
    
    double GPSTime;// = 0;      // unique identifier of each pulse
    double SensorX;// = 0;      // Sensor Position in X coordinates TODO: Decide what to do, when no Sensor position is known
    double SensorY;// = 0;      // Sensor Position in Y coordinates
    double SensorZ;// = 0;      // Sensor Position in Z coordinates
    int NumberOfReturns;// = 0;// Total number of returns in this pulse
    double ScanAngle;// = 0;
    int flightNr;// = 0;

    map<int,boost::shared_ptr<Echo> > Echoes; //A map containing all laser echoes with the return number of each echo as the key

    int missingReturns;// = 99;

    bool wasSensorPos_shifted = false;

    double SensorShift_X;
    double SensorShift_Y;
    double SensorShift_Z;

    // indicator for wether pulse is empty, i.e. has no returns.
    // different then pulse with no echoes! this is to indicate intentional empty pulse
    bool empty;
    // direction is usually calculated from sensorpos to last return, must be given for empty pulse
    double DirectionX;
    double DirectionY;
    double DirectionZ;

    //bool wasTraversed;

    
public:
    bool hasSensorPosition;// = false;

    Pulse (double,double,double,double,int,double,int,map<int,boost::shared_ptr<Echo> >);
    ~Pulse();
    
    Pulse ();
    Pulse(int); //Constructor only taking the number of returns inside this pulse

    Pulse(double, double, double, double, double, double, double) ; 
    //contructor for empty pulse, takes gpstime, sensor position and direction

    void addSensorXPosition (double);
    void addSensorYPosition (double);
    void addSensorZPosition (double);
    void addNumberOfReturns (int);
    void addGPSTime (double);
    void addScanAngle (double);
    void addFlightNumber(int);
    void addEcho(boost::shared_ptr<Echo>);
    void addMultipleEchoes(map<int,boost::shared_ptr<Echo> >);
    
    double getGPSTime(){return this->GPSTime;}
    double getSensorX(){return this->SensorX;}
    double getSensorY(){return this->SensorY;}
    double getSensorZ(){return this->SensorZ;}
    int getNumberOfReturns(){return this->NumberOfReturns;}
    double getScanAngle(){return this->ScanAngle;}
    int getFlightNr(){return this->getFlightNr();}

    void moveSensorPosToLine();

    double getSensorShiftX();
    double getSensorShiftY();
    double getSensorShiftZ();

    double getDirectionX();
    double getDirectionY();
    double getDirectionZ();

    map<int,boost::shared_ptr<Echo> > getEchoes(){return this->Echoes;}
    
    boost::shared_ptr<Echo> getEcho(int);
    
    bool iscomplete();

    void cleanupPulse();

};


#endif /* Pulse_hpp */
