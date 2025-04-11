//
//  Pulse.cpp
//  AmanatidesWooAlgorithm_OO
//
//  Created by Daniel Kükenbrink on 30.10.15.
//  Copyright © 2015 Daniel Kükenbrink. All rights reserved.
//

#include "Pulse.hpp"

Pulse::Pulse (){
    
    this->hasSensorPosition = false;
//    this->Echoes = map<int,boost::shared_ptr<Echo> >;
    //this->wasTraversed = false;
    
    
}

Pulse::Pulse (double gps, double sensX, double sensY, double sensZ, int numRet, double scAngle, int fNum, map<int,boost::shared_ptr<Echo> > echoes) {
    // This consturctor assumes that only complete pulses (without missing returns) are inserted!!
    this->GPSTime = gps;
    this->SensorX = sensX;
    this->SensorY = sensY;
    this->SensorZ = sensZ;
    this->NumberOfReturns = numRet;
    this->ScanAngle = scAngle;
    this->flightNr = fNum;
    
    this->Echoes = echoes;
    
    
    this->missingReturns = 0;
    
    this->hasSensorPosition = true;

    // note: empty does not indicate a normal pulse with missing returns, but an intentionally empty pulse
    // never set this based on echo map
    this->empty = false;

    //this->wasTraversed = false;

    
    
}

Pulse::Pulse(int numRets) {
    
    this->missingReturns = numRets;
    
//    this->Echoes = new map<int,Echo*>();
    
    this->hasSensorPosition = false;

    //this->wasTraversed = false;

    this->empty = false;
    
}

Pulse::Pulse(double gps, double sensX, double sensY, double sensZ, double directionX, double directionY, double directionZ){
    this->GPSTime = gps;
    this->SensorX = sensX;
    this->SensorY = sensY;
    this->SensorZ = sensZ;
    this->DirectionX = directionX;
    this->DirectionY = directionY;
    this->DirectionZ = directionZ;

    this->empty = true;
}


Pulse::~Pulse() {

    if (!this->empty) {
        map<int,boost::shared_ptr<Echo> >::iterator it = this->Echoes.begin();
        while ( it != this->Echoes.end() ){
            //Iterate thorugh all Echoes and delete them
            this->Echoes.erase(it++);
        }
    }
    
}

void Pulse::cleanupPulse(){

    map<int,boost::shared_ptr<Echo> > tmp;

    map<int,boost::shared_ptr<Echo> >::iterator it = this->Echoes.begin();

    int retNum = 1;
    while ( it != this->Echoes.end() ){
        //Iterate through all echoes and change the return number in order to get a clear order from 1 to NumberOfReturns
        boost::shared_ptr<Echo> e = it->second;
        e->setReturnNumber(retNum);
        tmp.insert( std::make_pair(e->getReturnNumber(),e));
        retNum++;
        it++;
    }


    this->Echoes = tmp;
    this->NumberOfReturns = this->Echoes.size();
    this->missingReturns = 0;

}

/* function to move selected point towards line */
void Pulse::moveSensorPosToLine(){

    // Calculate the line direction vector
    vector<double> lineDirection(3);
    lineDirection[0] = this->Echoes.at(this->NumberOfReturns)->getX() - this->Echoes.at(1)->getX();
    lineDirection[1] = this->Echoes.at(this->NumberOfReturns)->getY() - this->Echoes.at(1)->getY();
    lineDirection[2] = this->Echoes.at(this->NumberOfReturns)->getZ() - this->Echoes.at(1)->getZ();

    // Calculate the point's projection onto the line
    double dotProduct = 0;
    double lineMagnitudeSquared = 0;
    dotProduct += (this->SensorX - this->Echoes.at(1)->getX()) * lineDirection[0];
    dotProduct += (this->SensorY - this->Echoes.at(1)->getY()) * lineDirection[1];
    dotProduct += (this->SensorZ - this->Echoes.at(1)->getZ()) * lineDirection[2];

    for (int i = 0; i < 3; i++){
        lineMagnitudeSquared += lineDirection[i] * lineDirection[i];
    }
    double t = dotProduct / lineMagnitudeSquared;

    double proj_sensor_x = this->Echoes.at(1)->getX() + t * lineDirection[0];
    double proj_sensor_y = this->Echoes.at(1)->getY() + t * lineDirection[1];
    double proj_sensor_z = this->Echoes.at(1)->getZ() + t * lineDirection[2];

    this->SensorShift_X = proj_sensor_x - this->SensorX;
    this->SensorShift_Y = proj_sensor_y - this->SensorY;
    this->SensorShift_Z = proj_sensor_z - this->SensorZ;
    this->SensorX = proj_sensor_x;
    this->SensorY = proj_sensor_y;
    this->SensorZ = proj_sensor_z;
    this->wasSensorPos_shifted = true;

}

void Pulse::addGPSTime(double gps) {
    this->GPSTime = gps;
}

void Pulse::addNumberOfReturns(int numRet) {
    this->NumberOfReturns = numRet;
}


void Pulse::addScanAngle(double scAngle) {
    this->ScanAngle = scAngle;
}

void Pulse::addSensorXPosition(double sensX) {
    this->SensorX = sensX;
    this->hasSensorPosition = true;
}

void Pulse::addSensorYPosition(double sensY) {
    this->SensorY = sensY;
    this->hasSensorPosition = true;
}

void Pulse::addSensorZPosition(double sensZ) {
    this->SensorZ = sensZ;
    this->hasSensorPosition = true;
}

void Pulse::addFlightNumber(int fNum) {
    this->flightNr = fNum;
}

void Pulse::addEcho(boost::shared_ptr<Echo> ret){
    
    this->Echoes.insert( std::make_pair(ret->getReturnNumber(),ret) );
    
    
    this->missingReturns--;
    
}

void Pulse::addMultipleEchoes(map<int,boost::shared_ptr<Echo> > echoes) {
    this->Echoes.insert(echoes.begin(),echoes.end());
    
    this->missingReturns = this->missingReturns - echoes.size();
}

boost::shared_ptr<Echo> Pulse::getEcho(int retNum) {
    return this->Echoes.at(retNum);
}

double Pulse::getSensorShiftX(){
    if (this->wasSensorPos_shifted){
        return this->SensorShift_X;
    } else {
        return 0;
    }
}

double Pulse::getSensorShiftY(){
    if (this->wasSensorPos_shifted){
        return this->SensorShift_Y;
    } else {
        return 0;
    }
}

double Pulse::getSensorShiftZ(){
    if (this->wasSensorPos_shifted){
        return this->SensorShift_Z;
    } else {
        return 0;
    }
}


double Pulse::getDirectionX(){
    if (this->empty) {
        return this->DirectionX;
    } else {
        //calculate based on last return
        return this->getEchoes().at(this->getNumberOfReturns())->getX() - this->getSensorX();
    }
}

double Pulse::getDirectionY(){
    if (this->empty) {
        return this->DirectionY;
    } else {
        //calculate based on last return
        return this->getEchoes().at(this->getNumberOfReturns())->getY() - this->getSensorY();
    }
}

double Pulse::getDirectionZ(){
    if (this->empty) {
        return this->DirectionZ;
    } else {
        //calculate based on last return
        return this->getEchoes().at(this->getNumberOfReturns())->getZ() - this->getSensorZ();
    }
}


bool Pulse::iscomplete(){
    
    
    unsigned long numEchoes = this->Echoes.size();
    bool flag = false;
    if (numEchoes==this->NumberOfReturns) {
        
        // now check if all needed keys exist
        for (int k=1; k<=this->NumberOfReturns; k++) {
            if (this->Echoes.find(k) == this->Echoes.end()) {
                flag = false;
                break;
            } else {
                if (k==this->NumberOfReturns) {
                    flag = true;    //Now every needed key is available, and the pulse should be complete now!
                } else {
                    flag = false;
                }
            }
        }
       
    } else {
        flag = false;
    }
    
    return flag;
        
}

