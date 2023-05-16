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

    
    
}

Pulse::Pulse(int numRets) {
    
    this->missingReturns = numRets;
    
//    this->Echoes = new map<int,Echo*>();
    
    this->hasSensorPosition = false;
    
}


Pulse::~Pulse() {
    
    
//    cout << "Deallocating Pulse object " << endl;
    
    map<int,boost::shared_ptr<Echo> >::iterator it = this->Echoes.begin();
    
    while ( it != this->Echoes.end() ){
        //Iterate thorugh all Echoes and delete them
        this->Echoes.erase(it++);
//        delete it->second;
//        this->Echoes->erase(it);
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

