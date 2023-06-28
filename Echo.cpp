//
//  Echo.cpp
//  AmanatidesWooAlgorithm_OO
//
//  Created by dkuekenb on 12/11/15.
//  Copyright © 2015 Daniel Kükenbrink. All rights reserved.
//

#include "Echo.hpp"

Echo::Echo() {
    this->x = 0;
    this->y = 0;
    this->z = 0;
    this->intensity = 0;
    this->returnNumber = 0;
    
}

Echo::Echo(double xEcho, double yEcho, double zEcho, int intEcho, int retNum){
    this->x = xEcho;
    this->y = yEcho;
    this->z = zEcho;
    this->intensity = intEcho;
    this->returnNumber = retNum;
}

Echo* Echo::getEcho(){
    
    return this;
}

void Echo::setX(double xEcho) {
    this->x = xEcho;
}

void Echo::setY(double yEcho) {
    this->y = yEcho;
}

void Echo::setZ(double zEcho) {
    this->z = zEcho;
}

void Echo::setIntensity(int intEcho) {
    this->intensity = intEcho;
}

void Echo::setReturnNumber(int retNum) {
    this->returnNumber = retNum;
}

double Echo::getX() {
    return this->x;
}

double Echo::getY() {
    return this->y;
}

double Echo::getZ() {
    return this->z;
}

int Echo::getIntensity() {
    return this->intensity;
}

int Echo::getReturnNumber() {
    return this->returnNumber;
}
