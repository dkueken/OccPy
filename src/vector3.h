//
//  vector3.h
//
//  Created by Daniel Kükenbrink on 10.05.18.
//  Copyright © 2018 Daniel Kükenbrink. All rights reserved.
//

#ifndef vector3_h
#define vector3_h

class vector3 {

    double x;
    double y;
    double z;

public:
    vector3(double newX, double newY, double newZ) {
        this->x = newX;
        this->y = newY;
        this->z = newZ;
    }

    vector3() {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    ~vector3() {}
    
    double dot(const vector3& other)
    {
        return this->x * other.x + this->y * other.y + this->z * other.z;
    }
    
    double getX(){
        return this->x;
    }
    double getY(){
        return this->y;
    }
    double getZ(){
        return this->z;
    }
    
    void setX(double newX){
        this->x = newX;
    }
    void setY(double newY){
        this->y = newY;
    }
    void setZ(double newZ){
        this->z = newZ;
    }
private:

};


#endif /* vector3_h*/
