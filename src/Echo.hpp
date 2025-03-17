//
//  Echo.hpp
//  AmanatidesWooAlgorithm_OO
//
//  Created by dkuekenb on 12/11/15.
//  Copyright © 2015 Daniel Kükenbrink. All rights reserved.
//

#ifndef Echo_hpp
#define Echo_hpp

#include <cstring>
#include <iostream>

using namespace std;


class Echo {

    double x;
    double y;
    double z;
    int intensity;
    int returnNumber;

public:

    Echo ();
    Echo (double, double, double, int, int);


    Echo* getEcho();

    void setX(double);
    double getX();
    void setY(double);
    double getY();
    void setZ(double);
    double getZ();
    void setIntensity(int);
    int getIntensity();
    void setReturnNumber(int);
    int getReturnNumber();

};
#endif /* Echo_hpp */
