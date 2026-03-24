//
// Raytracer.cpp
//
// Created by Daniel Kükenbrink on 27.04.2023
// Copyright © 2023 Daniel Kükenbrink. All rights reserved.
//

#include "Raytracer.hpp"

Raytracer::Raytracer()
{
    this->hasDTM = false;
}

Raytracer::~Raytracer()
{
    //TODO: add a verbose flag that outputs these kinds of statements!
    //cout << "Cleaning up and destroy Raytracer object" << endl;

    vector<vector<vector<int > > >().swap(this->Nhit);
    vector<vector<vector<int > > >().swap(this->Nmiss);
    vector<vector<vector<int > > >().swap(this->Nocc);
    /*
    vector<double >().swap(this->X);
    vector<double >().swap(this->Y);
    vector<double >().swap(this->Z);
    vector<double >().swap(this->SensorX);
    vector<double >().swap(this->SensorY);
    vector<double >().swap(this->SensorZ);
    vector<double >().swap(this->GPSTime);
    */

}

void Raytracer::defineGrid(vector<double> minBound, vector<double> maxBound, int nx, int ny, int nz, float voxSize){
    //Test without adapting the grid origin: Meaning, we will take the minBound vector as the grid origin
    this->gridDim.origin = minBound;

    this->gridDim.minBound = minBound;
    this->gridDim.maxBound = maxBound;
    this->gridDim.voxSize = voxSize;

    this->gridDim.nx = nx;
    this->gridDim.ny = ny;
    this->gridDim.nz = nz;

    this->Nhit.resize(this->gridDim.nx, vector<vector<int > >(this->gridDim.ny, vector<int>(this->gridDim.nz,0)));
    this->Nmiss.resize(this->gridDim.nx, vector<vector<int > >(this->gridDim.ny, vector<int>(this->gridDim.nz,0)));
    this->Nocc.resize(this->gridDim.nx, vector<vector<int > >(this->gridDim.ny, vector<int>(this->gridDim.nz,0)));

}

static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
    if ( (x != n) && (x % (n/100+1) != 0) ) return;

    float ratio  =  x/(float)n;
    int   c      =  ratio * w;

    cout << std::setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
}

void Raytracer::addPointData(vector<double> X, vector<double> Y, vector<double> Z,
                             vector<double> sensor_x, vector<double> sensor_y, vector<double> sensor_z,
                             vector<double> gps_time, vector<int> return_number, vector<int> number_of_returns){


    for (int i = 0; i < gps_time.size(); i++) {

        map<double,boost::shared_ptr<Pulse> >::iterator it;

        this->storedPoints++; //TODO: this variable should probably be changed and maybe added as a parameter of the Pulse dataset?

        double gps = gps_time.at(i);

        it = this->incompletePulses.find(gps);

        if (it != this->incompletePulses.end()) {
            // pulse is already in pulsedataset and echo can be added. -> extract laser return information
            boost::shared_ptr<Echo> e = boost::make_shared<Echo>();

            e->setX(X.at(i));
            e->setY(Y.at(i));
            e->setZ(Z.at(i));
            e->setReturnNumber(return_number.at(i));

            // Check if an echoe with this return number already exists for this pulse (There are cases where there are duplicate returns (the same ones)
            if (it->second->getEchoes().count(e->getReturnNumber())==1){
                this->duplicateReturns++;
            }

            it->second->addEcho( e );

            if (it->second->iscomplete()) {

                this->pulsedataset.insert( std::make_pair(it->second->getGPSTime(), it->second) );
                this->incompletePulses.erase(it);

            }
        } else {
            // Pulse is not yet within the pulsedataset and should be added, before adding the echo

            // TODO: check if pulse already in completePulses set 
            // (so we can warn user that gps_time has duplicates which overwrites the pulses that are already in the set)


            // extract general pulse information
            boost::shared_ptr<Pulse> p = boost::make_shared<Pulse>();

            p->addGPSTime(gps_time.at(i));
            p->addNumberOfReturns(number_of_returns.at(i));
            p->addSensorXPosition(sensor_x.at(i));
            p->addSensorYPosition(sensor_y.at(i));
            p->addSensorZPosition(sensor_z.at(i));

            // extract laser return information
            boost::shared_ptr<Echo> e = boost::make_shared<Echo>();

            e->setX(X.at(i));
            e->setY(Y.at(i));
            e->setZ(Z.at(i));
            e->setReturnNumber(return_number.at(i));

            p->addEcho( e );

            if (p->iscomplete()) {
                this->pulsedataset.insert(std::make_pair(p->getGPSTime(),p));
            } else {
                this->incompletePulses.insert( std::make_pair(p->getGPSTime(), p) );
            }

        }


    }

}

void Raytracer::moveSensorPos2Collinearity(){

    this->SensorShifts.resize((int)this->pulsedataset.size(), vector<double>(3, 0) );

    // iterate over pulsedataset
    int p = 0;
    for (map<double,boost::shared_ptr<Pulse> >::iterator it = this->pulsedataset.begin(); it != this->pulsedataset.end(); ++it){
        // only move sensor position if we have at least 2 returns
        if (it->second->getNumberOfReturns() >= 2){
            it->second->moveSensorPosToLine();
            this->SensorShifts.at(p).at(0) = it->second->getSensorShiftX();
            this->SensorShifts.at(p).at(1) = it->second->getSensorShiftY();
            this->SensorShifts.at(p).at(2) = it->second->getSensorShiftZ();
            // cout << "### Shifted Sensor by X: " << sensor_shift.at(p).at(0) << " Y: " << sensor_shift.at(p).at(1) << " Z: " << sensor_shift.at(p).at(2) << " m" << endl;
        } else {
            this->SensorShifts.at(p).at(0) = 0;
            this->SensorShifts.at(p).at(1) = 0;
            this->SensorShifts.at(p).at(2) = 0;
        }

        p++;
    }
}

void Raytracer::getPulseDatasetReport(){
    cout << "#### Pulse Dataset Report ####" << endl;
    cout << "### " << this->storedPoints << " Returns have been read and stored! " << endl;
    cout << "### " << this->duplicateReturns << " duplicate returns have been recognized! " << endl;
    cout << "### " << this->pulsedataset.size() << " pulses are stored in the dataset! " << endl;
    cout << "### " << this->incompletePulses.size() << " incomplete pulses are not going to be analysed! " << endl;
}

/*
void Raytracer::getIncompletePulseDatasetReport(){

    //TODO: implement this function!
    cout << "#### Incomplete Pulse Dataset Report ####" << endl;

    int numMissingFirst = 0;
    int numMissingSecond = 0;
    int numMissingThird = 0;
    int numMissingFourth = 0;

    for (map<double,boost::shared_ptr<Pulse> >::iterator it = this->incompletePulses.begin(); it != this->incompletePulses.end(); ++it) {


    }
}
*/

void Raytracer::cleanUpPulseDataset(){
    //This function should clean up the incomplete pulses in order to include them within the voxel traversal

    // Iterate over incomplete pulses and clean up
    for (auto it = this->incompletePulses.begin(); it != this->incompletePulses.cend();){
        it->second->cleanupPulse();

        // check if pulse is now complete
        if (it->second->iscomplete()){
            this->pulsedataset.insert(std::make_pair(it->second->getGPSTime(), it->second));
            it = this->incompletePulses.erase(it);
        } else {
            ++it;
        }
    }
}

void Raytracer::doRaytracing()
{
    int pulsecount = 0;
    int traversedPulses = 0;
    int echoesOutside = 0;
    int numMissingReturns = 0;
    int numNoGridIntersection = 0;

    int total_pulses2iterate = (int)this->pulsedataset.size();

    //TODO: implement functionality for single return pulses with no known sensor positions

    int regHit = 0; //Number of registered hit voxels

    //iterating through each pulse, initialise traversing of the pulse and traverse the pulse through the voxel grid
    map<double,boost::shared_ptr<Pulse> >::iterator it = this->pulsedataset.begin();
    //while (it != this->pulsedataset.end()){
    for (map<double,boost::shared_ptr<Pulse> >::iterator it = this->pulsedataset.begin(); it != this->pulsedataset.end(); ++it) {
        pulsecount++;

        bool isfound = false;

        //update progressbar
        loadbar(pulsecount, total_pulses2iterate, 20);

        //init several necessary variables befor starting
        int flag = 0;               //flag whether pulse crosses voxel grid or not
        double tmin = 0;             //how far along the pulse direction we have to travel until we reach the voxel grid. based on the vector representation of a line: y = u + tmin*v (y,u and v are vectors)

        double tVoxelX = 0;
        double tVoxelY = 0;
        double tVoxelZ = 0;
        int stepX = 0;
        int stepY = 0;
        int stepZ = 0;

        double voxelMaxX = 0;
        double voxelMaxY = 0;
        double voxelMaxZ = 0;

        double voxelSizeX = 0;
        double voxelSizeY = 0;
        double voxelSizeZ = 0;

        double tDeltaX = 0;
        double tDeltaY = 0;
        double tDeltaZ = 0;

        double tMaxX = 0;
        double tMaxY = 0;
        double tMaxZ = 0;

        vector<double> origin (3,0);
        vector<double> direction (3,0);

        if (it->second->hasSensorPosition) {
            origin.at(0) = it->second->getSensorX();
            origin.at(1) = it->second->getSensorY();
            origin.at(2) = it->second->getSensorZ();

            //TODO: check if we should have to subtract 1 from the Number of Returns, because the vector indices are also 0 based
            //TODO: this should work, as the getEchoes() function returns a map with a shared ptr.
            direction.at(0) = it->second->getEchoes().at(it->second->getNumberOfReturns())->getX() - origin.at(0);
            direction.at(1) = it->second->getEchoes().at(it->second->getNumberOfReturns())->getY() - origin.at(1);
            direction.at(2) = it->second->getEchoes().at(it->second->getNumberOfReturns())->getZ() - origin.at(2);

            /*
            //TODO: Check if this is realy necessary. One small test showed, that not many pulses are actually affected by this effect!!!
            // This algorithm works on the assumption of an infinitesimally small pulse. It is therefore crucial that all laser returns lie exactly on top of the laser pulse line, otherwise it could happen that the traversal algorithm will not traverse a laser echo because it lies in another voxel which was not traversed. We therefore project all echoes for pulses with more than 1 laser return onto the vector which represents the laser pulse
            if (it->second->getNumberOfReturns()>1) {

                for (int r=1; r<it->second->getNumberOfReturns(); r++) {
                    double pX = it->second->getEcho(r)->getX();
                    double pY = it->second->getEcho(r)->getY();
                    double pZ = it->second->getEcho(r)->getZ();

                    //difference of point with line origin
                    double dy = pY - origin.at(0);
                    double dx = pX - origin.at(1);
                    double dz = pZ - origin.at(2);

                    // Position of projection on line using dot product
                    double delta = direction.at(1)*direction.at(1) + direction.at(0)*direction.at(0) + direction.at(2)*direction.at(2);
                    double tp = (dx*direction.at(1) + dy*direction.at(0) + dz*direction.at(2))/delta;

                    // convert position on line to cartesian coordinates
                    double pYnew = origin.at(0)+tp*direction.at(0);
                    double pXnew = origin.at(1)+tp*direction.at(1);
                    double pZnew = origin.at(2)+tp*direction.at(2);

                    //change coordinates of the Echo
                    it->second->getEcho(r)->setX(pXnew);
                    it->second->getEcho(r)->setY(pYnew);
                    it->second->getEcho(r)->setZ(pZnew);

                }
            }
            */
        } else {
            //TODO: implement traversal if sensor position is not known -> this should not be promoted as we should advocate
            // that the sensor position should always be known!
        }

        // if direction vector is == 0 assign a very small number to overcome problems with division through 0 in later steps. TODO: Check if this value is small enough. This could be a problem especially with TLS data!
        if (direction.at(0)==0) { direction.at(0) = numeric_limits<double>::min(); }
        if (direction.at(1)==0) { direction.at(1) = numeric_limits<double>::min(); }
        if (direction.at(2)==0) { direction.at(2) = numeric_limits<double>::min(); }

        //check if pulse intersects voxel grid
        rayBoxIntersection(origin, direction, this->gridDim.minBound, this->gridDim.maxBound, flag, tmin);

        if (flag==0) {
            numNoGridIntersection++;
            // pulse not relevant remove from pulse dataset
            //it = this->pulsedataset.erase(it);
            //cout << "The ray does not intersect the grid" << endl;
        } else {
            //cout << "The ray does intersect the grid" << endl;

            traversedPulses++;

            int echoesOfPulseOutside = 0;

            if (tmin<0) {
                tmin = 0;
            }

            vector<double> start (3,0);
            start.at(0) = origin.at(0) + (tmin*direction.at(0)); //+ tmin*direction.at(1) + tmin*direction.at(2);
            start.at(1) = origin.at(1) + (tmin*direction.at(1)); //+ tmin*direction.at(1) + tmin*direction.at(2);
            start.at(2) = origin.at(2) + (tmin*direction.at(2)); //+ tmin*direction.at(1) + tmin*direction.at(2);

            vector<double> boxSize (3,0);
            boxSize.at(0) = (double)this->gridDim.maxBound.at(0) - (double)this->gridDim.minBound.at(0);
            boxSize.at(1) = (double)this->gridDim.maxBound.at(1) - (double)this->gridDim.minBound.at(1);
            boxSize.at(2) = (double)this->gridDim.maxBound.at(2) - (double)this->gridDim.minBound.at(2);

            // Get the voxel indices of where the laser return lie in the grid
            vector<int> XInd;
            XInd.resize(it->second->getNumberOfReturns());
            vector<int> YInd;
            YInd.resize(it->second->getNumberOfReturns());
            vector<int> ZInd;
            ZInd.resize(it->second->getNumberOfReturns());

            for (int j=1; j<=it->second->getNumberOfReturns(); j++){
                XInd.at(j-1) = floor( ((it->second->getEchoes().at(j)->getX() - (double)this->gridDim.minBound.at(0))/boxSize.at(0))*(double)this->gridDim.nx );
                YInd.at(j-1) = floor( ((it->second->getEchoes().at(j)->getY() - (double)this->gridDim.minBound.at(1))/boxSize.at(1))*(double)this->gridDim.ny );
                ZInd.at(j-1) = floor( ((it->second->getEchoes().at(j)->getZ() - (double)this->gridDim.minBound.at(2))/boxSize.at(2))*(double)this->gridDim.nz );

                if (YInd.at(j-1)>=this->gridDim.ny || YInd.at(j-1)<0 || XInd.at(j-1)>=this->gridDim.nx || XInd.at(j-1)<0 ||
                    ZInd.at(j-1)>=this->gridDim.nz || ZInd.at(j-1)<0){

                    echoesOutside++;
                    echoesOfPulseOutside++;

                }
            }

            int x;
            int y;
            int z;

            x = floor( ((start.at(0)-(double)this->gridDim.minBound.at(0))/boxSize.at(0))*(double)this->gridDim.nx );
            y = floor( ((start.at(1)-(double)this->gridDim.minBound.at(1))/boxSize.at(1))*(double)this->gridDim.ny );
            z = floor( ((start.at(2)-(double)this->gridDim.minBound.at(2))/boxSize.at(2))*(double)this->gridDim.nz );

            if (x==(this->gridDim.nx)) {
                x = x-1;
            }
            if (y==(this->gridDim.ny)){
                y = y-1;
            }
            if (z==(this->gridDim.nz)){
                z = z-1;
            }

            if (direction.at(0)>=0) {
                tVoxelX = (double)(x+1)/(double)this->gridDim.nx;
                stepX = 1;
            } else {
                tVoxelX = ((double)(x+1)-1)/(double)this->gridDim.nx;
                stepX = -1;
            }
            if (direction.at(1)>=0) {
                tVoxelY = (double)(y+1)/(double)this->gridDim.ny;
                stepY = 1;
            } else {
                tVoxelY = ((double)(y+1)-1)/(double)this->gridDim.ny;
                stepY = -1;
            }
            if (direction.at(2)>=0) {
                tVoxelZ = (double)(z+1)/(double)this->gridDim.nz;
                stepZ = 1;
            } else {
                tVoxelZ = ((double)(z+1)-1)/(double)this->gridDim.nz;
                stepZ = -1;
            }

            voxelMaxX = (double)this->gridDim.minBound.at(0) + tVoxelX*boxSize.at(0);
            voxelMaxY = (double)this->gridDim.minBound.at(1) + tVoxelY*boxSize.at(1);
            voxelMaxZ = (double)this->gridDim.minBound.at(2) + tVoxelZ*boxSize.at(2);

            tMaxX = tmin + (voxelMaxX-start.at(0))/direction.at(0);
            tMaxY = tmin + (voxelMaxY-start.at(1))/direction.at(1);
            tMaxZ = tmin + (voxelMaxZ-start.at(2))/direction.at(2);

            voxelSizeX = boxSize.at(0)/(double)this->gridDim.nx;
            voxelSizeY = boxSize.at(1)/(double)this->gridDim.ny;
            voxelSizeZ = boxSize.at(2)/(double)this->gridDim.nz;

            vector<double> absDirection (3,0);
            absDirection.at(0) = (direction.at(0)<0) ? -1*direction.at(0) : direction.at(0);
            absDirection.at(1) = (direction.at(1)<0) ? -1*direction.at(1) : direction.at(1);
            absDirection.at(2) = (direction.at(2)<0) ? -1*direction.at(2) : direction.at(2);

            tDeltaX = voxelSizeX/absDirection.at(0);
            tDeltaY = voxelSizeY/absDirection.at(1);
            tDeltaZ = voxelSizeZ/absDirection.at(2);

            // First check, if the return is before the voxel grid

            int retNum = 1;

            for (int j=1; j<=it->second->getNumberOfReturns(); j++){

                //Distance Between origin and laser return
                double dP = sqrt( pow( (it->second->getEchoes().at(j)->getX())-origin.at(0), 2 ) + pow( (it->second->getEchoes().at(j)->getY())-origin.at(1), 2 ) + pow( (it->second->getEchoes().at(j)->getZ())-origin.at(2), 2) );
                //Distance between origin and the first traversed voxel (if origin == first return and origin inside voxel grid origin==start
                double dS = sqrt( pow( (start.at(0) - origin.at(0)), 2) + pow( (start.at(1) - origin.at(1)), 2) + pow( (start.at(2) - origin.at(2)), 2) );

                if (dP < dS) {
                    retNum++;
                }
            }

            /*
            vector<int> traversedX;
            vector<int> traversedY;
            vector<int> traversedZ;
            */
            // TODO: add support for DTM!
            while ( (x<this->gridDim.nx)&&(x>=0) && (y<this->gridDim.ny)&&(y>=0) && (z<this->gridDim.nz)&&(z>=0) && (!this->hasDTM || (this->hasDTM && z>=this->DTMind.at(x).at(y)))) {

                //TODO: Implement path length estimation inside voxel
                //TODO: Implement the estimation of the average incidence anlge output grid
                /*
                traversedX.push_back(x);
                traversedY.push_back(y);
                traversedZ.push_back(z);
                */

                if ( retNum <= it->second->getNumberOfReturns() && x==XInd.at(retNum-1) && y==YInd.at(retNum-1) && z==ZInd.at(retNum-1)) {
                    //the return is inside this voxel, feed information from return into this voxel

                    //cover the case where multiple returns are inside the same voxel:
                    while (retNum <= it->second->getNumberOfReturns() && x==XInd.at(retNum-1) && y==YInd.at(retNum-1) && z==ZInd.at(retNum-1)) {

                        this->Nhit.at(x).at(y).at(z)++;

                        //TODO: implement NhitWeigh and NhitWeighPathL!

                        retNum++;
                        regHit++;

                    }

                } else { //if there are no hits

                    // if there are still hits that were not yet traversed:
                    if (retNum <= it->second->getNumberOfReturns()) {
                        this->Nmiss.at(x).at(y).at(z)++;

                        //TODO: Implement NmissWeigh and NmissWeighPathL!

                    } else { //If all returns were traversed, the voxel is occluded for this pulse
                        this->Nocc.at(x).at(y).at(z)++;
                    }
                }

                //TODO: implement pulse density calculation when pulse reaches DTM! In general implement DTM support!

                // make step into next voxel:
                if (tMaxX < tMaxY) {
                    if (tMaxX < tMaxZ) {
                        x = x + stepX;
                        tMaxX = tMaxX + tDeltaX;
                    } else if (tMaxX > tMaxZ) {
                        z = z + stepZ;
                        tMaxZ = tMaxZ + tDeltaZ;
                    } else { // if they are equal
                        x = x + stepX;
                        z = z + stepZ;
                        tMaxX = tMaxX + tDeltaX;
                        tMaxZ = tMaxZ + tDeltaZ;
                    }

                } else if (tMaxX > tMaxY){
                    if (tMaxY < tMaxZ) {
                        y = y + stepY;
                        tMaxY = tMaxY + tDeltaY;
                    } else if ( tMaxY > tMaxZ) {
                        z = z + stepZ;
                        tMaxZ = tMaxZ + tDeltaZ;
                    } else { //if they are equal
                        y = y + stepY;
                        z = z + stepZ;
                        tMaxY = tMaxY + tDeltaY;
                        tMaxZ = tMaxZ + tDeltaZ;
//
                    }
                } else { // if they are equal
                    if (tMaxX < tMaxZ) {
                        x = x + stepX;
                        tMaxX = tMaxX + tDeltaX;
                        y = y + stepY;
                        tMaxY = tMaxY + tDeltaY;
                    } else if (tMaxX > tMaxZ) {
                        z = z + stepZ;
                        tMaxZ = tMaxZ + tDeltaZ;
                    } else if (tMaxX == tMaxZ) { // if all tMax are equal
                        x = x + stepX;
                        z = z + stepZ;
                        y = y + stepY;
                        tMaxX = tMaxX + tDeltaX;
                        tMaxZ = tMaxZ + tDeltaZ;
                        tMaxY = tMaxY + tDeltaY;
                    }
                }
            }
            if (retNum <= (it->second->getNumberOfReturns()-echoesOfPulseOutside)) {

                /*
                cout << "#### Not all echoes of the pulse with GPSTime " << it->second->getGPSTime() << " have been traversed!!! Number of missing echoes: " << it->second->getNumberOfReturns()-retNum+1 << endl;
                cout << " Missing return infos: " << endl;
                */
                for (int r=retNum; r<=it->second->getNumberOfReturns(); r++) {
                    numMissingReturns++;

                    /*
                    double xtmp = it->second->getEcho(r)->getX();
                    double ytmp = it->second->getEcho(r)->getY();
                    double ztmp = it->second->getEcho(r)->getZ();

                    cout << "Ret Num: " << it->second->getEcho(r)->getReturnNumber() << " || X: " << it->second->getEcho(r)->getX() << " || Y: " << it->second->getEcho(r)->getY() << " || Z: " << it->second->getEcho(r)->getZ() << " || XInd at: " << XInd.at(r-1) << " || YInd at: " << YInd.at(r-1) << " || ZInd at: " << ZInd.at(r-1) <<     endl;
                    */

                }
                /*
                cout << "###############################" << endl;
                */

                //Pulse tmpPls = *it->second;
                //map<int,shared_ptr<Echo> > tmpEch =  it->second->getEchoes();


            }

        }
        // TODO: fix this here similar to incomplete pulse fix - use auto iterator and re assign
        // The pulse has been traversed and can now be dropped from the map - !!! It seems that this causes problems!
        // - deleting an entry while iterating over the map. TODO: figure out how we could handle this.
        //it = this->pulsedataset.erase(it);
    }

    // update feedback on voxel traversl
    this->traversedPulses = this->traversedPulses + traversedPulses;
    this->totalPulsesInDataset = this->totalPulsesInDataset + pulsecount;
    this->regHit = this->regHit + regHit;
    this->echoesOutside = this->echoesOutside + echoesOutside;
    this->numMissingReturns = this->numMissingReturns + numMissingReturns;
    this->numNoGridIntersection = this->numNoGridIntersection + numNoGridIntersection;

    /*
    cout << "#### " << traversedPulses << " Pulses were traversed of possible " << pulsecount << " Pulses" << endl;
    cout << "#### " << regHit << " Returns have been registered by the algorithm " << endl;
    cout << "#### " << echoesOutside << " Returns were found outside the voxel grid " << endl;
    cout << "#### " << numMissingReturns << " Returns were missed during the traversal!" << endl;
    cout << "#### " << numNoGridIntersection << " Pulses did not intersect voxel grid!" << endl;
    */
}

string Raytracer::reportOnTraversal(){
    /*
    cout << "#### " << this->traversedPulses << " Pulses were traversed of possible " << this->totalPulsesInDataset << " Pulses" << endl;
    cout << "#### " << this->regHit << " Returns have been registered by the algorithm " << endl;
    cout << "#### " << this->echoesOutside << " Returns were found outside the voxel grid " << endl;
    cout << "#### " << this->numMissingReturns << " Returns were missed during the traversal!" << endl;
    cout << "#### " << this->numNoGridIntersection << " Pulses did not intersect voxel grid!" << endl;
    */
    stringstream multilineString;

    multilineString << "#### " << this->traversedPulses << " Pulses were traversed of possible " << this->totalPulsesInDataset << " Pulses" << "\n";
    multilineString << "#### " << this->regHit << " Returns have been registered by the algorithm " << "\n";
    multilineString << "#### " << this->echoesOutside << " Returns were found outside the voxel grid " << "\n";
    multilineString << "#### " << this->numMissingReturns << " Returns were missed during the traversal!" << "\n";
    multilineString << "#### " << this->numNoGridIntersection << " Pulses did not intersect voxel grid!" << "\n";

    // convert the stringstream to a string
    string out_str = multilineString.str();

    cout << out_str;

    return out_str;
}

void Raytracer::clearPulseDataset()
{
    this->pulsedataset.clear();
    this->emptypulsedataset.clear();
}

void Raytracer::doRaytracing_singleReturnPulses(vector<double> X, vector<double> Y, vector<double> Z, vector<double> sensor_x, vector<double> sensor_y, vector<double> sensor_z, vector<double> gps_time)
{
    //TODO: Add functionality for complete voxel traversal reporting as done in the doRaytracing() function!
    int pulsecount = 0;
    int traversedPulses = 0;
    int echoesOutside = 0;
    int numMissingReturns = 0;
    int numNoGridIntersection = 0;

    int regHit = 0; //Number of registered hit voxels

    // #pragma omp parallel for
    for (int i = 0; i < gps_time.size(); i++) {
        pulsecount++;

        bool isfound = false;

        //update progressbar
        // loadbar(pulsecount, (int)gps_time.size(), 20);

        //init several necessary variables befor starting
        int flag = 0;               //flag whether pulse crosses voxel grid or not
        double tmin = 0;             //how far along the pulse direction we have to travel until we reach the voxel grid. based on the vector representation of a line: y = u + tmin*v (y,u and v are vectors)

        double tVoxelX = 0;
        double tVoxelY = 0;
        double tVoxelZ = 0;
        int stepX = 0;
        int stepY = 0;
        int stepZ = 0;

        double voxelMaxX = 0;
        double voxelMaxY = 0;
        double voxelMaxZ = 0;

        double voxelSizeX = 0;
        double voxelSizeY = 0;
        double voxelSizeZ = 0;

        double tDeltaX = 0;
        double tDeltaY = 0;
        double tDeltaZ = 0;

        double tMaxX = 0;
        double tMaxY = 0;
        double tMaxZ = 0;

        vector<double> origin (3,0);
        vector<double> direction (3,0);

        // We are assuming that we do know the exact sensor position
        origin.at(0) = sensor_x.at(i);
        origin.at(1) = sensor_y.at(i);
        origin.at(2) = sensor_z.at(i);

        direction.at(0) = X.at(i) - origin.at(0);
        direction.at(1) = Y.at(i) - origin.at(1);
        direction.at(2) = Z.at(i) - origin.at(2);

        // if direction vector is == 0 assign a very small number to overcome problems with division through 0 in later steps. TODO: Check if this value is small enough. This could be a problem especially with TLS data!
        if (direction.at(0)==0) { direction.at(0) = numeric_limits<double>::min(); }
        if (direction.at(1)==0) { direction.at(1) = numeric_limits<double>::min(); }
        if (direction.at(2)==0) { direction.at(2) = numeric_limits<double>::min(); }

        //check if pulse intersects voxel grid
        rayBoxIntersection(origin, direction, this->gridDim.minBound, this->gridDim.maxBound, flag, tmin);

        if (flag==0) {
            numNoGridIntersection++;
            //cout << "The ray does not intersect the grid" << endl;
        } else {
            //cout << "The ray does intersect the grid" << endl;

            traversedPulses++;

            // int echoesOfPulseOutside = 0; // not needed as Horizon only has one return

            if (tmin<0) {
                tmin = 0;
            }

            vector<double> start (3,0);
            start.at(0) = origin.at(0) + (tmin*direction.at(0)); //+ tmin*direction.at(1) + tmin*direction.at(2);
            start.at(1) = origin.at(1) + (tmin*direction.at(1)); //+ tmin*direction.at(1) + tmin*direction.at(2);
            start.at(2) = origin.at(2) + (tmin*direction.at(2)); //+ tmin*direction.at(1) + tmin*direction.at(2);

            vector<double> boxSize (3,0);
            boxSize.at(0) = (double)this->gridDim.maxBound.at(0) - (double)this->gridDim.minBound.at(0);
            boxSize.at(1) = (double)this->gridDim.maxBound.at(1) - (double)this->gridDim.minBound.at(1);
            boxSize.at(2) = (double)this->gridDim.maxBound.at(2) - (double)this->gridDim.minBound.at(2);

            // Get the voxel indices of where the laser return lie in the grid
            int XInd = 0; // There are only one return for the used MLS sensor TODO: Check on that!
            int YInd = 0;
            int ZInd = 0;

            XInd = floor( ((X.at(i) - (double)this->gridDim.minBound.at(0))/boxSize.at(0))*(double)this->gridDim.nx );
            YInd = floor( ((Y.at(i) - (double)this->gridDim.minBound.at(1))/boxSize.at(1))*(double)this->gridDim.ny );
            ZInd = floor( ((Z.at(i) - (double)this->gridDim.minBound.at(2))/boxSize.at(2))*(double)this->gridDim.nz );

            if ( YInd >= this->gridDim.ny || YInd < 0 || XInd >= this->gridDim.nx || XInd < 0 || ZInd >= this->gridDim.nz || ZInd < 0 ){
                echoesOutside++;
                // echoesOfPulseOutside++; // Not needed as Horizon only has one return!
            }

            int x;
            int y;
            int z;

            x = floor( ((start.at(0)-(double)this->gridDim.minBound.at(0))/boxSize.at(0))*(double)this->gridDim.nx );
            y = floor( ((start.at(1)-(double)this->gridDim.minBound.at(1))/boxSize.at(1))*(double)this->gridDim.ny );
            z = floor( ((start.at(2)-(double)this->gridDim.minBound.at(2))/boxSize.at(2))*(double)this->gridDim.nz );

            if (x==(this->gridDim.nx)) {
                x = x-1;
            }
            if (y==(this->gridDim.ny)){
                y = y-1;
            }
            if (z==(this->gridDim.nz)){
                z = z-1;
            }

            if (direction.at(0)>=0) {
                tVoxelX = (double)(x+1)/(double)this->gridDim.nx;
                stepX = 1;
            } else {
                tVoxelX = ((double)(x+1)-1)/(double)this->gridDim.nx;
                stepX = -1;
            }
            if (direction.at(1)>=0) {
                tVoxelY = (double)(y+1)/(double)this->gridDim.ny;
                stepY = 1;
            } else {
                tVoxelY = ((double)(y+1)-1)/(double)this->gridDim.ny;
                stepY = -1;
            }
            if (direction.at(2)>=0) {
                tVoxelZ = (double)(z+1)/(double)this->gridDim.nz;
                stepZ = 1;
            } else {
                tVoxelZ = ((double)(z+1)-1)/(double)this->gridDim.nz;
                stepZ = -1;
            }

            voxelMaxX = (double)this->gridDim.minBound.at(0) + tVoxelX*boxSize.at(0);
            voxelMaxY = (double)this->gridDim.minBound.at(1) + tVoxelY*boxSize.at(1);
            voxelMaxZ = (double)this->gridDim.minBound.at(2) + tVoxelZ*boxSize.at(2);

            tMaxX = tmin + (voxelMaxX-start.at(0))/direction.at(0);
            tMaxY = tmin + (voxelMaxY-start.at(1))/direction.at(1);
            tMaxZ = tmin + (voxelMaxZ-start.at(2))/direction.at(2);

            voxelSizeX = boxSize.at(0)/(double)this->gridDim.nx;
            voxelSizeY = boxSize.at(1)/(double)this->gridDim.ny;
            voxelSizeZ = boxSize.at(2)/(double)this->gridDim.nz;

            vector<double> absDirection (3,0);
            absDirection.at(0) = (direction.at(0)<0) ? -1*direction.at(0) : direction.at(0);
            absDirection.at(1) = (direction.at(1)<0) ? -1*direction.at(1) : direction.at(1);
            absDirection.at(2) = (direction.at(2)<0) ? -1*direction.at(2) : direction.at(2);

            tDeltaX = voxelSizeX/absDirection.at(0);
            tDeltaY = voxelSizeY/absDirection.at(1);
            tDeltaZ = voxelSizeZ/absDirection.at(2);

            // First check, if the return is before the voxel grid

            int retNum = 1;

            //Distance Between origin and laser return
            double dP = sqrt( pow(X.at(i)-origin.at(0), 2) + pow(Y.at(i)-origin.at(1), 2) + pow(Z.at(i)-origin.at(2), 2));
            //Distance between origin and the first traversed voxel (if origin == first return and origin inside voxel grid origin==start
            double dS = sqrt( pow( (start.at(0) - origin.at(0)), 2) + pow( (start.at(1) - origin.at(1)), 2) + pow( (start.at(2) - origin.at(2)), 2) );

            if ( dP < dS ) {
                retNum++;
            }

            while ( (x<this->gridDim.nx)&&(x>=0) && (y<this->gridDim.ny)&&(y>=0) && (z<this->gridDim.nz)&&(z>=0) && (!this->hasDTM || (this->hasDTM && z>=this->DTMind.at(x).at(y)))) {

                if ( retNum <= 1 && x==XInd && y==YInd && z==ZInd) {
                    //the return is inside this voxel, feed information from return into this voxel

                    //cover the case where multiple returns are inside the same voxel:
                    while (retNum <= 1 && x==XInd && y==YInd && z==ZInd) {

                        this->Nhit.at(x).at(y).at(z)++;

                        //TODO: implement NhitWeigh and NhitWeighPathL!

                        retNum++;
                        regHit++;

                    }

                } else { //if there are no hits

                    // if there are still hits that were not yet traversed:
                    if (retNum <= 1) {
                        this->Nmiss.at(x).at(y).at(z)++;

                        //TODO: Implement NmissWeigh and NmissWeighPathL!

                    } else { //If all returns were traversed, the voxel is occluded for this pulse
                        this->Nocc.at(x).at(y).at(z)++;
                    }
                }

                //TODO: implement pulse density calculation when pulse reaches DTM! In general implement DTM support!

                // make step into next voxel:
                if (tMaxX < tMaxY) {
                    if (tMaxX < tMaxZ) {
                        x = x + stepX;
                        tMaxX = tMaxX + tDeltaX;
                    } else if (tMaxX > tMaxZ) {
                        z = z + stepZ;
                        tMaxZ = tMaxZ + tDeltaZ;
                    } else { // if they are equal
                        x = x + stepX;
                        z = z + stepZ;
                        tMaxX = tMaxX + tDeltaX;
                        tMaxZ = tMaxZ + tDeltaZ;
                    }

                } else if (tMaxX > tMaxY){
                    if (tMaxY < tMaxZ) {
                        y = y + stepY;
                        tMaxY = tMaxY + tDeltaY;
                    } else if ( tMaxY > tMaxZ) {
                        z = z + stepZ;
                        tMaxZ = tMaxZ + tDeltaZ;
                    } else { //if they are equal
                        y = y + stepY;
                        z = z + stepZ;
                        tMaxY = tMaxY + tDeltaY;
                        tMaxZ = tMaxZ + tDeltaZ;
//
                    }
                } else { // if they are equal
                    if (tMaxX < tMaxZ) {
                        x = x + stepX;
                        tMaxX = tMaxX + tDeltaX;
                        y = y + stepY;
                        tMaxY = tMaxY + tDeltaY;
                    } else if (tMaxX > tMaxZ) {
                        z = z + stepZ;
                        tMaxZ = tMaxZ + tDeltaZ;
                    } else if (tMaxX == tMaxZ) { // if all tMax are equal
                        x = x + stepX;
                        z = z + stepZ;
                        y = y + stepY;
                        tMaxX = tMaxX + tDeltaX;
                        tMaxZ = tMaxZ + tDeltaZ;
                        tMaxY = tMaxY + tDeltaY;
                    }
                }
            }
        }
    }

    // update feedback on voxel traversl
    this->traversedPulses = this->traversedPulses + traversedPulses;
    this->totalPulsesInDataset = this->totalPulsesInDataset + pulsecount;
    this->regHit = this->regHit + regHit;
    this->echoesOutside = this->echoesOutside + echoesOutside;
    this->numMissingReturns = this->numMissingReturns + numMissingReturns;
    this->numNoGridIntersection = this->numNoGridIntersection + numNoGridIntersection;

    /*
    cout << "#### " << traversedPulses << " Pulses were traversed of possible " << pulsecount << " Pulses" << endl;
    cout << "#### " << regHit << " Returns have been registered by the algorithm " << endl;
    cout << "#### " << echoesOutside << " Returns were found outside the voxel grid " << endl;
    cout << "#### " << numMissingReturns << " Returns were missed during the traversal!" << endl;
    */
}

void Raytracer::rayBoxIntersection (vector<double> origin, vector<double> direction, vector<double> vmin, vector<double> vmax, int & flag, double & tmin){
    // Ray/box intersection using the Smits' algorithm
    // Input:
    //  origin
    //  direction
    //  box = (vmin,vmax)
    // Output:
    //  result: result[0] = flag: (0) Reject, (1) Intersect
    //          result[1] = tmin: distance from the ray origin


    double tmax = 0;
    double tymin = 0;
    double tymax = 0;
    double tzmin = 0;
    double tzmax = 0;


    if (direction[0] >= 0) {
        tmin = (vmin[0] - origin[0]) / direction[0];
        tmax = (vmax[0] - origin[0]) / direction[0];
    } else {
        tmin = (vmax[0] - origin[0]) / direction[0];
        tmax = (vmin[0] - origin[0]) / direction[0];
    }

    if (direction[1] >= 0) {
        tymin = (vmin[1] - origin[1]) / direction[1];
        tymax = (vmax[1] - origin[1]) / direction[1];
    } else {
        tymin = (vmax[1] - origin[1]) / direction[1];
        tymax = (vmin[1] - origin[1]) / direction[1];
    }

    if ( (tmin > tymax) || (tymin > tmax) ) {
        flag = 0;
        tmin = -1;
        return;
    }

    if (tymin > tmin) {
        tmin = tymin;
    }

    if (tymax < tmax) {
        tmax = tymax;
    }

    if (direction[2] >= 0) {
        tzmin = (vmin[2] - origin[2]) / direction[2];
        tzmax = (vmax[2] - origin[2]) / direction[2];
    } else {
        tzmin = (vmax[2] - origin[2]) / direction[2];
        tzmax = (vmin[2] - origin[2]) / direction[2];
    }

    if ((tmin > tzmax) || (tzmin > tmax)) {
        flag = 0;
        tmin = -1;
        return;
    }

    if (tzmin > tmin) {
        tmin = tzmin;
    }

    if (tzmax < tmax) {
        tmax = tzmax;
    }

    flag = 1;

}

vector<vector<vector<int > > >& Raytracer::getNhit(){
    return this->Nhit;
}

vector<vector<vector<int > > >& Raytracer::getNmiss(){
    return this->Nmiss;
}

vector<vector<vector<int > > >& Raytracer::getNocc(){
    return this->Nocc;
}

vector<vector<double > >& Raytracer::reportSensorShifts(){
    return this->SensorShifts;
}

vector<int > Raytracer::getGridDimensions(){
    vector<int> grid_dim (3,0);
    grid_dim.at(0) = this->gridDim.nx;
    grid_dim.at(1) = this->gridDim.ny;
    grid_dim.at(2) = this->gridDim.nz;

    return grid_dim;
}

vector<double > Raytracer::getGridOrigin(){
    return this->gridDim.origin;
}

vector<double> Raytracer::getPulsesIntersectingBox(vector<double> x, vector<double> y, vector<double> z, vector<double> sensor_x, vector<double> sensor_y, vector<double> sensor_z, vector<double> gps_time, vector<double> vmin, vector<double> vmax){

    int pulsecount = 0;
    int numNoGridIntersection = 0;
    int numGridIntersection = 0;

    vector<double> pulses_intersecting_box;

    int total_pulses2iterate = gps_time.size();

    // cout << "Total number of pulses to check: " << total_pulses2iterate << endl;
    // cout << "minbounds: " << vmin.at(0) << " " << vmin.at(1) << " " << vmin.at(2) << endl;
    // cout << "maxbounds: " << vmax.at(0) << " " << vmax.at(1) << " " << vmax.at(2) << endl;

    //#pragma omp parallel for
    for (int i = 0; i < gps_time.size(); i++) {
        //vector<int> local_intersecting_pulses; // needed when using omp parallel as pulses_intersecting_box is shared between threads
        pulsecount++;

        // init several necessary variables before starting
        int flag = 0;               //flag whether pulse crosses voxel grid or not
        double tmin = 0;             //how far along the pulse direction we have to travel until we reach the voxel grid. based on the vector representation of a line: y = u + tmin*v (y,u and v are vectors)

        vector<double> origin (3,0);
        vector<double> direction (3,0);

        // We are assuming that we do know the exact sensor position
        origin.at(0) = sensor_x.at(i);
        origin.at(1) = sensor_y.at(i);
        origin.at(2) = sensor_z.at(i);

        direction.at(0) = x.at(i) - origin.at(0);
        direction.at(1) = y.at(i) - origin.at(1);
        direction.at(2) = z.at(i) - origin.at(2);

        // if direction vector is == 0 assign a very small number to overcome problems with division through 0 in later steps. TODO: Check if this value is small enough. This could be a problem especially with TLS data!
        if (direction.at(0)==0) { direction.at(0) = numeric_limits<double>::min(); }
        if (direction.at(1)==0) { direction.at(1) = numeric_limits<double>::min(); }
        if (direction.at(2)==0) { direction.at(2) = numeric_limits<double>::min(); }

        //check if pulse intersects voxel grid
        rayBoxIntersection(origin, direction, vmin, vmax, flag, tmin);

        if (flag==0) {
            numNoGridIntersection++;
            // pulse not relevant remove from pulse dataset
            //it = this->pulsedataset.erase(it);
            //cout << "The ray does not intersect the grid" << endl;
        } else {
            //cout << "The ray does intersect the grid" << endl;

            numGridIntersection++;

            pulses_intersecting_box.push_back(gps_time.at(i));
        }

        /*
        #pragma omp critical -> TODO: this just does not work right! Find out why and decide whether we want to go down that route!
        {
            pulses_intersecting_box.insert(pulses_intersecting_box.end(), local_intersecting_pulses.begin(), local_intersecting_pulses.end());
        }
        */

    }

    // TODO: Add verbose flag!
    /*
    cout << "Number of pulses intersecting the grid: " << numGridIntersection << endl;
    cout << "Number of pulses NOT intersecting the grid: " << numNoGridIntersection << endl;
    */

    // return pulses_intersecting_box vector
    return pulses_intersecting_box; // this vector just includes the GPSTimes of the pulses intersecting the box.
}




void Raytracer::addEmptyPulseData(vector<double> sensor_x, vector<double> sensor_y, vector<double> sensor_z, vector<double> direction_x, vector<double> direction_y, vector<double> direction_z, vector<double> gps_time){
    // add empty pulses
    int duplicate_empty_pulses = 0;
    for (int i = 0; i < gps_time.size(); i++) {
        map<double,boost::shared_ptr<Pulse> >::iterator it;

        double gps = gps_time.at(i);

        it = this->emptypulsedataset.find(gps);

        if (it != this->emptypulsedataset.end()) {
            duplicate_empty_pulses++;
        } else {
            boost::shared_ptr<Pulse> p = boost::make_shared<Pulse>(gps, sensor_x.at(i), sensor_y.at(i), sensor_z.at(i), direction_x.at(i), direction_y.at(i), direction_z.at(i));
            this->emptypulsedataset.insert(std::make_pair(p->getGPSTime(),p));
        }
    }
    // TODO: verbose flag
    if (duplicate_empty_pulses > 0) {
        cout << "Found " << duplicate_empty_pulses << "number of duplicate pulses (based on gps time), only first occurence is kept" << endl;
    }

}

void Raytracer::doRaytracingEmptyPulses(){
    // this function works very similar to doRaytracing, except that these pulses have no returns and should always move to edge of voxelgrid

    // the code has a lot of overlap, so could potentially be merged, but this would complicate the doRaytracing function
    // and currently this code is only used for RIEGL scanners in an experimental setup anyways
    // sooo for now seperate functions and some duplicate code

    int pulsecount = 0;
    int traversedPulses = 0;
    int numNoGridIntersection = 0;

    int total_pulses2iterate = (int)this->emptypulsedataset.size();
    // init boxsize and voxelsize, remains constant for all pulses
    vector<double> boxSize (3,0);
    boxSize.at(0) = (double)this->gridDim.maxBound.at(0) - (double)this->gridDim.minBound.at(0);
    boxSize.at(1) = (double)this->gridDim.maxBound.at(1) - (double)this->gridDim.minBound.at(1);
    boxSize.at(2) = (double)this->gridDim.maxBound.at(2) - (double)this->gridDim.minBound.at(2);
    vector<double> voxelSize (3,0);
    voxelSize.at(0) = boxSize.at(0)/(double)this->gridDim.nx;
    voxelSize.at(1) = boxSize.at(1)/(double)this->gridDim.ny;
    voxelSize.at(2) = boxSize.at(2)/(double)this->gridDim.nx;

    for (map<double,boost::shared_ptr<Pulse> >::iterator it = this->emptypulsedataset.begin(); it != this->emptypulsedataset.end(); ++it) {
        pulsecount++;

        //update progressbar
        loadbar(pulsecount, total_pulses2iterate, 20);

        //init several necessary variables befor starting
        int flag = 0;               //flag whether pulse crosses voxel grid or not
        double tmin = 0;             //how far along the pulse direction we have to travel until we reach the voxel grid. based on the vector representation of a line: y = u + tmin*v (y,u and v are vectors)

        vector<double> origin (3,0);
        vector<double> direction (3,0);

        origin.at(0) = it->second->getSensorX();
        origin.at(1) = it->second->getSensorY();
        origin.at(2) = it->second->getSensorZ();

        // get direction directly from pulse, as we have no returns
        direction.at(0) = it->second->getDirectionX();
        direction.at(1) = it->second->getDirectionY();
        direction.at(2) = it->second->getDirectionZ();

        // if direction vector is == 0 assign a very small number to overcome problems with division through 0 in later steps.
        if (direction.at(0)==0) { direction.at(0) = numeric_limits<double>::min(); }
        if (direction.at(1)==0) { direction.at(1) = numeric_limits<double>::min(); }
        if (direction.at(2)==0) { direction.at(2) = numeric_limits<double>::min(); }
        
        //check if pulse intersects voxel grid
        rayBoxIntersection(origin, direction, this->gridDim.minBound, this->gridDim.maxBound, flag, tmin);
        // if no intersection, go to next pulse
        if (flag==0) {
            numNoGridIntersection++;
            continue;
        }

        traversedPulses++;
        if (tmin<0) {
            tmin = 0;
        }

        vector<double> start (3,0);
        start.at(0) = origin.at(0) + (tmin*direction.at(0)); //+ tmin*direction.at(1) + tmin*direction.at(2);
        start.at(1) = origin.at(1) + (tmin*direction.at(1)); //+ tmin*direction.at(1) + tmin*direction.at(2);
        start.at(2) = origin.at(2) + (tmin*direction.at(2)); //+ tmin*direction.at(1) + tmin*direction.at(2);
        
        // get voxel indices of start
        int x = floor( ((start.at(0)-(double)this->gridDim.minBound.at(0))/boxSize.at(0))*(double)this->gridDim.nx );
        int y = floor( ((start.at(1)-(double)this->gridDim.minBound.at(1))/boxSize.at(1))*(double)this->gridDim.ny );
        int z = floor( ((start.at(2)-(double)this->gridDim.minBound.at(2))/boxSize.at(2))*(double)this->gridDim.nz );

        if (x==(this->gridDim.nx)) {
            x = x-1;
        }
        if (y==(this->gridDim.ny)){
            y = y-1;
        }
        if (z==(this->gridDim.nz)){
            z = z-1;
        }

        // initialize step direction for all three dimensions, 
        // also get tVoxel, denotes the fractional voxel index of the next voxel along the ray in each dimension
        double tVoxelX = 0;
        double tVoxelY = 0;
        double tVoxelZ = 0;
        int stepX = 0;
        int stepY = 0;
        int stepZ = 0;

        if (direction.at(0)>=0) {
            tVoxelX = (double)(x+1)/(double)this->gridDim.nx;
            stepX = 1;
        } else {
            tVoxelX = ((double)(x+1)-1)/(double)this->gridDim.nx;
            stepX = -1;
        }
        if (direction.at(1)>=0) {
            tVoxelY = (double)(y+1)/(double)this->gridDim.ny;
            stepY = 1;
        } else {
            tVoxelY = ((double)(y+1)-1)/(double)this->gridDim.ny;
            stepY = -1;
        }
        if (direction.at(2)>=0) {
            tVoxelZ = (double)(z+1)/(double)this->gridDim.nz;
            stepZ = 1;
        } else {
            tVoxelZ = ((double)(z+1)-1)/(double)this->gridDim.nz;
            stepZ = -1;
        }

        // voxelmax denotes the coordinate of the next voxel start along the ray
        double voxelMaxX = (double)this->gridDim.minBound.at(0) + tVoxelX*boxSize.at(0);
        double voxelMaxY = (double)this->gridDim.minBound.at(1) + tVoxelY*boxSize.at(1);
        double voxelMaxZ = (double)this->gridDim.minBound.at(2) + tVoxelZ*boxSize.at(2);

        // tMax denotes the distance we have to travel along the ray to reach the next voxel, in each dimension
        // will be used to determine which voxel change occurs first, and updated when we move along ray
        double tMaxX = tmin + (voxelMaxX-start.at(0))/direction.at(0);
        double tMaxY = tmin + (voxelMaxY-start.at(1))/direction.at(1);
        double tMaxZ = tmin + (voxelMaxZ-start.at(2))/direction.at(2);

        vector<double> absDirection (3,0);
        absDirection.at(0) = (direction.at(0)<0) ? -1*direction.at(0) : direction.at(0);
        absDirection.at(1) = (direction.at(1)<0) ? -1*direction.at(1) : direction.at(1);
        absDirection.at(2) = (direction.at(2)<0) ? -1*direction.at(2) : direction.at(2);
        // tdelta denotes how far along ray we have to move so the component equals a voxel length in that direction
        // if we change voxel in a dimension, we have to add tDelta to tMax in that dimension
        double tDeltaX = voxelSize.at(0)/absDirection.at(0);
        double tDeltaY = voxelSize.at(1)/absDirection.at(1);
        double tDeltaZ = voxelSize.at(2)/absDirection.at(2);

        // travel along ray as long as we are in voxelgrid
        while ( (x<this->gridDim.nx)&&(x>=0) && (y<this->gridDim.ny)&&(y>=0) && (z<this->gridDim.nz)&&(z>=0) && (!this->hasDTM || (this->hasDTM && z>=this->DTMind.at(x).at(y)))) {
            
            // as we are modelling an empty pulse, every traversal here is a miss
            this->Nmiss.at(x).at(y).at(z)++;

            // make next step
            if (tMaxX < tMaxY) {
                if (tMaxX < tMaxZ) {
                    x = x + stepX;
                    tMaxX = tMaxX + tDeltaX;
                } else if (tMaxX > tMaxZ) {
                    z = z + stepZ;
                    tMaxZ = tMaxZ + tDeltaZ;
                } else { // if they are equal
                    x = x + stepX;
                    z = z + stepZ;
                    tMaxX = tMaxX + tDeltaX;
                    tMaxZ = tMaxZ + tDeltaZ;
                }

            } else if (tMaxX > tMaxY){
                if (tMaxY < tMaxZ) {
                    y = y + stepY;
                    tMaxY = tMaxY + tDeltaY;
                } else if ( tMaxY > tMaxZ) {
                    z = z + stepZ;
                    tMaxZ = tMaxZ + tDeltaZ;
                } else { //if they are equal
                    y = y + stepY;
                    z = z + stepZ;
                    tMaxY = tMaxY + tDeltaY;
                    tMaxZ = tMaxZ + tDeltaZ;
                }
            } else { // if they are equal
                if (tMaxX < tMaxZ) {
                    x = x + stepX;
                    tMaxX = tMaxX + tDeltaX;
                    y = y + stepY;
                    tMaxY = tMaxY + tDeltaY;
                } else if (tMaxX > tMaxZ) {
                    z = z + stepZ;
                    tMaxZ = tMaxZ + tDeltaZ;
                } else if (tMaxX == tMaxZ) { // if all tMax are equal
                    x = x + stepX;
                    z = z + stepZ;
                    y = y + stepY;
                    tMaxX = tMaxX + tDeltaX;
                    tMaxZ = tMaxZ + tDeltaZ;
                    tMaxY = tMaxY + tDeltaY;
                }
            }
        }
    }

    // update feedback on voxel traversl
    this->traversedPulses = this->traversedPulses + traversedPulses;
    this->totalPulsesInDataset = this->totalPulsesInDataset + pulsecount;
    this->numNoGridIntersection = this->numNoGridIntersection + numNoGridIntersection;
}


// TODO: change tMax, step and tDelta to vectors
// then move to seperate function we can use in both raytracing functions
// def Raytracer::stepVoxel(tMax, step, tDelta){

// Get functions for status reporting of raytracing

int Raytracer::get_num_traversed_pulses() const{
    return this->traversedPulses;
}
int Raytracer::get_total_pulses_in_dataset() const{
    return this->totalPulsesInDataset;
}
int Raytracer::get_num_registered_hits() const{
    return this->regHit;
}
int Raytracer::get_num_echoes_outside() const{
    return this->echoesOutside;
}
int Raytracer::get_num_missing_returns() const{
    return this->numMissingReturns;
}
int Raytracer::get_num_pulses_no_intersection() const{
    return this->numNoGridIntersection;
}




/*
void Raytracer_MLS::readDTM(vector<vector< double > > DTM){
    this->hasDTM = true;

    //extract DTM indices

    this->DTMind.resize(this->gridDim.ny, vector< double >(this->gridDim.nx,0));


    for (int y=0; y<this->DTMind.size(); y++) {
        for (int x=0; x<this->DTMind.at(y).size(); x++) {

            this->DTMind.at(y).at(x) = floor( ((DTM.at(y).at(x) - this->gridDim.minBound.at(2)) / (this->gridDim.maxBound.at(2) - this->gridDim.minBound.at(2))) * this->gridDim.nz ) - 4; //The index of where the DTM is located for each (y,x) coordinate pair. Additionally we subtract four voxels from this index in order to assure that we traverse far enough so that we have a thorough Classification grid until a few meters below the DTM.

        }
    }
}
*/