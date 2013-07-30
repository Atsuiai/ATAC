

//Header Guard
#ifndef GEOMETRYFUNCTIONSATAC_H
#define	GEOMETRYFUNCTIONSATAC_H

//#include <stdio.h>//For chdir

#include "armadillo"

using namespace arma;
using namespace std;



//class GeometryFunctionsATAC
//{
//public:
//};
    
    //this function takes in positions, and deletes all of the ones that are outside a certain cylindrical radius.
//this is so we can seperate the inner and outer tubes of a DWNT

void tubeShaver(mat& pos);

//This function fills out 8 arrays with the atoms of 4 electrodes (8 terminals). xp1 are the
//electrodes in the +x extreme, xp2 are the electrodes in the terminal next to xp1,
// xn1 is the -x extreme, xn2 in the terminal next to xn1, etc.

void electrodeDeterminer(mat pos, imat& xp1Electrodes, imat& xn1Electrodes, imat& yp1Electrodes, imat& yn1Electrodes, imat& xp2Electrodes, imat& xn2Electrodes, imat& yp2Electrodes, imat& yn2Electrodes, int& numberOfElectrodes, double& xpElecWallMod, double& xnElecWallMod, double& ypElecWallMod, double& ynElecWallMod);


//Determines the x, y and z lattice parameters (these values fed to it are blank,
//to be assigned in the algorithm). Parameters can either be be statically given,
//or can be relative to the system given.

void latticeParameterDeterminer(mat pos, bool& isStaticPeriod, int& periodicityNum, double& xLatticeParameter, double& xPeriodFactor, double& yLatticeParameter, double& yPeriodFactor, double& zLatticeParameter, double& zPeriodFactor);


//determines the shortest path in a lattice cell between two points.

void shortestLineDeterminer(mat pos, int& i, int& j, double& distance, int& periodicityNum, double& xLatticeParameter, double& yLatticeParameter, double& zLatticeParameter, vec latticeVec1Norm, vec latticeVec1, vec latticeVec2Norm, vec latticeVec2, vec latticeVec3Norm, vec latticeVec3);


#endif	/* GEOMETRYFUNCTIONSATAC_H */
