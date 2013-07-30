
//Header Guard
#ifndef PHYSICSATAC_H
#define	PHYSICSATAC_H

//#include <stdio.h> //For chdir
#include <iomanip> // for setprecision()
#include <stdio.h>//For chdir()
#include <unistd.h> //For chdir())
#include "armadillo"
#include "ReadWriteFunctionsATAC.h"


using namespace std;
using namespace arma;


void forceParallelWalls(mat positions, mat& forceMatrix, double minWallPos, double maxWallPos, double force, int wallAxis);

void forceTube(mat positions, mat& forceMatrix, double force, double tubeRadius);

void forceVice(mat positions, mat forceMatrix, double minViceWall, double maxViceWall, double viceForce);

//simple Newtonian cooling

double newtonianCooling(double& tempS, double& tempF, int& acceptedRunNum, double& coolR) ;

//simple linear cooling. If we reach tempFinal, do not change temp.

double linearCooling(double& tempS, double& tempF, int& acceptedRunNum, double& coolR) ;




//generates three random percentages, as a vector of doubles

vector<double> threeRandomPercentages() ;


// oldEnergy is the energy of the previous structure.
// If the new allotrope of carbon is energetically favorable, it passes the energy check.
// If it does not pass the check, there is still a chance (arrhenius) that
// it will take on the metastable structure.
// Updates oldEnergy with NewEnergy if passes.
// Energy units are kcal/mol
// CHECK IF POSITIONS IS NEEDED, KCAL PER MOL?

bool energyCheck(double& oldEnergy, double& temperature, int spPossibility, double randomPercent, double ratioRotate, double ratioDimerDisplacement, string masterBatchName, int& namdRunNum, bool& relaxAfterSp2Sp3Failure, double& boltzmann, double& energyFailureTolerance, bool& acceptAnyMutation, int& relaxationMethod) ;



#endif	/* PHYSICSATAC_H */
