//Header Guard
#ifndef READWRITEFUNCTIONSATAC_H
#define	READWRITEFUNCTIONSATAC_H

#include <stdio.h>//For chdir()
#include <unistd.h> //For chdir())
#include "armadillo"
#include "GeometryFunctionsATAC.h"
#include "DataStructureFunctionsATAC.h"


using namespace arma;
using namespace std;

//given a string, tries to match that string in a given variable file, and returns the
//value immediately next to the keyword.

string variableReader(string keyWord);

//given a string, tries to match that string in a given configuration file, and returns the
//value immediately next to the keyword. Currently reads in the second token as the key, which means it only works for reading temperature.

string confReader(string keyWord, string masterBatchName, int keyWordPosition) ;



//This subfunction writes generates a datafile to be used in final movie
//showing the evolution of the chiral angles for all 6 rings.
//assumes tube axis is in x direction.

void recordChirality(vector<pair<imat, int> > rings, mat pos, bool findIndexes) ;




//This subfunction writes generates a datafile to be used in a final movie
//showing the evolution of the sizes of the carbon rings.

void recordRingSizes(vector<pair<imat, int> > rings, int& globalRingSizeRunNum, bool& generateRingHistogramImages);



//Writes an int in the C000 format used by vmd. 1000 atom limit.

string writeLabels(int i) ;


void writeTopFile(mat positions, imat bondList, string masterBatchName) ;

void writePgnFile(mat positions, string masterBatchName) ;




//Note about pdb files: positions are aloted 8 characters, with the 5th character being a decimal.
//for instance: -088.880
// we must round forces so that they fit this format, including adding white space
//before and after if the force is not already exactly 8 characters.

void writeForcePDB(mat positions, mat forceMatrix) ;

void writeConfFile(mat positions, double temperature, string masterBatchName, string& proteinDataBaseFile, int& forceNum, int& periodicityNum, double& xLatticeParameter, double& yLatticeParameter, double& zLatticeParameter, vec latticeVec1, vec latticeVec2, vec latticeVec3) ;

void writePdbPsfFiles(string masterBatchName) ;

//strips a coor xyz file to make it a matrix.

void coorStripper(mat& positions, string masterBatchName) ;

//strips xyz file to make it a matrix.

void xyzStripper(mat& positions, string masterBatchName) ;



//Reads the coor file of the last run, and writes to an .xyz file. This is then
//appended to movie.xyz.

void movieMaker(mat positions, string masterBatchName, int& numberOfElectrodes, imat& xp1Electrodes, imat& xp2Electrodes, imat& xn1Electrodes, imat& xn2Electrodes, imat& yp1Electrodes, imat& yp2Electrodes, imat& yn1Electrodes, imat& yn2Electrodes) ;

//this subfunction reads in an energy from a NAMD .log file.

double energyReader(double& newEnergy, string masterBatchName, int& relaxationMethod) ;


//This function updates the run number by reading the last inputted value of the energy history file.
//If the run type is 0, it will find the number of total runs.
//If 1, will find the number of accepted runs.

int runReader(int& runs, int runType, string batchName) ;


#endif	/* READWRITEFUNCTIONSATAC_H */
