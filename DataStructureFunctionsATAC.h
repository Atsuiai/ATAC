
//Header Guard
#ifndef DATASTRUCTUREFUNCTIONSATAC_H
#define	DATASTRUCTUREFUNCTIONSATAC_H

//forward declare
//class GeometryFunctionsATAC;

//#include <stdio.h>//For chdir


#include "armadillo"

#include "GeometryFunctionsATAC.h"
#include "ReadWriteFunctionsATAC.h"
#include "PhysicsATAC.h"
#include "GraphTheoryATAC.h"
#include "atac_interface.h"

using namespace std;
//using namespace lemon;
using namespace arma;

struct DataStructureFunctionsATAC{
    
public:
    
    //MEMBER VARIABLES
    //The following variables describe the environment the system is tested in,
    //and generally do not change from run to run.

    int relaxationMethod;

    bool isRestart;

    bool generateRingHistogramImages;

    bool superfluousSp3Sp2;
    bool relaxAfterSp2Sp3Failure = 0;
    string numProcessors;
    bool haveNonInteracting;
    bool periodicOnlyOnNonInteracting;
    bool findChiralAngles = 0;
    bool findRingSizes = 0;
    bool acceptAnyMutation;
    int numOfNonInteracting;
    int numOfNonInteractingBonds;
    int numOfNonInteractingNeighbors;
    int acceptedRunNum = 0;
    int globalChiralAngleRunNum = 0; //this is the same as the number of accepted steps
    int globalRingSizeRunNum = 0; //this is the same as the number of accepted steps
    int namdRunNum = 0; //this is the number of times NAMD has run.
    int forceNum = 0;
    int spDenialPathSize;
    mat nonInteracting;
    imat nonInteractingBonds;
    imat nonInteractingNeighbors;
    double sp2tosp3Tolerance; // in angstroms. Experimentally derived : 2 A.
    int periodicStretching;
    double periodicStretchRate;
    bool isStaticPeriod;
    double xPeriodFactor;
    double yPeriodFactor;
    double zPeriodFactor;
    double xLatticeParameter; //NOTE: these parameters are used only if the periodic number is 1, 2, or 3.
    double yLatticeParameter;
    double zLatticeParameter;
    int numberOfElectrodes;
    imat xp1Electrodes; //positive x electrodes.
    imat xn1Electrodes; //negative x electrodes.
    imat yp1Electrodes; //etc
    imat yn1Electrodes;
    imat xp2Electrodes; //positive x electrodes adjacent to the xp1 electrodes.
    imat xn2Electrodes; //negative x electrodes adjacent to the xn1 electrodes.
    imat yp2Electrodes; //etc
    imat yn2Electrodes;
    double xpElecWallMod; //positive x electrode wall modifier (how they are selected). should be negative
    double xnElecWallMod; // should be positive
    double ypElecWallMod; // should be negative
    double ynElecWallMod; // should be positive
    int rotationsPerRun;
    double noFluxRatioShift;
    double ratioShiftRate;
    double ratioRotateFinal;
    double ratioSp2ToSp3Final;
    double reversableSpRunNum;
    double vecRoundingErrorCutOff;
    int periodicityNum;
    int minChiralIndex;
    int maxChiralIndex;
    string proteinDataBaseFile;
    int ringSizeMin;
    int ringSizeMax;


    bool haveBeam; //presence of a beam (electron)
    int beamType;
    //Type 0: we will only mutate dimers in the beam path, and system temperature is used.

    //Type 1: keep a seperate temperature for thermal mutations which can happen system-wide (temp),
    //and beam-assisted mutations which can only happen in the beam path (beam temp). Must also have a ratio of attempted thermal to beam mutations (perhaps a few happening at each step?)

    //Type 2, we do the same as in type 1. However, we submerge the beamed atoms in a heat bath, and update the temperature for the mutations dynamically.

    double beamTemperature;
    int beamAxis; //axis of beam: 1 is x, 2 is y, 3 is z.
    double beamRadius; //diameter, in angstroms
    double beamAxisIntercept1; // Where the beam's central axis intersects the coordinate plane perpendicular to it.
    double beamAxisIntercept2; // For shifting the beam's position. First value should be lower axis intercept, second higher.
    double beamRadiusGrowthRate;

    double tempS; //starting temperature, in kelvin.
    double tempF; // finishing/asymptotic temperature, in kelvin.
    int coolMethod; //cooling method: 0 is none (starting temp does not change), 1 is linear, 2 is Newtonian
    double coolR; //cooling rate constant in K/run for linear, and 1/runs for Newtonian.
    vec latticeVec1 = vec(3); //these lattice vectors are used for arbitrary lattices
    vec latticeVec2 = vec(3);
    vec latticeVec3 = vec(3);
    vec latticeVec1Norm = vec(3);
    vec latticeVec2Norm = vec(3);
    vec latticeVec3Norm = vec(3);

    //double mol = 6.02214e23;
    double boltzmann = 0.0019872041; //the boltzmann constant, in units of kilocalories/Kelvin/mol
    double lowerNNLimit = 2.3; //In angstroms. These are acceptable bounds, a range youll find the NNN atom in (of the same ring);
    double higherNNLimit = 2.5;
    double energyFailureTolerance = 1; //for NAMD.
    double bondDistanceToleranceMax;
    double bondDistanceToleranceMin;





// Builds a list of Neighbors.

imat buildNeigh(mat pos);

//declaring prototype function for mutation for mutationBody, so C++ can compile.

void mutation(int mutationNum, imat& neighList, imat& bonds, mat& positions, imat& remainingBonds, string masterBatchName, int& inBeam);

//mutationBody contains the actual sp2 topological manipulation algorithms for our mutation alphabet.

void mutationBody(int mutationNum, imat& neighList, imat& bonds, mat& positions, int& chosen, int& neighbor, int& NN1, int& NN2, int& NNN1, int& NNN2, imat& remainingBonds, int currentRemainingBond, string masterBatchName, int& inBeam) ;

//INCLUDE DOUBLE BOND CHECK
//checks if atoms i and j can go through an sp2 to sp3 condition to bond with
//each other.

bool sp2Tosp3Conditions(int i, int j, mat positions, imat& neighbors) ;

// the sp3 to sp2 transition occurs via two even-numbered rings (including size 2) that are connected face to face.
// same could be said for 2 to 3, but in reverse.
// Only rings of size 2 are considered, as the 3 to 2 barrier is significantly lower in such a case.

void sp3Tosp2(imat& neighbors, imat& bonds, mat& positions);

//This algorithm converts sp2 bonds to sp3 bonds, via the mechanism of
//forcing two double bonded pairs of atoms within close vicinity of each other.

//First, looks for two atoms that are within 2 angstroms, and not otherwise bonded.
//if these atoms are a part of mutually exclusive six rings,
//then look for neighbors that also fulfill the previous conditions.

//if so, bond atoms to atoms and neighbors to neighbors

//Might later modify slightly so it can bond edge atoms as well.

void sp2Tosp3(mat positions, imat& neighbors, imat& bonds, int& spPossibility) ;

//adds non-interacting atoms to our position and bond lists.

void nonInteractingAppender(mat& positions, imat& bondList, imat& neighborList) ;

//removes the non-interacting atoms from our position and bond lists,
//while updating the non-interacting component
//the non interacting bond list never changes.

void nonInteractingDeAppender(mat& positions, imat& bondList, imat& neighList) ;

//finds chirality based on diameter. Assumes nanotube is x aligned.

void findChiralityBasedOnDiameter(mat positions, imat bonds, imat neighList, bool haveNonInteracting) ;

void NAMDrun(int mutationGo, double& oldEnergy, double& temperature, mat& positions, imat& bondList, imat& neighborList, double& ratioRotate, double& ratioDimerDisplacement, double& ratioSp2ToSp3, int& spPossibility, string masterBatchName);

//the flux refers to mass. Only needs the rotate ratio, as the DV and add dimer operations occur on a 1:1 ratio.

void noFluxMutation(int& spPossibility, double rotateRatio, double ratioDimerDisplacement, double ratioSp2ToSp3, imat& neighList, imat& bonds, mat& positions, double& randomPercent, double& oldEnergy, double& temperature, imat& neighborList, string masterBatchName, int& inBeam) ;

//This simulates chemical vapor depostion, which at the time is simply add dimer.

void CVD(imat& neighList, imat& bonds, mat& positions, double& oldEnergy, double& temperature, imat& neighborList, string masterBatchName, int& inBeam) ;

//simple linear cooling. If we reach tempFinal, do not change temp.

double ratioShift(double& ratioRotate, double& ratioDimerDisplacment, double& ratioSp2ToSp3);

void recordChirality(vector<pair<imat, int> > rings, mat pos, bool findIndexes);


};


#endif	/* DATASTRUCTUREFUNCTIONSATAC_H */
