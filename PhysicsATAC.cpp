#include "PhysicsATAC.h"


//walls are perpendicular to the wall access direction.
//A positive force will be pushing inward, a negative one outward.

void forceParallelWalls(mat positions, mat& forceMatrix, double minWallPos, double maxWallPos, double force, int wallAxis) {
    for (int i = 0; i < forceMatrix.n_rows; i++) {
        for (int j = 0; j < forceMatrix.n_cols; ++j) {
            double Pos = positions(i, j);
            if (wallAxis == 1) {
                //determines force of minimum wall cap. IN X DIRECTION
                if (Pos < minWallPos && j == 0) {
                    forceMatrix(i, j) = forceMatrix(i, j) + force;
                }
                //determines force of maximum wall cap. IN X DIRECTION
                if (Pos > maxWallPos && j == 0) {
                    forceMatrix(i, j) = forceMatrix(i, j) - force;
                }
            }
            if (wallAxis == 2) {
                //determines force of minimum wall cap. IN Y DIRECTION
                if (Pos < minWallPos && j == 1) {
                    forceMatrix(i, j) = forceMatrix(i, j) + force;
                }
                //determines force of maximum wall cap. IN Y DIRECTION
                if (Pos > maxWallPos && j == 1) {
                    forceMatrix(i, j) = forceMatrix(i, j) - force;
                }
            }
            if (wallAxis == 3) {
                //determines force of minimum wall cap. IN Y DIRECTION
                if (Pos < minWallPos && j == 2) {
                    forceMatrix(i, j) = forceMatrix(i, j) + force;
                }
                //determines force of maximum wall cap. IN Y DIRECTION
                if (Pos > maxWallPos && j == 2) {
                    forceMatrix(i, j) = forceMatrix(i, j) - force;
                }
            }
        }
    }
}

//force of a tube with axis along x direction. walls currently only repulse like a lennard-jones potential, so slightly compressive.

void forceTube(mat positions, mat& forceMatrix, double force, double tubeRadius) {
    for (int i = 0; i < forceMatrix.n_rows; i++) {
        for (int j = 0; j < forceMatrix.n_cols; ++j) {
            //force in y direction
            forceMatrix(i, 1) += -force * (pow((tubeRadius - positions(i, 1)), -13) - pow((tubeRadius + positions(i, 1)), -13));
            //force in z direction
            forceMatrix(i, 2) += -force * (pow((tubeRadius - positions(i, 2)), -13) - pow((tubeRadius + positions(i, 2)), -13));

        }
    }
}

//Force of two flat walls perpendicular to the z direction.

void forceVice(mat positions, mat forceMatrix, double minViceWall, double maxViceWall, double viceForce) {
    for (int i = 0; i < forceMatrix.n_rows; i++) {
        for (int j = 0; j < forceMatrix.n_cols; ++j) {
            double Pos = positions(i, j);
            //determines force of minimum wall cap. IN Z DIRECTION
            if (Pos < minViceWall && j == 2) {
                forceMatrix(i, j) = forceMatrix(i, j) + viceForce;
            }
            //determines force of maximum wall cap. IN Z DIRECTION
            if (Pos > maxViceWall && j == 2) {
                forceMatrix(i, j) = forceMatrix(i, j) - viceForce;
            }
        }
    }
}

//simple Newtonian cooling

double newtonianCooling(double& tempS, double& tempF, int& acceptedRunNum, double& coolR) {
    return tempF + (tempS - tempF) * exp(-acceptedRunNum * coolR);
}

//simple linear cooling. If we reach tempFinal, do not change temp.

double linearCooling(double& tempS, double& tempF, int& acceptedRunNum, double& coolR) {
    if (tempS - acceptedRunNum * coolR < tempF) {
        return tempF;
    } else {
        return tempS - acceptedRunNum*coolR;
    }
}





//generates three random percentages, as a vector of doubles

vector<double> threeRandomPercentages() {

    double randsmallest;
    double randlargest;
    //generates two random numbers between 0 and 1
    //srand(time(NULL));
    double random1 = (double) rand() / (double) (RAND_MAX);
    double random2 = (double) rand() / (double) (RAND_MAX);

    if (random1 > random2) {
        randlargest = random1;
        randsmallest = random2;
    } else {
        randlargest = random2;
        randsmallest = random1;
    }

    vector<double> ratios(3);
    ratios[0] = randsmallest;
    ratios[1] = (randlargest - randsmallest);
    ratios[2] = (1 - randlargest);
    return ratios;
}


// oldEnergy is the energy of the previous structure.
// If the new allotrope of carbon is energetically favorable, it passes the energy check.
// If it does not pass the check, there is still a chance (arrhenius) that
// it will take on the metastable structure.
// Updates oldEnergy with NewEnergy if passes.
// Energy units are kcal/mol
// CHECK IF POSITIONS IS NEEDED, KCAL PER MOL?

bool energyCheck(double& oldEnergy, double& temperature, int spPossibility, double randomPercent, double ratioRotate, double ratioDimerDisplacement, string masterBatchName, int& namdRunNum, bool& relaxAfterSp2Sp3Failure, double& boltzmann, double& energyFailureTolerance, bool& acceptAnyMutation, int& relaxationMethod) {

    bool pass;
    string runNumS;
    stringstream ssRun;
    ssRun << namdRunNum;
    runNumS = ssRun.str();


    if (spPossibility == -1 && relaxAfterSp2Sp3Failure == false) {
        pass = false;
        return pass;
    }

    double newEnergy;
    string oldEnergyS;
    string newEnergyS;
    string arrheniusS;
    string randomPercentS;
    string numeratorS;
    string denominatorS;

    //initilizes new energy
    energyReader(newEnergy, masterBatchName, relaxationMethod);

    //next block of code converts doubles to strings
    stringstream ssOld;
    ssOld << oldEnergy;
    oldEnergyS = ssOld.str();
    stringstream ssNew;
    ssNew << newEnergy;
    newEnergyS = ssNew.str();
    double numerator = ((-newEnergy + oldEnergy));
    double denominator = (boltzmann * temperature);
    //find the temperature the simulation is opperating at. newEnergy units are kcal/mol.
    double arrhenius = exp(numerator / denominator);
    stringstream ssArr;
    ssArr << (arrhenius * 100);
    arrheniusS = ssArr.str();
    stringstream ssNumer;
    ssNumer << numerator;
    numeratorS = ssNumer.str();
    stringstream ssDem;
    ssDem << denominator;
    denominatorS = ssDem.str();
    //current temp:
    stringstream ssTemp1;
    ssTemp1 << temperature;
    string tempStr1 = ssTemp1.str();

    //FOLLOWING IS REDUNDANT SECTION THAT CAN BE REMOVED

    //check for NAMD failure, which generates energies less than one (usually from .tcl file)
    if (newEnergy < energyFailureTolerance) {

        cout << "NAMD FAILED" << endl;
        cout << "OLD ENERGY; " << oldEnergyS << endl;
        cout << "NEW ENERGY: " << newEnergyS << endl;
        cout << "ARRHENIUS PERCENT THRESHOLD (arrhenius * 100): " << arrheniusS << endl;
        cout << "RANDOM PERCENT ATTEMPTED: " << randomPercent * 100 << endl;
        cout << "NUMERATOR: " << numeratorS << endl;
        cout << "DENOMINATOR " << denominatorS << endl;
        cout << "TEMPERATURE " << tempStr1 << endl;

        //system("read  -p \"namd failed, paused\" key");
        pass = false;
        return pass;
    }
    //PREVIOUS IS REDUNDANT SECTION THAT CAN BE REMOVED.

    // cout << "ENERGIES " << newEnergy << " " << oldEnergy << endl;
    //check if energy is lower,
    if ((newEnergy <= oldEnergy)) {



        cout << "NEW ALLOTROPE ENERGETICALLY FAVORABLE" << endl;
        cout << "OLD ENERGY; " << oldEnergyS << endl;
        cout << "NEW ENERGY: " << newEnergyS << endl;
        cout << "ARRHENIUS PERCENT THRESHOLD (arrhenius * 100): " << arrheniusS << endl;
        cout << "RANDOM PERCENT ATTEMPTED: " << 100 << endl;
        cout << "NUMERATOR: " << numeratorS << endl;
        cout << "DENOMINATOR " << denominatorS << endl;
        cout << "TEMPERATURE " << tempStr1 << endl;


        oldEnergy = newEnergy;
        pass = true;
    } else {
        double randomPercent = (double) rand() / (double) (RAND_MAX);




        if (acceptAnyMutation == true) {
            //the immediate next line will always allow success, no matter how energetic the result.
            randomPercent = arrhenius * 0.5;
        }

        //sees if the meta stable form is taken. Lower energy = higher arhenius = higher chance of random percent being less than that = higher chance of a pass.
        stringstream ssRand;
        ssRand << (randomPercent * 100);
        randomPercentS = ssRand.str();
        //if the random percent is lower than the arrhenius, or if a successful sp2 to sp3 transition...
        if (randomPercent <= arrhenius || (spPossibility == 1)) {
            cout << "NEW ALLOTROPE UNFAVORABLE, ARRHENIUS PASS" << endl;
            cout << "OLD ENERGY; " << oldEnergyS << endl;
            cout << "NEW ENERGY: " << newEnergyS << endl;
            cout << "ARRHENIUS PERCENT THRESHOLD (arrhenius * 100): " << arrheniusS << endl;
            cout << "RANDOM PERCENT ATTEMPTED: " << randomPercent * 100 << endl;
            cout << "NUMERATOR: " << numeratorS << endl;
            cout << "DENOMINATOR " << denominatorS << endl;
            cout << "TEMPERATURE " << temperature << endl;


            oldEnergy = newEnergy;
            pass = true;
        } else {
            cout << "NEW ALLOTROPE UNFAVORABLE, ARRHENIUS FAIL" << endl;
            cout << "OLD ENERGY; " << oldEnergyS << endl;
            cout << "NEW ENERGY: " << newEnergyS << endl;
            cout << "ARRHENIUS PERCENT THRESHOLD (arrhenius * 100): " << arrheniusS << endl;
            cout << "RANDOM PERCENT ATTEMPTED: " << randomPercent * 100 << endl;
            cout << "NUMERATOR: " << numeratorS << endl;
            cout << "DENOMINATOR " << denominatorS << endl;
            cout << "TEMPERATURE " << temperature << endl;


            pass = false;
        }
    }
    return pass;
}






