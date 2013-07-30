#ifndef SYSTEMATAC_H
#define	SYSTEMATAC_H


#include <sstream>
#include <math.h>
#include <iostream>
#include <stdio.h>//For chdir()
#include <unistd.h> //For chdir())
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <fstream>
#include "armadillo"
#include <time.h>// needed to timestamp logs
#include <vector> // needed for index sort
#include <utility>// needed for index sort
#include <algorithm>// needed for index sort
#include <lemon/list_graph.h>//lemon
#include <lemon/matching.h>//lemon
#include <iomanip> // for setprecision()
#include "DataStructureFunctionsATAC.h"
#include "ReadWriteFunctionsATAC.h"
#include "GeometryFunctionsATAC.h"
#include "PhysicsATAC.h"


using namespace arma;
using namespace lemon;
using namespace std;


//The system is a environment (the structure is a system is more appropriate),
//so system inherits from the environment.

struct systemATAC: public DataStructureFunctionsATAC {
public:


    //MEMBER VARIABLES

    //These variables are particular to the structure of the system at a given frame.


    mat positions;
    imat neighborList;
    imat bondList;
    int spPossibility = 0;
    double oldEnergy = 0;
    double newEnergy;
    double temperature;
    double ratioRotate;
    double ratioDimerDisplacement;
    double ratioSp2ToSp3;
    string batchName;

    int mutationGo;
    







    //MEMBER FUNCTIONS    


    void initialization(int argc, char* argv[]);


    int runSimulation();

};




#endif	/* ATAC_H */

