/*
 ATAC Beta 0.01 

 This program is the Accelerated Topological Annealing of Carbon (or ATAC).
 It is used to see the various transformation processes of carbon under
 irradiation and heat treatment. Please see the user manual on the group wiki at:
 
 http://icmp.phys.rpi.edu/index.php/ATAC_User_Manual

 ***Dependencies:*** 
 Requires the armadillo and LEMON libraries for base functionality.
 A relaxation method must be accessible as well. Currently programs include:
 * NAMD
 * Tersoff-Brenner (in house ICMP code).
 Gnuplot is required for some optional abilities.
 
 ***Files required in working directory:***
 * An xyz file of your system that MUST be named "pos.xyz," containing cartesian coordinates.
 * A file containing all the appropriately named variables that MUST be named "variables.txt"
 * If you are using NAMD, a protein database parameter file, which can be named as anything.

 One command line argument is required: a batch name identifier.
 The batch name identifier is not the name of the system; 
 this is used to differentiate different jobs running in parallel.

 ***General program structure*** 
 The core data structures in this calculation:
 * Is a cartesian x y z postion matrix, which is 3 by number of atoms in system
 * A bond matrix, which is 2 by number of bonds. Atoms are identified by their row number
   the position list.
 * A neighbor matrix, which is 4 by number of atoms. Row number corresponds to the atom
   of the same row on the position matrix, and column values are the atom indexes
   of the atoms its bonded too. A -1 index indicates lack of a neighbor.

 A mutation is then induced in the bond table.
 
 If NAMD is selected:
 With this information, a .top and a .pgn file are created.
 These files are feed into psfgen which to creates a .psf and .pdb file.
 NAMD uses the psf and pdb files to relax the structure, and generate an energy.
 
 If Tersoff-Brenner is selected:  
 
 Based on the generated energy, the mutation is either accepted or rejected.

 .xyz movie is generated: the first frame is the given xyz file, and the rest
 are minimizations/mutations.

 ***Atomic System Input Limitations:***
  
 Requires at least 5 atoms. 
  
 Code will not handle the removal of atoms if the system is composed of just
 two abutting rings, with one of them a three member triangle. This situation is
 so energetically unfavorable, however, it is very unlikely to occur.

 ***Coding Style Notes*** 
 I've tried to adhere to a specific style throughout this code for
 easier readability.
 * Capitalization is in camel caps.
 * Matrixes with 1 dimension are row major (i.e. many rows but 1 column).
 * Matrixes with positions are row major wrt the atom, and column major wrt coordiantes. i.e.

  1.132 2.123 -1.232 \\coordinates of atom 1
  -1.321 2.333 6.323 \\coordinates of atom 2
  9.234 1.323 -9.232 \\ etc.

 ***Contact:***
 Any questions or comments can be sent to: 
 bullaz@rpi.edu 
  
 ~Zac Bullard 7-29-2014
 */

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
#include "systemATAC.h"

using namespace arma;
using namespace lemon;
using namespace std;

//USES -std=c++11 flag 

int main(int argc, char* argv[]) {

    //create system.
    systemATAC system;




    // initialize and set system environment.
    system.initialization(argc, argv);


    //start simulation

    system.runSimulation();



    return 0;
}
