

#include "DataStructureFunctionsATAC.h"





// Builds a list of Neighbors.

imat DataStructureFunctionsATAC::buildNeigh(mat pos) {

    int posn_rows = pos.n_rows;
    int maxNumberOfNeighbors = 4;

    imat neigh(posn_rows, maxNumberOfNeighbors);
    /* row is atom number on original position list, columns are 4
    nearist neighbors. Negative -1 indicates no neighbor, i.e. an edge atom
     -2 indicates the proper values have not been assigned yet.*/
    neigh.fill(-2);

    //keeping this variable here in case of more generalization later.

    //initialize with junk value
    // imat bondedAtoms;
    //bondedAtoms << -1 << endr;

    for (int i = 0; i < posn_rows; ++i) // Iterates position atom.
    {

        // mat dist;
        // dist << 0 << -1 << endr;
        //first row is closest, etc. intialize with junk values.
        // mat sortedDist;


        //first column is distance, second column is neighbor number.
        mat dist(maxNumberOfNeighbors, 2);
        dist << 0 << -1 << endr
                << 0 << -1 << endr
                << 0 << -1 << endr
                << 0 << -1 << endr;
        int numberOfAcceptedNeighbors = 0;


        if (i % 100 == 0) {
            cout << "Finding neighbors for atoms: " << i << "-" << i + 99 << endl;
        }
        //cout << "Finding neighbors for atom " << i << endl;

        for (int j = i + 1; j < posn_rows; ++j) {
            // Iterates atom to check distance wrt.
            //The reason that j is always greater than i, is that we will find and force all neighbors for i as we go along
            //i.e., even if the jth atom has neighbors that are a better fit than i, too bad.
            //essentially a greedy algorithm that can be improved upon.
            // j - (i + 1) = zero indexed number of atoms we have tested for j so far.
            //  if(j%100 == 0){
            //    cout << "and atoms " << j << "-" << j+99 << endl;
            // }

            double distance;


            shortestLineDeterminer(pos, i, j, distance, periodicityNum, xLatticeParameter, yLatticeParameter, zLatticeParameter, latticeVec1Norm, latticeVec1, latticeVec2Norm, latticeVec2, latticeVec3Norm, latticeVec3);
            //cout << "AFTER FINDING DIST " << endl;


            //accept as neighbor if distance is between the ranges given, and we have not already bonded this atom before.
            if (distance < bondDistanceToleranceMax && distance > bondDistanceToleranceMin /*&& contains(bondedAtoms,j) == 0*/) {
                //       cout << "before dist setting " << endl;

                dist(numberOfAcceptedNeighbors, 0) = distance;
                dist(numberOfAcceptedNeighbors, 1) = j;
                //  cout << "after setting " << endl;
                //iterate accepted Neighbors, and add atom to bonded atom list.


                imat toBeAdded;
                toBeAdded << j;
                //  bondedAtoms = join_cols(bondedAtoms,toBeAdded);

                numberOfAcceptedNeighbors++;


                //If we have enough satisfactory neighbors, we can move onto the next atom.
                if (numberOfAcceptedNeighbors == maxNumberOfNeighbors) {
                    break;
                }

            }

            //fill out our distance matrix with the first initial guesses.
            //    if (j < maxNumberOfNeighbors){

            //    dist((j-i-1), 0) = distance;
            //    dist((j-i-1), 1) = j;

            //  continue;
            //  }

            //we can sort our list of guesses now.
            //                             if (j == maxNumberOfNeighbors){
            //    sortedDist = betterSortRows(dist, 0);
            //                           }

        }


        /*
        //loop through satisfying atoms.
                for (int s = 0; s < maxNumberOfNeighbors; s++) {
                    //loop over neighbor atoms for atom i
                          for (int n = 0; n < neigh.n_cols; n++) {
                    //dist(k,1) is a satisfing atom.

                              //
                             imat neighRow = neigh.row(i);
                    if (contains(neighRow,dist(s, 1))) {
                        continue

                        neigh(i, s) = dist(s, 1);
                    }

                    if (neigh(dist(s, 1), maxNumberOfNeighbors - 1 - s) == -2) {
                        neigh(dist(s, 1), maxNumberOfNeighbors - 1 - s) = i;
                    }
                          }
                }
         */


        ///*
        //now add our found neighbors to the neighbor list.
        //since -1 means no neighbor, we don't need to worry about clearing junk values from dist, as it would only have junk values
        //if there weren't enough atoms to satisfying all the possible neighbors.
        // k +1 should equal maxNumberOfNeighbors equals dist.n_rows.
        //loop through satisfying atoms.

        bool foundAllNeighs = false;
        int satifying = 0;
        //loop over neighbor atoms for atom i
        for (int n = 0; n < maxNumberOfNeighbors; n++) {
            //break if we have found the max number of neighbors, or if the satifying atom does not exist.
            if (foundAllNeighs == true || dist(satifying, 1) == -1) {
                break;
            }

            //need to make sure we aren't overwriting a pre-existing bond/
            if (neigh(i, n) == -2) {
                //now to make sure that the satisfying atom can accept atom i
                for (int n2 = 0; n2 < maxNumberOfNeighbors; n2++) {
                    if (neigh(dist(satifying, 1), n2) == -2) {
                        //dist(k,1) is a satisfing atom.
                        neigh(i, n) = dist(satifying, 1);

                        neigh(dist(satifying, 1), n2) = i;

                        //cout <<"YESS HAAA i, n, and dist(satifying, 1), n2: " << i << " " << n << " " << dist(satifying, 1) << " " << n2 << endl;
                        //cout << neigh << endl;
                        //system("read -p \"sdklfjdk \" key");


                        satifying++;

                        if (satifying == maxNumberOfNeighbors) {
                            foundAllNeighs = true;
                        }

                        break;

                    }
                }
            }

        }

        // */

    }

    //need to turn all undetermined values into -1 now.
    for (int i = 0; i < posn_rows; i++) {
        for (int j = 0; j < maxNumberOfNeighbors; j++) {
            if (neigh(i, j) == -2) {
                neigh(i, j) = -1;
            }
        }
    }


    return (neigh); // first column has closest atom, second has second, etc.
}





//these prototype functions are for debugging; They allow for an intermediate view of the structure,esp after adding a dimer.
//void writeTopFile(mat positions, imat bondList, string masterBatchName);
//void writePgnFile(mat positions, string masterBatchName);
//void writePdbPsfFiles(string masterBatchName);

//declaring prototype function for mutation for mutationBody, so C++ can compile.
//void DataStructureFunctionsATAC::mutation(int mutationNum, imat& neighList, imat& bonds, mat& positions, imat& remainingBonds, string masterBatchName, int& inBeam);


//mutationBody contains the actual sp2 topological manipulation algorithms for our mutation alphabet.

void DataStructureFunctionsATAC::mutationBody(int mutationNum, imat& neighList, imat& bonds, mat& positions, int& chosen, int& neighbor, int& NN1, int& NN2, int& NNN1, int& NNN2, imat& remainingBonds, int currentRemainingBond, string masterBatchName, int& inBeam) {
    imat toBeAppended;
    int atomAddedToNNRingAfterRotation;
    int atomAddedToNNNRingAfterRotation;
    int NNnotFound = -1;
    int NNNFound = -1;
    int NNNnotFound = -1;
    int NNFound = -1;
    int beforeEnd;
    int afterEnd;
    int prospectiveNeighbor;
    int randomAdvancement = 0;

    //mutation for rotate
    if (mutationNum == 1) {
        cout << "ENTERING MUTATION BODY ROTATE " << endl;

        //Immediately, we must check if the chosen and neighbor can be in a double bond. if no double bond, call the mutation again.
        if (isDoubleBond(neighList, chosen, neighbor) == false) {

            remainingBonds.shed_row(currentRemainingBond);
            if (remainingBonds.is_empty()) {
                cout << "NO MUTATION POSSIBLE, FAILED AT TRYING TO FIND DOUBLE BOND FOR ROTATION " << endl;
                return;
            }

            mutation(mutationNum, neighList, bonds, positions, remainingBonds, masterBatchName, inBeam);
            return;
        }

        //FOLLOWING CODE USES RINGS TO FIND NEIGHBORS. AT THE MOMENT, COMPUTATIONALLY EXPENSIVE.

        //first tries to find a ring with the chosen and neighbor atoms as one side.
        //if this fails, finds atoms to switch by checking the positions.


        /*
        //initilize our stub of a ring, which is the chosen atom and neighbor
        imat ring;
        ring << chosen << neighbor << endr;
        imat winningRing;


        bool victory = false;
        int ringSizeMax = atoi(variableReader("ringSizeMax").c_str());


        // must start depth at 0, as 1 is added at the beginning.
        findSmallestRingForAtom(neighbor, ring, winningRing, neighList, 0, victory, ringSizeMax);

        //if a ring is found, find the found and not found atoms
        if (winningRing.is_empty() == false) {

            //now we must determine which ring is the Found Ring.
            uvec NN1Index = find(winningRing == NN1);
            uvec NN2Index = find(winningRing == NN2);
            uvec NNN1Index = find(winningRing == NNN1);
            uvec NNN2Index = find(winningRing == NNN2);

            if (NN1Index.is_empty() == 0) {
                NNnotFound = NN2;
                NNFound = NN1;
            } else {
                NNnotFound = NN1;
                NNFound = NN2;
            }

            if (NNN1Index.is_empty() == 0) {
                NNNFound = NNN1;
                NNNnotFound = NNN2;
            } else {
                NNNFound = NNN2;
                NNNnotFound = NNN1;
            }
        } else */

        {//if a ring is not found, we must use spatial coordinates.
            //The distance between NNFound and NNNFound in a normal graphene sheet is about 3.07 A...
            //The distance between NNFound and NNNnotFound in a normal graphene sheet is about 3.77 A...
            double distance1;
            shortestLineDeterminer(positions, NN1, NNN1, distance1, periodicityNum, xLatticeParameter, yLatticeParameter, zLatticeParameter, latticeVec1Norm, latticeVec1, latticeVec2Norm, latticeVec2, latticeVec3Norm, latticeVec3);

            double distance2;
            shortestLineDeterminer(positions, NN1, NNN2, distance2, periodicityNum, xLatticeParameter, yLatticeParameter, zLatticeParameter, latticeVec1Norm, latticeVec1, latticeVec2Norm, latticeVec2, latticeVec3Norm, latticeVec3);

            //if the distance from NN1 to NNN1 is shortest, then NNN1 would be considered a part of the (imaginary) found ring.
            // likewise for NN1 and NNN2
            if (distance1 <= distance2) {
                NNFound = NN1;
                NNnotFound = NN2;
                NNNFound = NNN1;
                NNNnotFound = NNN2;
            } else {
                NNFound = NN1;
                NNnotFound = NN2;
                NNNFound = NNN2;
                NNNnotFound = NNN1;
            }

        }

        //first index is row in bonds, second is column in bonds.
        vec NNnotFoundIndex(2);
        vec NNNFoundIndex(2);
        vec NNFoundIndex(2);
        vec NNNnotFoundIndex(2);

        //loops over rows in bonds to find indexes of those to be switched.
        for (int j = 0; j < bonds.n_rows; ++j) {

            if ((bonds(j, 0) == NNFound) && (bonds(j, 1) == chosen)) {
                NNFoundIndex(0) = j;
                NNFoundIndex(1) = 0;
            }
            if ((bonds(j, 1) == NNFound) && (bonds(j, 0) == chosen)) {
                NNFoundIndex(0) = j;
                NNFoundIndex(1) = 1;
            }
            if ((bonds(j, 0) == NNnotFound) && (bonds(j, 1) == chosen)) {
                NNnotFoundIndex(0) = j;
                NNnotFoundIndex(1) = 0;
            }
            if ((bonds(j, 1) == NNnotFound) && (bonds(j, 0) == chosen)) {
                NNnotFoundIndex(0) = j;
                NNnotFoundIndex(1) = 1;
            }
            if ((bonds(j, 1) == NNNFound) && (bonds(j, 0) == neighbor)) {
                NNNFoundIndex(0) = j;
                NNNFoundIndex(1) = 1;
            }
            if ((bonds(j, 0) == NNNFound) && (bonds(j, 1) == neighbor)) {
                NNNFoundIndex(0) = j;
                NNNFoundIndex(1) = 0;
            }
            if ((bonds(j, 1) == NNNnotFound) && (bonds(j, 0) == neighbor)) {
                NNNnotFoundIndex(0) = j;
                NNNnotFoundIndex(1) = 1;
            }
            if ((bonds(j, 0) == NNNnotFound) && (bonds(j, 1) == neighbor)) {
                NNNnotFoundIndex(0) = j;
                NNNnotFoundIndex(1) = 0;
            }
        }


        //imat bondsBefore = bonds;
        //switches bonds. Randomly determines which to switch (randomizes rotation direction)
        int randomDirection = rand() % 2;
        //NOT RANDOM, CHANGE LATER
        randomDirection = 1;
        if (randomDirection == 0) {
            bonds(NNnotFoundIndex(0), NNnotFoundIndex(1)) = NNNFound;
            bonds(NNNFoundIndex(0), NNNFoundIndex(1)) = NNnotFound;

        } else {
            bonds(NNNnotFoundIndex(0), NNNnotFoundIndex(1)) = NNFound;
            bonds(NNFoundIndex(0), NNFoundIndex(1)) = NNNnotFound;
        }

        // These values might not be nessesarily used. The renaming is arbitrary
        //SHOULD BE RANDOM BETWEEN NEIGHBOR AND CHOSEN, IS NOT, CHANGE LATER
        atomAddedToNNRingAfterRotation = neighbor;
        atomAddedToNNNRingAfterRotation = chosen;

        //expensive routene, replace if speed benefits are really needed.
        buildNeighFromBonds(bonds, neighList);

        //does not need to rebuild bonds
    }

    //mutation for remove dimer
    if (mutationNum == 2) {

        for (int i = 0; i < bonds.n_rows; i++) {
            for (int j = i; j < bonds.n_rows; j++) {
                if (j != i) {
                    if ((bonds(i, 0) == bonds(j, 0) && bonds(i, 1) == bonds(j, 1))
                            || (bonds(i, 0) == bonds(j, 1) && bonds(i, 1) == bonds(j, 0))) {
                        cout << "OFFENDING BONDS: THE FIRST " << bonds.row(i) << "THE SECOND " << bonds.row(j) << endl;
                        cout << "OFFENDING ROW NUM: " << i << " " << j << endl;
                        system("read  -p \"paused, same bond counted twice BEFORE REMOVAL.\" key");
                    }
                }
            }
        }


        //Immediately, we must check if the chosen and neighbor can be in a double bond. if no double bond, call the mutation again.
        if (isDoubleBond(neighList, chosen, neighbor) == false) {
            remainingBonds.shed_row(currentRemainingBond);
            if (remainingBonds.is_empty()) {
                cout << "NO MUTATION POSSIBLE, FAILED AT TRYING TO FIND DOUBLE BOND FOR REMOVE DIMER (GRAPH THEORY) " << endl;
                return;
            }
            mutation(mutationNum, neighList, bonds, positions, remainingBonds, masterBatchName, inBeam);
            return;
        }

        //adjusts neighbor indexes to account for a shorter neighbor matrix
        int greaterIndex;
        int lesserIndex;
        if (chosen > neighbor) {
            greaterIndex = chosen;
            lesserIndex = neighbor;
        } else {
            greaterIndex = neighbor;
            lesserIndex = chosen;
        }

        //handy code
        /*
    writeTopFile(positions, bonds, batchName);
    writePgnFile(positions, batchName);
    writePdbPsfFiles(batchName);

    //Handy debugging code colors the ring in question; this pdb file is kept seperate, in batchnameModified.pdb
    system(("rm " + batchName + "Modified.pdb").c_str());
    string line;
    ifstream original((batchName + ".pdb").c_str());
    ofstream modified((batchName + "Modified.pdb").c_str());
    int lineNumber = 0;
    while (getline(original, line)) {

        //does not modify the first or last line.
        //looks for only atoms in the winning ring.
        if ((lineNumber > 0) && (line.size() > 4) && (chosen == lineNumber - 1 || neighbor == lineNumber - 1)) {
            line.erase(77, 79);
            line = line + "N";

        }
        modified << line << "\n";
        ++lineNumber;
    }
    modified.close();
    original.close();
    /// */



        //links NNs to each other, and also the NNNs
        for (int i = 0; i < neighList.n_cols; ++i) {

            if (neighList(NN1, i) == chosen) {
                neighList(NN1, i) = NN2;
            }
            if (neighList(NN2, i) == chosen) {
                neighList(NN2, i) = NN1;
            }

            if (neighList(NNN1, i) == neighbor) {
                neighList(NNN1, i) = NNN2;
            }
            if (neighList(NNN2, i) == neighbor) {
                neighList(NNN2, i) = NNN1;
            }
        }

        //loops over neighbor list, adjusting for the soon to be missing rows
        for (int j = 0; j < neighList.n_rows; ++j) {
            for (int i = 0; i < neighList.n_cols; ++i) {
                if (neighList(j, i) > greaterIndex) {
                    neighList(j, i) = neighList(j, i) - 2;
                    continue;
                }
                if ((neighList(j, i) > lesserIndex) && (neighList(j, i) < greaterIndex)) {
                    neighList(j, i) = neighList(j, i) - 1;
                    continue;
                }
            }
        }
        //loops over neighbor list, adjusting indexes

        // sheds the rows corrosponding to the divacancy atoms. Need to shed rows from larger to smaller, or indexes will be off (falling stack of boxes or sand analogy);
        positions.shed_row(greaterIndex);
        positions.shed_row(lesserIndex);


        //sheds atoms in neighbor list.
        neighList.shed_row(greaterIndex);
        neighList.shed_row(lesserIndex);

        /*
        //following block of code should work. replace the build bonds line
         * on. possible speed benefits.
        //update bonds
        int bondsRowMax = bonds.n_rows;
        cout << "chosend neighb" << chosen << " " << neighbor << endl;
        for(int i = 0; i < bondsRowMax; ++i){
            // delete all bonds containing the removed atoms
            if ((bonds(i,0)== chosen)||(bonds(i,0)== neighbor)||
                    (bonds(i,1)== chosen)||(bonds(i,1)== neighbor)){
                bonds.shed_row(i);
                --bondsRowMax;
                --i;
            }
        }
        //add two new chasm-crossing bonds
        imat toBeAdded;
        toBeAdded << NN1 << NN2 << endr << NNN1 << NNN2 << endr;
        bonds = join_cols(bonds,toBeAdded);
                                //loops over bond list, adjusting for the missing rows
        for (int j = 0; j < bonds.n_rows; ++j) {
            for (int i = 0; i < 2; ++i) {
                if (bonds(j, i) > greaterIndex) {
                    bonds(j, i) = bonds(j, i) - 2;

                }
                if ((bonds(j, i) > lesserIndex) && (bonds(j, i) < greaterIndex)) {
                    bonds(j, i) = bonds(j, i) - 1;

                }
            }
        }
         */
        //if bonds are adjusted directly with previous bit of code, might result
        //in speed benefits.
        bonds = buildBond(neighList);


        for (int i = 0; i < bonds.n_rows; i++) {
            for (int j = i; j < bonds.n_rows; j++) {
                if (j != i) {
                    if ((bonds(i, 0) == bonds(j, 0) && bonds(i, 1) == bonds(j, 1))
                            || (bonds(i, 0) == bonds(j, 1) && bonds(i, 1) == bonds(j, 0))) {
                        cout << "OFFENDING BONDS: THE FIRST " << bonds.row(i) << "THE SECOND " << bonds.row(j) << endl;
                        cout << "OFFENDING ROW NUM: " << i << " " << j << endl;
                        system("read  -p \"paused, same bond counted twice AFTER REMOVAL.\" key");
                    }
                }
            }
        }



    }

    //This mutation for add dimer. The two bonds that the dimer intersects must NOT be double bonds.
    if (mutationNum == 3) {


        //initilize our stub of a ring, which is the chosen atom and neighbor
        //imat ring;
        //ring << chosen << neighbor << endr;

        imat winningRing;

        vector<pair<imat, int> > ringsOfSizes = findUniqueRingsForDimer(neighList, chosen, neighbor, ringSizeMin, ringSizeMax);


        //dimer might not have any rings.
        if (ringsOfSizes.empty()) {
            remainingBonds.shed_row(currentRemainingBond);
            mutation(mutationNum, neighList, bonds, positions, remainingBonds, masterBatchName, inBeam);
            return;
        }

        // ADD CHECK THAT AT LEAST ONE RING EXISTS; ELSE THIS PROGRAM WILL RUN INFINITELY
        //This algoritm DOES require at least one ring to be found. if it is not found, recursively call mutation
        if (ringsOfSizes.empty() == 1) {
            //cout << "DIDN'T FIND RING, FAILED" << endl;
            //calls mutation again to attempt another chosen and neighbor
            remainingBonds.shed_row(currentRemainingBond);
            if (remainingBonds.is_empty()) {
                cout << "NO MUTATION POSSIBLE, FAILED AT TRYING TO FIND RING FOR ADD DIMER " << endl;
                return;
            }

            mutation(mutationNum, neighList, bonds, positions, remainingBonds, masterBatchName, inBeam);
            return;
        }

        //  cout << "FLAG " << endl;

        //we must choose a winning ring out of ringsOfSizes. We first must find total number of rings.
        int totalNumberOfRings = 0;
        for (int i = 0; i < ringsOfSizes.size(); ++i) {

            //      cout << "Size of rings: " << ringsOfSizes[i].second << ". Number of rings: " << ringsOfSizes[i].first.n_rows << endl << "endl" << endl;
            // if (ringsOfSizes[i].first.n_rows != 0) {
            // cout << ringsOfSizes[i].first;
            // }
            //no not count empty rings and rings of 3.
            if (ringsOfSizes[i].first.n_rows != 0 && ringsOfSizes[i].first.n_cols != 3) {
                //cout << "TOTAL NUMBER OF RINGS " << totalNumberOfRings << endl;

                totalNumberOfRings = totalNumberOfRings + ringsOfSizes[i].first.n_rows;
            }
        }

        //cout << "TOTAL NUM OF RINGS " << totalNumberOfRings << endl;

        for (int i = 0; i < ringsOfSizes.size(); i++) {
            for (int j = 0; j < ringsOfSizes[i].first.n_rows; j++) {
                //      cout << "RINGS FOR DIMER: " << ringsOfSizes[i].first.row(j) << endl;
            }
        }

        //choose a random, acceptable ring to be the winner.
        int randomRing = rand() % totalNumberOfRings;

        // keep track of how many rings you've seen.

        int runningTotal = 0;
        bool victory = false;
        //advance through rings until you reach the random ring.
        //advance through ring sizes first

        for (int i = 0; i < ringsOfSizes.size(); ++i) {
            //then advance through rings in that size.
            //cout << "RINGS SIZE SIZE " << ringsOfSizes.size() << "num of rings in this size " << ringsOfSizes[i].first.n_rows << endl;
            //cout << "Running total and random ring  " << runningTotal << " " << randomRing << endl;

            //we will not consider addind dimers to triangular formations
            if (ringsOfSizes[i].second != 3) {
                //cout << ringsOfSizes[i].second << "RING SIZE THE SECOND " << endl;
                for (int j = 0; j < ringsOfSizes[i].first.n_rows; ++j) {

                    //             cout << "I J " << i << " " << j << endl;
                    //           cout << "Running total and random ring  " << runningTotal << " " << randomRing << endl;
                    if (runningTotal == randomRing) {
                        winningRing = ringsOfSizes[i].first.row(j);
                        //             cout << "Victory" << endl;
                        victory = true;
                        break;
                    } else {
                        ++runningTotal;
                    }
                }

                //victory should always eventually be true, do to check that
                //ringsOfSizes was not empty earlier on.
                if (victory == true) {
                    break;
                }
            } else {
                //   cout << "RING OF SIZE THREE CHOSEN " << endl;
            }
        }



        //   cout << "FLAGyes " << endl;


        //advance around the ring by a random amount. Must eliminate the possibility of advancing 0.
        // 3 rings forbidden...



        int numOfRandomAdvancements = winningRing.n_cols - 3;
        imat advancementList;
        //initilize junk value
        advancementList << -1 << endr;


        //  cout << "numOfRandomAdvancements " << numOfRandomAdvancements << endl;
        //    cout << "winning ring " << winningRing << endl;
        for (int i = 0; i < numOfRandomAdvancements; ++i) {
            imat toBeAdded;
            toBeAdded << i + 2 << endr;
            advancementList = join_cols(advancementList, toBeAdded);
        }
        //shed junk value
        advancementList.shed_row(0);
        // cout << "FLAGmider " << endl;
        //randomAdvancement = (rand() % (winner.n_cols - 3)) + 2;

        int addDimerSuccess = 0;
        do {
            //randomly pick an advancement left on the list.
            int advancementListRandIndex = rand() % (advancementList.n_rows);
            int randomAdvancement = advancementList(advancementListRandIndex);
            //path variable actually not used; for debugging and possible use elsewhere.
            imat path;
            //initialize with starting dimer

            path << chosen << neighbor << endr;

            //   cout << "FLAGmid " << endl;

            beforeEnd = chosen;
            afterEnd = neighbor;

            // starting at one end of a ring, travels to the opposite end (for the two adatoms).
            for (int i = 0; i < randomAdvancement; ++i) {
                //loops over neighbors
                for (int j = 0; j < neighList.n_cols; ++j) {
                    // condition for advancing is the neighbor is not the previous atom,
                    //that the neighbor is on the ring, and ,
                    //and the neighbor exists
                    if ((neighList(afterEnd, j) != beforeEnd) &&
                            (contains(winningRing, neighList(afterEnd, j))) &&
                            (neighList(afterEnd, j) != -1)) {
                        beforeEnd = afterEnd;
                        afterEnd = neighList(afterEnd, j);
                        imat toBeAdded;
                        toBeAdded << afterEnd << endr;

                        path = join_rows(path, toBeAdded);
                        break;
                    }
                }
            }
            //destroy junk value
            path.shed_col(0);

            // cout << "FLAGno " << endl;

            //We must check if the beforeEnd and afterEnd are NOT in a double bond. if condtion fails, call the mutation again.
            if (isNotDoubleBond(neighList, chosen, neighbor, NN1, NN2, NNN1, NNN2, beforeEnd, afterEnd) == 1) {
                addDimerSuccess = 1;
            } else {
                //shed failed end,
                advancementList.shed_row(advancementListRandIndex);
                //if advancementList is empty, shed the picked bond and try again.
                if (advancementList.is_empty() == true) {
                    remainingBonds.shed_row(currentRemainingBond);

                    if (remainingBonds.is_empty()) {
                        cout << "NO MUTATION POSSIBLE, FAILED AT TRYING TO FIND SINGLE BOND FOR ADD DIMER " << endl;
                        return;
                    }
                    mutation(mutationNum, neighList, bonds, positions, remainingBonds, masterBatchName, inBeam);
                    return;
                }
            }
            // loop continues as long as the advancement list is NOT empty
        } while (addDimerSuccess == 0);



        //    cout << "FLAG " << endl;


        // note, originally this coded created a new atom, with a position that is the average of the first chosen and neighbor atoms.
        //However, NAMD cannot handle atoms that are too close to one another;
        //thus this new code generates a circle of a given radius whose plane is perpendicular to the chosen-neighbor vector through the midpoint,
        //and creates a new atom that is on a random angle on this circle.
        //repeats with the end atoms.



        double addDimerCircleRadius = 0.4;
        double randomAngle = (rand() % 359) * 3.14159265 / 180;
        vec xVec;
        xVec << 1 << 0 << 0 << endr;
        rowvec yVec;
        yVec << 0 << 1 << 0 << endr;
        vec u;
        vec v;
        vec circleCenter;
        vec circleNormal;
        double parallelCrossingTolerance = 0.1;
        //Formula for a plane is A(x - x0) + B(y - y0) + C(z - z0) = 0,
        //where A,B,and C are the slopes
        //(derived from the chosen and neighbor coordinates), and the x0, y0, and z0 values are
        //the midpoint positions.

        //Parametric formula for a circle is R(t)=Center + r*cos(t)u + r*sin(t)v,
        //where Center is the coordinates of the center, r is radius, and u and v are orthogonal unit vectors in the circle.

        circleCenter << (positions(chosen, 0) + positions(neighbor, 0)) / 2 << (positions(chosen, 1) + positions(neighbor, 1)) / 2 << (positions(chosen, 2) + positions(neighbor, 2)) / 2 << endr;

        circleNormal << positions(chosen, 0) - positions(neighbor, 0) << positions(chosen, 1) - positions(neighbor, 1) << positions(chosen, 2) - positions(neighbor, 2) << endr;

        // we must now find two orthogonal directions to the normal vector.
        //we can do this by crossing the normal with an arbitrary vector, and crossing this result with the norm again, normalizing as we go along.
        //We will arbitrarily choose [1 0 0], however if this is coincidentally parallel to the normal, then we will choose [0 1 0]

        u = cross(circleNormal, xVec);
        mat absu = abs(u);
        vec vecConvertedFromMat = absu;
        if (sum(vecConvertedFromMat) < parallelCrossingTolerance) {
            //the vectors are too close to parallel, choose a new arbitrary.
            u = cross(circleNormal, yVec);
            //    cout << "TOO CLOSE THE FIRST TIME" << endl;
        }
        //normalize
        u = u / norm(u, "inf");

        //second vector.
        v = cross(circleNormal, u);
        v = v / norm(v, "inf");

        mat newAtom1Positions;
        //now use the paramentric formula discussed from before.
        newAtom1Positions << circleCenter(0) + addDimerCircleRadius * (cos(randomAngle) * u(0) + sin(randomAngle) * v(0))
                << circleCenter(1) + addDimerCircleRadius * (cos(randomAngle) * u(1) + sin(randomAngle) * v(1))
                << circleCenter(2) + addDimerCircleRadius * (cos(randomAngle) * u(2) + sin(randomAngle) * v(2))
                << endr;



        //repeat for second atom to be added
        circleCenter.clear();
        circleCenter << (positions(beforeEnd, 0) + positions(afterEnd, 0)) / 2 << (positions(beforeEnd, 1) + positions(afterEnd, 1)) / 2 << (positions(beforeEnd, 2) + positions(afterEnd, 2)) / 2 << endr;

        circleNormal.clear();
        circleNormal << positions(beforeEnd, 0) - positions(afterEnd, 0) << positions(beforeEnd, 1) - positions(afterEnd, 1) << positions(beforeEnd, 2) - positions(afterEnd, 2) << endr;

        // we must now find two orthogonal directions to the normal vector.
        //we can do this by crossing the normal with an arbitrary vector, and crossing this result with the norm again, normalizing as we go along.
        //We will arbitrarily choose [1 0 0], however if this is coincidentally parallel to the normal, then we will choose [0 1 0]
        u.clear();

        u = cross(circleNormal, xVec);
        absu = abs(u);
        vecConvertedFromMat = absu;
        if (sum(vecConvertedFromMat) < parallelCrossingTolerance) {
            //the vectors are too close to parallel, choose a new arbitrary.
            u = cross(circleNormal, yVec);

        }
        //normalize
        u = u / norm(u, "inf");

        //second vector.
        v.clear();
        v = cross(circleNormal, u);
        v = v / norm(v, "inf");

        mat newAtom2Positions;
        //now use the paramentric formula discussed from before.
        newAtom2Positions << circleCenter(0) + addDimerCircleRadius * (cos(randomAngle) * u(0) + sin(randomAngle) * v(0))
                << circleCenter(1) + addDimerCircleRadius * (cos(randomAngle) * u(1) + sin(randomAngle) * v(1))
                << circleCenter(2) + addDimerCircleRadius * (cos(randomAngle) * u(2) + sin(randomAngle) * v(2))
                << endr;




        //handy code
        /*
        writeTopFile(positions, bonds, batchName);
        writePgnFile(positions, batchName);
        writePdbPsfFiles(batchName);

        //Handy debugging code colors the ring in question; this pdb file is kept seperate, in batchnameModified.pdb
        system(("rm " + batchName + "Modified.pdb").c_str());
        string line;
        ifstream original((batchName + ".pdb").c_str());
        ofstream modified((batchName + "Modified.pdb").c_str());
        int lineNumber = 0;
        while (getline(original, line)) {

            //does not modify the first or last line.
            //looks for only atoms in the winning ring.
            if ((lineNumber > 0) && (line.size() > 4) && contains(winningRing, lineNumber - 1)) {
                line.erase(77, 79);
                line = line + "N";

            }
            modified << line << "\n";
            ++lineNumber;
        }
        modified.close();
        original.close();
         */


        //adds atoms
        positions = join_cols(positions, newAtom1Positions);
        positions = join_cols(positions, newAtom2Positions);


        //updates neighbor list by finding references to chosen atom pairs, and replacing them with the indexes
        //of the inserted atoms (which will be appended to the end of the neighbor list)

        for (int j = 0; j < neighList.n_cols; ++j) {
            if (neighList(chosen, j) == neighbor) {
                neighList(chosen, j) = (neighList.n_rows);
                break;
            }
        }
        for (int j = 0; j < neighList.n_cols; ++j) {
            if (neighList(neighbor, j) == chosen) {
                neighList(neighbor, j) = (neighList.n_rows);
                break;
            }
        }
        for (int j = 0; j < neighList.n_cols; ++j) {
            if (neighList(beforeEnd, j) == afterEnd) {
                neighList(beforeEnd, j) = (neighList.n_rows + 1);
                break;
            }
        }
        for (int j = 0; j < neighList.n_cols; ++j) {
            if (neighList(afterEnd, j) == beforeEnd) {
                neighList(afterEnd, j) = (neighList.n_rows + 1);
                break;
            }
        }
        //new atom1 has new atom2 as a neighbor, and vice versa
        imat newAtom1Neighbors;
        newAtom1Neighbors << chosen << neighbor << neighList.n_rows + 1 << -1 << endr;
        imat newAtom2Neighbors;
        newAtom2Neighbors << beforeEnd << afterEnd << neighList.n_rows << -1 << endr;


        /*
                for (int i = 0; i < bonds.n_rows; i++) {
                    for (int j = i; j < bonds.n_rows; j++) {
                        if (j != i) {
                            if ((bonds(i, 0) == bonds(j, 0) && bonds(i, 1) == bonds(j, 1))
                                    || (bonds(i, 0) == bonds(j, 1) && bonds(i, 1) == bonds(j, 0))) {
                                cout << "OFFENDING BONDS: THE FIRST " << bonds.row(i) << "THE SECOND " << bonds.row(j) << endl;
                                cout << "OFFENDING ROW NUM: " << i << " " << j << endl;
                                system("read  -p \"paused, same bond counted twice BEFORE ADDITION.\" key");
                            }
                        }
                    }
                }
         */



        neighList = join_cols(neighList, newAtom1Neighbors);
        neighList = join_cols(neighList, newAtom2Neighbors);

        bonds = buildBond(neighList);



    }
}

//the main point of this subroutine is find a dimer and its neighbors around which feed to
//mutationBody.
//Either rotates a bond, creates a divacancy, or adds a dimer, depending on
//the mutation number given (1,2,and 3 respectively
//The mutations all have the same set-up: choose an atom, find it's Nearist Neighbors NN, and it's
//Next Nearest Neighbors NNN. Modify the rings, delete the outdated versions, and then add the new ones.
//remainingBonds starts off as a full bond list. As the bonds fail various conditions and this subfunction is recursively called again, this list is chopped down.
//This ensures the code will not run forever if a specific mutation is impossible, and will improve run times.

void DataStructureFunctionsATAC::mutation(int mutationNum, imat& neighList, imat& bonds, mat& positions, imat& remainingBonds, string masterBatchName, int& inBeam) {

    int chosen;
    int neighbor;
    int NN1;
    int NN2;
    int NNN1;
    int NNN2;
    int randomBond;
    uvec NNIndexes;
    uvec NNNIndexes;
    uvec chosenUnbonded;
    uvec neighborUnbonded;
    bool dimerFound = false;
    do {

        //check that we still have bonds to Test...
        if (remainingBonds.is_empty()) {
            //    cout << "BEAM RADIUS " << beamRadius << endl;
            cout << "NO MUTATION POSSIBLE, FAILED TRYING TO FIND APPROPRIATE DIMER " << endl;
            return;
        }
        randomBond = rand() % remainingBonds.n_rows;
        //    cout << "RANDOM BOND" << endl;
        chosen = remainingBonds(randomBond, 0);
        //  cout << " CHO " << chosen << endl;
        neighbor = remainingBonds(randomBond, 1);
        //    cout << "NEIGH " << neighbor << endl;
        NNIndexes = find(neighList.row(chosen) != neighbor);
        //            cout <<"NNIndexes " << NNIndexes << endl;
        NN1 = neighList(chosen, NNIndexes[0]);
        //        cout <<"NN1 " << NN1 << endl;
        NN2 = neighList(chosen, NNIndexes[1]);
        //      cout <<"NN2 " << NN2<< endl;
        NNNIndexes = find(neighList.row(neighbor) != chosen);
        //                    cout <<"NNNIndexes " << NNNIndexes << endl;
        NNN1 = neighList(neighbor, NNNIndexes[0]); //
        //       cout <<"NNN1 " << NNN1 << endl;
        NNN2 = neighList(neighbor, NNNIndexes[1]);
        //         cout <<"NNN2 " << NNN2 << endl;
        chosenUnbonded = find(neighList.row(chosen) == -1);
        neighborUnbonded = find(neighList.row(neighbor) == -1);

        bool isElectrode = false;
        //if this system has electrodes, no atoms may participate in the mutation.
        if (numberOfElectrodes != 0) {

            //NOTE: TO UPDATE ELECTRODES AT EVERY INSTANCE, UNCOMMENT FOLLOWING LINE. DO NOT DELETE-ISH.
            //electrodeDeterminer(positions, xp1Electrodes, xn1Electrodes, yp1Electrodes, yn1Electrodes, xp2Electrodes, xn2Electrodes, yp2Electrodes, yn2Electrodes);

            //since ordering is not important, we can just stack the electrodes into one matrix.
            imat electrodes;
            if (numberOfElectrodes == 2) {
                imat xpElectrodes = join_cols(xp1Electrodes, xp2Electrodes);
                imat xnElectrodes = join_cols(xn1Electrodes, xn2Electrodes);
                electrodes = join_cols(xpElectrodes, xnElectrodes);
            }
            if (numberOfElectrodes == 4) {
                imat xpElectrodes = join_cols(xp1Electrodes, xp2Electrodes);
                imat xnElectrodes = join_cols(xn1Electrodes, xn2Electrodes);
                imat intermediate1;
                intermediate1 = join_cols(xpElectrodes, xnElectrodes);
                imat ypElectrodes = join_cols(yp1Electrodes, yp2Electrodes);
                imat ynElectrodes = join_cols(yn1Electrodes, yn2Electrodes);
                imat intermediate2;
                intermediate2 = join_cols(ypElectrodes, ynElectrodes);

                electrodes = join_cols(intermediate1, intermediate2);
            }

            //loop through our electrode atoms.
            for (int i = 0; i < electrodes.n_rows; i++) {

                //if the electrode atom is any of the atoms we are currently looking at,
                //disreguard this dimer.
                if ((electrodes(i) == chosen) || (electrodes(i) == neighbor)
                        || (electrodes(i) == NN1) || (electrodes(i) == NN2)
                        || (electrodes(i) == NNN1) || (electrodes(i) == NNN2)) {
                    isElectrode = true;
                }
            }

        }

        //test makes sure the atoms we chose are within the beam (if the beam is exclusionary like that).
        //assume the system passes the bool, until proven wrong.
        bool passBeamTest = true;

        //We check if the system has a irradation beam.
        if (haveBeam == 1) {
            //find axis intercept dimentions (0, 1, 2, for x, y, z).
            int intercept1, intercept2;
            if (beamAxis == 1) {
                intercept1 = 1;
                intercept2 = 2;
            } else if (beamAxis == 2) {
                intercept1 = 0;
                intercept2 = 2;
            } else if (beamAxis == 3) {
                intercept1 = 0;
                intercept2 = 1;
            }

            //check beam type.
            if (beamType == 1) {
                //This is a stationary beam that does not allow mutations outside the beam.
                //the only way we can fail this test is if either neighbor or chosen are outside the beam.

                //check that both chosen and neighbor atoms are in the beam.
                if (pow(positions(chosen, intercept1) - beamAxisIntercept1, 2) + pow(positions(chosen, intercept2) - beamAxisIntercept2, 2) <= pow(beamRadius, 2) &&
                        pow(positions(neighbor, intercept1) - beamAxisIntercept1, 2) + pow(positions(neighbor, intercept2) - beamAxisIntercept2, 2) <= pow(beamRadius, 2)) {
                    //passed
                    inBeam == 1;
                } else {
                    //failed
                    passBeamTest = false;
                    inBeam == 0;
                }

            }

            if (beamType == 2) {
                //This is a stationary beam that DOES allow mutations outside the beam.

                //always passes test technically, as mutaitons are allowed everywhere.
                //check that both chosen and neighbor atoms are in the beam.
                if (pow(positions(chosen, intercept1) - beamAxisIntercept1, 2) + pow(positions(chosen, intercept2) - beamAxisIntercept2, 2) <= pow(beamRadius, 2) &&
                        pow(positions(neighbor, intercept1) - beamAxisIntercept1, 2) + pow(positions(neighbor, intercept2) - beamAxisIntercept2, 2) <= pow(beamRadius, 2)) {
                    //passed
                    inBeam == 1;
                } else {
                    //failed
                    inBeam == 0;
                }
            }
        }



        // if neighbor or the NNs or the NNNs don't exist,
        // or  if the NN and NNN have the same values (a three-ring triangle situation)
        // or if the chosen or neighbor are sp3,
        // or if the NN's or NNN's are neighbors of each other (triangle at end of dimer situation),
        // or if any of the atoms are electrodes,
        // delete failed bond and continue loop.

        if ((neighbor < 0) || (NN1 < 0) || (NN2 < 0) || (NNN1 < 0) || (NNN2 < 0)
                || (NN1 == NN2) || (NN1 == NNN1) || (NN1 == NNN2) || (NN2 == NNN1) || (NN2 == NNN2) || (NNN1 == NNN2)
                || (chosenUnbonded.is_empty() == true) || (neighborUnbonded.is_empty() == true)
                || contains(neighList.row(NN1), NN2) || contains(neighList.row(NN2), NN1)
                || contains(neighList.row(NNN1), NNN2) || contains(neighList.row(NNN2), NNN1)
                || (isElectrode == true)
                || (passBeamTest == false)) {
            //          cout << "BOND FAILED " << endl;
            remainingBonds.shed_row(randomBond);
        } else {

            //cout << "chosen neighbor NN1 NN2 NNN1 NNN2 " << chosen <<" " << neighbor << " " << NN1 << " " << NN2 << " " << " " << NNN1 << " " << NNN2 << endl;
            dimerFound = true;
        }
    } while (dimerFound == false);
    //cout << " TRUE BONDS AFTER" << bonds << endl;
    mutationBody(mutationNum, neighList, bonds, positions, chosen, neighbor, NN1, NN2, NNN1, NNN2, remainingBonds, randomBond, masterBatchName, inBeam);
}

//INCLUDE DOUBLE BOND CHECK
//checks if atoms i and j can go through an sp2 to sp3 condition to bond with
//each other.

bool DataStructureFunctionsATAC::sp2Tosp3Conditions(int i, int j, mat positions, imat& neighbors) {
    bool pass = false;

    //check that i and j exist...
    if ((i < 0) || (j < 0)) {
        return pass;
    }

    uvec jNeighborOfiIndex = find(neighbors.row(i) == j);
    uvec emptySpotsi = find(neighbors.row(i) == -1);
    uvec emptySpotsj = find(neighbors.row(j) == -1);

    // fails if i and j are the same atom, or
    // are neighbors
    // or are not sp2 (i.e. bonded to three atoms, so one possible neighbor is missing...)
    if ((j == i)
            || (jNeighborOfiIndex.is_empty() == false)
            || (emptySpotsi.n_elem != 1)
            || (emptySpotsj.n_elem != 1)) {
        return pass;
    }

    double distance;

    /*
    int periodicityNumBackup;
    if (periodicOnlyOnNonInteracting == true) {
        //IMPORTANT: CURRENTLY THE DYNAMIC PERIODIC CONDITION WILL NOT LOOK AT THE NORMAL ATOMS, ONLY THE NON-INTERACTING ONES.
        periodicityNumBackup = periodicityNum;
        periodicityNum = 0;
    }*/

    shortestLineDeterminer(positions, i, j, distance, periodicityNum, xLatticeParameter, yLatticeParameter, zLatticeParameter, latticeVec1Norm, latticeVec1, latticeVec2Norm, latticeVec2, latticeVec3Norm, latticeVec3);

    /*
    if (periodicOnlyOnNonInteracting == true) {
        periodicityNum = periodicityNumBackup;
    }*/

    //atoms must be within a certain distance of each other.

    if (distance > sp2tosp3Tolerance) {
        //cout << "ATOMS TOO FAR FOR SP2 TO SP3" << endl;
        return false;
    }

    //"spDenialPathSize"
    //the next set of lines will not allow the sp2 to sp3 transition to occur
    //if a path of a given length can be found between the dimers in question.
    //essentially it is to prevent the sp2 to sp3 transition to occur within
    //the same structure, for instance, imploding a bucky ball. This is more
    //likely to occur for large sp2tosp3Tolerance values.

    bool victory = false;
    imat path;
    path << -1 << endr;
    imat winner;
    int layersDeep = 0 - 1;
    //ring size must be four, as we are looking for a 4 ring.
    int pathSizeMax = spDenialPathSize + 1; //a value of 4+1 will effectively not allow for sp2 to sp3 to be undone.

    findPathForAtoms(i, j, path, winner, neighbors, layersDeep, victory, pathSizeMax);

    // if a path can be found between the two atoms, then they are part of the immediate
    // surrounding structure, not seperate ones.
    if (victory) {
        return false;
    }

    //congrats, if the code as made it this far, it has passed all tests.
    pass = true;
    /*
    //now must see if the atoms are a part of two six rings, each!
    //CHANGE TO DOUBLE BOND FINDING ALGORITM!!!
    for (int k = 0; k < rings.size(); ++k) {
        if ((rings[k].second == 6) && (rings[k].first.n_rows != 0)) {
            uvec iInSixIndexes = find(i == rings[k].first);
            uvec jInSixIndexes = find(j == rings[k].first);
            if ((iInSixIndexes.n_elem >= 2) && (jInSixIndexes.n_elem >= 2)) {
                //congrats, all conditions have been passed for the chosen pair.
                pass = true;
                return pass;
            }else{
                return pass;
            }
        }else{return pass;}
    }
    return pass;
}else{return pass;}
    }*/
    return pass;
}

// the sp3 to sp2 transition occurs via two even-numbered rings (including size 2) that are connected face to face.
// same could be said for 2 to 3, but in reverse.
// Only rings of size 2 are considered, as the 3 to 2 barrier is significantly lower in such a case.

void DataStructureFunctionsATAC::sp3Tosp2(imat& neighbors, imat& bonds, mat& positions) {

    //SP3 to SP2: break the bonds between chosen and chosen's neighbor, other and other's neighbor

    //simply creates a list numbered 0 through # of atoms
    imat atomsLeftToBeTested(neighbors.n_rows, 1);
    for (
            int i = 0; i < neighbors.n_rows; ++i) {
        atomsLeftToBeTested(i) = i;
    }

    //must find a four ring in sp3...
    //make a neighlist of just sp3
    imat sp3Neighbors = neighbors;

    ///first we loop over neighbors, and eliminate all sp2 carbon atoms.
    for (int i = 0; i < sp3Neighbors.n_rows; ++i) {
        uvec negIndex = find(sp3Neighbors.row(i) == -1);
        if (negIndex.is_empty() != true) {
            //at least one negative index, thus less than 4 neighbor atoms, thus sp2, thus must be eliminated.
            sp3Neighbors(i, 0) = -1;
            sp3Neighbors(i, 1) = -1;
            sp3Neighbors(i, 2) = -1;
            sp3Neighbors(i, 3) = -1;

            atomsLeftToBeTested(i) = -1;
        }
    }

    int iterLimit = atomsLeftToBeTested.n_rows;
    for (int i = 0; i < iterLimit; ++i) {

        if (atomsLeftToBeTested(i) == -1) {
            atomsLeftToBeTested.shed_row(i);
            --i;
            --iterLimit;
        }
    }

    if (atomsLeftToBeTested.is_empty()) {
        cout << " NOT ENOUGH SP3 TO CONVERT TO SP2 " << endl;
        return;
    }

    //inefficient, change later.
    imat sp3Bonds = buildBond(sp3Neighbors);
    buildNeighFromBonds(sp3Bonds, sp3Neighbors);
    imat winningRing;
    int chosen;
    bool victory = false;

    do {
        //cout << "ATOMS LEFT" << atomsLeftToBeTested.n_rows << endl;

        int randomIndex = rand() % atomsLeftToBeTested.n_rows;
        chosen = atomsLeftToBeTested(randomIndex);

        //ring size must be 4, in this transition.
        int ringSizeMax = 4;
        imat ring;

        for (int i = 0; i < neighbors.n_cols; ++i) {
            ring << chosen << neighbors(chosen, i) << endr;

            // must start depth at 0, as 1 is added at the beginning.
            findSmallestRingForAtom(neighbors(chosen, i), ring, winningRing, neighbors, 0, victory, ringSizeMax);
            if (victory == true) {
                break;
            }
        }

        if (victory == false) {
            atomsLeftToBeTested.shed_row(randomIndex);
            if (atomsLeftToBeTested.is_empty()) {
                cout << " COULDN'T FIND 4 RING " << endl;
                return;
            }
        }
    } while (victory == false);

    //only working with  4 rings currently. So, these atoms were chosen, chosenNeighbor,
    //other, and otherNeighbor, from the sp2 to sp3 transition.
    //We must randomize which atom in the ring is chosen, etc.
    int chosenNeighbor;
    int otherNeighbor;
    int other;

    int randNum = rand() % 4;


    if (superfluousSp3Sp2 == false) {
        //Disallows superfluous sp3 to sp2 transitions, i.e. ones that reverse the effects of the previous sp2 to sp3 reaction.
        randNum = 2;
    }

    chosen = winningRing(randNum);
    randNum++;
    if (randNum >= 4) {
        randNum = randNum - 4;
    }
    chosenNeighbor = winningRing(randNum);
    randNum++;
    if (randNum >= 4) {
        randNum = randNum - 4;
    }
    otherNeighbor = winningRing(randNum);
    randNum++;
    if (randNum >= 4) {
        randNum = randNum - 4;
    }
    other = winningRing(randNum);

    //loops over rows in bonds to find indexes of those to be deleted.
    int iterMax = bonds.n_rows;
    for (int i = 0; i < iterMax; ++i) {

        if ((bonds(i, 0) == chosen) && (bonds(i, 1) == chosenNeighbor)) {
            bonds.shed_row(i);
            --i;
            --iterMax;
        }
        if ((bonds(i, 0) == chosenNeighbor) && (bonds(i, 1) == chosen)) {
            bonds.shed_row(i);
            --i;
            --iterMax;
        }
        if ((bonds(i, 0) == other) && (bonds(i, 1) == otherNeighbor)) {
            bonds.shed_row(i);
            --i;
            --iterMax;
        }
        if ((bonds(i, 0) == otherNeighbor) && (bonds(i, 1) == other)) {
            bonds.shed_row(i);
            --i;
            --iterMax;
        }
    }

    //do the same for neighbors
    uvec chosenNeighborIndex1 = find(neighbors.row(chosen) == chosenNeighbor);
    neighbors(chosen, chosenNeighborIndex1(0)) = -1;
    uvec chosenIndex1 = find(neighbors.row(chosenNeighbor) == chosen);
    neighbors(chosenNeighbor, chosenIndex1(0)) = -1;
    uvec otherNeighborIndex1 = find(neighbors.row(other) == otherNeighbor);
    neighbors(other, otherNeighborIndex1(0)) = -1;
    uvec otherIndex1 = find(neighbors.row(otherNeighbor) == other);
    neighbors(otherNeighbor, otherIndex1(0)) = -1;

    //we now rotate the neighbor other and neighbor neighbor bonds.
    //but first, must find the NN and NNN of each respectively.
    int chosenNN1;
    int chosenNN2;
    uvec chosenNNIndexes = find(neighbors.row(chosen) != other);
    uvec chosenNNPosIndexes = find(neighbors.row(chosen) != -1);
    int countNN = 0;
    for (int i = 0; i < chosenNNIndexes.n_elem; ++i) {
        for (int j = 0; j < chosenNNPosIndexes.n_elem; ++j) {
            if (chosenNNIndexes(i) == chosenNNPosIndexes(j) && (countNN == 0)) {
                chosenNN1 = neighbors(chosen, chosenNNIndexes(i));
                ++countNN;
                continue;
            }
            if (chosenNNIndexes(i) == chosenNNPosIndexes(j) && (countNN == 1)) {
                chosenNN2 = neighbors(chosen, chosenNNIndexes(i));
                ++countNN;
                continue;
            }
        }
    }

    int otherNN1;
    int otherNN2;
    uvec otherNNIndexes = find(neighbors.row(other) != chosen);
    uvec otherNNPosIndexes = find(neighbors.row(other) != -1);
    int countNNN = 0;
    for (int i = 0; i < otherNNIndexes.n_elem; ++i) {
        for (int j = 0; j < otherNNPosIndexes.n_elem; ++j) {
            if (otherNNIndexes(i) == otherNNPosIndexes(j) && (countNNN == 0)) {
                otherNN1 = neighbors(other, otherNNIndexes(i));
                ++countNNN;
                continue;
            }
            if (otherNNIndexes(i) == otherNNPosIndexes(j) && (countNNN == 1)) {
                otherNN2 = neighbors(other, otherNNIndexes(i));
                ++countNNN;
                continue;
            }
        }
    }

    //these are just dummy values, and should never be called in mutation body...
    imat remainingBonds;
    int currentBond;

    int chosenNeighborNN1;
    int chosenNeighborNN2;
    uvec chosenNeighborNNIndexes = find(neighbors.row(chosenNeighbor) != otherNeighbor);
    uvec chosenNeighborNNPosIndexes = find(neighbors.row(chosenNeighbor) != -1);
    int countNeighborNN = 0;
    for (int i = 0; i < chosenNeighborNNIndexes.n_elem; ++i) {
        for (int j = 0; j < chosenNeighborNNPosIndexes.n_elem; ++j) {
            if (chosenNeighborNNIndexes(i) == chosenNeighborNNPosIndexes(j) && (countNeighborNN == 0)) {
                chosenNeighborNN1 = neighbors(chosenNeighbor, chosenNeighborNNIndexes(i));
                ++countNeighborNN;
                continue;
            }
            if (chosenNeighborNNIndexes(i) == chosenNeighborNNPosIndexes(j) && (countNeighborNN == 1)) {
                chosenNeighborNN2 = neighbors(chosenNeighbor, chosenNeighborNNIndexes(i));
                ++countNN;
                continue;
            }
        }
    }

    int otherNeighborNN1;
    int otherNeighborNN2;
    uvec otherNeighborNNIndexes = find(neighbors.row(otherNeighbor) != chosenNeighbor);
    uvec otherNeighborNNPosIndexes = find(neighbors.row(otherNeighbor) != -1);
    int countNeighborNNN = 0;
    for (int i = 0; i < otherNeighborNNIndexes.n_elem; ++i) {
        for (int j = 0; j < otherNeighborNNPosIndexes.n_elem; ++j) {
            if (otherNeighborNNIndexes(i) == otherNeighborNNPosIndexes(j) && (countNeighborNNN == 0)) {
                otherNeighborNN1 = neighbors(otherNeighbor, otherNeighborNNIndexes(i));
                ++countNeighborNNN;
                continue;
            }
            if (otherNeighborNNIndexes(i) == otherNeighborNNPosIndexes(j) && (countNeighborNNN == 1)) {
                otherNeighborNN2 = neighbors(otherNeighbor, otherNeighborNNIndexes(i));
                ++countNeighborNNN;
                continue;
            }
        }
    }
}

//This algorithm converts sp2 bonds to sp3 bonds, via the mechanism of
//forcing two double bonded pairs of atoms within close vicinity of each other.

//First, looks for two atoms that are within 2 angstroms, and not otherwise bonded.
//if these atoms are a part of mutually exclusive six rings,
//then look for neighbors that also fulfill the previous conditions.

//if so, bond atoms to atoms and neighbors to neighbors

//Might later modify slightly so it can bond edge atoms as well.

void DataStructureFunctionsATAC::sp2Tosp3(mat positions, imat& neighbors, imat& bonds, int& spPossibility) {

    spPossibility = 0; //reset the spPossibility number, as we have performed the operations needed.

    bool chosenPairMeetConditions = false;
    bool neighborPairMeetConditions = false;
    int chosen;
    int other;
    int randomIndex1;
    int randomIndex2;
    //simply creates a list numbered 0 through # of atoms
    ivec atomsLeftToBeTested(positions.n_rows);
    for (
            int i = 0; i < positions.n_rows; ++i) {
        atomsLeftToBeTested(i) = i;
    }
    //we will randomly choose atoms to find a sp2 to sp3 candidate,
    //until all possibilities have been exhausted.

    // we keep testing until we succeed or have nothing left to test
    do {
        randomIndex1 = rand() % atomsLeftToBeTested.n_rows;
        chosen = atomsLeftToBeTested(randomIndex1);

        //simply creates a list numbered 0 through # of atoms
        ivec atomsToCheckAgainst(positions.n_rows);
        for (int i = 0; i < positions.n_rows; ++i) {
            atomsToCheckAgainst(i) = i;
        }

        do {
            randomIndex2 = rand() % atomsToCheckAgainst.n_rows;
            //other is the atom we are comparing chosen to
            other = atomsToCheckAgainst(randomIndex2);
            // and we check the conditions...
            // DataStructureFunctionsATAC molecule; //DANGER DANGER DANGER
            chosenPairMeetConditions = sp2Tosp3Conditions(chosen, other, positions, neighbors);
            //If the chosen pair meets the conditions, then sees if there are neighbors that do the same
            if (chosenPairMeetConditions == true) {
                //chosenNeighbor is the chosen neighbor atom we are currently looking at
                for (int chosenNeighborIter = 0; chosenNeighborIter < neighbors.n_cols; ++chosenNeighborIter) {
                    //otherNeighbor is the other neighbor atom we are comparing chosenNeighbor to
                    for (int otherNeighborIter = 0; otherNeighborIter < neighbors.n_cols; ++otherNeighborIter) {
                        int chosenNeighbor = neighbors(chosen, chosenNeighborIter);
                        int otherNeighbor = neighbors(other, otherNeighborIter);
                        // must also check that the neighbors we are picking out are not one of the already chosen atoms
                        if ((chosen != otherNeighbor) && (other != chosenNeighbor)) {

                            //we must find neighbors of our chosen pair, that also fufil the condition.
                            neighborPairMeetConditions = sp2Tosp3Conditions(chosenNeighbor, otherNeighbor, positions, neighbors);
                            //now see if we can execute the mutation... all sp2 conditions must be fufilled, and the double bonds confirmed.
                            if ((chosenPairMeetConditions == true) && (neighborPairMeetConditions == true) && isDoubleBond(neighbors, chosen, chosenNeighbor) && isDoubleBond(neighbors, other, otherNeighbor)) {
                                //success
                                //set the chosen pair to each other...
                                uvec emptySpoti = find(-1 == neighbors.row(chosen));
                                uvec emptySpotj = find(-1 == neighbors.row(other));
                                neighbors(chosen, emptySpoti(0)) = other;
                                neighbors(other, emptySpotj(0)) = chosen;
                                //and the neighbor pair
                                uvec emptySpotl = find(-1 == neighbors.row(chosenNeighbor));
                                uvec emptySpotm = find(-1 == neighbors.row(otherNeighbor));
                                neighbors(chosenNeighbor, emptySpotl(0)) = otherNeighbor;
                                neighbors(otherNeighbor, emptySpotm(0)) = chosenNeighbor;
                                //add to bonds
                                imat toBeAppended;
                                toBeAppended
                                        << chosen << other << endr
                                        << chosenNeighbor << otherNeighbor << endr;

                                bonds = join_cols(bonds, toBeAppended);
                                spPossibility = 1;
                                cout << "SP2 TO SP3 TRANSITION OK" << endl;
                                return;
                            }
                        }
                    }
                }
            }
            //if the code ever reaches this far, then it failed to fufil all the conditions
            // chosen and other failed, get new other atom
            atomsToCheckAgainst.shed_row(randomIndex2);
        } while (atomsToCheckAgainst.n_elem > 0);
        // chosen failed, get new chosen atom
        atomsLeftToBeTested.shed_row(randomIndex1);

    } while (atomsLeftToBeTested.n_elem > 0);
    cout << "SP2 TO SP3 TRANSITION CURRENTLY IMPOSSIBLE" << endl;
    spPossibility = -1;
}

//adds non-interacting atoms to our position and bond lists.

void DataStructureFunctionsATAC::nonInteractingAppender(mat& positions, imat& bondList, imat& neighborList) {




    // for appending the bonds and the neighbors, we must keep in mind that the non-interacting
    // indexes have to be shifted over
    // by the number of normal interacing atoms. Also, -1 must stay -1.
    imat bondsToBeAppended = nonInteractingBonds + positions.n_rows;
    bondList = join_cols(bondList, bondsToBeAppended);

    // Also, -1 must stay -1. We can find what values used to be -1, buy seeing which values equal added amount (positions.n_rows)-1.
    imat neighborsToBeAppended = nonInteractingNeighbors + positions.n_rows;
    for (int i = 0; i < neighborsToBeAppended.n_rows; ++i) {
        for (int j = 0; j < neighborsToBeAppended.n_cols; ++j) {
            //if the neighbor to be appended used to not exist (be -1)....
            if (neighborsToBeAppended(i, j) == positions.n_rows - 1) {
                neighborsToBeAppended(i, j) = -1;
            }
        }
    }

    neighborList = join_cols(neighborList, neighborsToBeAppended);

    //for positions we can simply append the non interacting on the back end of the normal atoms
    positions = join_cols(positions, nonInteracting);


}

//removes the non-interacting atoms from our position and bond lists,
//while updating the non-interacting component
//the non interacting bond list never changes.

void DataStructureFunctionsATAC::nonInteractingDeAppender(mat& positions, imat& bondList, imat& neighList) {



    //non interacting particles are appended at the end.
    nonInteracting = positions.rows(positions.n_rows - numOfNonInteracting, positions.n_rows - 1);
    positions = positions.rows(0, positions.n_rows - 1 - numOfNonInteracting);
    bondList = bondList.rows(0, bondList.n_rows - 1 - numOfNonInteractingBonds);
    neighList = neighList.rows(0, neighList.n_rows - 1 - numOfNonInteractingNeighbors);


}

//finds chirality based on diameter. Assumes nanotube is x aligned.

void DataStructureFunctionsATAC::findChiralityBasedOnDiameter(mat positions, imat bonds, imat neighList, bool haveNonInteracting) {

    if (haveNonInteracting == true) {
        nonInteractingDeAppender(positions, bonds, neighList);
    }

    //assume diamter is an average of the x and y max and min differences.
    vec zCol = positions.col(2);
    vec yCol = positions.col(1);
    double zDiff = abs(zCol.max() - zCol.min());
    double yDiff = abs(yCol.max() - yCol.min());
    double diameter = (yDiff + zDiff) / 2;

    //now we must generate a set of known diameter values to compare against.





    double a = 2.46;
    //double pi = 3.14159265;
    mat diameterDifferences;
    diameterDifferences << -1 << -1 << -1 << endr;

    for (int n = minChiralIndex; n <= maxChiralIndex; ++n) {
        for (int m = minChiralIndex - 1; m <= n; m++) {
            mat toBeAppended;
            toBeAppended << abs(diameter - a / 3.14159265 * sqrt(n * n + n * m + m * m)) << n << m << endr;
            diameterDifferences = join_cols(diameterDifferences, toBeAppended);
            //       cout << "N, M, and Diameter, sqrt " <<n<< " " << m<<" " <<  sqrt(n*n + n * m + m*m)*a/3.14159265 << " " << sqrt(n*n + n * m + m*m)*a/3.14159265 << endl;
        }
    }
    diameterDifferences.shed_row(0);
    //sort allong the diameter differences, min to max.
    mat sorted = betterSortRows(diameterDifferences, 0);

    //save

    // cout << "SORTED " << sorted << endl;
    // cout << "POS " << positions << endl;
    // cout << "z max min " <<zCol.max()<<" " <<zCol.min() << endl;
    //   cout << "y max min " <<yCol.max()<<" " <<yCol.min() << endl;
    //cout << "DIAMTER " << diameter << endl;
    for (int i = 0; i < diameterDifferences.n_rows; i++) {
        string xS;
        stringstream ssX;
        ssX << sorted(i, 0) << "    " << sorted(i, 1) << "    " << sorted(i, 2);

        xS = ssX.str();
        system(("echo '" + xS + "' >> 0ApproximateDiameterIndexes.txt").c_str());
        // system(("echo '" + sorted(i,0) + " " + sorted(i,1) + " "+sorted(i,2) + "' >> ApproximateDiameterIndexes.txt").c_str());
    }
    //sorted.save("ApproximateDiameterIndexes.txt", arma_ascii);


    //add non-interacting atoms to our system
    if (haveNonInteracting == true) {
        nonInteractingAppender(positions, bonds, neighList);
    }
}

//void DataStructureFunctionsATAC::NAMDrun(int mutationGo, double& oldEnergy, double& temperature, mat& positions, imat& bondList, imat& neighborList, double& ratioRotate, double& ratioDimerDisplacement, double& ratioSp2ToSp3, int& spPossibility, string masterBatchName);

//void movieMaker(mat positions, string batchname);
//the flux refers to mass. Only needs the rotate ratio, as the DV and add dimer operations occur on a 1:1 ratio.

void DataStructureFunctionsATAC::noFluxMutation(int& spPossibility, double rotateRatio, double ratioDimerDisplacement, double ratioSp2ToSp3, imat& neighList, imat& bonds, mat& positions, double& randomPercent, double& oldEnergy, double& temperature, imat& neighborList, string masterBatchName, int& inBeam) {



    //remove non-interacting atoms to our system
    if (haveNonInteracting == true) {

        nonInteractingDeAppender(positions, bonds, neighList);

    }



    imat remainingBonds = bonds;

    if (randomPercent < rotateRatio) {
        for (int i = 0; i < rotationsPerRun; i++) {
            cout << "ROTATE" << endl;
            mutation(1, neighList, bonds, positions, remainingBonds, masterBatchName, inBeam);
            cout << "ROTATED" << endl;
        }

    } else if ((randomPercent > rotateRatio) && (randomPercent < rotateRatio + ratioDimerDisplacement)) {
        //dimer displacement.

        cout << "REMOVE DIMER" << endl;
        mutation(2, neighList, bonds, positions, remainingBonds, masterBatchName, inBeam);
        cout << "REMOVED DIMER" << endl;


        int inBeamDummy = -1;
        cout << "ADD DIMER" << endl;
        imat remainingBonds = bonds; //must refresh the remaining bonds.
        mutation(3, neighList, bonds, positions, remainingBonds, masterBatchName, inBeamDummy);
        //note, if we have a beam, and a dimer was removed from the beam, we are considering the add dimer reaction to always have the energy
        // of the beam, reguardless if the dimer is added in the beam or not. Reason being: the displaced
        // dimer would carry the energy from the beam to where ever it goes.
        // This is also why we need to make a dummy inBeam varible, as mutation() expects a pointer to an int (int&),
        // not an int value.

        cout << "ADDED DIMER" << endl;


    } else {

        //As a reminder, spPossibility -1 means a sp2 to sp3 transition is impossible,
        //0 means that we are not considering that transition, and
        //1 means that the transition was taken.

        //int spPossibility = 0;

        cout << "SP2 to SP3 START" << endl;
        sp2Tosp3(positions, neighList, bonds, spPossibility);
        cout << "SP2 to SP3 END" << endl;


        if (relaxAfterSp2Sp3Failure == true || spPossibility == 1) {

            cout << "RELAXING STRUCTURE AFTER SP2 TO SP3" << endl;

            //add non-interacting atoms to our system
            if (haveNonInteracting == true) {
                nonInteractingAppender(positions, bonds, neighList);
            }


            NAMDrun(0, oldEnergy, temperature, positions, bonds, neighborList, rotateRatio, ratioDimerDisplacement, ratioSp2ToSp3, spPossibility, masterBatchName);

            //removenon-interacting atoms to our system
            if (haveNonInteracting == true) {
                nonInteractingDeAppender(positions, bonds, neighList);
            }

        }






        //if the transition was taken, then we can also take the sp3 to sp2 transition.
        //cout << "SP POSSIBILITY " << spPossibility << endl;
        if (spPossibility == 1) {



            // take sp3 to sp2 reaction.
            cout << "SP3 to SP2 START" << endl;
            sp3Tosp2(neighList, bonds, positions);
            cout << "SP3 to SP2 END" << endl;
        }

        //spPossibility = 0; //reset the spPossibility number, as we have performed the operations needed.
    }
    //cout << "Before adding noninteracting atoms" << endl;
    //add non-interacting atoms to our system
    if (haveNonInteracting == true) {
        nonInteractingAppender(positions, bonds, neighList);
    }
    //cout << "After adding noninteracting atoms" << endl;
}

//This simulates chemical vapor depostion, which at the time is simply add dimer.

void DataStructureFunctionsATAC::CVD(imat& neighList, imat& bonds, mat& positions, double& oldEnergy, double& temperature, imat& neighborList, string masterBatchName, int& inBeam) {
    //remove non-interacting atoms to our system
    if (haveNonInteracting == true) {

        nonInteractingDeAppender(positions, bonds, neighList);

    }

    cout << "ADD DIMER" << endl;
    imat remainingBonds = bonds; //must refresh the remaining bonds.
    mutation(3, neighList, bonds, positions, remainingBonds, masterBatchName, inBeam);
    cout << "ADDED DIMER" << endl;

    //add non-interacting atoms to our system
    if (haveNonInteracting == true) {
        nonInteractingAppender(positions, bonds, neighList);
    }
}

double DataStructureFunctionsATAC::ratioShift(double& ratioRotate, double& ratioDimerDisplacment, double& ratioSp2ToSp3) {
    //if noFluxRatioShift type is 1, increase rotation at the expense of sp2 to sp3.
    if (noFluxRatioShift == 1) {
        //if we go beyond the bounds of the ratios we have set, take the cap values and disable ratio shifting.
        if (ratioRotate + ratioShiftRate > ratioRotateFinal || ratioSp2ToSp3 - ratioShiftRate < ratioSp2ToSp3Final) {
            //cout << "RATIOooooooooooos " << ratioRotate << " " << ratioRotateFinal << " " << ratioSp2ToSp3 << " " << ratioSp2ToSp3Final << " " << ratioShiftRate << endl;
            ratioRotate = ratioRotateFinal;
            ratioSp2ToSp3 = ratioSp2ToSp3Final;
            noFluxRatioShift = 0;
        } else {
            ratioRotate += ratioShiftRate;
            ratioSp2ToSp3 -= ratioShiftRate;
        }
    }
}


//bit of code that will run namd, and will always accept configuration (if NAMD doesn't fail)

void DataStructureFunctionsATAC::NAMDrun(int mutationGo, double& oldEnergy, double& temperature, mat& positions, imat& bondList, imat& neighborList, double& ratioRotate, double& ratioDimerDisplacement, double& ratioSp2ToSp3, int& spPossibility, string batchName) {


    cout << "NAMD RUN NUMBER: " << namdRunNum << endl;
    cout << batchName << endl;
    //We need to back up the current information, in case the introduced defect is reduced.
    mat positionsBackup = positions; //do I really need this line?
    imat neighborListBackup = neighborList;
    imat bondListBackup = bondList;
    // mat nonInteractingPositionsBackup = nonInteracting;
    // imat nonInteractingBondListBackup = nonInteractingBonds;

    int inBeam = -1;
    //inBeam tells us if the current mutation occured inside or outside the irradiation beam, if it is present.
    //0 means the mutation is outside, 1 means the mutation is inside, and -1 means that it is unknown, or unimportant.

    double randomPercent;
    //pure relaxation.
    if (mutationGo == 0) {
        cout << "RELAXING WITH NO MUTATION" << endl;
    } else
        //no flux mutation
        if (mutationGo == 1) {
        // cout << "MUT GO INSIDE" << mutationGo << endl;
        spPossibility = 0;
        randomPercent = (double) rand() / (double) (RAND_MAX);
        noFluxMutation(spPossibility, ratioRotate, ratioDimerDisplacement, ratioSp2ToSp3, neighborList, bondList, positions, randomPercent, oldEnergy, temperature, neighborList, batchName, inBeam);


    } else
        //Chemical Vapor Deposition
        if (mutationGo == 2) {
        CVD(neighborList, bondList, positions, oldEnergy, temperature, neighborList, batchName, inBeam);
    }

    // the only case we can ignore the following code is if we tried an sp2 to sp3 run that failed.
    // and we don't want to record a mutationless frame.

//BAD CODING AHEAD
    //need to hide these variables if namd is selected, but still keep them in scope if namd is not selected...
      //  min_input_t input;
      //  min_output_t output;

    

    if (spPossibility != -1 || relaxAfterSp2Sp3Failure == true) {

        if (forceNum == 1) {

            mat xValues = positions.col(0);
            double minWallPos = xValues.min() + atof(variableReader("minWallPosModifier").c_str());
            double maxWallPos = xValues.max() + atof(variableReader("maxWallPosModifier").c_str());
            double wallForce = atof(variableReader("wallForce").c_str());
            double tubeRadius = atof(variableReader("tubeRadius").c_str());
            double tubeForce = pow(atof(variableReader("tubeForceBase").c_str()), atof(variableReader("tubeForceExponent").c_str())); //force in kcal/mol*A. Should be 7^12.
            //tubeForce = 138412;

            mat forceMatrix;
            forceMatrix.copy_size(positions);
            //fills with zeros, as we are initilizing.
            forceMatrix.zeros();
            forceParallelWalls(positions, forceMatrix, minWallPos, maxWallPos, wallForce, 1);
            forceTube(positions, forceMatrix, tubeForce, tubeRadius);
            writeForcePDB(positions, forceMatrix);

        } else if (forceNum == 2) {

            double minViceWall = atof(variableReader("minViceWall").c_str());
            double maxViceWall = atof(variableReader("maxViceWall").c_str());
            double viceForce = atof(variableReader("viceForce").c_str());

            mat forceMatrix;
            forceMatrix.copy_size(positions);
            //fills with zeros, as we are initilizing.
            forceMatrix.zeros();
            forceVice(positions, forceMatrix, minViceWall, maxViceWall, viceForce);
            writeForcePDB(positions, forceMatrix);
        }

        cout << "Writing Files...." << endl;
        writeTopFile(positions, bondList, batchName);
        writePgnFile(positions, batchName);
        writePdbPsfFiles(batchName);
        cout << "Writing Files Done " << endl;
        // cout << "POSITON " << positions << endl << "NEIGH " << neighborList << endl << "BONDs " << bondList << endl;


        //remove old log file, as we are about to make a new one.
        system(("rm " + batchName + ".log").c_str());

        cout << "RUNNING NAMD" << endl;

        if (relaxationMethod == 1) {
            //NAMD selected
            system(("namd2 +p" + numProcessors + " " + batchName + ".conf >> " + batchName + ".log ").c_str());
        } else {
            //TB selected
            // minimization method (enumeration in header file)
 //           input.method = CT_OPT;

            // inputs
   //         input.bondList = bondList;
//            input.positions = positions;

            // exit conditions, set either to 0 to disable
  //          input.interation_limit = 250;
    //        input.max_force_cutoff = 1e-3;

            // periodicity options, 0 to disable
      //      input.x_period = xLatticeParameter;
        //    input.y_period = yLatticeParameter;
//            input.z_period = zLatticeParameter;

            // pass to the minimize function and get an output structure
  //          output = minimize(input);

        }
    }

    //if we have a beam, and the mutation takes place in the beam region, we must use the beam temperature and not the ambient one.
    double ambientTempBackUp = temperature;
    if (haveBeam == true && inBeam == true) {
        temperature = beamTemperature;
    }


    //cout << "SP FSFDFS " << spPossibility << endl;
    if (energyCheck(oldEnergy, temperature, spPossibility, randomPercent, ratioRotate, ratioDimerDisplacement, batchName, namdRunNum, relaxAfterSp2Sp3Failure, boltzmann, energyFailureTolerance, acceptAnyMutation, relaxationMethod)) {
        // cout << "ENERGY WIN " << endl;
        //if energy check is passes, updates positions to be used for the next round
        //of defects, and adds the frame to the xyz movie.

        //copy over sucessful log file to movies directory.
        system(("cp " + batchName + ".log movies").c_str());

        //place "if coor not found" error check here.
        if (relaxationMethod == 1) {
            //Get positions from Namd Run
            coorStripper(positions, batchName);
        } else {
            //Get positions from TB
    //        positions = output.positions;
        }

        movieMaker(positions, batchName, numberOfElectrodes, xp1Electrodes, xp2Electrodes, xn1Electrodes, xn2Electrodes, yp1Electrodes, yp2Electrodes, yn1Electrodes, yn2Electrodes);

        ++acceptedRunNum;
        ++namdRunNum;

        string runNumS;
        stringstream ssRun;
        ssRun << namdRunNum;
        runNumS = ssRun.str();

        string successfulRunNumS;
        stringstream ssSRun;
        ssSRun << acceptedRunNum;
        successfulRunNumS = ssSRun.str();

        string oldEnergyS;
        stringstream ssEner;
        ssEner << oldEnergy;
        oldEnergyS = ssEner.str();



        //record energy
        chdir("energies");
        system(("echo '" + runNumS + " " + oldEnergyS + "' >> energyHist" + batchName + ".dat").c_str());
        system(("echo '" + successfulRunNumS + " " + oldEnergyS + "' >> energyHistMovie" + batchName + ".dat").c_str());
        chdir("..");




        //record chirality
        if (findChiralAngles == true) {

            //remove non-interacting atoms to our system
            if (haveNonInteracting == true) {
                nonInteractingDeAppender(positions, bondList, neighborList);
            }

            //generate rings.
            vector<pair<imat, int> > rings;
            rings = findRings(neighborList, ringSizeMin, ringSizeMax);
            //using this just to generate the angle.dat files, not for the final chiral to index conversion.
            recordChirality(rings, positions, false);
            //record ring sizes
            if (findRingSizes == true) {
                recordRingSizes(rings, globalRingSizeRunNum, generateRingHistogramImages);
            }

            //add non-interactig atoms back.
            if (haveNonInteracting == true) {
                nonInteractingAppender(positions, bondList, neighborList);
            }

        }

        //record ring sizes, if haven't done so already in the chiral angle finder.
        if ((findRingSizes == true) && (findChiralAngles != true)) {

            //remove non-interacting atoms to our system
            if (haveNonInteracting == true) {
                nonInteractingDeAppender(positions, bondList, neighborList);
            }
            //generate rings.
            vector<pair<imat, int> > rings;
            rings = findRings(neighborList, ringSizeMin, ringSizeMax);
            recordRingSizes(rings, globalRingSizeRunNum, generateRingHistogramImages);
            //add non-interactig atoms back.
            if (haveNonInteracting == true) {
                nonInteractingAppender(positions, bondList, neighborList);
            }

        }

        //adjust beam radius.
        beamRadius += beamRadiusGrowthRate;


        //restore ambient temperature

        if (haveBeam == true && inBeam == true) {
            temperature = ambientTempBackUp;
        }


        if (coolMethod != 0) {
            //adjust temperature

            if (coolMethod == 1) {
                temperature = linearCooling(tempS, tempF, acceptedRunNum, coolR);
            } else if (coolMethod == 2) {
                temperature = newtonianCooling(tempS, tempF, acceptedRunNum, coolR);
            }

            //string with updated energy
            stringstream ssTemp;
            ssTemp << temperature;
            string tempStr = ssTemp.str();

            system(("sed -i 's/set temperature.*/set temperature " + tempStr + "/' " + batchName + ".conf").c_str()); //" + batchName + ".conf"
        }

        if (noFluxRatioShift != 0) {
            ratioShift(ratioRotate, ratioDimerDisplacement, ratioSp2ToSp3);
            //cout << "rrrrrrrrrrrrrrratioSP2ToSp3 " << ratioSp2ToSp3 << endl;
        }

        if (periodicStretching != 0) {
            //adjust periodicity

            //see if we are updating in an orthogonal or arbitrary periodicity scheme.
            if (periodicityNum < 4) {
                if (periodicStretching == 1) {
                    //                   cout << "X LAT PARAM IS " << xLatticeParameter << endl;
                    //   system("read -p \"paused\" key");

                    xLatticeParameter += periodicStretchRate;
                    stringstream ssX;
                    ssX << setprecision(15) << xLatticeParameter;
                    string XS = ssX.str();

                    //assumes that system is periodic in at least x
                    system(("sed -i 's/cellBasisVector1.*/cellBasisVector1   " + XS + "    0    0/' " + batchName + ".conf").c_str());


                }
                //assumes system is at least perodic in x and y directions.
                if (periodicStretching == 2) {

                    xLatticeParameter += periodicStretchRate;

                    //cout << "Y LAT PARAM before" << yLatticeParameter <<  " rate stretch "<< periodicStretchRate <<endl;
                    yLatticeParameter = yLatticeParameter + periodicStretchRate;
                    //cout << "Y LAT PARAM after" << yLatticeParameter << endl;



                    stringstream ssX;
                    ssX << setprecision(15) << xLatticeParameter;
                    string XS = ssX.str();



                    stringstream ssY;
                    ssY << setprecision(15) << yLatticeParameter;
                    string YS = ssY.str();
                    //cout <<"STRING " << YS << endl;

                    system(("sed -i 's/cellBasisVector1.*/cellBasisVector1   " + XS + "    0    0/' " + batchName + ".conf").c_str());
                    system(("sed -i 's/cellBasisVector2.*/cellBasisVector2    0   " + YS + "    0/' " + batchName + ".conf").c_str());
                }

                if (periodicStretching == 3) {

                    xLatticeParameter += periodicStretchRate;

                    //cout << "Y LAT PARAM before" << yLatticeParameter <<  " rate stretch "<< periodicStretchRate <<endl;
                    yLatticeParameter = yLatticeParameter + periodicStretchRate;
                    //cout << "Y LAT PARAM after" << yLatticeParameter << endl;
                    zLatticeParameter = yLatticeParameter + periodicStretchRate;


                    stringstream ssX;
                    ssX << setprecision(15) << xLatticeParameter;
                    string XS = ssX.str();



                    stringstream ssY;
                    ssY << setprecision(15) << yLatticeParameter;
                    string YS = ssY.str();

                    stringstream ssZ;
                    ssZ << setprecision(15) << zLatticeParameter;
                    string ZS = ssZ.str();
                    //cout <<"STRING " << YS << endl;

                    system(("sed -i 's/cellBasisVector1.*/cellBasisVector1   " + XS + "    0    0/' " + batchName + ".conf").c_str());
                    system(("sed -i 's/cellBasisVector2.*/cellBasisVector2    0   " + YS + "    0/' " + batchName + ".conf").c_str());
                    system(("sed -i 's/cellBasisVector3.*/cellBasisVector3    0   " + ZS + "    0/' " + batchName + ".conf").c_str());
                }
            }

            if (periodicityNum == 4) {

                if (periodicStretching == 1) {
                    latticeVec1(0) += periodicStretchRate; //lattice vector 1, dimension 1.
                    //cout << "LATTICE VEC 1 " << latticeVec1(0) << endl;
                    //find the normalized values.
                    double normalizationFactor;
                    normalizationFactor = abs(latticeVec1(0)) + abs(latticeVec1(1)) + abs(latticeVec1(2));
                    for (int i = 0; i < 3; i++) {
                        latticeVec1Norm(i) = latticeVec1(i) / normalizationFactor;
                    }


                    stringstream commandSS;
                    commandSS << setprecision(15) << "sed -i 's/cellBasisVector1.*/cellBasisVector1   " << latticeVec1(0) << "    " << latticeVec1(1) << "    " << latticeVec1(2) << "/' " << batchName << ".conf";
                    string Command;
                    Command.append(commandSS.str());
                    system(Command.c_str());

                    // system(("sed -i 's/cellBasisVector1.*/cellBasisVector1   " + latticeVec1(0) + "    " + latticeVec1(1) + "    " + latticeVec1(2) + "/' " + batchName + ".conf").c_str());
                }

                if (periodicStretching == 2) {
                    latticeVec1(0) += periodicStretchRate; //lattice vector 1, dimension 1.

                    latticeVec2(1) += periodicStretchRate;

                    //find the normalized values.
                    double normalizationFactor;
                    normalizationFactor = abs(latticeVec1(0)) + abs(latticeVec1(1)) + abs(latticeVec1(2));
                    for (int i = 0; i < 3; i++) {
                        latticeVec1Norm(i) = latticeVec1(i) / normalizationFactor;
                    }
                    normalizationFactor = latticeVec2(0) + latticeVec2(1) + latticeVec2(2);
                    for (int i = 0; i < 3; i++) {
                        latticeVec2Norm(i) = latticeVec2(i) / normalizationFactor;
                    }

                    stringstream command1SS;
                    command1SS << setprecision(15) << "sed -i 's/cellBasisVector1.*/cellBasisVector1   " << latticeVec1(0) << "    " << latticeVec1(1) << "    " << latticeVec1(2) << "/' " << batchName << ".conf";
                    string Command1;
                    Command1.append(command1SS.str());
                    system(Command1.c_str());

                    stringstream command2SS;
                    command2SS << setprecision(15) << "sed -i 's/cellBasisVector2.*/cellBasisVector2   " << latticeVec2(0) << "    " << latticeVec2(1) << "    " << latticeVec2(2) << "/' " << batchName << ".conf";
                    string Command2;
                    Command2.append(command2SS.str());
                    system(Command2.c_str());
                }




                if (periodicStretching == 3) {
                    latticeVec1(0) += periodicStretchRate; //lattice vector 1, dimension 1.
                    latticeVec2(1) += periodicStretchRate;
                    latticeVec3(2) += periodicStretchRate;
                    //find the normalized values.
                    double normalizationFactor;
                    normalizationFactor = abs(latticeVec1(0)) + abs(latticeVec1(1)) + abs(latticeVec1(2));
                    for (int i = 0; i < 3; i++) {
                        latticeVec1Norm(i) = latticeVec1(i) / normalizationFactor;
                    }
                    normalizationFactor = latticeVec2(0) + latticeVec2(1) + latticeVec2(2);
                    for (int i = 0; i < 3; i++) {
                        latticeVec2Norm(i) = latticeVec2(i) / normalizationFactor;
                    }
                    normalizationFactor = latticeVec3(0) + latticeVec3(1) + latticeVec3(2);
                    for (int i = 0; i < 3; i++) {
                        latticeVec3Norm(i) = latticeVec3(i) / normalizationFactor;
                    }

                    stringstream command1SS;
                    command1SS << setprecision(15) << "sed -i 's/cellBasisVector1.*/cellBasisVector1   " << latticeVec1(0) << "    " << latticeVec1(1) << "    " << latticeVec1(2) << "/' " << batchName << ".conf";
                    string Command1;
                    Command1.append(command1SS.str());
                    system(Command1.c_str());

                    stringstream command2SS;
                    command2SS << setprecision(15) << "sed -i 's/cellBasisVector2.*/cellBasisVector2   " << latticeVec2(0) << "    " << latticeVec2(1) << "    " << latticeVec2(2) << "/' " << batchName << ".conf";
                    string Command2;
                    Command2.append(command2SS.str());
                    system(Command2.c_str());

                    stringstream command3SS;
                    command3SS << setprecision(15) << "sed -i 's/cellBasisVector3.*/cellBasisVector3   " << latticeVec3(0) << "    " << latticeVec3(1) << "    " << latticeVec3(2) << "/' " << batchName << ".conf";
                    string Command3;
                    Command3.append(command3SS.str());
                    system(Command3.c_str());
                }
            }
        }


    } else {
        //if energy check is not tripped, restores values

        ++namdRunNum;

        string runNumS;
        stringstream ssRun;
        ssRun << namdRunNum;
        runNumS = ssRun.str();

        string successfulRunNumS;
        stringstream ssSRun;
        ssSRun << acceptedRunNum;
        successfulRunNumS = ssSRun.str();

        string oldEnergyS;
        stringstream ssEner;
        ssEner << oldEnergy;
        oldEnergyS = ssEner.str();



        //record energy
        chdir("energies");
        system(("echo '" + runNumS + " " + oldEnergyS + "' >> energyHist" + batchName + ".dat").c_str());
        system(("echo '" + successfulRunNumS + " " + oldEnergyS + "' >> energyHistMovie" + batchName + ".dat").c_str());
        chdir("..");



        neighborList = neighborListBackup;
        bondList = bondListBackup;
        positions = positionsBackup;
    }


    //restore ambient temperature

    if (haveBeam == true && inBeam == true) {
        temperature = ambientTempBackUp;
    }

    cout << "---------------------------------------------------------" << endl;

}

//This subfunction writes generates a datafile to be used in final movie
//showing the evolution of the chiral angles for all 6 rings.
//assumes tube axis is in x direction.

void DataStructureFunctionsATAC::recordChirality(vector<pair<imat, int> > rings, mat pos, bool findIndexes) {

    //find if rings of six exist
    for (int i = 0; i < rings.size(); ++i) {
        if (rings[i].second == 6) {
            chdir("chiralAngles");
            imat sixRings = rings[i].first;
            //start file with chiral angle points.
            ofstream writeFile;

            string globalRunNumS;
            stringstream ssglobalRunNum;
            ssglobalRunNum << globalChiralAngleRunNum;
            globalRunNumS = ssglobalRunNum.str();

            string globalRunNum1S;
            stringstream ssglobalRunNum1;
            ssglobalRunNum1 << globalChiralAngleRunNum + 1;
            globalRunNum1S = ssglobalRunNum1.str();

            mat chiralList;
            chiralList << -1 << endr;

            writeFile.open((globalRunNum1S + ".dat").c_str());
            // loop through each ring...
            for (int j = 0; j < sixRings.n_rows /*&& j==0*/; ++j) {


                mat NNVec;
                //junk initilization value
                NNVec << -1 << -1 << -1 << endr;

                //loop through each atom in the ring...
                for (int k = 0; k < sixRings.n_cols; ++k) {
                    //loop through each atom again
                    for (int l = k; l < sixRings.n_cols; ++l) {
                        if (k != l) {
                            double distance;
                            //compare the two atoms, and see if they are within NNN distance of each other.
                            shortestLineDeterminer(pos, sixRings(j, k), sixRings(j, l), distance, periodicityNum, xLatticeParameter, yLatticeParameter, zLatticeParameter, latticeVec1Norm, latticeVec1, latticeVec2Norm, latticeVec2, latticeVec3Norm, latticeVec3);
                            if ((lowerNNLimit < distance) && (distance < higherNNLimit)) {
                                //generate "zig zag" direction vector
                                mat vectorToBeAppended;
                                vectorToBeAppended << pos(sixRings(j, k), 0) - pos(sixRings(j, l), 0) << pos(sixRings(j, k), 1) - pos(sixRings(j, l), 1) << pos(sixRings(j, k), 2) - pos(sixRings(j, l), 2) << endr;
                                NNVec = join_cols(NNVec, vectorToBeAppended);
                                // cout << "KAY EL " << k << " " << l << " "<< endl;
                                // cout <<pos(sixRings(j, k), 0) << " "<<pos(sixRings(j, k), 1) << " "<< pos(sixRings(j, k), 2)<<" "<< endl;
                                //cout <<pos(sixRings(j, l), 0) << " "<<pos(sixRings(j, l), 1) << " "<< pos(sixRings(j, l), 2)<<" "<< endl;
                                //cout << pos(sixRings(j, k), 0) - pos(sixRings(j, l), 0) <<" "<< pos(sixRings(j, k), 1) - pos(sixRings(j, l), 1) <<" "<< pos(sixRings(j, k), 2) - pos(sixRings(j, l), 2) << endl;


                            }
                        }
                    }
                }
                //shed junk initilization value
                NNVec.shed_row(0);

                //check we got at least one vector. This warning message should never trigger, since
                //we are looking at a list of rings. should be 3 vectors.

                if (NNVec.is_empty()) {

                    cout << "CANNOT FIND ANY NEXT NEARIST NEIGHBORS, THUS FINDING CHIRAL ANGLES IMPOSSIBLE" << endl;
                    exit;
                }
                mat NNVecAveraged;

                int useAveragedVectors = atoi(variableReader("useAveragedVectors").c_str());

                if (useAveragedVectors == 0) {
                    NNVecAveraged = NNVec;


                } else {
                    //we have 6 vectors, one for each zig-zag edge of the hexagon.
                    //However, it's really 3 pairs of parrallel vectors.
                    //But these vectors are not parrallel in most real cases, due to
                    //warping of the hexagon. So, we take an average of the pairs that
                    //are closest. This also makes our results more accurate, as
                    //often the warping can be canceled by such averaging.



                    NNVecAveraged << -1 << -1 << -1 << endr;

                    //iter over vector 1...
                    for (int vecIter1 = 0; vecIter1 < NNVec.n_rows; vecIter1++) {
                        mat averagedVector;
                        //first index is normalized dot product, or similarity. Second and third are the vec Iters.
                        //last 3 are the averaged vector.
                        averagedVector << -1 << -1 << -1 << -1 << -1 << -1 << endr;

                        //iter over vector 2...
                        for (int vecIter2 = vecIter1; vecIter2 < NNVec.n_rows; vecIter2++) {
                            // cout << "Vect Iter 1, 2, and num rows " << vecIter1 << " " << vecIter2 << " " << NNVec.n_rows << endl;
                            if (vecIter1 != vecIter2) {
                                double similarity;
                                similarity = abs(dot(NNVec.row(vecIter1) /
                                        sqrt(NNVec(vecIter1, 0) * NNVec(vecIter1, 0) + NNVec(vecIter1, 1) * NNVec(vecIter1, 1) + NNVec(vecIter1, 2) * NNVec(vecIter1, 2)),
                                        NNVec.row(vecIter2) /
                                        sqrt(NNVec(vecIter2, 0) * NNVec(vecIter2, 0) + NNVec(vecIter2, 1) * NNVec(vecIter2, 1) + NNVec(vecIter2, 2) * NNVec(vecIter2, 2))
                                        ));
                                //    cout << "SIMILARITY " << similarity << endl << endl;
                                //we record the vector with the greatest similiarity
                                if (similarity > averagedVector(0)) {
                                    //before we take the average, keep in mind that one of the vectors might be the
                                    // negative of the other
                                    //If the two vector components multiplied against each other are negative,
                                    //and are above a threshhold size (because of rounding errors near zero), then muliply one vector by -1

                                    if (((NNVec(vecIter1, 0) * NNVec(vecIter2, 0) < 0) && abs(NNVec(vecIter1, 0)) >= vecRoundingErrorCutOff && abs(NNVec(vecIter2, 0)) >= vecRoundingErrorCutOff) ||
                                            ((NNVec(vecIter1, 1) * NNVec(vecIter2, 1) < 0) && abs(NNVec(vecIter1, 1)) >= vecRoundingErrorCutOff && abs(NNVec(vecIter2, 1)) >= vecRoundingErrorCutOff) ||
                                            ((NNVec(vecIter1, 2) * NNVec(vecIter2, 2) < 0) && abs(NNVec(vecIter1, 2)) >= vecRoundingErrorCutOff && abs(NNVec(vecIter2, 2)) >= vecRoundingErrorCutOff)) {
                                        //vectors are antiparallel

                                        averagedVector << similarity << vecIter1 << vecIter2
                                                << (-1 * NNVec(vecIter1, 0) + NNVec(vecIter2, 0)) / 2
                                                << (-1 * NNVec(vecIter1, 1) + NNVec(vecIter2, 1)) / 2
                                                << (-1 * NNVec(vecIter1, 2) + NNVec(vecIter2, 2)) / 2 << endr;
                                    } else {
                                        //vectors are parallel.
                                        averagedVector << similarity << vecIter1 << vecIter2
                                                << (NNVec(vecIter1, 0) + NNVec(vecIter2, 0)) / 2
                                                << (NNVec(vecIter1, 1) + NNVec(vecIter2, 1)) / 2
                                                << (NNVec(vecIter1, 2) + NNVec(vecIter2, 2)) / 2 << endr;
                                        cout << "VEC 1 " << NNVec.row(vecIter1) << endl << "VEC 2 " << NNVec.row(vecIter2) << endl;
                                    }
                                }

                            }
                        }
                        //             cout << "AVERAGED VECTOR " << averagedVector << endl;
                        NNVec.shed_row(averagedVector(2));
                        NNVec.shed_row(averagedVector(1));
                        vecIter1--;



                        NNVecAveraged = join_cols(NNVecAveraged, averagedVector.cols(3, 5));


                    }
                    NNVecAveraged.shed_row(0);
                    //cout << "AVERAGED  " << NNVecAveraged << endl;
                }
                //*/

                //assume tube is x aligned:
                rowvec tubeAxis;
                tubeAxis << 1 << 0 << 0 << endr;

                //find smallest chiral angle (largest angle possible should be 30 deg,
                // so start with the assumption that it's 90 deg, and we'll update it
                //with with smaller values as we find them).
                double chiralAngle = 3.14159265 / 2;

                //          cout << "NN VEC " << NNVecAveraged << endl;
                for (int vecIter = 0; vecIter < NNVecAveraged.n_rows; vecIter++) {
                    //phi is angle between the NNN vector and the tube axis
                    double phi;
                    rowvec NNVecRow = NNVecAveraged.row(vecIter);
                    //phi is angle between the NNN vector and the tube axis
                    phi = acos(dot(tubeAxis, NNVecRow) / (1 * sqrt(NNVecRow(0) * NNVecRow(0) + NNVecRow(1) * NNVecRow(1) + NNVecRow(2) * NNVecRow(2))));

                    //we need the smallest angle of the intersection of these two vectors,
                    //so find the acute complimentary of phi if it is oblique.
                    // if (phi < 3.14259265 / 2) {
                    //       phi = 3.14259265 - phi;
                    // }
                    //  cout <<"PHI "<< abs(phi)*180/3.141 << endl;
                    //theta is the chiral angle. We must remember that the
                    //chiral angle is that between the zig-zag direction (NNN vector)
                    //and the tube diameter (the plane normal to the tube axis).
                    //to convert phi into theta, do the simple conversion:
                    double theta;
                    if (phi < 3.14159265 / 2) {
                        theta = 3.14159265 / 2 - phi;
                    } else {
                        theta = phi - 3.14159265 / 2;
                    }
                    //              cout <<"THATA "<< abs(theta)*180/3.141 << endl;
                    //find smallest chiral angle out of all zig-zag directions.
                    if (abs(theta) < chiralAngle) {

                        //theta should always be positive, but due to rounding, small angles could be negative.
                        chiralAngle = abs(theta);

                    }
                }

                mat toBeAppended;
                toBeAppended << chiralAngle * 180 / 3.1412 << endr;
                chiralList = join_cols(chiralList, toBeAppended);

                writeFile << chiralAngle * 180 / 3.14159265 << endl;
            }
            writeFile.close();
            chiralList.shed_row(0);
            //    cout << "CHIRAL LIST " << chiralList << endl;
            if (findIndexes == true) {


                double a = 2.46; //in angstroms.
                //double pi = 3.14159265;

                mat anglesAndIndexes;
                anglesAndIndexes << -1 << -1 << -1 << -1 << endr;

                //creat a matrix recording how many rings of a given n and m are seen.
                mat counts;
                counts << -1 << -1 << -1 << endr;
                for (int n = minChiralIndex; n <= maxChiralIndex; ++n) {
                    for (int m = minChiralIndex - 1; m <= n; m++) {
                        mat toBeAppended;
                        toBeAppended << 0 << n << m << endr;
                        counts = join_cols(counts, toBeAppended);

                    }
                }
                counts.shed_row(0);

                // for each angle...
                for (int angleIter = 0; angleIter < chiralList.n_rows; angleIter++) {
                    mat angleDifferences;
                    angleDifferences << -1 << -1 << -1 << -1 << endr;
                    for (int n = minChiralIndex; n <= maxChiralIndex; ++n) {
                        for (int m = minChiralIndex - 1; m <= n; m++) {
                            //must remember that N must always be greater than M in the nanotube indexed scheme.
                            if (n >= m) {
                                mat toBeAppended;
                                toBeAppended << abs(chiralList(angleIter) - atan(sqrt(3) * m / (2 * n + m))*180 / 3.14159265) << chiralList(angleIter) << n << m << endr;
                                angleDifferences = join_cols(angleDifferences, toBeAppended);
                                //cout << "N, M, ANGLE " << n << " " << m << " " << atan(sqrt(3) * m / (2 * n + m))*180 / 3.14159265 << " "<< endl;
                            }
                        }
                    }

                    angleDifferences.shed_row(0);
                    //sort to find the best match of a given angle to the index list.
                    mat sortedForAnAngle = betterSortRows(angleDifferences, 0);
                    // cout <<"SORTEDFORANANGLE " << sortedForAnAngle << endl;

                    bool findChiralSpread = atoi(variableReader("findChiralSpread").c_str());

                    if (findChiralSpread == 0) {
                        //best chiral fit found. will have the least difference
                        //add to our final list of indexes
                        mat toBeAppendedIndexes;
                        toBeAppendedIndexes << sortedForAnAngle(0, 0) << sortedForAnAngle(0, 1) << sortedForAnAngle(0, 2) << sortedForAnAngle(0, 3) << endr;
                        anglesAndIndexes = join_cols(anglesAndIndexes, toBeAppendedIndexes);

                        //have to find matching indexes on count, then add to it.
                        for (int i = 0; i < counts.n_rows; i++) {
                            if (counts(i, 1) == sortedForAnAngle(0, 2) && counts(i, 2) == sortedForAnAngle(0, 3)) {
                                counts(i, 0) = counts(i, 0) + 1;
                                break;
                            }
                        }

                    } else {
                        //same as immediatly above, except we record multiple possible theoretical indexes for the experimental angle.
                        double maxChiralAngleDifference = atof(variableReader("maxChiralAngleDifference").c_str());


                        for (int sortedIter = 0; sortedIter < sortedForAnAngle.n_rows; sortedIter++) {


                            if (sortedForAnAngle(sortedIter, 0) < maxChiralAngleDifference) {
                                mat toBeAppendedIndexes;
                                toBeAppendedIndexes << sortedForAnAngle(sortedIter, 0) << sortedForAnAngle(sortedIter, 1) << sortedForAnAngle(sortedIter, 2) << sortedForAnAngle(sortedIter, 3) << endr;
                                anglesAndIndexes = join_cols(anglesAndIndexes, toBeAppendedIndexes);

                                //have to find matching indexes on count, then add to it.
                                for (int i = 0; i < counts.n_rows; i++) {
                                    if (counts(i, 1) == sortedForAnAngle(sortedIter, 2) && counts(i, 2) == sortedForAnAngle(sortedIter, 3)) {
                                        counts(i, 0) = counts(i, 0) + 1;
                                        break;
                                    }
                                }
                            }
                        }

                    }

                }

                anglesAndIndexes.shed_row(0);

                //cout << "ANGLDS " << anglesAndIndexes << "FIN " << endl;

                // counts.shed_row(0);

                mat sortedForCount = betterSortRows(counts, 0);


                //save count
                //reverse iter for count.
                for (int i = sortedForCount.n_rows - 1; i >= 0; i--) {
                    string xS;
                    stringstream ssX;
                    ssX << sortedForCount(i, 0) << "    " << sortedForCount(i, 1) << "    " << sortedForCount(i, 2);

                    xS = ssX.str();
                    system(("echo '" + xS + "' >> 0IndexCountAngleBased.txt").c_str());
                    // system(("echo '" + sorted(i,0) + " " + sorted(i,1) + " "+sorted(i,2) + "' >> ApproximateDiameterIndexes.txt").c_str());
                }

                string xS1;
                stringstream ssX1;
                ssX1 << "Median Difference From Ideal: " << median(anglesAndIndexes.col(0)) << " Median Angle: " << median(anglesAndIndexes.col(1)) << " Average Difference From Ideal: " << mean(anglesAndIndexes.col(0)) << " Average Angle: " << mean(anglesAndIndexes.col(1));

                xS1 = ssX1.str();
                system(("echo '" + xS1 + "' >> 0EachRingsAnglesAndIndexes.txt").c_str());

                //save angles and indexes.
                for (int i = 0; i < anglesAndIndexes.n_rows; i++) {
                    string xS;
                    stringstream ssX;
                    ssX << anglesAndIndexes(i, 0) << "    " << anglesAndIndexes(i, 1) << "    " << anglesAndIndexes(i, 2) << "    " << anglesAndIndexes(i, 3);


                    xS = ssX.str();
                    system(("echo '" + xS + "' >> 0EachRingsAnglesAndIndexes.txt").c_str());
                    // system(("echo '" + sorted(i,0) + " " + sorted(i,1) + " "+sorted(i,2) + "' >> ApproximateDiameterIndexes.txt").c_str());
                }
            }

            //we must label the jpegs of the plots in an appropriate manner to be
            //converted into an mpeg later.
            string label;
            stringstream number;
            if ((globalChiralAngleRunNum + 1) < 10) {
                label = "0000";
            } else if ((globalChiralAngleRunNum + 1) < 100) {
                label = "000" + number.str();
            } else if ((globalChiralAngleRunNum + 1) < 1000) {
                label = "00" + number.str();
            } else if ((globalChiralAngleRunNum + 1) < 10000) {
                label = "0" + number.str();
            } else if ((globalChiralAngleRunNum + 1) < 100000) {
                label = "" + number.str();
            } else {
                cerr << "MORE THAN 99,999 frames" << endl;
            }
            number << (globalChiralAngleRunNum + 1);
            label.append(number.str());

            if (generateRingHistogramImages == 1) {
                //updates the histogram script with the current run data, executes it, and generates a jpeg.
                system(("sed \"s/0.dat/" + globalRunNum1S + ".dat/g\" histogramScript.plt | gnuplot > pic" + label + ".jpg").c_str());
            }

            globalChiralAngleRunNum++;
            chdir("..");

        }
    }
}


