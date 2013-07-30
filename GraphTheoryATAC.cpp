#include "GraphTheoryATAC.h"

// This subfunction takes in a 1 by N vector of doubles, and returns a vector
// of the same size with the low to high index of numbers.

vec betterIndexSort(vec toBeSorted) {
    static int index = 0;
    vector<pair<double, int> > myVector(toBeSorted.n_elem);

    for (int i = 0; i < toBeSorted.n_elem; i++) {
        pair<double, int> valueAndIndex(toBeSorted(i), i);
        myVector[i] = valueAndIndex;
    }

    sort(myVector.begin(), myVector.end());
    vec lowToHighIndex(toBeSorted.n_elem);

    for (int i = 0; i < toBeSorted.n_elem; i++) {
        lowToHighIndex(i) = myVector[i].second;
    }

    return lowToHighIndex;
}



// Sorts rows based on values in the column colNum, low to high.

mat betterSortRows(mat matrix, int colNum) {
    vec lowHighIndex = betterIndexSort(matrix.col(colNum)); //sorts in ascending order.
    mat sortedMatrix;
    sortedMatrix.copy_size(matrix);

    for (int i = 0; i < lowHighIndex.n_elem; i++) {
        for (int j = 0; j < matrix.n_cols; j++) {
            sortedMatrix(i, j) = matrix(lowHighIndex(i), j);
        }
    }
    return (sortedMatrix);
}


//This function is the reverse of the previous function; it takes a bond list
//and generates a neighbor list.

void buildNeighFromBonds(imat bonds, imat& neighToBeRelplaced) {
    // this neighbor list has one extra column in the front;
    //it is the atom for which the rest in the columns are neighbors.
    // It is used to sort the matrix and is shed later.

    imat neigh(1, 5);
    neigh << -1 << -1 << -1 << -1 << -1 << endr; //initializes with junk values to be shed later.

    //loops over bonds (the rows)
    for (int i = 0; i < bonds.n_rows; ++i) {
        // loops over bond atoms (the columns) should be 2
        for (int j = 0; j < bonds.n_cols; ++j) {

            //checks for atom already used as a base (.is_empty == 0). skips if so.
            uvec alreadyCountedAtoms = find(bonds(i, j) == neigh.col(0));
            if (alreadyCountedAtoms.is_empty() == 0) {
                continue;
            }
            // set base atom. To be appended to neigh later.
            imat baseAndNeigh(1, 1);
            baseAndNeigh << bonds(i, j) << endr;

            //finds all occurances of base atom, adds neighbors to base atom, and appends onto neigh.
            imat baseAtomFinder(bonds.n_rows, 1);
            baseAtomFinder.fill(bonds(i, j));

            // looks through each column, finds the base atom,
            // and adds the atom that it's bonded to (on other side of column)
            //to neigh list.

            uvec atomIndexColumn0 = find(baseAtomFinder == bonds.col(0));
            // loops over found atoms index.
            for (int k = 0; k < atomIndexColumn0.n_rows; ++k) {
                imat toBeAppended;
                toBeAppended << bonds(atomIndexColumn0[k], 1) << endr;
                baseAndNeigh = join_rows(baseAndNeigh, toBeAppended);
            }

            //for opposite column
            uvec atomIndexColumn1 = find(baseAtomFinder == bonds.col(1));
            for (int l = 0; l < atomIndexColumn1.n_rows; ++l) {
                imat toBeAppended;
                toBeAppended << bonds(atomIndexColumn1[l], 0) << endr;
                baseAndNeigh = join_rows(baseAndNeigh, toBeAppended);
            }

            //If an atom does not have as many neighbors as the rest of the
            //neighbor matrix atoms, fills with a negative value
            //(non-existant neighbor).
            while (baseAndNeigh.n_cols < neigh.n_cols) {
                imat nonExistantNeighbor;
                nonExistantNeighbor << -1 << endr;
                baseAndNeigh = join_rows(baseAndNeigh, nonExistantNeighbor);
            }
            //add to current neighbor list.
            neigh = join_cols(neigh, baseAndNeigh);
        }
    }
    neigh.shed_row(0); //sheds initiation values.
    //sorts neigh based on atoms
    neigh = conv_to<imat>::from(betterSortRows(conv_to<mat>::from(neigh), 0));
    neigh.shed_col(0);

    neighToBeRelplaced = neigh;
}


// This subfunction returns true if the vector contains at least one instance of the int.
// toBeSearched is a Nx1 or 1xN matrix.

bool contains(imat toBeSearched, int tryToFind) {
    uvec tryToFindVec;
    tryToFindVec << tryToFind << endr;
    uvec Index = find(toBeSearched == tryToFind);
    //Sees if the found indexes exist or not, then flips the bool (essentially .is_present)
    bool found = Index.is_empty();
    found = !found;
    return found;
}

//bit of code to eliminate duplicate rings in in list
//works for arbitrarly sized rings.

void eliminateDuplicates(imat& rings) {

    int ringListSize = rings.n_rows;

    // i iterates over individual rings, or rows to start check with
    for (int i = 0; i < ringListSize; ++i) {

        // k iterates over individual rings, or rows to start check against
        for (int k = i + 1; k < ringListSize; ++k) {
            int numOfDupes = 0; // number of duplicates found.

            // j iterates over the atoms in an individual ring, or columns
            for (int j = 0; j < rings.n_cols; ++j) {
                imat atomToCheckVector;
                imat sizeReference;
                sizeReference = rings.row(i);
                atomToCheckVector.copy_size(sizeReference);
                atomToCheckVector.fill(rings(i, j));

                umat matchBoolList = (atomToCheckVector == rings.row(k));

                numOfDupes += accu(matchBoolList); // adds match bool to number of duplicates.
            }

            //if number of duplicates equals the size of the ring, then the two rings are identical.
            if (numOfDupes == rings.n_cols) {

                rings.shed_row(k);
                --ringListSize;
                --k; // have to set "to check against" iterater back, since the row to be removed is actually replaced.
            }

        }
    }


}

// generates a matrix with the row as atom number, and the values as rings that the atom is a part of.

vector<imat> atomsAndTheirRings(imat neighbors, vector<imat> ringsWithoutPairs) {

    //generic variable.
    imat toBeAdded;

    //This matrix's size is as long as there are numbers of atoms
    vector<imat> atomRings(neighbors.n_rows);

    //Loops over rings
    for (int i = 0; i < ringsWithoutPairs.size(); ++i) {

        //Loops over atoms
        for (int j = 0; j < ringsWithoutPairs[i].n_cols; ++j) {

            //to be Added is the ring number to be added
            toBeAdded << i << endr;

            imat ring = ringsWithoutPairs[i];

            //checks if there are any rings recored for the current atom.
            //If not, initilizes the value.
            if (atomRings[ring(j)].is_empty()) {

                atomRings[ring(j)] = toBeAdded;



            } else {//Else, checks to see if this neighbor is already recorded

                if (contains(atomRings[ring(j)], i)) {
                    // THIS CHECK SUPERFLOUOUS?
                    //If already recored, do nothing and continue

                    continue;
                } else {
                    //if not recorded, will record.

                    atomRings[ring(j)] = join_rows(atomRings[ring(j)], toBeAdded);
                }
            }
        }
    }
    return atomRings;
}

//Takes in three 1xN matrices. Returns true if you can completely reconstruct
//the third matrix (C) using any of the available atoms in the first two (A, B).

bool ringsContainRing(imat A, imat B, imat C) {
    bool containsRing;

    //checks that all rings exist. If not, bool is false.
    if (A.is_empty() || B.is_empty() || C.is_empty()) {
        containsRing = 0;
        return containsRing;
    }


    //loops through atoms of C
    int i;
    for (i = 0; i < C.n_cols; ++i) {

        //checks if atom exists in either A or B, if so cotinues.
        //If atom is not found, entire subfunction returns false.
        if ((contains(A, C(i))) || (contains(B, C(i)))) {
            continue;
        } else {
            containsRing = 0;
            return containsRing;
        }

    }
    // if subfunction continues to run past the for loop, insinuates that all atoms were found.
    containsRing = 1;
    return containsRing;
}


//converts the format of ringsOfAllSizes from pairs to just imats

vector<imat> stripPairs(vector<pair<imat, int> > ringsOfAllSizes) {
    vector<imat> ringsWithoutPairs;
    for (int i = 0; i < ringsOfAllSizes.size(); ++i) {
        for (int j = 0; j < ringsOfAllSizes[i].first.n_rows; ++j) {

            ringsWithoutPairs.push_back(ringsOfAllSizes[i].first.row(j));
        }
    }
    return ringsWithoutPairs;
}



//Next bit of code eliminates rings that contain other rings.
//modifies ringsOfAllSizes

void eliminateRingsInRings(vector<pair<imat, int> >& ringsOfAllSizes, imat neighborList) {
    cout << "ELIMINATING RINGS" << endl;
    bool containsRing;
    //kill List names the rings that can be constructed from other rings, to be deleted later.
    //initiation with junk value.
    imat killList;
    killList << -1 << endr;

    //We must find an atom that is a part of 4 or more rings, which is not allowed.
    // All the rings of this atom become "suspects."

    //strippedPairs are a reformatted version of rings of All Sizes.
    vector<imat> strippedPairs = stripPairs(ringsOfAllSizes);
    /*
 cout << "STRIPPED RINGS" << endl;
 for (int i = 0; i < strippedPairs.size(); ++i) {
     cout << strippedPairs[i];
     cout << "RING END " << endl;
 }
   //  */
    vector<imat> atomRings = atomsAndTheirRings(neighborList, strippedPairs);

    /*
        cout << "ATOMS AND THEIR RINGS"<<endl;
        for(int i = 0; i < atomRings.size(); ++i){
          cout << atomRings[i] << endl;
    }
    // */
    int maxNumberOfOverlappingRings;
    uvec negValueIndices;

    //loops over atoms
    for (int i = 0; i < atomRings.size(); ++i) {

        //Checks if atom is on the edge.
        negValueIndices = find(neighborList.row(i) == -1);
        if (negValueIndices.is_empty()) {
            // no negative values, atom has 4 neighbors, and can NOT be a part of any rings. (sp3 bonded atom)
            maxNumberOfOverlappingRings = 0;
        } else if (negValueIndices.n_rows == 1) {
            //sp2 atom, can be a part of at most 3 ring.
            maxNumberOfOverlappingRings = 3;
        } else if (negValueIndices.n_rows == 2) {
            //edge atom, can be a part of at most 1 rings.
            maxNumberOfOverlappingRings = 1;
        } else if (negValueIndices.n_rows == 3) {
            //dangling atom, can be a part of at most 0 rings.
            maxNumberOfOverlappingRings = 0;
        }


        // if atom is not allowed to be a part of any rings, add all rings it is a part of (if they exist)
        // to kill list

        if ((maxNumberOfOverlappingRings == 0) && (atomRings[i].is_empty() == false)) {
            for (int j = 0; j < atomRings[i].n_cols; ++j) {
                imat toBeAppended;
                toBeAppended = atomRings[i](j);
                killList = join_rows(killList, toBeAppended);
            }
        }

        //checks for suspects
        if (atomRings[i].n_cols > maxNumberOfOverlappingRings) {

            imat suspects = atomRings[i];
            //cout << "SUSPECTS! " << suspects << endl;
            //if found, tries to find if any combonation of suspects equals any other ring.
            //the following indexes all represent the name of a ring, which is simply a number
            //from 0 to number of rings (the size of strippedPairs

            //loops over the first suspect, A
            for (int j = 0; j < suspects.n_cols; ++j) {
                imat A = strippedPairs[suspects(j)];
                //loops over the second suspect, B
                for (int k = 0; k < suspects.n_cols; ++k) {
                    imat B = strippedPairs[suspects(k)];

                    //cout << "I " << i << endl;
                    //cout << " J " <<j << endl;
                    //cout << "K " << k << endl;

                    //If we are testing an edge atom, then we need to check just two rings.
                    if (maxNumberOfOverlappingRings == 1) {
                        //... and all rings must be unique.

                        if (j != k) {

                            // if three rings can be found that can construct each other,
                            // determines which is largest.
                            //cout << "A: " << A << endl;
                            //cout << "B: " << B << endl;
                            //cout << "suspect rings " << suspects(j) << " " << suspects(k) << endl;

                            if (strippedPairs[suspects(j)].n_cols > strippedPairs[suspects(k)].n_cols) {
                                //cout << "A is greater, kill" << endl;
                                imat toBeAppended;
                                toBeAppended << suspects(j);
                                killList = join_rows(killList, toBeAppended);
                                continue;

                            } else if (strippedPairs[suspects(j)].n_cols < strippedPairs[suspects(k)].n_cols) {
                                //cout << "B is greater, kill" << endl;
                                imat toBeAppended;
                                toBeAppended << suspects(k);
                                killList = join_rows(killList, toBeAppended);
                                continue;
                            } else if (strippedPairs[suspects(j)].n_cols == strippedPairs[suspects(k)].n_cols) {
                                //cout << "A and B are equal in size... shouldn't happen " << endl;
                                continue;
                            }



                        }

                    }

                    // If testing an interior atom, then we need to check three rings at a time.
                    if (maxNumberOfOverlappingRings == 3) {

                        //loops over third suscpect, C
                        for (int l = 0; l < suspects.n_cols; ++l) {

                            imat C = strippedPairs[suspects(l)];

                            //... and all rings must be unique.
                            if ((j != l) && (k != l) && (j != k)) {
                                containsRing = ringsContainRing(A, B, C);
                                // if three rings can be found that can construct each other,
                                // determines which is largest.
                                // cout << "A: " << A << endl;
                                // cout << "B: " << B << endl;
                                //cout << "C: " << C << endl;
                                //  cout << "suspect rings " << suspects(j) << " " << suspects(k) << " " << suspects(l) << endl;

                                if (containsRing) {
                                    //cout << "CONTAINS!" << endl;
                                    if ((C.n_cols > A.n_cols) && (C.n_cols > B.n_cols)) {
                                        //C Largest check to see if it is not already on the kill List.
                                        if (contains(killList, suspects(l))) {
                                            //cout << "ALREADY RECORDED" << endl;
                                            continue;
                                        } else {

                                            imat toBeAppended;
                                            toBeAppended << suspects(l);
                                            killList = join_rows(killList, toBeAppended);
                                            //  cout << "C LARGEST, KILL C" << endl;
                                            continue;
                                        }
                                    } else {
                                        //cout << "C NOT LARGEST " << endl;
                                        continue;
                                    }

                                } else {
                                    // cout << "DOES NOT CONTAIN" << endl;
                                    continue;
                                }
                            }
                        }
                    }
                }
            }
        }

    }

    //Once we have a complete Kill list, delete those rings and reconstruct ringsAndSizes.
    killList.shed_col(0); //shed initiation value.
    //only need to delete rings if there are rings on the kill list.
    if (killList.is_empty() == false) {
        //cout << "KILLLIST" << killList << endl;
        ivec killed;
        killed = -1;
        //iterate over rings to kill
        for (int i = 0; i < killList.n_cols; ++i) {
            strippedPairs[killList(i)] = killed;
        }
        //erases the killed rings
        int strippedPairSize = strippedPairs.size();
        //iterates over the rings
        int killIter = 0;
        while (killIter < strippedPairSize) {
            uvec foundNegIndexes = find(strippedPairs[killIter] == -1);
            if (foundNegIndexes.is_empty()) {
                //this ring was not killed
                ++killIter;
                continue;
            } else {
                //ring was killed, remove
                strippedPairs.erase(strippedPairs.begin() + killIter);
                --strippedPairSize;
            }
        }
    }

    /*
        cout << "STRIPPED RINGS" << endl;
        for (int i = 0; i < strippedPairs.size(); ++i) {
            cout << strippedPairs[i];
            cout << "RING END " <<endl;
        }
     //*/
    //converting stripped pairs back to pairs. Yes I know this is inefficient, need to change later.

    //initialize ringsOfAllSizesNew
    pair<imat, int > initilization(strippedPairs[0], strippedPairs[0].n_cols);
    vector<pair<imat, int > > ringsOfAllSizesNew;
    ringsOfAllSizesNew.push_back(initilization);

    // continues to add other rings to the converted data format.
    //i iterates over rings. Start at 1 and not zero do to the first initilization step

    for (int i = 1; i < strippedPairs.size(); ++i) {
        bool ringCategoryFound = false;
        //loops through our ring size catagories to see if we have one of the size
        //of the current selected ring. If so, appends it to that collection of rings.
        //if not, creates the catagory.
        //j iterates over the ring ringsize pair vector that is being built

        for (int j = 0; j < ringsOfAllSizesNew.size(); ++j) {

            //strippedPairs[i].n_cols is the ring size catagory.
            if (ringsOfAllSizesNew[j].second == strippedPairs[i].n_cols) {
                ringsOfAllSizesNew[j].first = join_cols(ringsOfAllSizesNew[j].first, strippedPairs[i]);
                ringCategoryFound = true;
                //cout << "found pre-existing ring size category " << endl;
                break;
            }
        }
        //adds new catagory if it wasn't found
        if (ringCategoryFound == false) {
            pair<imat, int> ringsAndRingSize(strippedPairs[i], strippedPairs[i].n_cols);
            ringsOfAllSizesNew.push_back(ringsAndRingSize);
            // cout << "make new ring size category " << endl;

        }

        /*             //prints out the rings and their sizes.
         for (int i = 0; i < ringsOfAllSizesNew.size(); ++i) {
             //prints out  the size of the rings, then the rings.
             cout << "Size of rings: " << ringsOfAllSizesNew[i].second << ". Number of rings: " << ringsOfAllSizesNew[i].first.n_rows << endl;
             if (ringsOfAllSizesNew[i].first.n_rows != 0) {
                 cout << ringsOfAllSizesNew[i].first << endl;
             }
         }
         //*/

    }
    //cout << "RINGS WITHIN RINGS DONE" << endl;
    ringsOfAllSizes = ringsOfAllSizesNew;
}

//Recursive algorithm for finding the first smallest possible ring given an atom.
//currently only looks for rings size max (in variable file) and below.

void findSmallestRingForAtom(int atomNum, imat ring, imat& winner, imat& neigh, int layersDeep, bool& victory, int& ringSizeMax) {

    // Iterates over all neighbors.
    for (int neighborIter = 0; neighborIter < neigh.n_cols; ++neighborIter) {

        // records depth of iteration, only triggers first time in loop.
        if (neighborIter == 0) {
            ++layersDeep;
        }

        // trims the ring back down to the layer depth, if it had additions from other prospective atoms.
        ring = ring.cols(0, layersDeep);

        //Check for loop completion, i.e. the end atom is the same as the first
        // and the ring is at least 3 large.
        if ((ring(0) == neigh(atomNum, neighborIter)) && (ring.n_cols <= ringSizeMax) && (ring.n_cols >= 3)) {

            winner = ring; //sets winning ring
            victory = true;
            return;
        }

        // breaks loop if prospective ring completing neighbors is exhausted and max ring size is reached

        if ((neighborIter == (neigh.n_cols - 1)) && (ring.n_cols == ringSizeMax)) {
            break;
        }

        // Creates a matrix filled with the next neighbor to be tested.
        // To be used in a Boolian later to check if the next neighbor already
        // exists in the ring.
        imat nextNeighbor;
        nextNeighbor.copy_size(ring);
        nextNeighbor.fill(neigh(atomNum, neighborIter));

        // Checks to make sure there is no duplicates,
        ///and neighbor exists (positive value),
        //and ring is not at the max size
        //  adds atom to ring
        if ((accu(nextNeighbor == ring) == 0)
                && (neigh(atomNum, neighborIter) >= 0)
                && (ring.n_elem < ringSizeMax)) {

            //Joins atom to ring
            imat toBeAppended;
            toBeAppended << neigh(atomNum, neighborIter) << endr;
            ring = join_rows(ring, toBeAppended);

        } else {
            continue; // if duplicates or non-existant, skips to the next neighbor
        }

        findSmallestRingForAtom(neigh(atomNum, neighborIter), ring, winner, neigh, layersDeep, victory, ringSizeMax);
        if (victory == true) {
            return;
        }
    }
}

//Assuming at max three neighbors.
//Modifies the address of "rings" instead of returning value.

void findRingsForAtom(int atomNum, imat ring, int& ringSize, imat& neigh, imat& rings, int layersDeep) {

    // Iterates over all neighbors.
    for (int neighborIter = 0; neighborIter < neigh.n_cols; ++neighborIter) {

        // records depth of iteration, only triggers first time in loop.
        if (neighborIter == 0) {
            ++layersDeep;
        }

        // sets the current atom to the end of the ring chain.
        ring = ring.cols(0, layersDeep);

        // Creates a matrix filled with the next neighbor to be tested.
        // To be used in a Boolian later to check if the next neighbor already
        // exists in the ring.
        imat nextNeighbor;
        nextNeighbor.copy_size(ring);
        nextNeighbor.fill(neigh(atomNum, neighborIter));

        //Check for loop completion, i.e. the end atom is the same as the first,
        //and the ring is the right size.
        if ((ring(0) == neigh(atomNum, neighborIter)) && (ring.n_cols == ringSize)) {
            rings = join_cols(rings, ring); //Adds the ring list.
            // cout << "ring reached! " << endl << "ring base " << ring(0) << " next atom " << neigh(atomNum, neighborIter)  << " ringSize:" << ringSize << " ring: " << endl << ring << endl;
        }

        // breaks loop if ring size reached and prospective ring completing neighbors is exhausted.
        if ((ring.n_cols == ringSize) && (neighborIter == (neigh.n_cols - 1))) {

            break;
        }

        //we must check that the next neighbor's neighbors are not already part
        // of the ring path, unless the atom in question is either at the end of the ring
        //where we just came from), or if we are just before the final ring size:
        // the beginning of the ring (loop completion).
        //This is important, as it prevents many false compound rings from being found
        //by this algorithm, speeding up the eliminatingRingsInRings subfunction greatly.

        bool noIntersection = true;
        if (nextNeighbor(0) != -1) {
            for (int nextNeighborIter = 0; nextNeighborIter < neigh.n_cols; nextNeighborIter++) {

                bool nextNeighborsNeighborIsStartingAtom = false;
                if (ring.size() == ringSize - 1
                        && neigh(nextNeighbor(0), nextNeighborIter) == ring(0)) {
                    nextNeighborsNeighborIsStartingAtom = true;
                }

                if (contains(ring, neigh(nextNeighbor(0), nextNeighborIter))
                        && neigh(nextNeighbor(0), nextNeighborIter) != ring(ring.size() - 1)
                        && nextNeighborsNeighborIsStartingAtom != true) {

                    noIntersection = false;
                }
            }
        }

        // Checks to make sure next neighbor isn't allready on path,
        ///and neighbor exists (positive value),
        //and ring isn't already full (just checking neighbors for completion),
        //and next neighbor will not intersect the ring mid-way,

        //adds atom to ring
        if ((accu(nextNeighbor == ring) != 1)
                && (neigh(atomNum, neighborIter) >= 0)
                && (ring.n_cols < ringSize)
                && noIntersection) {

            //Joins atom to ring
            imat toBeAppended;
            toBeAppended << neigh(atomNum, neighborIter) << endr;
            ring = join_rows(ring, toBeAppended);


        } else {

            continue; // if duplicates or non-existant or crossing over path again, skips to the next neighbor
        }
        //go next layer down now that we have a new atom added to the ring...
        findRingsForAtom(neigh(atomNum, neighborIter), ring, ringSize, neigh, rings, layersDeep);
    }
}

// Finds all the rings of a given ring size range for a single atom. Returns ring matrix.

vector<pair<imat, int> > findUniqueRingsForDimer(imat neighbors, int chosen, int neighbor, int& ringSizeMin, int& ringSizeMax) {

    vector<pair<imat, int> > ringsOfAllSizes;

    for (int ringSize = ringSizeMin; ringSize <= ringSizeMax; ++ringSize) {

        imat rings(1, ringSize); //ring list, initiation values will be shed later

        imat ring; // starting the ring. This matrix is not modified, the address of rings is.

        //initilize our stub of a ring, which is the chosen atom and neighbor
        ring << chosen << neighbor << endr;


        //starting at index 0 since one is added when entering the first time

        findRingsForAtom(neighbor, ring, ringSize, neighbors, rings, 0);

        //if the rings matrix contains more than the initial dummy values,
        //then rings were found and can be added to master vector.
        if (rings.n_rows != 1) {
            rings.shed_row(0); //shed initiation value.

            eliminateDuplicates(rings);

            //pairs the rings matrix with it's size, pushes onto the existing vector of pairs.
            pair<imat, int> ringsAndRingSize(rings, ringSize);
            ringsOfAllSizes.push_back(ringsAndRingSize);

        }

    }

    if (ringsOfAllSizes.empty()) {
        cerr << "selected dimer is not part of a ring" << endl;
    } else {
        eliminateRingsInRings(ringsOfAllSizes, neighbors);
    }
    return ringsOfAllSizes;
}

// Finds all the rings of a given ring size range. Returns ring matrix.

vector<pair<imat, int> > findRings(imat neighbors, int& ringSizeMin, int& ringSizeMax) {
    cout << "FINDING RINGS" << endl;




    vector<pair<imat, int> > ringsOfAllSizes;

    //Iterates over rings between sizes min and max.
    for (int ringSize = ringSizeMin; ringSize <= ringSizeMax; ++ringSize) {

        imat rings(1, ringSize); //ring list, initiation values will be shed later

        imat ring(1, 1); // starting the ring. This matrix is not modified, the address of rings is.

        //Iterates over all atoms, finding the rings for each.
        for (int i = 0; i < neighbors.n_rows; ++i) {

            ring(0) = i; // sets base atom. This is a 1xN matrix of atoms in a ring for atom i.


            // -1 must be initial depth, as 1 is added at the beginning of the loop.
            findRingsForAtom(i, ring, ringSize, neighbors, rings, -1);
        }

        rings.shed_row(0); //shed initiation value.

        eliminateDuplicates(rings);



        //pairs the rings matrix with it's size, pushes onto the existing vector of pairs.
        pair<imat, int> ringsAndRingSize(rings, ringSize);
        ringsOfAllSizes.push_back(ringsAndRingSize);
    }

    /*
        for (int i = 0; i < ringsOfAllSizes.size(); ++i) {
        //prints out  the size of the rings, then the rings.
        cout << "Size of rings: " << ringsOfAllSizes[i].second << ". Number of rings: " << ringsOfAllSizes[i].first.n_rows << endl;
        if (ringsOfAllSizes[i].first.n_rows != 0) {
            cout << ringsOfAllSizes[i].first;
        }
    }
     */

    eliminateRingsInRings(ringsOfAllSizes, neighbors);
    return ringsOfAllSizes;
}

//Recursive algorithm for finding a path from atom A to atom B.
//Any path is considered, not just the shortest.

void findPathForAtoms(int& atomNum, int B, imat path, imat& winner, imat& neigh, int layersDeep, bool& victory, int& pathSizeMax) {

    // Iterates over all neighbors.
    for (int neighborIter = 0; neighborIter < neigh.n_cols; ++neighborIter) {
        // records depth of iteration, only triggers first time in loop.
        if (neighborIter == 0) {
            ++layersDeep;
        }

        // trims the ring back down to the layer depth, if it had additions from other prospective atoms.
        path = path.cols(0, layersDeep);

        //Check for loop completion, i.e. the end atom is B
        if ((B == neigh(atomNum, neighborIter)) && (path.n_cols <= pathSizeMax)) {
            winner = path; //sets winning ring
            victory = true;
            return;
        }

        // breaks loop if prospective ring completing neighbors is exhausted and max ring size is reached
        if ((neighborIter == (neigh.n_cols - 1)) && (path.n_cols == pathSizeMax)) {
            break;
        }

        // Creates a matrix filled with the next neighbor to be tested.
        // To be used in a Boolian later to check if the next neighbor already
        // exists in the ring.
        imat nextNeighbor;
        nextNeighbor.copy_size(path);
        nextNeighbor.fill(neigh(atomNum, neighborIter));

        // Checks to make sure there is no duplicates,
        ///and neighbor exists (positive value),
        //and ring is not at the max size
        //  adds atom to ring
        if ((accu(nextNeighbor == path) == 0)
                && (neigh(atomNum, neighborIter) >= 0) && (path.n_elem < pathSizeMax)) {

            //Joins atom to ring
            imat toBeAppended;
            toBeAppended << neigh(atomNum, neighborIter) << endr;
            path = join_rows(path, toBeAppended);

        } else {
            continue; // if duplicates or non-existant, skips to the next neighbor
        }

        findPathForAtoms(neigh(atomNum, neighborIter), B, path, winner, neigh, layersDeep, victory, pathSizeMax);
        if (victory == true) {
            return;
        }
    }
}

// Using a neighbor list, generates 2 by N bond list.

imat buildBond(imat neigh) {
    imat bonds(1, 2); //bonds list, initiation values will be shed later

    for (int i = 0; i < neigh.n_rows; i++) // loops through first atoms ie rows
    {
        for (int j = 0; j < neigh.n_cols; j++) // loops through neighbor atoms of the first atom, ie columns
        {
            // Condition for bond is that an atom is not missing (positive
            // value), and that each atom sees each other in their respective
            // neighbor lists.

            if ((neigh(i, j) >= 0) &&
                    ((i == neigh(neigh(i, j), 0)) ||
                    (i == neigh(neigh(i, j), 1)) ||
                    (i == neigh(neigh(i, j), 2)) ||
                    (i == neigh(neigh(i, j), 3)))) {

                imat bond;
                bond << i << neigh(i, j) << endr; // creates bond from 2 atoms
                bonds = join_cols(bonds, bond); //adds to bond list
                neigh(i, j) = -1; // "crosses off" atom from neighbors,
                // to prevent double counting.
            }
        }
    }

    bonds.shed_row(0); // sheds junk initiation values

    return (bonds);
}

// sees if, given two neighboring atoms, that they can be rotated.
// the condition is that the atoms are sp2,
// and form an edge of a perfect graph (i.e. have a double bond between the two).
// ADD CHECK LATER TO CONFIRM THAT ATOM 0 and 1 ARE REALLY NEIGHBORS.
// ALSO ADD CHECK THAT AT LEAST ONE APPROPRIATE DOUBLE BOND EXISTS FOR ADD DIMER.
// atoms 2 and 3 are for the special case of AddDimer, where the bonds between
// 0-1 and 2-3 must be single. Otherwise 2 and 3 are ignored.
// The bool findDouble is to tell the algorithm if you want to find a double bond
// (true), or a single bond (false).

bool isDoubleBond(imat neighborList, int atom0, int atom1) {

    ListGraph graph;
    imat atomIDs;
    atomIDs << -1 << endr;

    ///first we loop over neighbors, and eliminate all sp3 carbon atoms, as
    //they do not obey the rules of perfect matching. We add nodes as we go along.
    for (int i = 0; i < neighborList.n_rows; ++i) {
        uvec negIndex = find(neighborList.row(i) == -1);
        if (negIndex.is_empty() == true) {
            //no negative index, thus 4 neighbor atoms, thus sp3, thus must be eliminated.
            neighborList(i, 0) = -1;
            neighborList(i, 1) = -1;
            neighborList(i, 2) = -1;
            neighborList(i, 3) = -1;
        } else {
            ListGraph::Node n = graph.addNode();
            //nodesVec.push_back(n);
            imat toBeAdded;
            toBeAdded << i << endr;
            atomIDs = join_cols(atomIDs, toBeAdded);
        }
    }

    atomIDs.shed_row(0);

    imat sp2Bonds = buildBond(neighborList);
    int chosenBond;
    int endBond;
    ListGraph::EdgeMap<int> edgeMap(graph);

    //now add edges, keeping in mind which is the chosen edge.
    for (int i = 0; i < sp2Bonds.n_rows; ++i) {

        uvec atomIDIndex0 = find(atomIDs == sp2Bonds(i, 0));
        uvec atomIDIndex1 = find(atomIDs == sp2Bonds(i, 1));

        graph.addEdge(graph.nodeFromId((atomIDIndex0(0))), graph.nodeFromId((atomIDIndex1(0))));


        if (((sp2Bonds(i, 0) == atom0) && (sp2Bonds(i, 1) == atom1)) || (((sp2Bonds(i, 0) == atom1) && (sp2Bonds(i, 1) == atom0)))) {
            chosenBond = i;
            // edgeMap[graph.edgeFromId(i)] = 1;
        }
        //else{
        //edgeMap[graph.edgeFromId(i)] = 0;
        //}

        //if ((findDouble == false) && (((sp2Bonds(i, 0) == atom2) && (sp2Bonds(i, 1) == atom3)) || (((sp2Bonds(i, 0) == atom3) && (sp2Bonds(i, 1) == atom2))))) {
        //    endBond = i;
        // }
    }

    //if the chosen bond wasn't found, then returns false.
    if (chosenBond < 0) {
        ///cout << "COULDNT FIND DOUBLE BONDS" << endl;
        return false;
    }


    //Attempt to seed resonance structure, depending if you are looking for a single or double bond.
    // if (findDouble == true) {
    edgeMap[graph.edgeFromId(chosenBond)] = 1;
    //} else {
    //   edgeMap[graph.edgeFromId(chosenBond)] = 0;
    // }

    //edgeMap[graph.edgeFromId(1)] = 1;

    //Attempt perfect matching
    MaxWeightedPerfectMatching<ListGraph> mat_test(graph, edgeMap);

    mat_test.run();


    //Quickly test if we have more than one double bond (meaning we are asssuming the system has 2 connected atoms at least...
    //If not, then the system has no resonance.
    bool hasDouble = false;
    for (int i = 0; i < sp2Bonds.n_rows; i++) {
        if (mat_test.matching(graph.edgeFromId(i)) == 1) {
            hasDouble = true;
        }
    }
    if (hasDouble == false) {
        system("read  -p \"ERROR, CANNOT FIND DOUBLE BOND ARRANGEMENT SATISFYING SP2 CONDITIONS. ILLEGAL STRUCTURE\" key");
    }

    //  for (int i = 0; i < sp2Bonds.n_rows; ++i) {
    //   cout << "EDGE VALUE: " << graph.id(graph.edgeFromId(i)) << " Bond: " << sp2Bonds(i, 0) << " " << sp2Bonds(i, 1) << " NODES " << graph.id(graph.nodeFromId(i*2)) << " " << graph.id(graph.nodeFromId(i*2+1))<<" MATCHING " << mat_test.matching(graph.edgeFromId(i)) << endl;
    //  }


    //   if (findDouble == true) {
    return mat_test.matching(graph.edgeFromId(chosenBond));
    // } else {
    //     return ( (mat_test.matching(graph.edgeFromId(chosenBond)) == 0) && (mat_test.matching(graph.edgeFromId(endBond)) == 0));
    // }
    // */ return true;
}

//This algorithm tells if the bond between chosen and other is NOT a double.
//If fed non-negative ints for beforeEnd and afterEnd, then adds the aditional condition that
//the bond between those two must also be single.
//You cannot directly seed the perfect matching algorithm with a NOT bond,
//so you must seed the bonds of the NN and NNNs and see if one of the four configurations will be accepted.

bool isNotDoubleBond(imat neighborList, int chosen, int neighbor, int NN1, int NN2, int NNN1, int NNN2, int beforeEnd, int afterEnd) {

    ListGraph graph;
    imat atomIDs;
    //junk initation value
    atomIDs << -1 << endr;


    ///first we loop over neighbors, and eliminate all sp3 carbon atoms, as
    //they do not obey the rules of perfect matching. We add nodes as we go along.
    for (int i = 0; i < neighborList.n_rows; ++i) {
        uvec negIndex = find(neighborList.row(i) == -1);
        if (negIndex.is_empty() == true) {
            //no negative index, thus 4 neighbor atoms, thus sp3, thus must be eliminated.
            neighborList(i, 0) = -1;
            neighborList(i, 1) = -1;
            neighborList(i, 2) = -1;
            neighborList(i, 3) = -1;
        } else {
            ListGraph::Node n = graph.addNode();
            imat toBeAdded;
            toBeAdded << i << endr;
            atomIDs = join_cols(atomIDs, toBeAdded);
        }
    }

    //shed junk initiation value
    atomIDs.shed_row(0);

    imat sp2Bonds = buildBond(neighborList);

    int chosenNeighborBond;
    int chosenNN1Bond;
    int chosenNN2Bond;
    int neighborNNN1Bond;
    int neighborNNN2Bond;
    int endBond = -1;

    ListGraph::EdgeMap<int> edgeMap(graph);

    //now add edges, keeping in mind the important player edges.
    for (int i = 0; i < sp2Bonds.n_rows; ++i) {

        uvec atomIDIndex0 = find(atomIDs == sp2Bonds(i, 0));
        uvec atomIDIndex1 = find(atomIDs == sp2Bonds(i, 1));

        graph.addEdge(graph.nodeFromId(atomIDIndex0(0)), graph.nodeFromId(atomIDIndex1(0)));

        if (((sp2Bonds(i, 0) == chosen) && (sp2Bonds(i, 1) == neighbor)) || (((sp2Bonds(i, 0) == neighbor) && (sp2Bonds(i, 1) == chosen)))) {
            chosenNeighborBond = i;
        }
        if (((sp2Bonds(i, 0) == chosen) && (sp2Bonds(i, 1) == NN1)) || (((sp2Bonds(i, 0) == NN1) && (sp2Bonds(i, 1) == chosen)))) {
            chosenNN1Bond = i;
        }
        if (((sp2Bonds(i, 0) == chosen) && (sp2Bonds(i, 1) == NN2)) || (((sp2Bonds(i, 0) == NN2) && (sp2Bonds(i, 1) == chosen)))) {
            chosenNN2Bond = i;
        }
        if (((sp2Bonds(i, 0) == neighbor) && (sp2Bonds(i, 1) == NNN1)) || (((sp2Bonds(i, 0) == NNN1) && (sp2Bonds(i, 1) == neighbor)))) {
            neighborNNN1Bond = i;
        }
        if (((sp2Bonds(i, 0) == neighbor) && (sp2Bonds(i, 1) == NNN2)) || (((sp2Bonds(i, 0) == NNN2) && (sp2Bonds(i, 1) == neighbor)))) {
            neighborNNN2Bond = i;
        }
        if (((sp2Bonds(i, 0) == afterEnd) && (sp2Bonds(i, 1) == beforeEnd)) || (((sp2Bonds(i, 0) == beforeEnd) && (sp2Bonds(i, 1) == afterEnd)))) {
            endBond = i;
        }
    }

    // We now must go through and test 4 times, to see if we can get a single bond from both
    // the chosen-neighbor and the beforeEnd-afterEnd

    //start with the bond between chosen and NN1 and setting it to a double...
    edgeMap[graph.edgeFromId(chosenNN1Bond)] = 1;
    MaxWeightedPerfectMatching<ListGraph> mat_test0(graph, edgeMap);
    mat_test0.run();
    // sees if the chosen bond is still single...
    if (mat_test0.matching(graph.edgeFromId(chosenNeighborBond)) == 0) {
        //check if the operation is add dimer or not (if the afterEnd and beforeEnd are positive or not)
        if ((afterEnd >= 0) && (beforeEnd >= 0)) {
            //if adding dimer, we have the additional constraint that the afterEnd-beforeEnd has to be single too.
            if (mat_test0.matching(graph.edgeFromId(endBond)) == 0) {
                //we add dimer, and also the endBond is single too.
                return true;
            } else {
                //if failed, then blank the seed from edge graph.
                edgeMap[graph.edgeFromId(chosenNN1Bond)] = 0;
            }

        } else {
            //if no add dimer, we are complete and successful.
            return true;
        }
    } else {
        //if failed, then blank the seed from edge graph.
        edgeMap[graph.edgeFromId(chosenNN1Bond)] = 0;
    }

    //now for NN2..
    edgeMap[graph.edgeFromId(chosenNN2Bond)] = 1;
    MaxWeightedPerfectMatching<ListGraph> mat_test1(graph, edgeMap);
    mat_test1.run();
    // sees if the chosen bond is still single...
    if (mat_test1.matching(graph.edgeFromId(chosenNeighborBond)) == 0) {
        //check if the operation is add dimer or not (if the afterEnd and beforeEnd are positive or not)
        if ((afterEnd >= 0) && (beforeEnd >= 0)) {
            //if adding dimer, we have the additional constraint that the afterEnd-beforeEnd has to be single too.
            if (mat_test1.matching(graph.edgeFromId(endBond)) == 0) {
                //we add dimer, and also the endBond is single too.
                return true;
            } else {
                //if failed, then blank the seed from edge graph.
                edgeMap[graph.edgeFromId(chosenNN2Bond)] = 0;
            }
        } else {
            //if no add dimer, we are complete and successful.
            return true;
        }
    } else {
        //if failed, then blank the seed from edge graph.
        edgeMap[graph.edgeFromId(chosenNN2Bond)] = 0;
    }

    //NNN1
    edgeMap[graph.edgeFromId(neighborNNN1Bond)] = 1;
    MaxWeightedPerfectMatching<ListGraph> mat_test2(graph, edgeMap);
    mat_test2.run();
    // sees if the chosen bond is still single...
    if (mat_test2.matching(graph.edgeFromId(chosenNeighborBond)) == 0) {
        //check if the operation is add dimer or not (if the afterEnd and beforeEnd are positive or not)
        if ((afterEnd >= 0) && (beforeEnd >= 0)) {
            //if adding dimer, we have the additional constraint that the afterEnd-beforeEnd has to be single too.
            if (mat_test2.matching(graph.edgeFromId(endBond)) == 0) {
                //we add dimer, and also the endBond is single too.
                return true;
            } else {
                //if failed, then blank the seed from edge graph.
                edgeMap[graph.edgeFromId(neighborNNN1Bond)] = 0;
            }
        } else {
            //if no add dimer, we are complete and successful.
            return true;
        }
    } else {
        //if failed, then blank the seed from edge graph.
        edgeMap[graph.edgeFromId(neighborNNN1Bond)] = 0;
    }

    //now for NNN2..
    edgeMap[graph.edgeFromId(neighborNNN2Bond)] = 1;
    MaxWeightedPerfectMatching<ListGraph> mat_test3(graph, edgeMap);
    mat_test3.run();
    // sees if the chosen bond is still single...
    if (mat_test3.matching(graph.edgeFromId(chosenNeighborBond)) == 0) {
        //check if the operation is add dimer or not (if the afterEnd and beforeEnd are positive or not)
        if ((afterEnd >= 0) && (beforeEnd >= 0)) {
            //if adding dimer, we have the additional constraint that the afterEnd-beforeEnd has to be single too.
            if (mat_test3.matching(graph.edgeFromId(endBond)) == 0) {
                //we add dimer, and also the endBond is single too.
                return true;
            } else {
                //if failed, then blank the seed from edge graph.
                edgeMap[graph.edgeFromId(neighborNNN2Bond)] = 0;
            }
        } else {
            //if no add dimer, we are complete and successful.
            return true;
        }
        //OPTIMIZE LATER SO THAT IF THE REASON FOR FAILURE WAS THAT THE CHOSEN BOND WAS NEVER SINGLE, THEN ADD DIMER FAILS WITHOUT CHECKING OTHER ENDS
    } else {
        //if failed, then blank the seed from edge graph.
        edgeMap[graph.edgeFromId(neighborNNN2Bond)] = 0;
    }

    /*
     for (int i = 0; i < sp2Bonds.n_rows; ++i) {
        cout << "EDGE VALUE: " << graph.id(graph.edgeFromId(i)) << " Bond: " << sp2Bonds(i, 0) << " " << sp2Bonds(i, 1) << " NODES " << graph.id(graph.nodeFromId(i*2)) << " " << graph.id(graph.nodeFromId(i*2+1))<<" MATCHING " << mat_test.matching(graph.edgeFromId(i)) << endl;
     }
 // */

    //after exhaustively testing all 4 possible bond configurations, and having them all fail, we cannot find the single bonds we want.
    return false;
}

