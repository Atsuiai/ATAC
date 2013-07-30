/* 
 * File:   GraphTheoryATAC.h
 * Author: zacpavilion
 *
 * Created on July 28, 2013, 9:25 PM
 */

#ifndef GRAPHTHEORYATAC_H
#define	GRAPHTHEORYATAC_H

#include "armadillo"
#include <lemon/list_graph.h>//lemon
#include <lemon/matching.h>//lemon
#include <lemon/bfs.h>//lemon
#include <lemon/concepts/graph.h>//lemon

using namespace std;
using namespace arma;
using namespace lemon;


// This subfunction takes in a 1 by N vector of doubles, and returns a vector
// of the same size with the low to high index of numbers.

vec betterIndexSort(vec toBeSorted);

// Sorts rows based on values in the column colNum, low to high.

mat betterSortRows(mat matrix, int colNum);




void buildNeighFromBonds(imat bonds, imat& neighToBeRelplaced);


// This subfunction returns true if the vector contains at least one instance of the int.
// toBeSearched is a Nx1 or 1xN matrix.
bool contains(imat toBeSearched, int tryToFind);

//bit of code to eliminate duplicate rings in in list
//works for arbitrarly sized rings.
void eliminateDuplicates(imat& rings);
// generates a matrix with the row as atom number, and the values as rings that the atom is a part of.

vector<imat> atomsAndTheirRings(imat neighbors, vector<imat> ringsWithoutPairs);

//Takes in three 1xN matrices. Returns true if you can completely reconstruct
//the third matrix (C) using any of the available atoms in the first two (A, B).
bool ringsContainRing(imat A, imat B, imat C);

//converts the format of ringsOfAllSizes from pairs to just imats
vector<imat> stripPairs(vector<pair<imat, int> > ringsOfAllSizes);

//Next bit of code eliminates rings that contain other rings.
//modifies ringsOfAllSizes
void eliminateRingsInRings(vector<pair<imat, int> >& ringsOfAllSizes, imat neighborList);
//Recursive algorithm for finding the first smallest possible ring given an atom.
//currently only looks for rings size max (in variable file) and below.
void findSmallestRingForAtom(int atomNum, imat ring, imat& winner, imat& neigh, int layersDeep, bool& victory, int& ringSizeMax);

//Assuming at max three neighbors.
//Modifies the address of "rings" instead of returning value.
void findRingsForAtom(int atomNum, imat ring, int& ringSize, imat& neigh, imat& rings, int layersDeep);
  
// Finds all the rings of a given ring size range for a single atom. Returns ring matrix.

vector<pair<imat, int> > findUniqueRingsForDimer(imat neighbors, int chosen, int neighbor, int& ringSizeMin, int& ringSizeMax);

// Finds all the rings of a given ring size range. Returns ring matrix.

vector<pair<imat, int> > findRings(imat neighbors, int& ringSizeMin, int& ringSizeMax);

//Recursive algorithm for finding a path from atom A to atom B.
//Any path is considered, not just the shortest.

void findPathForAtoms(int& atomNum, int B, imat path, imat& winner, imat& neigh, int layersDeep, bool& victory, int& pathSizeMax);

// Using a neighbor list, generates 2 by N bond list.

imat buildBond(imat neigh);
// sees if, given two neighboring atoms, that they can be rotated.
// the condition is that the atoms are sp2,
// and form an edge of a perfect graph (i.e. have a double bond between the two).
// ADD CHECK LATER TO CONFIRM THAT ATOM 0 and 1 ARE REALLY NEIGHBORS.
// ALSO ADD CHECK THAT AT LEAST ONE APPROPRIATE DOUBLE BOND EXISTS FOR ADD DIMER.
// atoms 2 and 3 are for the special case of AddDimer, where the bonds between
// 0-1 and 2-3 must be single. Otherwise 2 and 3 are ignored.
// The bool findDouble is to tell the algorithm if you want to find a double bond
// (true), or a single bond (false).

bool isDoubleBond(imat neighborList, int atom0, int atom1);

//This algorithm tells if the bond between chosen and other is NOT a double.
//If fed non-negative ints for beforeEnd and afterEnd, then adds the aditional condition that
//the bond between those two must also be single.
//You cannot directly seed the perfect matching algorithm with a NOT bond,
//so you must seed the bonds of the NN and NNNs and see if one of the four configurations will be accepted.

bool isNotDoubleBond(imat neighborList, int chosen, int neighbor, int NN1, int NN2, int NNN1, int NNN2, int beforeEnd, int afterEnd);


#endif	/* GRAPHTHEORYATAC_H */

