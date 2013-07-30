
#include "GeometryFunctionsATAC.h"


//this function takes in positions, and deletes all of the ones that are outside a certain cylindrical radius.
//this is so we can seperate the inner and outer tubes of a DWNT

void tubeShaver(mat& pos) {
    //find central axis:

    vec zCol = pos.col(2);
    vec yCol = pos.col(1);
    double zAvg = (zCol.max() + zCol.min()) / 2;
    double yAvg = (yCol.max() + yCol.min()) / 2;
    // cout << "Y MAX MIN, Z MAX MIN " << yCol.min() << " " << yCol.max() << " " << zCol.min() << " " << zCol.max() << endl;
    // cout << "y average, z average " << yAvg << " " << zAvg << endl;

    double shaverCutOff = 6; //radius, in angstroms

    for (int i = 0; i < pos.n_rows; i++) {
        double y2Pos = (pos(i, 1) - yAvg) * (pos(i, 1) - yAvg);
        double z2Pos = (pos(i, 2) - zAvg) * (pos(i, 2) - zAvg);
        if (y2Pos + z2Pos > shaverCutOff * shaverCutOff) {
            pos.shed_row(i);
            i--;
        }
    }
}

//This function fills out 8 arrays with the atoms of 4 electrodes (8 terminals). xp1 are the
//electrodes in the +x extreme, xp2 are the electrodes in the terminal next to xp1,
// xn1 is the -x extreme, xn2 in the terminal next to xn1, etc.

void electrodeDeterminer(mat pos, imat& xp1Electrodes, imat& xn1Electrodes, imat& yp1Electrodes, imat& yn1Electrodes, imat& xp2Electrodes, imat& xn2Electrodes, imat& yp2Electrodes, imat& yn2Electrodes, int& numberOfElectrodes, double& xpElecWallMod, double& xnElecWallMod, double& ypElecWallMod, double& ynElecWallMod) {
    //find our cutoffs for what is and and isn't an electrode.
    double xp1CutOff;
    double xn1CutOff;
    double xp2CutOff;
    double xn2CutOff;
    if (numberOfElectrodes >= 2) {
        mat xValues = pos.col(0);
        xp1CutOff = xValues.max() + xpElecWallMod;
        xp2CutOff = xValues.max() + xpElecWallMod * 2;
        xn1CutOff = xValues.min() + xnElecWallMod;
        xn2CutOff = xValues.min() + xnElecWallMod * 2;
    }

    double yp1CutOff;
    double yn1CutOff;
    double yp2CutOff;
    double yn2CutOff;
    if (numberOfElectrodes == 4) {
        mat yValues = pos.col(1);
        yp1CutOff = yValues.max() + ypElecWallMod;
        yp2CutOff = yValues.max() + ypElecWallMod * 2;
        yn1CutOff = yValues.min() + ynElecWallMod;
        yn2CutOff = yValues.min() + ynElecWallMod * 2;
    }

    //cout << "CUTOFFS XP XN YP YN " << xpCutOff << " " << xnCutOff << " " << ypCutOff << " " << ynCutOff << endl;
    ///system("read  -p \"Number of electrodes changed! Push enter to unpause.\" key");


    //initialize matrixes with junk values.
    imat newxp1;
    newxp1 << -1 << endr;
    imat newxn1;
    newxn1 << -1 << endr;
    imat newyp1;
    newyp1 << -1 << endr;
    imat newyn1;
    newyn1 << -1 << endr;
    imat newxp2;
    newxp2 << -1 << endr;
    imat newxn2;
    newxn2 << -1 << endr;
    imat newyp2;
    newyp2 << -1 << endr;
    imat newyn2;
    newyn2 << -1 << endr;

    // loop over atoms to see if they can be electrodes.
    for (int i = 0; i < pos.n_rows; i++) {

        if (numberOfElectrodes >= 2) {
            if (pos(i, 0) > xp1CutOff) {
                imat toBeAdded;
                toBeAdded << i << endr;
                newxp1 = join_cols(newxp1, toBeAdded);
            }
            if (pos(i, 0) > xp2CutOff && pos(i, 0) < xp1CutOff) {
                imat toBeAdded;
                toBeAdded << i << endr;
                newxp2 = join_cols(newxp2, toBeAdded);
            }
            if (pos(i, 0) < xn1CutOff) {
                imat toBeAdded;
                toBeAdded << i << endr;
                newxn1 = join_cols(newxn1, toBeAdded);
            }
            if (pos(i, 0) < xn2CutOff && pos(i, 0) > xn1CutOff) {
                imat toBeAdded;
                toBeAdded << i << endr;
                newxn2 = join_cols(newxn2, toBeAdded);
            }
        }
        if (numberOfElectrodes == 4) {
            if (pos(i, 1) > yp1CutOff) {
                imat toBeAdded;
                toBeAdded << i << endr;
                newyp1 = join_cols(newyp1, toBeAdded);
                //cout << "YP1";
            }
            if (pos(i, 1) > yp2CutOff && pos(i, 1) < yp1CutOff) {
                imat toBeAdded;
                toBeAdded << i << endr;
                newyp2 = join_cols(newyp2, toBeAdded);
                //              cout << "YP2";
            }
            if (pos(i, 1) < yn1CutOff) {
                imat toBeAdded;
                toBeAdded << i << endr;
                newyn1 = join_cols(newyn1, toBeAdded);
                //            cout << "YN1";
            }
            if (pos(i, 1) < yn2CutOff && pos(i, 1) > yn1CutOff) {
                imat toBeAdded;
                toBeAdded << i << endr;
                newyn2 = join_cols(newyn2, toBeAdded);
                //          cout << "YN2";
            }
        }
    }
    //shave junk values.
    newxp1.shed_row(0);
    newxn1.shed_row(0);
    newyp1.shed_row(0);
    newyn1.shed_row(0);
    newxp2.shed_row(0);
    newxn2.shed_row(0);
    newyp2.shed_row(0);
    newyn2.shed_row(0);

    //We need to compare our newly derived electrode atoms with the old ones;
    //they should not have changed.

    //check if this is the first run, i.e. the arrays are null
    if (xp1Electrodes.is_empty()) {
        //simple initilization.
        if (numberOfElectrodes >= 2) {
            xp1Electrodes = newxp1;
            xn1Electrodes = newxn1;
            xp2Electrodes = newxp2;
            xn2Electrodes = newxn2;
        }
        if (numberOfElectrodes == 4) {
            yp1Electrodes = newyp1;
            yn1Electrodes = newyn1;
            yp2Electrodes = newyp2;
            yn2Electrodes = newyn2;
        }
    } else {
        //check that there are the same number of new and the old electrode atoms
        if (numberOfElectrodes >= 2) {
            if ((newxp1.size() != xp1Electrodes.size()) || (newxn1.size() != xn1Electrodes.size())
                    ) {
                //cout << "NUMBER OF X ELECTRODES CHANGED!!!!!!!!!!!!!!!!!! " << endl;
                system("read  -p \"Number x1 of electrodes changed! Push enter to unpause.\" key");
            }
            if ((newxp2.size() != xp2Electrodes.size()) || (newxn2.size() != xn2Electrodes.size())
                    ) {
                //cout << "NUMBER OF X ELECTRODES CHANGED!!!!!!!!!!!!!!!!!! " << endl;
                system("read  -p \"Number x2 of electrodes changed! Push enter to unpause.\" key");
            }

        }
        if (numberOfElectrodes == 4) {
            if ((newyp1.size() != yp1Electrodes.size()) || (newyn1.size() != yn1Electrodes.size())) {
                //cout << "NUMBER OF Y ELECTRODES CHANGED!!!!!!!!!!!!!!!!!! " << endl;
                system("read  -p \"Number of y1 electrodes changed! Push enter to unpause.\" key");
            }
            if ((newyp2.size() != yp2Electrodes.size()) || (newyn2.size() != yn2Electrodes.size())) {
                //cout << "NUMBER OF Y ELECTRODES CHANGED!!!!!!!!!!!!!!!!!! " << endl;
                system("read  -p \"Number of y2 electrodes changed! Push enter to unpause.\" key");
            }
        }
        //number of electrodes the same, fairly ok to assume they did not change.
        if (numberOfElectrodes >= 2) {
            xp1Electrodes = newxp1;
            xp2Electrodes = newxp2;
            xn1Electrodes = newxn1;
            xn2Electrodes = newxn2;
        }
        if (numberOfElectrodes == 4) {
            yp1Electrodes = newyp1;
            yp2Electrodes = newyp2;
            yn1Electrodes = newyn1;
            yn2Electrodes = newyn2;
        }
    }


}


//Determines the x, y and z lattice parameters (these values fed to it are blank,
//to be assigned in the algorithm). Parameters can either be be statically given,
//or can be relative to the system given.

void latticeParameterDeterminer(mat pos, bool& isStaticPeriod, int& periodicityNum, double& xLatticeParameter, double& xPeriodFactor, double& yLatticeParameter, double& yPeriodFactor, double& zLatticeParameter, double& zPeriodFactor) {

    mat xValues;
    mat yValues;
    mat zValues;


    //These period factors are context dependent. With a static period, these
    //factors will be the lattice parameters. Without a static period, these
    //factors will be added to the difference of atomic position axis extremes
    //for the lattice parameters.

    if (isStaticPeriod == true) {
        xLatticeParameter = xPeriodFactor;

        if (periodicityNum == 1) {
            xLatticeParameter = xPeriodFactor;
        }

        if (periodicityNum == 2) {
            xLatticeParameter = xPeriodFactor;
            yLatticeParameter = yPeriodFactor;
        }

        if (periodicityNum == 3) {
            xLatticeParameter = xPeriodFactor;
            yLatticeParameter = yPeriodFactor;
            zLatticeParameter = zPeriodFactor;
        }

    } else {

        if (periodicityNum == 1) {
            xValues = pos.col(0);
            xLatticeParameter = abs(xValues.max() - xValues.min()) + xPeriodFactor;
        }

        if (periodicityNum == 2) {
            xValues = pos.col(0);
            yValues = pos.col(1);

            xLatticeParameter = abs(xValues.max() - xValues.min()) + xPeriodFactor;
            yLatticeParameter = abs(yValues.max() - yValues.min()) + yPeriodFactor;

        }
        if (periodicityNum == 3) {
            xValues = pos.col(0);
            yValues = pos.col(1);
            zValues = pos.col(2);
            xLatticeParameter = abs(xValues.max() - xValues.min()) + xPeriodFactor;
            yLatticeParameter = abs(yValues.max() - yValues.min()) + yPeriodFactor;
            zLatticeParameter = abs(zValues.max() - zValues.min()) + zPeriodFactor;
        }
    }
}


//determines the shortest path in a lattice cell between two points.

void shortestLineDeterminer(mat pos, int& i, int& j, double& distance, int& periodicityNum, double& xLatticeParameter, double& yLatticeParameter, double& zLatticeParameter, vec latticeVec1Norm, vec latticeVec1, vec latticeVec2Norm, vec latticeVec2, vec latticeVec3Norm, vec latticeVec3) {

    //difference in coordinates between the chosen atom and a prospective
    //neighbor test atom.
    double xDiff = (pos(i, 0) - pos(j, 0));
    double yDiff = (pos(i, 1) - pos(j, 1));
    double zDiff = (pos(i, 2) - pos(j, 2));


    //Follow check sees if the direct distance between two atoms should be used,
    //or if wrapping around via lattice paramater is shorter.
    ///*
    if (periodicityNum == 1) {

        while (xDiff > xLatticeParameter / 2) {
            xDiff = xDiff - xLatticeParameter;
        }

        while (xDiff<-xLatticeParameter / 2) {
            xDiff = xDiff + xLatticeParameter;
        }
    }

    if (periodicityNum == 2) {

        while (xDiff > xLatticeParameter / 2) {
            xDiff = xDiff - xLatticeParameter;
        }

        while (xDiff<-xLatticeParameter / 2) {
            xDiff = xDiff + xLatticeParameter;
        }

        while (yDiff > yLatticeParameter / 2) {
            yDiff = yDiff - yLatticeParameter;
        }

        while (yDiff<-yLatticeParameter / 2) {
            yDiff = yDiff + yLatticeParameter;
        }

    }

    if (periodicityNum == 3) {

        while (xDiff > xLatticeParameter / 2) {
            xDiff = xDiff - xLatticeParameter;
        }

        while (xDiff<-xLatticeParameter / 2) {
            xDiff = xDiff + xLatticeParameter;
        }

        while (yDiff > yLatticeParameter / 2) {
            yDiff = yDiff - yLatticeParameter;
        }

        while (yDiff<-yLatticeParameter / 2) {
            yDiff = yDiff + yLatticeParameter;
        }

        while (zDiff > zLatticeParameter / 2) {
            zDiff = zDiff - zLatticeParameter;
        }

        while (zDiff<-zLatticeParameter / 2) {
            zDiff = zDiff + zLatticeParameter;
        }

    }

    //*/

    /*FOR 3D periodic conditions
    if ((pos(i, 2) < pos(j, 2) &&
            (abs(pos(i, 2)-(pos(j, 2) - zLatticeParameter)) <
            abs(pos(i, 2) - pos(j, 2))))) {
        zDiff = pos(i, 2)-(pos(j, 2) - zLatticeParameter);
    }
    if ((pos(i, 2) > pos(j, 2) &&
            (abs(pos(i, 2)-(pos(j, 2) + zLatticeParameter)) <
            abs(pos(i, 2) - pos(j, 2))))) {
        zDiff = pos(i, 2)-(pos(j, 2) + zLatticeParameter);
    }
     */

    if (periodicityNum == 4) {

        //the difference in distances in the direction of the first lattice vector.
        double Vec1Diff1 = (pos(i, 0) - pos(j, 0)) * latticeVec1Norm(0);
        double Vec1Diff2 = (pos(i, 1) - pos(j, 1)) * latticeVec1Norm(1);
        double Vec1Diff3 = (pos(i, 2) - pos(j, 2)) * latticeVec1Norm(2);

        //cout << "I< J " << i << " " << j<< endl;

        double magVec1Diff = sqrt(
                pow(Vec1Diff1, 2) +
                pow(Vec1Diff2, 2) +
                pow(Vec1Diff3, 2));
        double magVec1Lattice = sqrt(
                pow(latticeVec1(0), 2) +
                pow(latticeVec1(1), 2) +
                pow(latticeVec1(2), 2));


        //if we are not within half a lattice vector length, then it is possible to translate some more.
        //while (magVec1Diff > magVec1Lattice / 2) {


        bool transFlag = false;
        //need to find out if we are adding or subtracting the lattice vector.
        //The one that reduces the magVec1Diff will be used.

        double AddVec1Diff1 = Vec1Diff1 + latticeVec1(0);
        double AddVec1Diff2 = Vec1Diff2 + latticeVec1(1);
        double AddVec1Diff3 = Vec1Diff3 + latticeVec1(2);
        double AddMagVec1Diff = sqrt(
                pow(AddVec1Diff1, 2) +
                pow(AddVec1Diff2, 2) +
                pow(AddVec1Diff3, 2));

        double SubVec1Diff1 = Vec1Diff1 - latticeVec1(0);
        double SubVec1Diff2 = Vec1Diff2 - latticeVec1(1);
        double SubVec1Diff3 = Vec1Diff3 - latticeVec1(2);
        double SubMagVec1Diff = sqrt(
                pow(SubVec1Diff1, 2) +
                pow(SubVec1Diff2, 2) +
                pow(SubVec1Diff3, 2));


        //           cout <<"magVec1Lattice / 2 " <<magVec1Lattice / 2 << endl;
        //     cout <<"magVec1Diff " << magVec1Diff << endl;
        //       cout <<"SubmagVec1Diff " << SubMagVec1Diff << endl;
        //         cout <<"AddmagVec1Diff " << AddMagVec1Diff << endl;
        //                   system("read  -p \"Paused.\" key");


        //choose if we are adding or subtracting the vector based on which is smaller in an absolute way.
        if (SubMagVec1Diff < magVec1Diff && SubMagVec1Diff < AddMagVec1Diff) {
            Vec1Diff1 = Vec1Diff1 - latticeVec1(0);
            Vec1Diff2 = Vec1Diff2 - latticeVec1(1);
            Vec1Diff3 = Vec1Diff3 - latticeVec1(2);
            magVec1Diff = sqrt(
                    pow(Vec1Diff1, 2) +
                    pow(Vec1Diff2, 2) +
                    pow(Vec1Diff3, 2));


            transFlag = true;
        } else if (AddMagVec1Diff < magVec1Diff && AddMagVec1Diff < SubMagVec1Diff) {
            Vec1Diff1 = Vec1Diff1 + latticeVec1(0);
            Vec1Diff2 = Vec1Diff2 + latticeVec1(1);
            Vec1Diff3 = Vec1Diff3 + latticeVec1(2);
            magVec1Diff = sqrt(
                    pow(Vec1Diff1, 2) +
                    pow(Vec1Diff2, 2) +
                    pow(Vec1Diff3, 2));

            transFlag = true;
        } //if neither of these conditions are met, then magVec1Diff is already teh smallest Diff.
        // }

        //not really the x difference, just eqivalently used for the last line of this subfunction.
        xDiff = magVec1Diff;


        //the difference in distances in the direction of the SECOND lattice vector.
        double Vec2Diff1 = (pos(i, 0) - pos(j, 0)) * latticeVec2Norm(0);
        double Vec2Diff2 = (pos(i, 1) - pos(j, 1)) * latticeVec2Norm(1);
        double Vec2Diff3 = (pos(i, 2) - pos(j, 2)) * latticeVec2Norm(2);
        double magVec2Diff = sqrt(
                pow(Vec2Diff1, 2) +
                pow(Vec2Diff2, 2) +
                pow(Vec2Diff3, 2));
        double magVec2Lattice = sqrt(
                pow(latticeVec2(0), 2) +
                pow(latticeVec2(1), 2) +
                pow(latticeVec2(2), 2));

        //if we are not within half a lattice vector length, then it is possible to translate some more.
        //while (magVec2Diff > magVec2Lattice / 2) {


        //need to find out if we are adding or subtracting the lattice vector.
        //The one that reduces the magVec2Diff will be used.

        double AddVec2Diff1 = Vec2Diff1 + latticeVec2(0);
        double AddVec2Diff2 = Vec2Diff2 + latticeVec2(1);
        double AddVec2Diff3 = Vec2Diff3 + latticeVec2(2);
        double AddMagVec2Diff = sqrt(
                pow(AddVec2Diff1, 2) +
                pow(AddVec2Diff2, 2) +
                pow(AddVec2Diff3, 2));

        double SubVec2Diff1 = Vec2Diff1 - latticeVec2(0);
        double SubVec2Diff2 = Vec2Diff2 - latticeVec2(1);
        double SubVec2Diff3 = Vec2Diff3 - latticeVec2(2);
        double SubMagVec2Diff = sqrt(
                pow(SubVec2Diff1, 2) +
                pow(SubVec2Diff2, 2) +
                pow(SubVec2Diff3, 2));




        //choose if we are adding or subtracting the vector based on which is smaller in an absolute way.
        //choose if we are adding or subtracting the vector based on which is smaller in an absolute way.
        if (SubMagVec2Diff < magVec2Diff && SubMagVec2Diff < AddMagVec2Diff) {
            Vec2Diff1 = Vec2Diff1 - latticeVec2(0);
            Vec2Diff2 = Vec2Diff2 - latticeVec2(1);
            Vec2Diff3 = Vec2Diff3 - latticeVec2(2);
            magVec2Diff = sqrt(
                    pow(Vec2Diff1, 2) +
                    pow(Vec2Diff2, 2) +
                    pow(Vec2Diff3, 2));

            transFlag = true;

        } else if (AddMagVec2Diff < magVec2Diff && AddMagVec2Diff < SubMagVec2Diff) {
            Vec2Diff1 = Vec2Diff1 + latticeVec2(0);
            Vec2Diff2 = Vec2Diff2 + latticeVec2(1);
            Vec2Diff3 = Vec2Diff3 + latticeVec2(2);
            magVec2Diff = sqrt(
                    pow(Vec2Diff1, 2) +
                    pow(Vec2Diff2, 2) +
                    pow(Vec2Diff3, 2));

            transFlag = true;
        }
        // }

        yDiff = magVec2Diff;

        //the difference in distances in the direction of the THIRD lattice vector.
        double Vec3Diff1 = (pos(i, 0) - pos(j, 0)) * latticeVec3Norm(0);
        double Vec3Diff2 = (pos(i, 1) - pos(j, 1)) * latticeVec3Norm(1);
        double Vec3Diff3 = (pos(i, 2) - pos(j, 2)) * latticeVec3Norm(2);
        double magVec3Diff = sqrt(
                pow(Vec3Diff1, 2) +
                pow(Vec3Diff2, 2) +
                pow(Vec3Diff3, 2));
        double magVec3Lattice = sqrt(
                pow(latticeVec3(0), 2) +
                pow(latticeVec3(1), 2) +
                pow(latticeVec3(2), 2));

        //if we are not within half a lattice vector length, then it is possible to translate some more.
        // while (magVec3Diff > magVec3Lattice / 2) {
        //need to find out if we are adding or subtracting the lattice vector.
        //The one that reduces the magVec1Diff will be used.

        double AddVec3Diff1 = Vec3Diff1 + latticeVec3(0);
        double AddVec3Diff2 = Vec3Diff2 + latticeVec3(1);
        double AddVec3Diff3 = Vec3Diff3 + latticeVec3(2);
        double AddMagVec3Diff = sqrt(
                pow(AddVec3Diff1, 2) +
                pow(AddVec3Diff2, 2) +
                pow(AddVec3Diff3, 2));

        double SubVec3Diff1 = Vec3Diff1 - latticeVec3(0);
        double SubVec3Diff2 = Vec3Diff2 - latticeVec3(1);
        double SubVec3Diff3 = Vec3Diff3 - latticeVec3(2);
        double SubMagVec3Diff = sqrt(
                pow(SubVec3Diff1, 2) +
                pow(SubVec3Diff2, 2) +
                pow(SubVec3Diff3, 2));

        //choose if we are adding or subtracting the vector based on which is smaller in an absolute way.
        //choose if we are adding or subtracting the vector based on which is smaller in an absolute way.
        if (SubMagVec3Diff < magVec3Diff && SubMagVec3Diff < AddMagVec3Diff) {
            Vec3Diff1 = Vec3Diff1 - latticeVec3(0);
            Vec3Diff2 = Vec3Diff2 - latticeVec3(1);
            Vec3Diff3 = Vec3Diff3 - latticeVec3(2);
            magVec3Diff = sqrt(
                    pow(Vec3Diff1, 2) +
                    pow(Vec3Diff2, 2) +
                    pow(Vec3Diff3, 2));


            transFlag = true;
        } else if (AddMagVec3Diff < magVec3Diff && AddMagVec3Diff < SubMagVec3Diff) {
            Vec3Diff1 = Vec3Diff1 + latticeVec3(0);
            Vec3Diff2 = Vec3Diff2 + latticeVec3(1);
            Vec3Diff3 = Vec3Diff3 + latticeVec3(2);
            magVec3Diff = sqrt(
                    pow(Vec3Diff1, 2) +
                    pow(Vec3Diff2, 2) +
                    pow(Vec3Diff3, 2));
            transFlag = true;

        }

        if (transFlag == true) {
            //cout << "PERIODIC TRANSLATION BETWEEN ATOMS i, j " << i << " " << j << endl;
        }
        // }

        zDiff = magVec3Diff;




    }

    distance = sqrt(
            pow(xDiff, 2) +
            pow(yDiff, 2) +
            pow(zDiff, 2));
    if (distance == 0) {
        cerr << "two atoms occupy the same position!" << endl;
        system("read  -p \"two atoms occupy the same position! Push enter to unpause.\" key");

    }

}
