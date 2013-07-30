#include "ReadWriteFunctionsATAC.h"
//given a string, tries to match that string in a given variable file, and returns the
//value immediately next to the keyword.

string variableReader(string keyWord) {

    ifstream inF("variables.txt");
    int lineNumber = 0;
    string line;

    while (getline(inF, line)) {
        istringstream source(line);
        string token;
        source >> token;

        //cout << keyWord << " TOKEN: " << token << endl;
        if (keyWord.compare(token.c_str()) == 0) {
            source >> token;
            string result;
            result = (token);
            //cout << "RESULTS" << result << endl;
            inF.close();
            return result;
        }
    }
    //++lineNumber;
    inF.close();
    cout << "CANNOT FIND GIVEN KEYWORD " << keyWord << " IN VARIABLE FILE " << endl;
    return 0;
}

//given a string, tries to match that string in a given configuration file, and returns the
//value immediately next to the keyword. Currently reads in the second token as the key, which means it only works for reading temperature.

string confReader(string keyWord, string masterBatchName, int keyWordPosition) {

    ifstream inF((masterBatchName + ".conf").c_str());
    int lineNumber = 0;
    string line;


    while (getline(inF, line)) {
        istringstream source(line);
        string token;

        //cout << keyWord << " keyWordPostition: " << keyWordPosition << endl;

        //   cout << keyWord << " TOKEN: " << token << endl;

        for (int i = 0; i < keyWordPosition; i++) {
            source >> token;
            // cout << keyWord << " TOKEN: " << token << endl;
        };

        //cout <<" TOKEN : " << token;



        if (keyWord.compare(token.c_str()) == 0) {
            source >> token;
            string result;
            result = (token);
            //cout << "RESULTS" << result << endl;
            inF.close();
            return result;
        }
    }
    //++lineNumber;
    inF.close();
    cout << "CANNOT FIND GIVEN KEYWORD " << keyWord << " IN CONF FILE " << endl;
    return 0;
}




//This subfunction writes generates a datafile to be used in a final movie
//showing the evolution of the sizes of the carbon rings.

void recordRingSizes(vector<pair<imat, int> > rings, int& globalRingSizeRunNum, bool& generateRingHistogramImages) {

    int ringSizeMin = atoi(variableReader("ringSizeMin").c_str());
    int ringSizeMax = atoi(variableReader("ringSizeMax").c_str());
    bool ringOfAnySizeFound = 0;




    //find if rings of six exist
    for (int i = 0; i < rings.size(); ++i) {
        //we only want rings that are in between the min and max ring size (inclusive).

        if ((rings[i].second >= ringSizeMin) && (rings[i].second <= ringSizeMax)) {
            ringOfAnySizeFound = 1;
            chdir("ringSizes");
            imat ringsOfASize = rings[i].first;
            //start file with chiral angle points.


            string globalRunNumS;
            stringstream ssglobalRunNum;
            ssglobalRunNum << globalRingSizeRunNum;
            globalRunNumS = ssglobalRunNum.str();

            string globalRunNum1S;
            stringstream ssglobalRunNum1;
            ssglobalRunNum1 << globalRingSizeRunNum + 1;
            globalRunNum1S = ssglobalRunNum1.str();

            string sizeS;
            stringstream sssize;
            sssize << ringsOfASize.n_cols;
            sizeS = sssize.str();



            //writes down the size of the ring for each ring.
            for (int j = 0; j < ringsOfASize.n_rows; j++) {
                system(("echo \"" + sizeS + "\" >> " + globalRunNum1S + ".dat").c_str());
            }


            //we must label the jpegs of the plots in an appropriate manner to be
            //converted into an mpeg later.
            string label;
            stringstream number;
            if ((globalRingSizeRunNum + 1) < 10) {
                label = "0000";
            } else if ((globalRingSizeRunNum + 1) < 100) {
                label = "000" + number.str();
            } else if ((globalRingSizeRunNum + 1) < 1000) {
                label = "00" + number.str();
            } else if ((globalRingSizeRunNum + 1) < 10000) {
                label = "0" + number.str();
            } else if ((globalRingSizeRunNum + 1) < 100000) {
                label = "" + number.str();
            } else {
                cerr << "MORE THAN 99,999 frames" << endl;
            }
            number << (globalRingSizeRunNum + 1);
            label.append(number.str());

            if (generateRingHistogramImages == 1) {
                //updates the histogram script with the current run data, executes it, and generates a jpeg.
                system(("sed \"s/0.dat/" + globalRunNum1S + ".dat/g\" ringSizePlotScript.plt | gnuplot > pic" + label + ".jpg").c_str());
            }


            chdir("..");

        }


    }
    if (ringOfAnySizeFound == true) {
        globalRingSizeRunNum++;

    }
}



//Writes an int in the C000 format used by vmd. 1000 atom limit.

string writeLabels(int i) {
    string label;
    stringstream number;


    if (i < 10) {
        label = "C00";
    } else if (i < 100) {
        label = "C0" + number.str();
    } else if (i < 1000) {
        label = "C" + number.str();
    } else {
        label = "C" + number.str();
        //cerr << "More than 999 atoms ";
    }

    number << i;
    label.append(number.str());

    // */
    return label;
}


void writeTopFile(mat positions, imat bondList, string masterBatchName) {
    //Writes topology file.



    //system("read  -p \"writingtop file\" key");





    ofstream writeFile;
    writeFile.open((masterBatchName + ".top").c_str());
    writeFile << "0 0" << endl;
    writeFile << "MASS     1 CA    12.01100 C" << endl;
    writeFile << "AUTO ANGLES DIHE" << endl;
    writeFile << "RESI TUB          0.00" << endl;
    for (int i = 0; i < positions.n_rows; i++) {

        writeFile << "ATOM " << writeLabels(i) << " CA" << "      0.00" << endl;
    }
    // cout << "BONDS AFST2" << bondList << endl;
    for (int i = 0; i < bondList.n_rows; i++) {
        writeFile << "BOND  " << writeLabels(bondList(i, 0)) << " "
                << writeLabels(bondList(i, 1)) << endl;
    }
    writeFile.close();


}

void writePgnFile(mat positions, string masterBatchName) {


    //Writes .pgn file.
    ofstream writeFile2;
    writeFile2.open((masterBatchName + ".pgn").c_str());
    writeFile2 << " topology " << masterBatchName << ".top" << endl;
    writeFile2 << " segment NT {residue 1 TUB}" << endl;
    for (int i = 0; i < positions.n_rows; i++) {
        writeFile2 << "coord NT 1 " << writeLabels(i) << " {  "
                << positions(i, 0) << "  "
                << positions(i, 1) << "  "
                << positions(i, 2) << " }" << endl;
    }
    writeFile2 << "writepdb " << masterBatchName << ".pdb" << endl;
    writeFile2 << "writepsf " << masterBatchName << ".psf" << endl;
    writeFile2 << "exit" << endl;
    writeFile2.close();



}




//Note about pdb files: positions are aloted 8 characters, with the 5th character being a decimal.
//for instance: -088.880
// we must round forces so that they fit this format, including adding white space
//before and after if the force is not already exactly 8 characters.

void writeForcePDB(mat positions, mat forceMatrix) {


    //must erase old force file.
    cout << "REMOVING OLD forces.pdb FILE" << endl;
    system("rm forces.pdb");

    ofstream forceFile;
    forceFile.open("forces.pdb"); //change later

    forceFile << "REMARK force file" << endl;
    //following code loops through x, y, and z, but variable naming is as if we took the x case.

    double force;

    for (int i = 0; i < forceMatrix.n_rows; i++) {
        string line;
        stringstream lineSS;

        for (int j = 0; j < forceMatrix.n_cols; ++j) {
            force = forceMatrix(i, j);
            //find number of digits before the period, INCLUDING the minus sign as a digit.
            int xDigits = 0;

            if (force >= 0) {

                //check for overflow.
                if (force > 999.999) {
                    //cout << j << " OVERFLOW, RESET TO 999.999" << endl;
                    force = 999.999;

                }
                if (force < 0.001 && force != 0) {
                    //cout << j << " UNDERFLOW, RESET TO 0" << endl;
                    force = 0;
                }
                //divides and counts digits as it loops
                int temp = floor(force);
                do {
                    ++xDigits;
                    temp /= 10;
                } while (temp > 0);

            } else {//this is the negative case.
                //check for overflow.
                if (force < -999.999) {
                    //cout << j << " OVERFLOW, RESET TO -999.999" << endl;
                    force = -999.999;
                }
                if (force > -0.001 && force != 0) {
                    //cout << j << " UNDERFLOW, RESET TO 0" << endl;
                    force = 0;
                }
                //divides and counts digits as it loops
                int temp = floor(-force);
                while (temp > 0) {
                    xDigits++;
                    temp /= 10;
                }
                //include minus sign
                ++xDigits;
            }

            //if starting string stream, start off with known values...
            if (j == 0) {
                //keeping in mind white space...
                //HOW MANY ATOMS CAN NAMD HANDLE? MY CODE HANDLE? CHECK LATER AND CHANGE THIS SECTION.
                if (i + 1 < 10) {
                    lineSS << "ATOM" << "     " << i + 1 << " " << writeLabels(i) << " TUB     1    ";
                } else {
                    if (i + 1 < 100) {
                        lineSS << "ATOM" << "    " << i + 1 << " " << writeLabels(i) << " TUB     1    ";
                    } else {
                        if (i + 1 < 1000) {
                            lineSS << "ATOM" << "   " << i + 1 << " " << writeLabels(i) << " TUB     1    ";
                        } else {
                            if (i + 1 < 10000) {
                                lineSS << "ATOM" << "  " << i + 1 << " " << writeLabels(i) << " TUB     1    ";
                            }
                        }
                    }
                }
            }
            string xS;
            stringstream ssX;
            ssX << force;
            xS = ssX.str();

            //number of decimals. The period IS considered a decimal for the sake of this algorithm.
            int xDecimals = xS.length() - xDigits;
            //checks that we have the decimal point (will not if given integer.)
            if (xDecimals <= 0) {
                xS.append(".");
                xDecimals = 1;
            }
            //now, we fill out the xString so that it 8 characters long, adding spaces where appropriate.
            if (xDigits < 4) {
                for (int i = 0; i < (4 - xDigits); ++i) {
                    xS.insert(0, " ");
                }
            }
            if (xDecimals < 4) {
                for (int i = 0; i < (4 - xDecimals); ++i) {
                    xS.append(" ");
                }
            }
            //must also look for having too many decimals.
            if (xDecimals > 4) {
                xS.resize(8);
            }
            //append
            lineSS << xS;
        }
        lineSS << "  1.00  0.00      NT   C";
        line = lineSS.str();

        forceFile << line << endl;
    }
    forceFile << "END" << endl;
    forceFile.close();



}

void writeConfFile(mat positions, double temperature, string masterBatchName, string& proteinDataBaseFile, int& forceNum, int& periodicityNum, double& xLatticeParameter, double& yLatticeParameter, double& zLatticeParameter, vec latticeVec1, vec latticeVec2, vec latticeVec3) {
    //Writes .conf file
    ofstream writeConfig;
    writeConfig.open((masterBatchName + ".conf").c_str()); //change later
    writeConfig << "structure          " << masterBatchName << ".psf" << endl;
    writeConfig << "coordinates        " << masterBatchName << ".pdb" << endl;
    writeConfig << endl;
    writeConfig << "set temperature    " << temperature << endl;
    writeConfig << "set outputname     " << masterBatchName << endl;
    writeConfig << "binaryoutput	no" << endl;
    writeConfig << endl;
    writeConfig << "firsttimestep      0" << endl;
    writeConfig << endl;
    writeConfig << "paraTypeCharmm	    on" << endl;
    writeConfig << "parameters          " << proteinDataBaseFile << endl; //par_all27_prot_lipid.inp
    writeConfig << "temperature         $temperature" << endl;
    writeConfig << endl;
    writeConfig << "exclude             scaled1-4" << endl;
    writeConfig << "1-4scaling          1.0" << endl;
    writeConfig << "cutoff              12.0" << endl;
    writeConfig << "switching           on" << endl;
    writeConfig << "switchdist          10.0" << endl;
    writeConfig << "pairlistdist        13.5" << endl;
    writeConfig << endl;
    writeConfig << "timestep            0.5" << endl;
    writeConfig << "rigidBonds          all" << endl;
    writeConfig << "nonbondedFreq       1" << endl;
    writeConfig << "fullElectFrequency  2" << endl;
    writeConfig << "stepspercycle       10" << endl;
    writeConfig << endl;
    writeConfig << "langevin            on" << endl;
    writeConfig << "langevinDamping     5" << endl;
    writeConfig << "langevinTemp        $temperature" << endl;
    writeConfig << "langevinHydrogen    off" << endl;
    writeConfig << endl;

    //if we have non-zero forces, declare force file.
    // constant external force
    if (forceNum != 0) {
        writeConfig << "constantForce yes" << endl;
        writeConfig << "consForceFile forces.pdb" << endl;
    }

    if (periodicityNum == 1) {

        mat xValues = positions.col(0);

        writeConfig << "cellBasisVector1    " << (xLatticeParameter) << "    0     0" << endl;
        writeConfig << "cellOrigin   " << xValues.min() << "    0" << "    0" << endl;
    }

    if (periodicityNum == 2) {

        mat xValues = positions.col(0);
        mat yValues = positions.col(1);

        writeConfig << "cellBasisVector1    " << (xLatticeParameter) << "    0     0" << endl;
        writeConfig << "cellBasisVector2    0    " << (yLatticeParameter) << "     0" << endl;
        writeConfig << "cellOrigin   " << xValues.min() << "    " << yValues.min() << "    0" << endl;
    }
    if (periodicityNum == 3) {

        mat xValues = positions.col(0);
        mat yValues = positions.col(1);

        writeConfig << "cellBasisVector1    " << (xLatticeParameter) << "    0     0" << endl;
        writeConfig << "cellBasisVector2    0    " << (yLatticeParameter) << "     0" << endl;
        writeConfig << "cellBasisVector3    0     0     " << (zLatticeParameter) << endl;
        writeConfig << "cellOrigin   " << xValues.min() << "    " << yValues.min() << "    0" << endl;
    }

    if (periodicityNum == 4) {

        mat xValues = positions.col(0);
        mat yValues = positions.col(1);
        mat zValues = positions.col(2);

        writeConfig << "cellBasisVector1    " << latticeVec1(0) << "    " << latticeVec1(1) << "    " << latticeVec1(2) << endl;
        writeConfig << "cellBasisVector2    " << latticeVec2(0) << "    " << latticeVec2(1) << "    " << latticeVec2(2) << endl;
        writeConfig << "cellBasisVector3    " << latticeVec3(0) << "    " << latticeVec3(1) << "    " << latticeVec3(2) << endl;

        writeConfig << "cellOrigin   " << xValues.min() << "    " << yValues.min() << "    " << zValues.min() << endl;
    }

    writeConfig << endl;
    writeConfig << "outputName          $outputname" << endl;
    writeConfig << "restartfreq         500 " << endl;
    writeConfig << "dcdfreq             100000000000" << endl;
    writeConfig << "outputEnergies      1000" << endl;
    writeConfig << "outputPressure      1000" << endl;
    writeConfig << endl;
    writeConfig << "minimize            250" << endl;
    writeConfig << "reinitvels          $temperature" << endl;
    writeConfig << "run 250" << endl;
    writeConfig << "minimize            250" << endl;

    writeConfig.close();
}

void writePdbPsfFiles(string masterBatchName) {

    // Executes vmd in the command line, generating the .pdb and .psf files.
    cout << ("PSFGEN RUNNING") << endl;
    //system("mesg n");
    system(("psfgen -dispdev text -eofexit < " + masterBatchName + ".pgn >/dev/null ").c_str());
    //system("read  -p \"namd failed, paused\" key");
    //  system("mesg y");

    //system(("vmd "+batchName+".psf "+batchName+".pdb").c_str());

}

//strips a coor xyz file to make it a matrix.

void coorStripper(mat& positions, string masterBatchName) {

    string line;
    ifstream coorFile((masterBatchName + ".coor").c_str());
    ofstream coorFileStripped((masterBatchName + "Stripped.coor").c_str());
    int lineNumber = 0;
    while (getline(coorFile, line)) {
        //skips first line
        ++lineNumber;
        if ((lineNumber > 1) && (line.size() > 4)) {
            line.erase(0, 26);
            line.erase(27, 51);
            coorFileStripped << line << "\n";
        }
    }
    coorFile.close();
    coorFileStripped.close();
    positions.load((masterBatchName + "Stripped.coor").c_str());
    system(("rm " + masterBatchName + "Stripped.coor").c_str());
}

//strips xyz file to make it a matrix.

void xyzStripper(mat& positions, string masterBatchName) {

    string line;
    ifstream coorFile((masterBatchName + ".xyz").c_str());
    ofstream coorFileStripped((masterBatchName + "Stripped.xyz").c_str());
    int lineNumber = 0;
    while (getline(coorFile, line)) {
        //skips first two lines
        ++lineNumber;
        if ((lineNumber > 2) && (line.size() > 4)) {
            line.erase(0, 2);
            coorFileStripped << line << "\n";
        }
    }
    coorFile.close();
    coorFileStripped.close();
    positions.load((masterBatchName + "Stripped.xyz").c_str());
    system(("rm " + masterBatchName + "Stripped.xyz").c_str());
}



//Reads the coor file of the last run, and writes to an .xyz file. This is then
//appended to movie.xyz.

void movieMaker(mat positions, string masterBatchName, int& numberOfElectrodes, imat& xp1Electrodes, imat& xp2Electrodes, imat& xn1Electrodes, imat& xn2Electrodes, imat& yp1Electrodes, imat& yp2Electrodes, imat& yn1Electrodes, imat& yn2Electrodes) {

    system(("cp " + masterBatchName + ".coor movies").c_str());
    chdir("movies");
    ifstream coorFile((masterBatchName + ".coor").c_str());
    ofstream lastRunFile("lastAcceptedRun.xyz");
    string line;
    //preps the .xyz file, number of atoms, and a blank line.
    lastRunFile << positions.n_rows << "\n";
    lastRunFile << "\n";
    int lineNumber = 0;

    while (getline(coorFile, line)) {
        ++lineNumber;
        //skips first line, and ignores the last END line.
        if ((lineNumber > 1) && (line.size() > 4)) {
            line.erase(0, 26);
            line.erase(27, 51);

            //keep in mind electrode labeling...
            if (numberOfElectrodes != 0) {
                bool electrodeFound = 0;

                if (numberOfElectrodes >= 2) {

                    if (electrodeFound == 0) {
                        //loop through our electrode atoms.
                        for (int j = 0; j < xp1Electrodes.n_rows; j++) {
                            //see if the line we are currently on is an electrode atom. Atoms start at line
                            //3 (starting our count from 1), so subtract 3 from lineNumber.
                            if ((xp1Electrodes(j) == lineNumber - 2)) {

                                lastRunFile << "O     " << line << "\n";
                                electrodeFound = 1;
                                break;
                            }
                        }
                    }

                    if (electrodeFound == 0) {
                        //loop through our electrode atoms.
                        for (int j = 0; j < xp2Electrodes.n_rows; j++) {
                            //see if the line we are currently on is an electrode atom. Atoms start at line
                            //3 (starting our count from 1), so subtract 3 from lineNumber.
                            if ((xp2Electrodes(j) == lineNumber - 2)) {
                                lastRunFile << "F     " << line << "\n";
                                electrodeFound = 1;
                                break;
                            }
                        }
                    }

                    if (electrodeFound == 0) {
                        //loop through our electrode atoms.
                        for (int j = 0; j < xn1Electrodes.n_rows; j++) {
                            //see if the line we are currently on is an electrode atom. Atoms start at line
                            //3 (starting our count from 1), so subtract 3 from lineNumber.
                            if ((xn1Electrodes(j) == lineNumber - 2)) {

                                lastRunFile << "N     " << line << "\n";
                                electrodeFound = 1;
                                break;
                            }
                        }
                    }
                    if (electrodeFound == 0) {
                        //loop through our electrode atoms.
                        for (int j = 0; j < xn2Electrodes.n_rows; j++) {
                            //see if the line we are currently on is an electrode atom. Atoms start at line
                            //3 (starting our count from 1), so subtract 3 from lineNumber.
                            if ((xn2Electrodes(j) == lineNumber - 2)) {


                                lastRunFile << "B     " << line << "\n";
                                electrodeFound = 1;
                                break;
                            }
                        }
                    }

                    if (numberOfElectrodes == 4) {

                        if (electrodeFound == 0) {
                            //loop through our electrode atoms.
                            for (int j = 0; j < yp1Electrodes.n_rows; j++) {
                                // if(lineNumber-2 > 450){
                                // cout <<" line number-2, yp1 electrodes: " << lineNumber - 2 << "  " << yp1Electrodes(j)<< endl;

                                //}
                                //see if the line we are currently on is an electrode atom. Atoms start at line
                                //3 (starting our count from 1), so subtract 3 from lineNumber.
                                if ((yp1Electrodes(j) == lineNumber - 2)) {
                                    lastRunFile << "S     " << line << "\n";
                                    electrodeFound = 1;
                                    break;
                                }
                            }
                        }

                        if (electrodeFound == 0) {
                            //loop through our electrode atoms.
                            for (int j = 0; j < yp2Electrodes.n_rows; j++) {
                                //see if the line we are currently on is an electrode atom. Atoms start at line
                                //3 (starting our count from 1), so subtract 3 from lineNumber.
                                if ((yp2Electrodes(j) == lineNumber - 2)) {
                                    lastRunFile << "Y     " << line << "\n";
                                    electrodeFound = 1;
                                    break;
                                }
                            }
                        }

                        if (electrodeFound == 0) {
                            //loop through our electrode atoms.
                            for (int j = 0; j < yn1Electrodes.n_rows; j++) {
                                //see if the line we are currently on is an electrode atom. Atoms start at line
                                //3 (starting our count from 1), so subtract 3 from lineNumber.
                                if ((yn1Electrodes(j) == lineNumber - 2)) {
                                    lastRunFile << "P     " << line << "\n";
                                    electrodeFound = 1;
                                    break;
                                }
                            }
                        }
                        if (electrodeFound == 0) {
                            //loop through our electrode atoms.
                            for (int j = 0; j < yn2Electrodes.n_rows; j++) {
                                //see if the line we are currently on is an electrode atom. Atoms start at line
                                //3 (starting our count from 1), so subtract 3 from lineNumber.
                                if ((yn2Electrodes(j) == lineNumber - 2)) {
                                    lastRunFile << "I     " << line << "\n";
                                    electrodeFound = 1;
                                    break;
                                }
                            }
                        }
                    }
                }
                //if current atom is not an electrode, make carbon.
                if (electrodeFound == 0) {
                    lastRunFile << "C     " << line << "\n";
                }

            } else {
                //not considering electrodes
                lastRunFile << "C     " << line << "\n";
            }

        }
    }
    coorFile.close();
    lastRunFile.close();
    //makes a copy of last run
    // system("cp lastRun.xyz lastAcceptedRun.xyz");

    //appends xyz to movie.
    system(("cat lastAcceptedRun.xyz >> movie" + masterBatchName + ".xyz").c_str());
    system(("mv " + masterBatchName + ".coor lastAcceptedRun.coor").c_str());
    chdir("..");
}

//this subfunction reads in an energy from a NAMD .log file.

double energyReader(double& newEnergy, string masterBatchName, int& relaxationMethod) {
    //namd selected
    if(relaxationMethod == 1){
    ifstream logFile((masterBatchName + ".log").c_str());
    string line2;
    //simply finds the number of lines in the log file
    int lineNumber2 = 0;
    while (getline(logFile, line2)) {
        ++lineNumber2;
    }
    logFile.close();

    //read file to get to line "last line" - 11, where the final energy data is.
    ifstream logFile2((masterBatchName + ".log").c_str());
    int lineNumber3 = 0;
    string line3;
    while (getline(logFile2, line3)) {
        if (lineNumber3 == (lineNumber2 - 12)) {
            // reads the fourteenth word in line3, which is the energy.
            int wordNum = 0;
            istringstream source(line3);
            string token;
            while (source.good()) {
                source >> token;
                if (wordNum == 11) {
                    newEnergy = strtod(token.c_str(), NULL);
                    break;
                }
                ++wordNum;
            }
        }
        ++lineNumber3;
    }
    logFile2.close();
    }
}


//This function updates the run number by reading the last inputted value of the energy history file.
//If the run type is 0, it will find the number of total runs.
//If 1, will find the number of accepted runs.

int runReader(int& runs, int runType, string batchName) {

    string line;
    string modifier;
    if (runType == 0) {
        modifier = "energyHist";
    } else if (runType == 1) {
        modifier = "energyHistMovie";
    } else {
        system("read -p \"runType of runReader is an unknown number! Push enter to unpause.\" key");
    }

    ifstream logFile((modifier + batchName + ".dat").c_str());
    string line2;
    //simply finds the number of lines in the log file
    int lineNumber2 = 0;
    while (getline(logFile, line2)) {
        ++lineNumber2;
    }
    logFile.close();

    //read file to get to line "last line", where the final run number is.


    ifstream logFile2((modifier + batchName + ".dat").c_str());


    int lineNumber3 = 0;
    string line3;
    while (getline(logFile2, line3)) {
        if (lineNumber3 == (lineNumber2 - 1)) {
            // reads the first word in line3, which is the run number.
            int wordNum = 0;
            istringstream source(line3);
            string token;
            while (source.good()) {
                source >> token;
                if (wordNum == 0) {
                    runs = strtod(token.c_str(), NULL);
                    ///cout << "WOOOOOOOOOOOOOOOOO " << runs << endl;
                    break;
                }
                ++wordNum;
            }
        }
        //cout << "LINENUMBERS 2 and 3 " << lineNumber3 << " " << lineNumber2 << endl;
        ++lineNumber3;
    }
    logFile2.close();

}
