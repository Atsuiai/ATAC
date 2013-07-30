#include "systemATAC.h"


    
    
void systemATAC::initialization(int argc, char* argv[]) {


        // system("pwd");
    // system("read -p \"paused\" key");

    isRestart = atof(variableReader("isRestart").c_str());
    cout << setprecision(15);
    //get batch name identifier.
    if (argc < 2) {
        cout << "Not enough arguments, ATAC needs 1: an a batch name identifier" << endl;
        exit(0);
    } else if (argc > 2) {
        cout << "Too many arguments, ATAC needs 1: an a batch name identifier." << endl;
        exit(0);
    }
    char* batchIdentifier = argv[1];

    //the batch identifier string is converted to a number
    //and added to the random seed (hopefully people will not run
    //parallel runs with the exact same name at the exact same time...)
    char* givenInt = argv[1];
    int givenIntConverted = atoi(givenInt);
    srand(time(NULL) + givenIntConverted);

    string masterBatchName = variableReader("masterBatchName").c_str();
    batchName = (masterBatchName + batchIdentifier).c_str();
    
    string proteinDataBaseFile = variableReader("proteinDataBaseFile").c_str();

    if (isRestart == false) {

        //this clears out the working directory, handy if you are debugging and
        //using the same master batch name again and again.

        cout << "REMOVING OLD FILES AND DIRECTORIES" << endl;
        system(("rm -r " + batchName + " >/dev/null").c_str());

        //create a directory of the masterBatchName. copy over neccessary files.
        system(("mkdir " + batchName).c_str());
        system(("cp pos.xyz " + batchName).c_str());
        system(("cp " + proteinDataBaseFile + " " + batchName).c_str());
        //system(("cp test.slurm "+masterBatchName).c_str());
        system(("cp variables.txt " + batchName).c_str());
        system(("cp nonInteractingAtoms.xyz " + batchName).c_str());
    }
    // cout << "BEFORE MASTER " << endl;
    // cout << system("pwd") << endl;
    // cout << "MASTER BATCH NAME " << masterBatchName << endl;
    system("pwd");
    cout << "MASTER BATCH NAME " << batchName << endl;
    chdir((batchName).c_str());
    // cout << "AFTER MASTER " << endl;
    // cout << system("pwd") << endl;
    system("pwd");





    globalChiralAngleRunNum = 0; //this is the same as the number of accepted steps
    globalRingSizeRunNum = 0; //this is the same as the number of accepted steps

    double randomPercent;
    //temperature in kelvin
    temperature;

 

    vector<double> ratiosF = threeRandomPercentages();
     mutationGo = atof(variableReader("mutationGo").c_str());

    ratioRotate;
    ratioDimerDisplacement;
    ratioSp2ToSp3;

    //if doing no-flux mutations, call ratios
    if (mutationGo == 1) {
        ratioRotate = atof(variableReader("ratioRotate").c_str());
        ratioDimerDisplacement = atof(variableReader("ratioDimerDisplacement").c_str());
        ratioSp2ToSp3 = atof(variableReader("ratioSp2ToSp3").c_str());
        //check that ratios add up to one

        bool inequalityBool = ratioRotate + ratioDimerDisplacement + ratioSp2ToSp3 != 1;

        if (inequalityBool != 1) {
            cout << "Rotate ratio: " << ratioRotate << " Dimer displacement ratio: " << ratioDimerDisplacement << " sp2 to sp3 parallel ratio: " << ratioSp2ToSp3 << endl;
            //    system("read -p \"ratios do NOT add up to unity! Push enter to unpause.\" key");

        }
    }
    
    
    
    
    
    
    
    
    
    relaxationMethod = atof(variableReader("relaxationMethod").c_str());
    noFluxRatioShift = atof(variableReader("noFluxRatioShift").c_str());
    ratioShiftRate = atof(variableReader("ratioShiftRate").c_str());
    ratioRotateFinal = atof(variableReader("ratioRotateFinal").c_str());
    ratioSp2ToSp3Final = atof(variableReader("ratioSp2ToSp3Final").c_str());
    reversableSpRunNum = atof(variableReader("reversableSpRunNum").c_str());
    superfluousSp3Sp2 = atof(variableReader("superfluousSp3Sp2").c_str());
    rotationsPerRun = atof(variableReader("rotationsPerRun").c_str());
    proteinDataBaseFile = variableReader("proteinDataBaseFile").c_str();
    bondDistanceToleranceMax = atof(variableReader("bondDistanceToleranceMax").c_str());
    bondDistanceToleranceMin = atof(variableReader("bondDistanceToleranceMin").c_str());
    ringSizeMin = atoi(variableReader("ringSizeMin").c_str());
    ringSizeMax = atoi(variableReader("ringSizeMax").c_str());

    //get number of processors.
    numProcessors = variableReader("numProcessors").c_str();

    forceNum = atof(variableReader("forceNum").c_str());
    generateRingHistogramImages = atof(variableReader("generateRingHistogramImages").c_str());
    tempS = atof(variableReader("tempS").c_str());

    if (isRestart == false) {
        temperature = tempS;
    } else {
        temperature = atof(confReader("temperature", batchName, 2).c_str());
    }

    coolMethod = atof(variableReader("coolMethod").c_str());
    if (coolMethod != 0) {
        tempF = atof(variableReader("tempF").c_str());

        coolR = atof(variableReader("coolR").c_str());
    }

    periodicityNum = atof(variableReader("periodicityNum").c_str());
    if (periodicityNum != 0) {
        periodicStretching = atof(variableReader("periodicStretching").c_str());
        if (periodicStretching != 0) {
            periodicStretchRate = atof(variableReader("periodicStretchRate").c_str());
        }
        isStaticPeriod = (atof(variableReader("isStaticPeriod").c_str()));
        xPeriodFactor = (atof(variableReader("xPeriodFactor").c_str()));
        yPeriodFactor = (atof(variableReader("yPeriodFactor").c_str()));
        zPeriodFactor = (atof(variableReader("zPeriodFactor").c_str()));

        //in this case arbitrary periodicity
        if (periodicityNum == 4) {
            double LVec11 = (atof(variableReader("LVec11").c_str()));
            double LVec12 = (atof(variableReader("LVec12").c_str()));
            double LVec13 = (atof(variableReader("LVec13").c_str()));
            double LVec21 = (atof(variableReader("LVec21").c_str()));
            double LVec22 = (atof(variableReader("LVec22").c_str()));
            double LVec23 = (atof(variableReader("LVec23").c_str()));
            double LVec31 = (atof(variableReader("LVec31").c_str()));
            double LVec32 = (atof(variableReader("LVec32").c_str()));
            double LVec33 = (atof(variableReader("LVec33").c_str()));

            latticeVec1(0) = LVec11; //lattice vector 1, dimension 1.
            latticeVec1(1) = LVec12;
            latticeVec1(2) = LVec13;
            latticeVec2(0) = LVec21; //lattice vector 2, dimension 1.
            latticeVec2(1) = LVec22;
            latticeVec2(2) = LVec23;
            latticeVec3(0) = LVec31;
            latticeVec3(1) = LVec32;
            latticeVec3(2) = LVec33;
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
        }
    }
    acceptAnyMutation = atof(variableReader("acceptAnyMutation").c_str());
    sp2tosp3Tolerance = atof(variableReader("sp2tosp3Tolerance").c_str()); // in angstroms. Experimentally derived : 2 A.
    spDenialPathSize = atof(variableReader("spDenialPathSize").c_str());

    numberOfElectrodes = atof(variableReader("numberOfElectrodes").c_str());
    if (numberOfElectrodes != 0) {
        if (numberOfElectrodes >= 2) {
            xpElecWallMod = atof(variableReader("xpElecWallMod").c_str());
            xnElecWallMod = atof(variableReader("xnElecWallMod").c_str());
        }

        if (numberOfElectrodes == 4) {
            ypElecWallMod = atof(variableReader("ypElecWallMod").c_str());
            ynElecWallMod = atof(variableReader("ynElecWallMod").c_str());
        }
    }


    //initilize beam variables.
    haveBeam = atoi(variableReader("haveBeam").c_str());
    if (haveBeam == 1) {

        beamType = atoi(variableReader("beamType").c_str());
        beamTemperature = atof(variableReader("beamTemperature").c_str());
        beamRadius = atof(variableReader("beamRadius").c_str());
        beamAxisIntercept1 = atof(variableReader("beamAxisIntercept1").c_str());
        beamAxisIntercept2 = atof(variableReader("beamAxisIntercept2").c_str());
        beamAxis = atoi(variableReader("beamAxis").c_str());
        beamRadiusGrowthRate = atof(variableReader("beamRadiusGrowthRate").c_str());
    }

    if (isRestart == false) {
        system("mkdir energies");
        chdir("energies");
        //writes file that will determine how energy plot is made
        system(("echo > plot movieEnergies.plt"));

        system(("echo 'set style data linespoints' >>  movieEnergies.plt"));
        system(("echo 'unset key' >>  movieEnergies.plt"));
        system(("echo 'set tmargin 3' >>  movieEnergies.plt"));
        system(("echo 'set bmargin 4.5' >>  movieEnergies.plt"));
        system(("echo 'set xtics border in scale 1,0.5 nomirror offset character 0, 0, 0' >>  movieEnergies.plt"));
        system(("echo 'set title \"Energy over Runs\" font \"Times-Roman,20\"' >>  movieEnergies.plt"));
        system(("echo 'set ylabel \"Total System Energy (kcal/mol)\" font \"Times-Roman,17\" offset -1,-1' >>  movieEnergies.plt"));
        system(("echo 'set xlabel \"# of Accepted Steps\" font \"Times-Roman,20\"' >>  movieEnergies.plt"));
        system(("echo 'plot \"energyHistMovie" + batchName + ".dat\" using 1:2' >>  movieEnergies.plt").c_str());


        system(("echo 'set style data linespoints' >>  attemptedEnergies.plt"));
        system(("echo 'unset key' >>  attemptedEnergies.plt"));
        system(("echo 'set tmargin 3' >>  attemptedEnergies.plt"));
        system(("echo 'set bmargin 4.5' >>  attemptedEnergies.plt"));
        system(("echo 'set xtics border in scale 1,0.5 nomirror offset character 0, 0, 0' >>  attemptedEnergies.plt"));
        system(("echo 'set title \"Energy over Runs\" font \"Times-Roman,20\"' >> attemptedEnergies.plt"));
        system(("echo 'set ylabel \"Total System Energy (kcal/mol)\" font \"Times-Roman,17\" offset -1,-1' >>  attemptedEnergies.plt"));
        system(("echo 'set xlabel \"# of Accepted Steps\" font \"Times-Roman,20\"' >>  attemptedEnergies.plt"));
        system(("echo 'plot \"energyHist" + batchName + ".dat\" using 1:2' >>  attemptedEnergies.plt").c_str());

        chdir("..");
    }

    findChiralAngles = atof(variableReader("findChiralAngles").c_str());



    //if we are finding the chiral angles, we must find the upper and lower distance
    //bounds which we use to find the next nearist neighbor atoms.
    //also, create directory to hold angle data files.
    if (findChiralAngles == true) {
        minChiralIndex = atoi(variableReader("minChiralIndex").c_str());
        maxChiralIndex = atoi(variableReader("maxChiralIndex").c_str());
        vecRoundingErrorCutOff = atof(variableReader("vecRoundingErrorCutOff").c_str());
        if (isRestart == false) {
            //make directory for chiral angles, set it up with a variables file
            system("mkdir chiralAngles");
            system("cp variables.txt chiralAngles");
            chdir("chiralAngles");
            //write file with chiral angle graph format
            ofstream chiralGraphFile;
            chiralGraphFile.open("histogramScript.plt");
            chiralGraphFile << "set terminal jpeg" << endl;
            chiralGraphFile << "unset key;" << endl;
            chiralGraphFile << "set border front;" << endl;
            chiralGraphFile << "red = \"#080000\";" << endl;
            chiralGraphFile << "set boxwidth 1 absolute" << endl;
            chiralGraphFile << "set xzeroaxis linetype 0 linewidth 1.000" << endl;
            chiralGraphFile << "set yzeroaxis linetype 0 linewidth 1.000" << endl;
            chiralGraphFile << "set zzeroaxis linetype 0 linewidth 1.000" << endl;
            chiralGraphFile << "set xlabel'Chiral Angle'" << endl;
            chiralGraphFile << "set ylabel'Hexagon Frequency'" << endl;
            chiralGraphFile << "set title \"Chiral Angle Distribution\" " << endl;
            chiralGraphFile << "set style fill solid border 3" << endl;
            chiralGraphFile << "set xtics border in scale 1,0.5" << endl;
            chiralGraphFile << "bin(x, s) = s*int(x/s)" << endl;
            chiralGraphFile << "GPFUN_bin = \"bin(x, s) = s*int(x/s)\"" << endl;
            chiralGraphFile << "plot [0:30] '0.dat'  u (bin($1,1)):(1)  lc rgb red s f t '' w boxes" << endl;
            //special note: this script takes in data of the name 0.dat. this file name is changed to all the other data sets (2,3, etc) using the unix command sed later on.
            //u 1:(0.25*rand(0)-.35) lc rgb red t '',      \"\"
            chiralGraphFile.close();

            chdir("..");
        }
    }

    findRingSizes = atof(variableReader("findRingSizes").c_str());
    //create directory to hold ring data files.
    if (findRingSizes == true) {
        int ringSizeMin = atoi(variableReader("ringSizeMin").c_str());
        int ringSizeMax = atoi(variableReader("ringSizeMax").c_str());
        if (isRestart == false) {
            //make directory for ring sizes, set it up with a variables file
            system("mkdir ringSizes");
            system("cp variables.txt ringSizes");
            chdir("ringSizes");
            //write file with ring size graph format
            ofstream chiralGraphFile;
            chiralGraphFile.open("ringSizePlotScript.plt");
            chiralGraphFile << "set terminal jpeg" << endl;
            chiralGraphFile << "unset key;" << endl;
            chiralGraphFile << "set border front;" << endl;
            chiralGraphFile << "red = \"#080000\";" << endl;
            chiralGraphFile << "set boxwidth 1 absolute" << endl;
            chiralGraphFile << "set xzeroaxis linetype 0 linewidth 1.000" << endl;
            chiralGraphFile << "set yzeroaxis linetype 0 linewidth 1.000" << endl;
            chiralGraphFile << "set zzeroaxis linetype 0 linewidth 1.000" << endl;
            chiralGraphFile << "set xlabel'Ring Size'" << endl;
            chiralGraphFile << "set ylabel'Frequency'" << endl;
            chiralGraphFile << "set title \"Ring Size Distribution\" " << endl;
            chiralGraphFile << "set style fill solid border 3" << endl;
            chiralGraphFile << "set xtics border in scale 1,0.5" << endl;
            chiralGraphFile << "bin(x, s) = s*int(x/s)" << endl;
            chiralGraphFile << "GPFUN_bin = \"bin(x, s) = s*int(x/s)\"" << endl;

            string ringSizeMinS;
            stringstream ssringSizeMin;
            ssringSizeMin << ringSizeMin;
            ringSizeMinS = ssringSizeMin.str();

            string ringSizeMaxS;
            stringstream ssringSizeMax;
            ssringSizeMax << ringSizeMax;
            ringSizeMaxS = ssringSizeMax.str();

            chiralGraphFile << ("plot [" + ringSizeMinS + ":" + ringSizeMaxS + "] '0.dat'  u (bin($1,1)):(1)  lc rgb red s f t '' w boxes").c_str() << endl;
            //special note: this script takes in data of the name 0.dat. this file name is changed to all the other data sets (2,3, etc) using the unix command sed later on.
            //u 1:(0.25*rand(0)-.35) lc rgb red t '',      \"\"
            chiralGraphFile.close();

            chdir("..");
        }
    }

    //checks if we have non-interacting components, and if so, if we will only apply periodic conditions to them.
    haveNonInteracting = atof(variableReader("haveNonInteracting").c_str());
    if (haveNonInteracting == true) {
        periodicOnlyOnNonInteracting = atof(variableReader("periodicOnlyOnNonInteracting").c_str());
    }

    if (isRestart == false) {
        //strips xyz file to make it a matrix.
        xyzStripper(positions, "pos");



        int periodicityNumBackup;
        if (periodicOnlyOnNonInteracting == true) {
            //IMPORTANT: CURRENTLY THE DYNAMIC PERIODIC CONDITION WILL NOT LOOK AT THE NORMAL ATOMS, ONLY THE NON-INTERACTING ONES.
            periodicityNumBackup = periodicityNum;
            periodicityNum = 0;
        }
        //       cout << "PSO " << positions << endl;



        if (periodicOnlyOnNonInteracting == true) {
            periodicityNum = periodicityNumBackup;
        }



        if (haveNonInteracting == true) {
            //get position and bond list matrixes for the non-interacting component.
            //strips xyz file to make it a matrix.
            string lineNon;
            ifstream xyzFileNon("nonInteractingAtoms.xyz");
            ofstream xyzFileStrippedNon(("nonStripped.xyz"));
            int lineNumberNon = 0;

            while (getline(xyzFileNon, lineNon)) {

                //skips first two lines
                ++lineNumberNon;
                if (lineNumberNon > 2) {
                    //erases first two characters, space and atom identifier.
                    xyzFileStrippedNon << lineNon.erase(0, 2) << "\n";
                }
            }

            xyzFileNon.close();
            xyzFileStrippedNon.close();
            nonInteracting.load(("nonStripped.xyz"));
            if (periodicityNum != 0) {
                latticeParameterDeterminer(nonInteracting, isStaticPeriod, periodicityNum, xLatticeParameter,  xPeriodFactor,  yLatticeParameter,  yPeriodFactor, zLatticeParameter, zPeriodFactor);
            }
            system(("rm nonStripped.xyz >/dev/null"));
            cout << "building non-interacting neighbors" << endl;
            nonInteractingNeighbors = buildNeigh(nonInteracting);
            //cout << "NONINTERACTING NEIGHBORS " << nonInteractingNeighbors;

            nonInteractingBonds = buildBond(nonInteractingNeighbors);


            numOfNonInteracting = nonInteracting.n_rows;
            numOfNonInteractingBonds = nonInteractingBonds.n_rows;
            numOfNonInteractingNeighbors = nonInteractingNeighbors.n_rows;



            //add non-interacting atoms to our system
            nonInteractingAppender(positions, bondList, neighborList);
        } else {
            //we are periodic but with only interacting atoms, so use the nomral positions.
            if (periodicityNum != 0) {
                latticeParameterDeterminer(positions, isStaticPeriod, periodicityNum, xLatticeParameter,  xPeriodFactor,  yLatticeParameter,  yPeriodFactor, zLatticeParameter, zPeriodFactor);
            }

        }
        neighborList = buildNeigh(positions);
        bondList = buildBond(neighborList);
        //cout << bondList << endl;


    } else {
        //This is a restart, so we must grab appropriate data from existing files.

        //Finds the positions of the last appropriate run.
        chdir("movies");
        //get postion matrix, and resulting neighbor and bond lists.
        xyzStripper(positions, "lastAcceptedRun");

        chdir("..");

        if (haveNonInteracting == true) {
            //If we have non-interacting atoms, we must seperate them from the
            //interacting ones. This will depend on the system in question,
            //but since I've only used the non-interacting feature with peapods,
            //I'm going to assume that the current system is a peapod.
            nonInteracting = positions;
            //shaves off the excess positions
            tubeShaver(positions);
            //subtract the normal atoms from the noninteracting matrix
            nonInteracting.shed_rows(0, positions.n_rows - 1);

        }

        int periodicityNumBackup;
        if (periodicityNum != 0) {
            if (periodicOnlyOnNonInteracting == true) {
                //IMPORTANT: CURRENTLY THE DYNAMIC PERIODIC CONDITION WILL NOT LOOK AT THE NORMAL ATOMS, ONLY THE NON-INTERACTING ONES.
                periodicityNumBackup = periodicityNum;
                periodicityNum = 0;

                latticeParameterDeterminer(nonInteracting, isStaticPeriod, periodicityNum, xLatticeParameter,  xPeriodFactor,  yLatticeParameter,  yPeriodFactor, zLatticeParameter, zPeriodFactor);
            } else {
                latticeParameterDeterminer(positions, isStaticPeriod, periodicityNum, xLatticeParameter,  xPeriodFactor,  yLatticeParameter,  yPeriodFactor, zLatticeParameter, zPeriodFactor);
            }
        }

        neighborList = buildNeigh(positions);
        //neighborList.save("neighborList", arma_ascii);
        //cout << neighborList << endl;

        if (periodicOnlyOnNonInteracting == true) {
            periodicityNum = periodicityNumBackup;
        }

        bondList = buildBond(neighborList);
        //cout << bondList << endl;

        if (haveNonInteracting == true) {
            //find non interacting atoms here.


            cout << "building non-interacting neighbors" << endl;
            nonInteractingNeighbors = buildNeigh(nonInteracting);
            //cout << "NONINTERACTING NEIGHBORS " << nonInteractingNeighbors;

            nonInteractingBonds = buildBond(nonInteractingNeighbors);

            numOfNonInteracting = nonInteracting.n_rows;
            numOfNonInteractingBonds = nonInteractingBonds.n_rows;
            numOfNonInteractingNeighbors = nonInteractingNeighbors.n_rows;

            //add non-interacting atoms to our system
            nonInteractingAppender(positions, bondList, neighborList);
        }
    }
    //write configuration file
   // writeConfFile(positions, temperature, batchName, periodicityNum,  xLatticeParameter,  yLatticeParameter, zLatticeParameter, latticeVec1, latticeVec2, latticeVec3);
 writeConfFile( positions, temperature, batchName, proteinDataBaseFile, forceNum, periodicityNum, xLatticeParameter, yLatticeParameter, zLatticeParameter, latticeVec1, latticeVec2, latticeVec3);
    
    
    //start electrodes
    if (numberOfElectrodes != 0) {
        electrodeDeterminer(positions, xp1Electrodes, xn1Electrodes, yp1Electrodes, yn1Electrodes, xp2Electrodes, xn2Electrodes, yp2Electrodes, yn2Electrodes, numberOfElectrodes,  xpElecWallMod,  xnElecWallMod,  ypElecWallMod, ynElecWallMod);
    }

    if (isRestart == false) {
        system("mkdir movies");
        chdir("movies");
        //Starts movie, top, pgn, pdb, and psf files.
        ofstream movieStart(("movie" + batchName + ".xyz").c_str());
        string lineM;
        //preps the .xyz file, number of atoms, and a blank line.

        movieStart << positions.n_rows << "\n";
        movieStart << "\n";
        int lineNumber = 2;
        for (int i = 0; i < positions.n_rows; ++i) {
            lineNumber++;

            //keep in mind electrode labeling...
            if (numberOfElectrodes != 0) {
                bool electrodeFound = 0;
                if (numberOfElectrodes >= 2) {



                    //loop through our electrode atoms.
                    for (int j = 0; j < xp1Electrodes.n_rows; j++) {
                        //see if the line we are currently on is an electrode atom. Atoms start at line
                        //3 (starting our count from 1), so subtract 3 from lineNumber.
                        if ((xp1Electrodes(j) == lineNumber - 3)) {
                            movieStart << "O     " << positions(i, 0) << " " << positions(i, 1) << " " << positions(i, 2) << "\n";
                            electrodeFound = 1;
                            break;
                        }
                    }

                    if (electrodeFound == 0) {

                        //loop through our electrode atoms.
                        for (int j = 0; j < xp2Electrodes.n_rows; j++) {
                            //see if the line we are currently on is an electrode atom. Atoms start at line
                            //3 (starting our count from 1), so subtract 3 from lineNumber.
                            if ((xp2Electrodes(j) == lineNumber - 3)) {
                                movieStart << "F     " << positions(i, 0) << " " << positions(i, 1) << " " << positions(i, 2) << "\n";
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
                            if ((xn1Electrodes(j) == lineNumber - 3)) {
                                movieStart << "N     " << positions(i, 0) << " " << positions(i, 1) << " " << positions(i, 2) << "\n";
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
                            if ((xn2Electrodes(j) == lineNumber - 3)) {
                                movieStart << "B     " << positions(i, 0) << " " << positions(i, 1) << " " << positions(i, 2) << "\n";
                                electrodeFound = 1;
                                break;
                            }
                        }
                    }

                    if (numberOfElectrodes == 4) {

                        if (electrodeFound == 0) {
                            //loop through our electrode atoms.
                            for (int j = 0; j < yp1Electrodes.n_rows; j++) {
                                //see if the line we are currently on is an electrode atom. Atoms start at line
                                //3 (starting our count from 1), so subtract 3 from lineNumber.
                                if ((yp1Electrodes(j) == lineNumber - 3)) {
                                    movieStart << "S     " << positions(i, 0) << " " << positions(i, 1) << " " << positions(i, 2) << "\n";

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
                                if ((yp2Electrodes(j) == lineNumber - 3)) {
                                    movieStart << "Y     " << positions(i, 0) << " " << positions(i, 1) << " " << positions(i, 2) << "\n";

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
                                if ((yn1Electrodes(j) == lineNumber - 3)) {
                                    movieStart << "P     " << positions(i, 0) << " " << positions(i, 1) << " " << positions(i, 2) << "\n";

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
                                if ((yn2Electrodes(j) == lineNumber - 3)) {
                                    movieStart << "I     " << positions(i, 0) << " " << positions(i, 1) << " " << positions(i, 2) << "\n";

                                    electrodeFound = 1;
                                    break;
                                }
                            }
                        }
                    }
                }
                //if current atom is not an electrode, make carbon.
                if (electrodeFound == 0) {
                    movieStart << "C     " << positions(i, 0) << " " << positions(i, 1) << " " << positions(i, 2) << "\n";

                }

            } else {
                //not considering electrodes
                movieStart << "C     " << positions(i, 0) << " " << positions(i, 1) << " " << positions(i, 2) << "\n";

            }

        }

        movieStart.close();
        chdir("..");

    }

    //if this is the first time running, we make the last accepted run the first given position.
    if (isRestart == false) {
        chdir("movies");
        //system("ls");
        system(("cp movie" + batchName + ".xyz lastAcceptedRun.xyz").c_str());
        chdir("..");
    }
    //cout << neighborList << endl;
    writeTopFile(positions, bondList, batchName);
    writePgnFile(positions, batchName);
    writePdbPsfFiles(batchName);

    if (forceNum == 1) {
        mat xValues = positions.col(0);
        double minWallPos = xValues.min() + atof(variableReader("minWallPosModifier").c_str());
        double maxWallPos = xValues.max() + atof(variableReader("maxWallPosModifier").c_str());
        double wallForce = atof(variableReader("wallForce").c_str());
        double tubeRadius = atof(variableReader("tubeRadius").c_str());
        double tubeForce = pow(atof(variableReader("tubeForceBase").c_str()), atof(variableReader("tubeForceExponent").c_str())); //force in kcal/mol*A. Should be 7^12.
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

    //start here if you want to start with sp2 to sp3 bond
    //Starts log file, minimizes starting structure, and updates movie file and histogram
    //initialize Arrhenius

    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    string timeS = asctime(timeinfo);

    stringstream ssTemp;
    ssTemp << temperature;
    string TempS = ssTemp.str();

    stringstream ssRatioRotate;
    ssRatioRotate << ratioRotate;
    string RRS = ssRatioRotate.str();
    stringstream ssratioDimerDisplacement;
    ssratioDimerDisplacement << ratioDimerDisplacement;
    string RFS = ssratioDimerDisplacement.str();
    stringstream ssRatioSp2ToSp3;
    ssRatioSp2ToSp3 << ratioSp2ToSp3;
    string RSPS = ssRatioSp2ToSp3.str();




    if (isRestart == false) {
        //"clean start," in this case, refers to if we can get the initial energy
        // value without crashing.

        bool cleanStart = false;
        while (cleanStart == false) {
            cout << "---------------------------------------------------------" << endl;
            cout << "RUNNING NAMD FOR INITIALIZATION VALUES, POSITIONS NOT KEPT" << endl;


            if(relaxationMethod == 1){
            //NAMD selected
                system(("namd2 +p" + numProcessors + " " + batchName + ".conf >> " + batchName + ".log ").c_str());
            }else{
                //TB selected
                
            }
            
            energyReader(newEnergy, batchName, relaxationMethod);



            if (newEnergy >= energyFailureTolerance) {
                // new and old energies are the same for the first step.
                oldEnergy = newEnergy;
                cleanStart = true;

                cout << "CLEAN START ACHIEVED " << endl;

                string runNumS;
                stringstream ssRun;
                ssRun << namdRunNum;
                runNumS = ssRun.str();

                string successfulRunNumS;
                stringstream ssSRun;
                ssSRun << acceptedRunNum;
                successfulRunNumS = ssSRun.str();
                stringstream ssOldEnergy;
                ssOldEnergy << oldEnergy;
                string oldEnergyS = ssOldEnergy.str();

                if (isRestart == false) {
                    //record energy
                    chdir("energies");
                    system(("echo '" + runNumS + " " + oldEnergyS + "' >> energyHist" + batchName + ".dat").c_str());
                    system(("echo '" + successfulRunNumS + " " + oldEnergyS + "' >> energyHistMovie" + batchName + ".dat").c_str());
                    chdir("..");
                }
                cleanStart = true;
                cout << "---------------------------------------------------------" << endl;
            } else {
                cout << "FAILED AT START" << endl;
            }

        }
    } else {

        //in this case it is a system restart, and we can read the last accepted energy from the last successful log file.
        chdir("movies");
        energyReader(oldEnergy, batchName, relaxationMethod);
        chdir("..");

        //we might as well find the number of accepted and total runs at this point too.
        chdir("energies");
        runReader(namdRunNum, 0, batchName);
        runReader(acceptedRunNum, 1, batchName);
        globalChiralAngleRunNum = acceptedRunNum; //this is the same as the number of accepted steps
        globalRingSizeRunNum = acceptedRunNum; //this is the same as the number of accepted steps
        chdir("..");

    }


    for (int relaxationSteps = 0; relaxationSteps < atoi(variableReader("relaxationSteps").c_str()); ++relaxationSteps) {

        //relax by number of initilization relaxations steps
        //rotateRatio = 0;
        ratioDimerDisplacement = 0;
         spPossibility = 0;
        NAMDrun(0, oldEnergy, temperature, positions, bondList, neighborList, ratioRotate, ratioDimerDisplacement, ratioSp2ToSp3, spPossibility, batchName);

    }



    /*
    string runNumS;
    stringstream ssRun;
    ssRun << namdRunNum;
    runNumS = ssRun.str();

    string successfulRunNumS;
    stringstream ssSRun;
    ssSRun << acceptedNum;
    successfulRunNumS = ssSRun.str();

    string oldEnergyS;
    stringstream ssEner;
    ssEner << oldEnergy;
    oldEnergyS = ssEner.str();


    chdir("energies");
    system(("echo '" + runNumS + " " + oldEnergyS + "' >> energyHist" + batchName + ".dat").c_str());
    system(("echo '" + successfulRunNumS + " " + oldEnergyS + "' >> energyHistMovie" + batchName + ".dat").c_str());
    chdir("..");
     */

    
    
    
    
    
    
     if (numberOfElectrodes != 0 && ratioDimerDisplacement != 0) {
        system("read -p \"Cannot have nonzero electrodes AND dimer displacement at same time in current code!\" key");

    }

}


int systemATAC::runSimulation() {

    for (int numberOfRuns = 0; numberOfRuns < atof(variableReader("numberOfRuns").c_str()); numberOfRuns++) {

        //Switches on the ability of sp2 to sp3 transitions to be reversable.
        if (reversableSpRunNum == numberOfRuns) {
            spDenialPathSize = 2;
        }

        NAMDrun(mutationGo, oldEnergy, temperature, positions, bondList, neighborList, ratioRotate, ratioDimerDisplacement, ratioSp2ToSp3, spPossibility, batchName);

    }

    bool findChiralIndexesBasedOnDiameter = atoi(variableReader("findChiralIndexesBasedOnDiameter").c_str());

    if (findChiralIndexesBasedOnDiameter == 1) {
        chdir("chiralAngles");

        findChiralityBasedOnDiameter(positions, bondList, neighborList, haveNonInteracting);
        chdir("..");
    }

    bool findChiralIndexesBasedOnAngles = atoi(variableReader("findChiralIndexesBasedOnAngles").c_str());

    if (findChiralIndexesBasedOnAngles == 1) {

        //removenon-interacting atoms to our system
        if (haveNonInteracting == true) {
            nonInteractingDeAppender(positions, bondList, neighborList);
        }

        vector<pair<imat, int> > rings;
        rings = findRings(neighborList, ringSizeMin, ringSizeMax);

        //if we don't have a single dat file, we
        ///   if (findRingSizes == 0) {
        //      recordRingSizes(rings);
        //}

        //using to generate the final angle.dat file and for the final chiral to index conversion.

        recordChirality(rings, positions, true);

        //add non-interacting atoms to our system
        if (haveNonInteracting == true) {
            nonInteractingAppender(positions, bondList, neighborList);
        }
    }
    //  cout << "GENERAGTE RING HISTOGRAM IMAGAE AND CHIRAL " <<generateRingHistogramImages <<" " <<findChiralAngles << endl;
    if (generateRingHistogramImages == 1) {

        //generates movie of histograms
        if (findChiralAngles == true) {
            chdir("chiralAngles");
            //ffmpeg note %05d splices together images with the format 00000, 00001, etc.:
            // system("ffmpeg -i pic%05d.jpg movie.mpeg");

            system("ffmpeg -f image2 -r 5 -i ./pic%05d.jpg -b 600k ./00out.mp4");
            //-r is the frames per second.
            chdir("..");

        }

        if (findRingSizes == true) {
            chdir("ringSizes");
            //ffmpeg note %05d splices together images with the format 00000, 00001, etc.:
            // system("ffmpeg -i pic%05d.jpg movie.mpeg");
            system("ffmpeg -f image2 -r 5 -i ./pic%05d.jpg -b 600k ./00out.mp4");
            //system("ffmpeg -r 0.1 -f image2 -i pic%05d.jpg movie.mpeg");
            chdir("..");
        }

        /*

        //Simulates movie
        chdir("movies");
        system(("vmd movie" + batchName + ".xyz").c_str());
        chdir("..");

        //Simulates final form
        writeTopFile(positions, bondList, batchName);
        writePgnFile(positions, batchName);
        writePdbPsfFiles(batchName);
        system(("vmd " + batchName + ".psf " + batchName + ".pdb").c_str());

        //executes energy plot code
        chdir("energies");
        system(("gnuplot -persist plot" + batchName + ".plt").c_str());
        chdir("..");



    //*/
    }


    cout << "BATCH IS FINISHED! " << endl;
    return 0;
}