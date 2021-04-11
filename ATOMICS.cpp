//Ahtr Thermal behaviOr Modeling and Iterative Criticality Suite (ATOMICS)
//Made for FHR-AHTR
//Features supporting multiphysics and control rod movement
//Written by Kyle Ramey, Georgia Tech

//Initializing libraries to be used
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;  //Indicates standard naming conventions within C++

//Functions to be used:
void startup (void);  // Splash function
void ReadUserInput (void); //Reads in user input parameters from options.txt
void TriMeshMapping (void); // Mapping of triangular mesh bins to script naming convention from TriMap.txt
void HexMapping (void); // Mapping from script hex to mesh hex bins (for symmetric results) from HexMap.txt for use in radial power profile (RadialPower.txt)
void TriMeshExtraction (void); //Extract results from SERPENT output for script (input.txt_det#.m)
void PrintRadialPower (void); //Print radial power output file RadialPower.txt
void PrintAxialPower (void); //Print axial power output file AxialPower.txt
void CRScheduleRead (void); //Read-in the CR schedule file CRSchedule.txt
void resmRead (void); //Obtain the eigenvalue from the previous Serpent simulation from input.txt_res.m
void depRead (void); //Read-in the depletion schedule user input file dep.txt
void iterationIn (void); //Read-in the iteration file iteration.txt
void iterationOut (void); //Write-out to the iteration file iteration.txt (by appending results)
void CRScheduleWrite (void); //Done only if running the CR Schedule Search Mode; finds highest CR positions
void CRSearchMode (void); //Find the next guess for CB movement or accept current insertion and continue to next iterative step
void RodWithdrawSchedule (void); // Order in which control rod banks are moved
void RodInsertionAmount (void); // Sets the proper insertion amount for CRs
void THTemps (void); // Assign temperatures to materials
void THTempsPlot (void); //Create discretized temperature profile THChannelProfile.txt for hottest location
void dims (void); // Thermally expanded dimensions
void XYLocs (void); //Locations of geometric stuff
void densities (void); // Thermally expanded densities
void materials (void); //Materials output file materials.txt
void geoOut (void);  // Outputs the geometry and tallies to geometry.txt

//Global variables to be used:
char userIn[] = "options.txt";
char trimapIn[] = "TriMap.txt";
char hexmapIn[] = "HexMap.txt";
char triresultsIn[] = "input.txt_det0.m";
char tridepIn[] = "input.txt_det1.m";
char CRScheduleFile[] = "CRSchedule.txt";
char withdrawFile[] = "CRWithdrawal.txt";
char resmFile[] = "input.txt_res.m";
char depFile[] = "dep.txt";
char iterationFile[] = "iteration.txt";
char radialpowerOut[] = "RadialPower.txt";
char axialpowerOut[] = "AxialPower.txt";
char serptoolOut[] = "ST_det0.m";
char THProfileOut[] = "THChannelProfile.txt";
char materialsOut[] = "materials.txt";
char geometryOut[] = "geometry.txt";

//Module modes
bool CRScheduleSearch = false; //Mode for finding where rods should be inserted (highest peaking)
bool CRWithdrawal = false; //Mode for withdrawal schedule
bool CRCritSearch = false; //Mode for criticality search
bool THMode = false; //Mode for thermal hydraulic iterations
bool depMode = false; //Mode for depletion

// Options
int userCodeMode = 0;
bool coldDim = true; //For Hot/Cold Dimensions
bool BPUsage = true; //Use burnable poison spheres or not
bool reflHole = true; //Radial reflector 2cm hole in the middle or not
bool printTH = true; //Whether the TH profile for highest power channel is printed
bool heteroFuelStripe = true; //Whether the heterogeneous fuel stripe temperature distribution is used (versus homogeneous)
bool gFluence = false; // Whether to use fluence-dependent thermal conductivity and thermal expansion coefficient for graphite
bool equilXenon = true; //Whether to use the equilibrium xenon feature for fissile burnable materials in a depletion sequence

//These are mostly for testing (faster)
bool cubic = true; //Cubic TRISO
bool singleFuel = true; //Same fuel material everywhere (faster for testing)
bool singleTrisoL = true;
bool singleGraph  = true;
bool singleCB = true;
bool singleBP = true; //Same BP material everywhere (faster for testing)
bool singleCL = true; //Just one axial core layer region but it spans the entire core height (faster for testing)
bool singleFlibe = true;

// Axial and expansion options
bool axialFlower = false; //False: prismatic by averaging top/bottom plate expansions. True: linearization between top and bottom independent free plate expansions.

// Functions for section mapping
int SectionMap[84][3][3][3];
int HexMap[19][19];

// Functions for power
double SectionPower[102][84][3][2];
int PeakAssemblyGroup; //assm group with the peak power
int PeakCRGroup; //for in case the peak is where a rod currently is
int PeakSectionGroup[2];
double PeakAssemblyPPF; //peak assm group PPF
double PeakCRPPF; //peak for non-CR group (only for CR search)
double PeakSectionPPF;
double PeakLASPPF;
int MinCRGroup;
int MinCRSectionGroup[2];
double MinCRGroupPPF = 0;
double MinCRSectionPPF = 0;
double RadialPower[19][19];
double WholeAssm[84][2];
double RPError[19][19];

//CR Placement
bool controlBladeIn [102][84] = {false};
int CRScheduleArray[84]; //Ordered list of CR insertion groups
int CRGroupNumb = 0; //Number of assm groups in the CRSchedule.txt file (first line of said file)
int CRActiveNumb = 0; //The ith CR group actively being moved (WRT CRScheduleArray)
int CRPartialNumb = 0; //Number of 1/16ths inserted from the top in the active group number
int CRActivePrev = 0; //Previous step active assembly
int CRPartialPrev = 0; //Previous step partial amount
int withdrawNumb = 0;
int withdrawArray[84];

//Initialize counters and global constants frequently used
int assm=0; //assembly #
int sect=0; //Assembly section #
int layr=0; //Vertical layer #
int plnk=0;
int colm=0;
int strp=0;
int radi=0; //Radial section number (for tally collapsing)
int nL = 3; //Vertical layer model # 16+2
int nLt = 18; //Vertical layer tally # 16+2
const int nA = 84; // Assembly # 84
const int nS = 3; //Assembly section # 3
int nTL = 4; //TRISO Layers # 4
int nTW = 202; //TRISO Width # 202
int nBPC = 5; //Number of BP Columns # 5
int nBPW = 40; //BP pitch separation (in integer number of BP axial pitches) # 40
double const pi=4*atan(1);
double const sqrt3=sqrt(3);

//Thermal expansion data/parameters
double alphag = 5E-6; //graphite
double alphacc = 5E-6; //Carbon-carbon composite
double alphasic = 5E-6; //Silicon carbide
double alphaf = 7.6E-6; //Fuel
double alphab = 5.5E-6; //Buffer
double alphaipyc = 5.5E-6; //Inner pyrolytic carbon
double alphaopyc = 5.5E-6; //Outer pyrolytic carbon
double alphacb = 4.8E-6; //Control blade (since mostly Mo, used Mo alpha)
double alphaeu = 7.5E-6; //Eu3O2
double alphabc = 5E-6; //Boron Carbide @ 2.37 g/cc https://www.azom.com/properties.aspx?ArticleID=75
double alphaalloyN = 13.6E-6; //http://www.haynesintl.com/alloys/alloy-portfolio_/Corrosion-resistant-Alloys/hastelloy-n-alloy/physical-properties
double alpha800H=17.3E-6; //http://www.specialmetals.com/assets/smc/documents/alloys/incoloy/incoloy-alloys-800h-800ht.pdf
double B[18][84][3]; //TE fuel matrix radial compensator for fixed axial expansion

//Temperatures
bool firstTH = true;
int THiters = 0; // Number of TH iterations
int THstep = 0; //Current TH iteration step
double corePowDens = 1.953376E-01; //Total Core Thermal Power [kW/g]
double coreMFR = 26750; //Total Core Mass Flow Rate [kg/s]
double coldT = 293; //Cold temperature, 20 C (293 K)
double inletT = 923; //650 C
double outletT = 973; //700 C
double averageT = 0; //Average of inlet and outlet temperature
double daverageT = 0; //Average expansion temperature
double Tf [102][84][3]; //Fuel
double Tb [102][84][3]; //Buffer
double Tipyc [102][84][3]; //Inner pyrolytic carbon
double Tsic [102][84][3]; //Sic
double Topyc [102][84][3]; //Outer pyrolytic carbon
double Tmat [102][84][3]; //Matrix graphite
double Tmatave [102]; //Matrix average
double Tmeat [102][84][3]; //Graphite meat
double Tsleeve [102][84][3];
double Tg [102][84][3];
double Tspacer [102][84][3];
double Tcc [102][84][3]; //Carbon-Carbon composite
double Tcb [102][84][3]; //Control Blade
double Teu [102][84][3]; //Eu3O2
double Tflibe [102][84][3]; //Flibe
double Tcoolantboundary [102][84][3]; //Temperature of FLiBe between axial sections (all other temps taken as axial section average)
double THetero[102][84][3][12][11];
double Tregion[102][84][3][12][6];
double Taxialave [102]; // Axial average coolant temperature (for reflector assemblies)
double hotSection [3]; //Location of the highest power section in model
double dTf [102][84][3]; //Fuel
double dTb [102][84][3]; //Buffer
double dTipyc [102][84][3]; //Inner pyrolytic carbon
double dTsic [102][84][3]; //Sic
double dTopyc [102][84][3]; //Outer pyrolytic carbon
double dTmat [102][84][3]; //Matrix graphite
double dTmatave [102]; //Matrix average
double dTmeat [102][84][3]; //Graphite meat
double dTsleeve [102][84][3];
double dTg [102][84][3]; //Graphite meat
double dTspacer [102][84][3];
double dTcc [102][84][3]; //Carbon-Carbon composite
double dTcb [102][84][3]; //Control Blade
double dTeu [102][84][3]; //Eu3O2
double dTaxialave [102];
double dTAxial [2];

//Cold Dimensions
double dimColdFuelR=0.02135;
double dimColdBufferR=0.03135;
double dimColdIpycR=0.03485;
double dimColdSicR=0.03835;
double dimColdOpycR=0.04235;
double dimColdZPitch=0.09266;
double dimColdXPitchAsym=0.09406;
double dimColdYPitchAsym=0.09128;
double dimColdEuR=0.035;
double dimColdZEuPitch=0.09936;
double dimColdPlankWidth = 2.55;
double dimColdSleeveWidth = 0.1;
double dimColdSpacer=0.7;
double dimColdSpacerSeparation=14;
double dimColdReflHoleRadius=2;
double dimColdReflHex=22.5;
double dimColdAxReflHeight=25;
double dimColdAxSuppPlateHeight=35;
double dimColdAssemblyApothem=22.5;
double dimColdAssemblyPitch=46.8;
double dimColdActiveCoreHeight=550;
double dimColdPermRadReflRadius=478;
double dimColdBoronCarbideRadius=479;
double dimColdCoreBarrelRadius=481;
double dimColdDowncomerRadius=519;
double dimColdAlloyNLinerRadius=520;
double dimColdPressureVesselRadius=525;

//Active Core Dimensions
double dimLayrHeight[102];
double dimActiveCoreHeight=550.02976;
double dimModelHeight=670;
double dimAxSuppPlateHeight[2];
double dimFuelRadius[102][84][3];
double dimBufferRadius[102][84][3];
double dimIpycRadius[102][84][3];
double dimSicRadius[102][84][3];
double dimOpycRadius[102][84][3];
double dimXpitch[102][84][3];
double dimYpitch[102][84][3];
double dimZpitch[102];
double dimZEupitch[102];
double dimSpacer[102][84][3];
double dimXLength[102][84][3];
double dimYWidth[102][84][3];
double dimEuRadius[102][84][3];
double dimWrapperInner[102][84][3];
double dimWrapperOuter[102][84][3];
double dimFlibeHex[102];
double dimAssemblyPitch[102];
double dimReflHoleRadius=2;
double dimReflHex=22.5;
double dimPermReflRadius=478;
double dimBoronCarbideRadius=479;
double dimCoreBarrelRadius=481;
double dimDowncomerRadius=519;
double dimAlloyNLinerRadius=520;
double dimPressureVesselRadius=525;

//Locations
double locSP[17][102][84][3];
double lociL[102][84][3];
double lociR[102][84][3];
double locBP[7][102][84][3];
double locCP[7][102][84][3];
double locPP[13][102][84][3];
double locSL[13][102][84][3];
double locFy[13][102][84][3];
double locFxL[13][102][84][3];
double locFxR[13][102][84][3];
double locFLx[13][102][84][3];
double locFLy[13][102][84][3];
double locSpL[8][102][84][3];
double locSpR[8][102][84][3];
double locBPC[7][230][102][84][3];

//Cold Densities
double denColdFuel=10.9;
double denColdBuffer=1.0;
double denColdIPyc=1.9;
double denColdSic=3.1;
double denColdOpyc=1.87;
double denColdGraphite=1.75;
double denColdCCC=1.95; //Carbon-Carbon composite
double denColdEu=1.25;
double denColdCB=10.28;
double denColdBC=2.37;
double denColdAlloyN=8.93;
double denColdH800=7.92;

//Densities
double denFuel[102][84][3];
double denBuffer[102][84][3];
double denIpyc[102][84][3];
double denSic[102][84][3];
double denOpyc[102][84][3];
double denMat[102][84][3];
double denEu[102][84][3];
double denMeat[102][84][3];
double denSleeve[102][84][3];
double denPlank[102][84][3];
double denSpacer[102][84][3];
double denCB[102][84][3];
double denCCC[102][84][3];
double denflibe[102][84][3];
double denaxialavemat [102];
double denaxialavegraph [102];
double denaxialaveflibe [102];
double denbarrel = 1.75;
double denBC = 2.37;
double denAlloyN = 8.93;
double denH800 = 7.92;
double denSupportPlate[2];

//Thermal Conductivities
double coldkgrph = 15;
double coldkmatrix = 15;
double coldkfuel  = 3.7;
double coldkbuff = 0.5;
double coldkipyro = 4;
double coldksic = 16;
double coldkopyro = 4;
double coldkflibe = 1;

//Coolant Properties
double flibeVisc = 0.0056; //Viscosity [Pa*s]
double flibePr = 0.0056; //Rrandtl Number
double flibecp = 2415; //Heat Capacity [J/(kg*K)]

//Burnup Related Parameters
double volumeFuel[102][84][3]; //Volume of fuel kernels (variable due to thermal expansion)
double volumeBP[102][84][3]; //Volume of europia BP particles
int BPBurnZones = 1; //Number of Europia BP burnable zones
int depType = 1; //Type of depletion feature from SERPENT
int depStepTotal = 1; //Total number of depletion steps for run
int depStep = 0; //Active depletion step
double depArray[100]; //Array for desired depletion history
int CRMoveStep = 0; //Number of CR movements in current depStep
bool nextDepStep = false; //whether to move to the next dep step or move CBs again
bool firstRun = false; //If the first time the code is running, certain files won't be available. Alternate methods are used.
int initialCRGroups = 0; //Initial guess at how many CB groups to insert.
double kprev = 1; //Previous iteration eigenvalue
double keff = 0; //Eigenvalue of previous step
double kerr = 1; //Eigenvalue error of previous step
double ktarget = 1; //Target value for eigenvalue during Control Blade movement
double kepsilon = 0.001; //Eigenvalue tolerance for converging on k
long int depActiveCycles = 20;
long int depInactiveCycles = 20;
long int depCycleSize = 1000;
long int kActiveCycles = 20;
long int kInactiveCycles = 20;
long int kCycleSize = 1000;

ifstream user_file; //Allows for the user-specified input file to be read
ifstream trimap_file;  //Reads in mapping from serpent-output hex notation to script indexing
ifstream hexmap_file; //Inverse mapping (script to hex)
ifstream triresults_file; //SERPENT output from using the triangular mesh tally
ifstream CRScheduleIn_file; //Reads schedule of which CB groups are inserted
ifstream iterationIn_file; //Read in entries from the power iteration file
ifstream resm_file; //.m Serpent output file reading (for k and other parameters)
ifstream dep_file; // List of burnup steps to be used for depletion sequence
ifstream withdrawIn_file; //Reads CB withdrawal schedule
ofstream CRScheduleOut_file; //Writes file for which CB groups are inserted (legacy)
ofstream withdrawOut_file; //Writes file for which CB groups are withdrawn
ofstream radialresults_file; //Write out 19x19 radial power profile
ofstream axialresults_file; //Write out axial power profile
ofstream iterationOut_file; //Write to the power iteration file
ofstream st_file; //Normalized results to look like Serpent output again for SerpentTools use
ofstream th_file; //Optional output of thermal hydraulic data
ofstream mat_file; //Materials output file
ofstream out_file;  //Geometry output file

int main (void) { //The entire program runs in this function
    startup(); //Splash function
    ReadUserInput(); //Read user input file
    TriMeshMapping(); //Map basis for triangular mesh tally region to basis for script
    HexMapping(); //Map basis for hex mesh tally to script assembly groups
    CRScheduleRead(); //REad CB schedule
    resmRead(); //Read in _res.m file (for k)
    if (depMode){
        depRead(); //Read in depletion schedule (dictated by user)
    }
    if ((CRCritSearch) || (THMode) ) { //If in critical search mode (also part of depletion)
        iterationIn(); //Read in iteration file for key parameters
        TriMeshExtraction(); //Extract triangular mesh tally results
        if (CRCritSearch){
            if (THstep == 0){CRSearchMode();} //Search for the predicted critical geometry
                else{nextDepStep = true;}}
        iterationOut(); //Write out the iteration summary
    }
        else{
            TriMeshExtraction(); //Extract triangular mesh tally results
        }
    if (CRScheduleSearch){
        CRScheduleWrite(); //If in CB insertion search mode, modify schedule file
    }
    if (CRWithdrawal){
        RodWithdrawSchedule(); //If in CB withdrawal search mode, modify this file
    }
    RodInsertionAmount(); //Modify CB insertion array based on CB requirements
    PrintRadialPower(); //Print radial power profile
    for(int THCounter = 0; THCounter < 10; THCounter++){
        dims(); //Geometric dimensions
        XYLocs(); //Locations for surfaces
        densities(); //Densities for materials
        THTemps(); //Thermal hydraulic temperatures
        firstTH = false;
    }
    if (printTH) {
        THTempsPlot();
    }
    materials(); //Materials output in SERPENT-specific format. Also includes some neutron physics options.
    geoOut(); //Geometry output in SERPENT-specific format. Also includes some tallies for post-processing and figures.
    PrintAxialPower(); //Print axial power profile
    return 0;
}

void startup (void){ //Splash text for when the program runs
    cout << "ATOMICS - Ahtr Thermal behaviOr Modeling and Iterative Criticality Suite\n";
    cout << "          =    =             =  =            =         =           =\n\n";
}

void ReadUserInput (void) {
    double fvalue = 0;
    int ivalue = 0;
    long int livalue = 0;
    int k = 0;
    int l = 0;
    char charnumb [20];
    string wholeline;
    string component;
    user_file.open(userIn);
    if (user_file.fail()) {cout << "\nCould not open options.txt file.\n";}
        else {
            cout << "\n=== Values Extracted from the User Input File options.txt ===\n";
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            userCodeMode = ivalue;
            cout << "\nScript Mode (line " << k << "): " << userCodeMode;
            if (userCodeMode == 0) {CRScheduleSearch = false;}
                else if (userCodeMode == 1) {CRScheduleSearch = true;}
                    else if (userCodeMode == 2) {CRWithdrawal = true;}
                        else if (userCodeMode == 3) {CRCritSearch = true;}
                            else if (userCodeMode == 4){THMode = true;}
                                else if (userCodeMode == 5){
                                    CRCritSearch = true;
                                    THMode = true;}
                                    else if (userCodeMode == 6) {depMode = true;}
                                        else if (userCodeMode == 7) {
                                            CRCritSearch = true;
                                            depMode = true;}
                                            else if(userCodeMode == 8){
                                                depMode = true;
                                                THMode = true;}
                                                else if (userCodeMode == 9){
                                                    depMode = true;
                                                    CRCritSearch = true;
                                                    THMode = true;}
                                                    else {cout << "\n   ***ERROR: invalid script mode. Code may not function as intended.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            coldDim = ivalue;
            cout << "\nUse of Cold Dimensions (line " << k << "): " << coldDim;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            livalue = atoi(charnumb);
            kCycleSize = livalue;
            cout << "\nStatepoint Particles per Cycle (line " << k << "): " << kCycleSize;
            if (livalue < 0) {cout << "\n   ***ERROR: cannot run negative particles.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            livalue = atoi(charnumb);
            kActiveCycles = livalue;
            cout << "\nStatepoint Active Cycles (line " << k << "): " << kActiveCycles;
            if (livalue < 0) {cout << "\n   ***ERROR: cannot run negative cycles.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            livalue = atoi(charnumb);
            kInactiveCycles = livalue;
            cout << "\nStatepoint Inactive Cycles (line " << k << "): " << kInactiveCycles;
            if (livalue < 0) {cout << "\n   ***ERROR: cannot run negative cycles.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdFuelR = fvalue;
            cout << "\n\nCold Fuel Kernel Radius (line " << k << "): " << dimColdFuelR << " cm";
            if (dimColdFuelR < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdBufferR = fvalue;
            cout << "\nCold Buffer Radius (line " << k << "): " << dimColdBufferR << " cm";
            if (dimColdBufferR < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdIpycR = fvalue;
            cout << "\nCold Inner Pyrolytic Radius (line " << k << "): " << dimColdIpycR << " cm";
            if (dimColdIpycR < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdSicR = fvalue;
            cout << "\nCold Silicon Carbide Radius (line " << k << "): " << dimColdSicR << " cm";
            if (dimColdSicR < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdOpycR = fvalue;
            cout << "\nCold Outer Pyrolytic Radius (line " << k << "): " << dimColdOpycR << " cm";
            if (dimColdOpycR < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            cubic = ivalue;
            cout << "\nCubic Lattice (line " << k << "): " << cubic;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            nTL = ivalue;
            cout << "\nWidth of Fuel Stripe in Layers (line " << k << "): " << nTL;
            if (nTL > 12) {cout << "\n   ***ERROR: Fuel Stripe is too thick for the size of the fuel plank.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            nTW = ivalue;
            cout << "\nLength of Fuel Stripe in Particles (line " << k << "): " << nTW;
            if (nTW > 225) {cout << "\n   ***WARNING: Fuel Stripe is too wide for cuboidal lattice.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdXPitchAsym = fvalue;
            cout << "\nCold X Lattice Pitch (line " << k << "): " << dimColdXPitchAsym << " cm";
            if (cubic) {cout << "\n   ***NOTE: since cubic, value will be overwritten with z value.***";}
            if (dimColdXPitchAsym < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdYPitchAsym = fvalue;
            cout << "\nCold Y Lattice Pitch (line " << k << "): " << dimColdYPitchAsym << " cm";
            if (dimColdYPitchAsym < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            if (cubic) {cout << "\n   ***NOTE: since cubic, value will be overwritten with z value.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdZPitch = fvalue;
            cout << "\nCold Z Lattice Pitch (line " << k << "): " << dimColdZPitch << " cm";
            if (dimColdZPitch < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            BPUsage = ivalue;
            cout << "\nBurnable Poison Sphere Usage (line " << k << "): " << BPUsage;
             if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdEuR = fvalue;
            cout << "\nCold Burnable Poison Sphere Radius (line " << k << "): " << dimColdEuR << " cm";
            if (dimColdEuR < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdZEuPitch = fvalue;
            cout << "\nCold Burnable Poison Sphere Axial Pitch (line " << k << "): " << dimColdZEuPitch << " cm";
            if (dimColdZEuPitch < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            nBPC = ivalue;
            cout << "\nNumber of Burnable Poison Columns (line " << k << "): " << nBPC;
            if (ivalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            if ((ivalue%2) == 0) {cout << "\n   ***ERROR: even columns not supported. Must be odd.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            nBPW = ivalue;
            cout << "\nBurnable Poison Spacing Along Plank (line " << k << "): " << nBPW;
            if (ivalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdPlankWidth = fvalue;
            cout << "\nCold Plank Width (line " << k << "): " << dimColdPlankWidth << " cm";
            if (dimColdSpacer < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdSleeveWidth = fvalue;
            cout << "\nCold Sleeve Width (line " << k << "): " << dimColdSleeveWidth << " cm";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdSpacerSeparation = fvalue;
            cout << "\nCold Sleeve Width (line " << k << "): " << dimColdSpacerSeparation << " cm";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdAssemblyApothem = fvalue;
            cout << "\nCold Assembly Apothem (line " << k << "): " << dimColdAssemblyApothem << " cm";
            if (dimColdAssemblyPitch < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdAssemblyPitch = fvalue;
            cout << "\nCold Assembly Pitch (line " << k << "): " << dimColdAssemblyPitch << " cm";
            if (dimColdAssemblyPitch < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            reflHole = ivalue;
            cout << "\nReflector Assembly Central Cooling Hole Usage (line " << k << "): " << reflHole;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdReflHoleRadius = fvalue;
            cout << "\nCold Reflector Assembly Central Cooling Hole Radius (line " << k << "): " << dimColdReflHoleRadius << " cm";
            if (dimColdReflHoleRadius < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdReflHex = fvalue;
            cout << "\nCold Reflector Assembly Apothem (line " << k << "): " << dimColdReflHex << " cm";
            if (dimColdReflHex < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            nL = ivalue + 2; // Need to increase the index because axial reflector regions are included in value.
            cout << "\nNumber of Modeled Axial Partitions in Active Core (line " << k << "): " << nL - 2;
            if ((ivalue < 1) || (ivalue > 16)) {cout << "\n   ***ERROR: Needs to be at least 1 and no greater than 16.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdActiveCoreHeight = fvalue;
            cout << "\nCold Active Core Height (line " << k << "): " << dimColdActiveCoreHeight << " cm";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdAxReflHeight = fvalue;
            cout << "\nCold Axial Reflector Height (line " << k << "): " << dimColdAxReflHeight << " cm";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdAxSuppPlateHeight = fvalue;
            cout << "\nCold Axial Support Plate Height (line " << k << "): " << dimColdAxSuppPlateHeight << " cm";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            axialFlower = ivalue;
            cout << "\nUse of Axial Flowering Interassembly Expansion (line " << k << "): " << axialFlower;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdPermRadReflRadius = fvalue;
            cout << "\nCold Permanent Radial Reflector Outer Radius (line " << k << "): " << dimColdPermRadReflRadius << " cm";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdBoronCarbideRadius = fvalue;
            cout << "\nCold Boron Carbide Outer Radius (line " << k << "): " << dimColdBoronCarbideRadius << " cm";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdCoreBarrelRadius = fvalue;
            cout << "\nCold Core Barrel Outer Radius (line " << k << "): " << dimColdCoreBarrelRadius << " cm";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdDowncomerRadius = fvalue;
            cout << "\nCold Downcomer Outer Radius (line " << k << "): " << dimColdDowncomerRadius << " cm";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdAlloyNLinerRadius = fvalue;
            cout << "\nCold Downcomer Outer Radius (line " << k << "): " << dimColdAlloyNLinerRadius << " cm";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            dimColdPressureVesselRadius = fvalue;
            cout << "\nCold Pressure Vessel Outer Radius (line " << k << "): " << dimColdPressureVesselRadius << " cm";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            singleFuel = ivalue;
            cout << "\n\nUniform Fuel Usage (line " << k << "): " << singleFuel;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            singleTrisoL = ivalue;
            cout << "\nUniform Other TRISO Layer Usage (line " << k << "): " << singleTrisoL;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            singleGraph = ivalue;
            cout << "\nUniform Structural Graphite Usage (line " << k << "): " << singleGraph;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            singleBP = ivalue;
            cout << "\nUniform Burnable Poison Usage (line " << k << "): " << singleBP;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            singleCB = ivalue;
            cout << "\nUniform Control Blade Usage (line " << k << "): " << singleCB;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            singleFlibe = ivalue;
            cout << "\nUniform FLiBe Coolant Usage (line " << k << "): " << singleFlibe;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdFuel = fvalue;
            cout << "\nCold Fuel Density (line " << k << "): " << denColdFuel << " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdBuffer = fvalue;
            cout << "\nCold Buffer Density (line " << k << "): " << denColdBuffer << " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdIPyc = fvalue;
            cout << "\nCold Inner Pyrolytic Carbon Density (line " << k << "): " << denColdIPyc << " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdSic = fvalue;
            cout << "\nCold Silicon Carbide Density (line " << k << "): " << denColdSic << " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdOpyc = fvalue;
            cout << "\nCold Outer Pyrolytic Carbon Density (line " << k << "): " << denColdOpyc << " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdGraphite = fvalue;
            cout << "\nCold Graphite Density (line " << k << "): " << denColdGraphite<< " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdCCC = fvalue;
            cout << "\nCold Carbon-Carbon Composite Density (line " << k << "): " << denColdCCC << " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdEu = fvalue;
            cout << "\nCold Europia (BP) Density (line " << k << "): " << denColdEu << " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdCB = fvalue;
            cout << "\nCold MHC (CB) Density (line " << k << "): " << denColdCB << " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdBC = fvalue;
            cout << "\nCold Boron Carbide Density (line " << k << "): " << denColdBC << " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdAlloyN = fvalue;
            cout << "\nCold Alloy N Density (line " << k << "): " << denColdAlloyN << " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            denColdH800 = fvalue;
            cout << "\nCold Hastelloy 800 Density (line " << k << "): " << denColdH800 << " g/cm3";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alphaf = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of Fuel (line " << k << "): " << alphaf;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alphab = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of Buffer (line " << k << "): " << alphab;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alphaipyc = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of Inner Pyrolytic Carbon (line " << k << "): " << alphaipyc;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alphasic = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of Silicon Carbide (line " << k << "): " << alphasic;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alphaopyc = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of Outer Pyrolytic Carbon (line " << k << "): " << alphaopyc;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alphag = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of Graphite (line " << k << "): " << alphag;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alphacc = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of Carbon-Carbon Composite (line " << k << "): " << alphacc;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alphaeu = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of Europia (BP) (line " << k << "): " << alphaeu;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alphacb = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of MHC (CB) (line " << k << "): " << alphacb;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alphabc = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of Boron Carbide (line " << k << "): " << alphabc;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alphaalloyN = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of Alloy N (line " << k << "): " << alphaalloyN;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            alpha800H = fvalue/1000000;
            cout << "\nThermal Expansion Coefficient of Hastelloy 800 (line " << k << "): " << alpha800H;
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            coldkfuel = fvalue;
            cout << "\nThermal Conductivity of Fuel (line " << k << "): " << coldkfuel << " W/(m*K)";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            coldkbuff = fvalue;
            cout << "\nThermal Conductivity of Buffer (line " << k << "): " << coldkbuff << " W/(m*K)";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            coldkipyro = fvalue;
            cout << "\nThermal Conductivity of Inner Pyrolytic Carbon (line " << k << "): " << coldkipyro << " W/(m*K)";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            coldksic = fvalue;
            cout << "\nThermal Conductivity of Silicon Carbide (line " << k << "): " << coldksic << " W/(m*K)";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            coldkopyro = fvalue;
            cout << "\nThermal Conductivity of Outer Pyrolytic Carbon (line " << k << "): " << coldkopyro << " W/(m*K)";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            coldkgrph = fvalue;
            cout << "\nThermal Conductivity of Unirradiated Graphite (line " << k << "): " << coldkgrph << " W/(m*K)";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            coldkmatrix = fvalue;
            cout << "\nThermal Conductivity of Unirradiated Matrix (line " << k << "): " << coldkmatrix << " W/(m*K)";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            THiters = ivalue;
            cout << "\nNumber of Thermal Hydraulic Iterations (line " << k << "): " << THiters;
            if (ivalue < 0) {cout << "\n   ***ERROR: value must be nonnegative.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            corePowDens = fvalue;
            cout << "\n\nCore Power Density (line " << k << "): " << corePowDens << " kW/g";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            coreMFR = fvalue;
            cout << "\nTotal Core Mass Flow Rate (line " << k << "): " << coreMFR << " kg/s";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            coldT = fvalue;
            cout << "\nCold Component Reference Temperature (line " << k << "): " << coldT << " K";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            inletT = fvalue;
            cout << "\nCore Inlet Temperature (line " << k << "): " << inletT << " K";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            printTH = ivalue;
            cout << "\nPrinting of TH Profile for Highest Power Zone (line " << k << "): " << printTH;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            heteroFuelStripe = ivalue;
            cout << "\nUse of Particle Reconstructed Fuel Stripe Temperature Profile (line " << k << "): " << heteroFuelStripe;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            flibeVisc = fvalue;
            cout << "\nFlibe Viscosity (line " << k << "): " << flibeVisc << " Pa*s";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            flibePr = fvalue;
            cout << "\nFlibe Prandtl (line " << k << "): " << flibePr;
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            flibecp = fvalue;
            cout << "\nFlibe Heat Capacity (line " << k << "): " << flibecp << " J/(kg*K)";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            coldkflibe = fvalue;
            cout << "\nThermal Conductivity Flibe (line " << k << "): " << coldkflibe << " W/(m*K)";
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            gFluence = ivalue;
            cout << "\n\nUse of Fluence-Dependent Graphite Properties (line " << k << "): " << gFluence;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            equilXenon = ivalue;
            cout << "\nUse of Equilibrium Xenon feature for Fuel Depletion (line " << k << "): " << equilXenon;
            if ((ivalue != 0) && (ivalue != 1)) {cout << "\n   ***ERROR: boolean needs to be either 0 or 1.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            BPBurnZones = ivalue;
            cout << "\nNumber of BP Burnable Zones (line " << k << "): " << BPBurnZones;
            if ((ivalue < 1) || (ivalue > 10)) {cout << "\n   ***ERROR: Serpent only supports divisions between 1 and 10.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            ktarget = fvalue;
            cout << "\nTarget eigenvalue for Control Blade Movement (line " << k << "): " << ktarget;
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            fvalue = atof(charnumb);
            kepsilon = fvalue/100000;
            cout << "\nEigenvalue Tolerance for Control Blade Movement (line " << k << "): " << kepsilon;
            if (fvalue < 0) {cout << "\n   ***ERROR: physical values must be positive.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            ivalue = atoi(charnumb);
            initialCRGroups = ivalue;
            cout << "\nInitial Number of CB Groups to Insert (line " << k << "): " << initialCRGroups;
            if (ivalue < 0) {cout << "\n   ***ERROR: cannot insert negative CBs.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            livalue = atoi(charnumb);
            depCycleSize = livalue;
            cout << "\nDepletion Particles per Cycle (line " << k << "): " << depCycleSize;
            if (livalue < 0) {cout << "\n   ***ERROR: cannot run negative particles.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            livalue = atoi(charnumb);
            depActiveCycles = livalue;
            cout << "\nDepletion Active Cycles (line " << k << "): " << depActiveCycles;
            if (livalue < 0) {cout << "\n   ***ERROR: cannot run negative cycles.***";}
            getline(user_file, wholeline);
            k++;
            component = wholeline.substr (0,18);
            for (l = 0; l < 18; l++){charnumb[l] = component[l];}
            livalue = atoi(charnumb);
            depInactiveCycles = livalue;
            cout << "\nDepletion Inactive Cycles (line " << k << "): " << depInactiveCycles;
            if (livalue < 0) {cout << "\n   ***ERROR: cannot run negative cycles.***";}
            cout << "\n";
            user_file.close();}
}

void TriMeshMapping (void) {
    bool doneReading = false;
    int tcounter=0;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    layr = 1;
    assm = 0;
    sect = 0;
    radi = 0;
    int linnumber = 1; //Line Number
    string wholeline;
    string component;
    char xynumb [2];
    char tnumb [1];
    double value = 0;
    int tvalue = 0;
    trimap_file.open(trimapIn);
    if (trimap_file.fail()) {
        cout << "Could not open TriMap.txt file.\n";}
        else {
            do {
                getline(trimap_file, wholeline);
                k++;
                if (k == linnumber){
                    component = wholeline.substr (0,2);
                    for (l = 0; l < 2; l++){
                        xynumb[l] = component[l];}
                    value = atoi(xynumb);
                    SectionMap[assm][sect][radi][tcounter] = value;
                    tcounter++;
                    component = wholeline.substr (3,2);
                    for (l = 0; l < 2; l++){
                        xynumb[l] = component[l];}
                    value = atoi(xynumb);
                    SectionMap[assm][sect][radi][tcounter] = value;
                    tcounter++;
                    component = wholeline.substr (5,2);
                    for (l = 0; l < 2; l++){
                        xynumb[l] = component[l];}
                    value = atoi(xynumb);
                    SectionMap[assm][sect][radi][tcounter] = value;
                    linnumber++;
                    tcounter=0;
                    radi++;
                    if(radi == 3){
                        radi = 0;
                        sect++;
                        if(sect == 3){
                            sect = 0;
                            assm++;
                            if (assm == 84){doneReading = true;}}}}
            }while (!doneReading);
            trimap_file.close();}
}

void HexMapping (void){
    bool doneReading = false;
    int x = 0;
    int y = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    assm = 0;
    int linnumber = 1; //Line Number
    string wholeline;
    string component;
    char numb [2];
    int value = 0;
    for (y=0;y<19;y++){for (x=0;x<19;x++){HexMap[y][x] = 0;}} // Initialize all to zero
    hexmap_file.open(hexmapIn);
    if (hexmap_file.fail()) {
        cout << "Could not open HexMap.txt file.\n";}
        else {
            do {
                getline(hexmap_file, wholeline );
                k++;
                if (k == linnumber){
                    component = wholeline.substr (0,2);
                    for (l = 0; l < 2; l++){
                        numb[l] = component[l];}
                    value = atoi(numb);
                    y = value;
                    component = wholeline.substr (3,2);
                    for (l = 0; l < 2; l++){
                        numb[l] = component[l];}
                    value = atoi(numb);
                    x = value;
                    component = wholeline.substr (6,2);
                    for (l = 0; l < 2; l++){
                        numb[l] = component[l];}
                    value = atoi(numb);
                    assm = value;
                    HexMap[y-1][x-1] = assm;
                    linnumber++;
                    if((y == 19) && (x == 8)){doneReading = true;}}
            }while (!doneReading);
            hexmap_file.close();}
}

void TriMeshExtraction (void){
    double inArrayT[34656];//Input file array
    double inErrorT[34656];//Associated errors with each result
    double resArrayT[18][19][19][6];
    double resErrorT[18][19][19][6];
    double resArrayT2[18][19][19][3];
    double resErrorT2[18][19][19][3];
    double radialArray[18][84][3][3];
    double radialError[18][84][3][3];
    layr=1;
    assm=0;
    radi=0;
    sect=0;
    int y=0;
    int x=0;
    int t=0;
    bool doneReading = false;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int linnumber = 3; //Line Number
    string wholeline;
    string component;
    string error;
    char charnumb [20];
    char errornumb [10];
    double value = 0;
    //Check if the files exist
    if ((depMode == true) && (depStep > 0) && (CRMoveStep == 1)){triresults_file.open(tridepIn);}
        else{triresults_file.open(triresultsIn);}
    if (triresults_file.fail()) {
        cout << "\nCould not open TriResults.txt file.\nNote: This is expected during the first iteration.\n";}
        else {
            do {
                getline(triresults_file, wholeline );
                k++;
                if (k == linnumber){
                    component = wholeline.substr (57,11);
                    for (l = 0; l < 11; l++){
                        charnumb[l] = component[l];}
                    value = atof(charnumb);
                    inArrayT[i] = value;
                    error = wholeline.substr (69,7);
                    for (l = 0; l < 7; l++){
                        errornumb[l] = error[l];}
                    value = atof(errornumb);
                    inErrorT[i] = value;
                    i++;
                    linnumber++;
                    if (i == 6*19*19*(nLt-2)) {doneReading = true;}}
            }while (!doneReading);
            triresults_file.close();
            //Place into corresponding locations
            i=0;
            for (layr=1;layr<(nLt-1);layr++){
                for (y=0;y<19;y++){
                    for (x=0;x<19;x++){
                        for (t=0;t<6;t++){
                            resArrayT[layr][y][x][t] = inArrayT[i];
                            resErrorT[layr][y][x][t] = inErrorT[i];
                            i++;}}}}
            //Combine results from triangles into sections
            for (layr=1;layr<(nLt-1);layr++){
                for (y=0;y<19;y++){
                    for (x=0;x<19;x++){
                        resArrayT2[layr][y][x][0] = resArrayT[layr][y][x][0] + resArrayT[layr][y][x][5];
                        resArrayT2[layr][y][x][1] = resArrayT[layr][y][x][3] + resArrayT[layr][y][x][4];
                        resArrayT2[layr][y][x][2] = resArrayT[layr][y][x][1] + resArrayT[layr][y][x][2];
                        if (resArrayT2[y][x][0] == 0){resErrorT2[layr][y][x][0] = 0;}
                            else{resErrorT2[layr][y][x][0] = sqrt(resArrayT[layr][y][x][0]*resErrorT[layr][y][x][0]*resArrayT[layr][y][x][0]*resErrorT[layr][y][x][0] + resArrayT[layr][y][x][5]*resErrorT[layr][y][x][5]*resArrayT[layr][y][x][5]*resErrorT[layr][y][x][5])/(resArrayT2[layr][y][x][0]);}
                        if (resArrayT2[y][x][1] == 0){resErrorT2[layr][y][x][1] = 0;}
                            else{resErrorT2[layr][y][x][1] = sqrt(resArrayT[layr][y][x][3]*resErrorT[layr][y][x][3]*resArrayT[layr][y][x][3]*resErrorT[layr][y][x][3] + resArrayT[layr][y][x][4]*resErrorT[layr][y][x][4]*resArrayT[layr][y][x][4]*resErrorT[layr][y][x][4])/(resArrayT2[layr][y][x][1]);}
                        if (resArrayT2[y][x][2] == 0){resErrorT2[layr][y][x][2] = 0;}
                            else{resErrorT2[layr][y][x][2] = sqrt(resArrayT[layr][y][x][1]*resErrorT[layr][y][x][1]*resArrayT[layr][y][x][1]*resErrorT[layr][y][x][1] + resArrayT[layr][y][x][2]*resErrorT[layr][y][x][2]*resArrayT[layr][y][x][2]*resErrorT[layr][y][x][2])/(resArrayT2[layr][y][x][2]);}}}}
            //Use mapping to go from tally grid to script locations
            for (layr=1;layr<(nLt-1);layr++){
                for (assm=0;assm<84;assm++){
                    for (sect=0;sect<3;sect++){
                        for (radi=0;radi<3;radi++){
                            radialArray[layr][assm][sect][radi]=resArrayT2[layr][SectionMap[assm][sect][radi][0]-1][SectionMap[assm][sect][radi][1]-1][SectionMap[assm][sect][radi][2]-1];
                            radialError[layr][assm][sect][radi]=resErrorT2[layr][SectionMap[assm][sect][radi][0]-1][SectionMap[assm][sect][radi][1]-1][SectionMap[assm][sect][radi][2]-1];
                            }}}}
            //Combine radial symmetric results
            for (layr=1;layr<(nLt-1);layr++){
                for (assm=0;assm<84;assm++){
                    for (sect=0;sect<3;sect++){
                        SectionPower[layr][assm][sect][0] = radialArray[layr][assm][sect][0] + radialArray[layr][assm][sect][1] + radialArray[layr][assm][sect][2];
                        if (SectionPower[layr][assm][sect][0] == 0){SectionPower[layr][assm][sect][1]=0;}
                            else{SectionPower[layr][assm][sect][1] = sqrt(radialArray[layr][assm][sect][0]*radialArray[layr][assm][sect][0]*radialError[layr][assm][sect][0]*radialError[layr][assm][sect][0] + radialArray[layr][assm][sect][1]*radialArray[layr][assm][sect][1]*radialError[layr][assm][sect][1]*radialError[layr][assm][sect][1] + radialArray[layr][assm][sect][2]*radialArray[layr][assm][sect][2]*radialError[layr][assm][sect][2]*radialError[layr][assm][sect][2])/SectionPower[layr][assm][sect][0];}}}}
            //Integral axial results
            double secttot = 0;
            double coretot = 0;
            double IntegralSection[84][3];
            for (assm=0;assm<84;assm++){
                for (sect=0;sect<3;sect++){
                    for (layr=1;layr<(nLt-1);layr++){
                        secttot += SectionPower[layr][assm][sect][0];}
                    IntegralSection[assm][sect] = secttot;
                    coretot += secttot;
                    secttot = 0;}}
            //Find Peak Axial Section
            int maxAsect[3];
            double maxAS = 0;
            for (layr=1;layr<(nLt-1);layr++){
                for (assm=0;assm<84;assm++){
                    for (sect=0;sect<3;sect++){
                        if (SectionPower[layr][assm][sect][0] > maxAS){
                            maxAS = SectionPower[layr][assm][sect][0];
                            maxAsect[0] = layr;
                            maxAsect[1] = assm+1;
                            maxAsect[2] = sect+1;}}}}
            PeakLASPPF = maxAS*252*(nL-2)/coretot;
            hotSection[0] = maxAsect[0];
            hotSection[1] = maxAsect[1];
            hotSection[2] = maxAsect[2];
            //Find Peak Integral Section
            int maxsect[2];
            int minCRsect[2];
            double maxval = 0;
            double minCRval = 1E38;
            for (assm=0;assm<84;assm++){
                for (sect=0;sect<3;sect++){
                    IntegralSection[assm][sect] *= 252/coretot;
                    if (IntegralSection[assm][sect] > maxval){
                        maxval = IntegralSection[assm][sect];
                        maxsect[0] = assm;
                        maxsect[1] = sect;}
                    if (IntegralSection[assm][sect] < minCRval){
                        for (i=0;i<CRGroupNumb;i++){
                            if (assm == CRScheduleArray[i]){
                                minCRval = IntegralSection[assm][sect];
                                minCRsect[0] = assm;
                                minCRsect[1] = sect;}}}}}
            PeakSectionGroup[0] = maxsect[0];
            PeakSectionGroup[1] = maxsect[1];
            PeakSectionPPF = maxval;
            MinCRSectionGroup[0] = minCRsect[0];
            MinCRSectionGroup[1] = minCRsect[1];
            MinCRSectionPPF = minCRval;
            //Integral assembly results
            double radtot = 0;
            double SPerror = 0;
            coretot = 0;
            for (assm=0;assm<84;assm++){
                for (layr=1;layr<(nLt-1);layr++){
                    for (sect=0;sect<3;sect++){
                        radtot += SectionPower[layr][assm][sect][0];
                        SPerror += SectionPower[layr][assm][sect][0]*SectionPower[layr][assm][sect][0]*SectionPower[layr][assm][sect][1]*SectionPower[layr][assm][sect][1];}}
                WholeAssm[assm][0] = radtot;
                coretot += radtot;
                radtot = 0;
                WholeAssm[assm][1] = sqrt(SPerror)/WholeAssm[assm][0];
                SPerror = 0;}
            //Find Peak Assembly
            int maxassm = 0;
            int minCRassm = 0;
            maxval = 0;
            minCRval = 1E38;
            int maxCRassm = 0;
            double CRval = 0;
            bool alreadyIn = false;
            for (assm=0;assm<84;assm++){
                WholeAssm[assm][0] *= 84/coretot;
                if (WholeAssm[assm][0] > maxval){
                    maxval = WholeAssm[assm][0];
                    maxassm = assm;}
                if (CRScheduleSearch){  //If CR schedule searching, cannot place rods in location where already present. Search for next highest location in this instance.
                    for (i=0;i<CRGroupNumb;i++){if (assm == CRScheduleArray[i]){alreadyIn = true;}} //Already in, do not select
                    if ((!alreadyIn) && (WholeAssm[assm][0] > CRval)){
                        CRval = WholeAssm[assm][0];
                        maxCRassm = assm;}
                    alreadyIn = false;}
                if (WholeAssm[assm][0] < minCRval){
                        for (i=0;i<CRGroupNumb;i++){
                            if (assm == CRScheduleArray[i]){
                                minCRval = WholeAssm[assm][0];
                                minCRassm = assm;}}}}
            PeakAssemblyGroup = maxassm;
            PeakCRGroup = maxCRassm;
            PeakAssemblyPPF = maxval;
            PeakCRPPF = CRval;
            if (!CRScheduleSearch){ //Overwrite CR values if not in search mode
                PeakCRGroup = PeakAssemblyGroup;
                PeakCRPPF = PeakAssemblyPPF;}
            MinCRGroup = minCRassm;
            MinCRGroupPPF = minCRval;
            //Unsymmetric radial results (unsure if will be used ever)
            double errorsum = 0;
            for (y=0;y<19;y++){
                for (x=0;x<19;x++){
                    RadialPower[y][x] = 0;
                    RPError[y][x] = 0;
                    for (layr=1;layr<(nLt-1);layr++){
                        for (t=0;t<6;t++){
                        RadialPower[y][x] += resArrayT[layr][y][x][t];
                        errorsum += resArrayT[layr][y][x][t]*resArrayT[layr][y][x][t]*resErrorT[layr][y][x][t]*resErrorT[layr][y][x][t];
                        }}
                    RPError[y][x] = sqrt(errorsum)/RadialPower[y][x];
                    }}
            double powersum = 0;
            for (y=0;y<19;y++){
                for (x=0;x<19;x++){
                    powersum += RadialPower[y][x];}}
            for (y=0;y<19;y++){
                for (x=0;x<19;x++){
                    RadialPower[y][x] *= 252/powersum;}}
        }
}

void PrintRadialPower (void){
    //Find Peak/Average/Min uncertainties in assemblies
    int minErrAss = 0;
    int maxErrAss = 0;
    double minEA = 1E38; //Arbitrarily large and will be overridden
    double maxEA = 0;
    double sumEA = 0;
    double avgEA = 0;
    double absErr = 0;
    int i = 0;
    for (i=0;i<nA;i++){
        absErr = WholeAssm[i][0]*WholeAssm[i][1];
        sumEA += absErr;
        if (absErr > maxEA){
            maxEA = absErr;
            maxErrAss = i;}
        if (absErr < minEA){
            minEA = absErr;
            minErrAss = i;}}
    avgEA = sumEA/nA;
    //This is for outputting human-readable PPF and relevant data
    radialresults_file.open(radialpowerOut);
    int x = 0;
    int y = 0;
    for (y=0;y<19;y++){ //Might want to invert these indices; look into this later (y-flip)
        for (x=0;x<19;x++){
            if (HexMap[y][x] == 0){radialresults_file << "0\t";}
                else {radialresults_file << WholeAssm[HexMap[y][x] - 1][0] << "\t";}}
        radialresults_file << endl;}
    radialresults_file << "\nPeak Assembly Group: " << PeakAssemblyGroup+1 << "\tPower Peaking Factor: " << PeakAssemblyPPF;
    if (PeakAssemblyGroup != PeakCRGroup){
        radialresults_file << "\nPeak Uncontrolled Group: " << PeakCRGroup+1 << "\tPower Peaking Factor: " << PeakCRPPF << " *Notice: CR already inserted in peak position. While likely persist for remainder of search.* ";}
    radialresults_file << "\nPeak Section Group: " << PeakSectionGroup[0]+1 << "-" << PeakSectionGroup[1]+1 << "\tPower Peaking Factor: " << PeakSectionPPF;
    radialresults_file << "\nMaximum Assembly Error: " << maxEA << "\tMinimum Assembly Error: " << minEA << "\tAverage Assembly Error: " << avgEA;
    radialresults_file.close();
    cout << "\nThe file " << radialpowerOut << " has been written with an assembly-wise power distribution and related parameters.\n";
    //This is for creating normalized SerpentTools plots
    int stc = 0; //st counter
    st_file.open(serptoolOut);
    st_file << "\nDEThex = [\n";
    for (y=1;y<20;y++){
        for(x=1;x<20;x++){
            stc++;
            if (stc < 10){st_file << "    ";}
                else if (stc < 100){st_file << "   ";}
                    else{st_file << "  ";}
            st_file << stc << "    1    1    1    1    1    1    1   ";
            if (y < 10){st_file << " ";}
            st_file << y << "   ";
            if (x < 10){st_file << " ";}
            st_file << x << "  " << fixed << setprecision(6) << WholeAssm[HexMap[y-1][x-1] - 1][0] << " " << WholeAssm[HexMap[y-1][x-1] - 1][1] << " \n";}} //Fix Errors later (error is just result again)
    st_file << "];\n\n";
    st_file.close();
    cout << "\nThe file " << serptoolOut << " has been written with detector data for plotting assembly-wise results.\n";
}

void PrintAxialPower (void){
    double inArray [112]; //Input file array
    double inError [112]; //Associated errors with each result
    double averageArray [112];
    int layers = 112;
    double axialOffset = 0;
    double AOE = 0;
    double maxError = 0;
    double aPPF = 0;
    bool doneReading = false;
    int i = 0;
    long int k = 0;
    int l = 0;
    long int linnumber = 37584; //Line Number
    string wholeline;
    string component;
    string errorA;
    char charnumb [20];
    char errornumb [10];
    double value = 0;
    triresults_file.open(triresultsIn);
    if (triresults_file.fail()) {i=0;}
        else{

            do {
                getline(triresults_file, wholeline );
                k++;
                if (k == linnumber){
                    component = wholeline.substr (52,11);
                    for (l = 0; l < 11; l++){
                        charnumb[l] = component[l];}
                    value = atof(charnumb);
                    inArray[i] = value;
                    errorA = wholeline.substr (64,8);
                    for (l = 0; l < 7; l++){
                        errornumb[l] = errorA[l];}
                    value = atof(errornumb);
                    inError[i] = value;
                    if (inError[i] > maxError) {maxError = inError[i];}
                    i++;
                    if (i == layers) {doneReading = true;}
                    linnumber++;}
            }while (!doneReading);
            triresults_file.close();
            double sumInput = 0;
            double sumInputBot = 0;
            double sumInputTop = 0;
            double sumErrorBot = 0;
            double sumErrorTop = 0;
            for (i=0;i<56;i++){
                sumInputBot += inArray[i];
                sumErrorBot += inArray[i]*inArray[i]*inError[i]*inError[i];}
            sumErrorBot = sqrt(sumErrorBot)/(sumInputBot);
            for (i=56;i<112;i++){
                sumInputTop += inArray[i];
                sumErrorTop += inArray[i]*inArray[i]*inError[i]*inError[i];}
            sumErrorTop = sqrt(sumErrorTop)/(sumInputTop);
            sumInput = sumInputBot + sumInputTop;
            double errorAO;
            axialOffset = (sumInputTop - sumInputBot)/sumInput;
            errorAO = (2/((sumInputBot+sumInputTop)*(sumInputBot+sumInputTop)))*sqrt(sumInputBot*sumInputBot*sumInputTop*sumErrorTop*sumInputTop*sumErrorTop + sumInputTop*sumInputTop*sumInputBot*sumErrorBot*sumInputBot*sumErrorBot)/axialOffset;
            AOE = axialOffset*errorAO;
            for (i=0;i<layers;i++){
                averageArray[i] = inArray[i]*112/sumInput;
                if (averageArray[i] > aPPF) {aPPF = averageArray[i];}}

            axialresults_file.open(axialpowerOut);
            axialresults_file << "Axial Index\tNormalized Axial Power\n";
            for (i = 0; i < layers; i++){
                axialresults_file << (i+1) << "\t" << averageArray[i] << "\n";}
            axialresults_file << "\nPPF:\t" << aPPF << "\nAxial Offset [%]:\t" << axialOffset*100 << "\nAxial Offset Error [%]:\t" << AOE*100 << endl;
            axialresults_file.close();
            cout << "\nThe file " << axialpowerOut << " has been written with the axial power profile.\n";}
}

void CRScheduleRead (void){
    int i = 0;
    int l = 0;
    string wholeline;
    char numb [3];
    int value = 0;
    CRScheduleIn_file.open(CRScheduleFile);
    if (CRScheduleIn_file.fail()) {
        cout << "Could not open CRSchedule.txt file.\n";}
        else {
            getline(CRScheduleIn_file, wholeline );
            for (l = 0; l < 2; l++){
                numb[l] = wholeline[l];}
            value = atoi(numb);
            CRGroupNumb = value;
            for(i=0;i<CRGroupNumb;i++){
                getline(CRScheduleIn_file, wholeline );
                for (l = 0; l < 2; l++){
                    numb[l] = wholeline[l];}
                value = atoi(numb);
                CRScheduleArray[i] = value-1;}
            CRScheduleIn_file.close();}
}

void CRScheduleWrite (void){
    int i = 0;
    CRScheduleOut_file.open(CRScheduleFile);
    CRScheduleArray[CRGroupNumb] = PeakCRGroup;
    CRGroupNumb++;
    CRActiveNumb = CRGroupNumb;
    CRScheduleOut_file << CRGroupNumb << " \n";
    for(i=0;i<CRGroupNumb;i++){CRScheduleOut_file << CRScheduleArray[i]+1 << " \n";} //If not in schedule search mode, should echo back what was originally read in
    CRScheduleOut_file.close();
    cout << "\nThe file " << CRScheduleFile << " has been written with updated Control Blade insertion information.\n";
}

void RodWithdrawSchedule (void){
    //Remove the lowest CR power from the schedule
    int i = 0;
    int l = 0;
    int removeGroup;
    for (i=0;i<CRGroupNumb;i++){
        if (CRScheduleArray[i] == MinCRGroup){removeGroup = i;}}
    CRActiveNumb--;
    if (CRActiveNumb < 0){cout << "***ERROR: Trying to withdraw control blades when none are inserted.***\n";}
    for (i=0;i<CRActiveNumb;i++){
        if (i >= removeGroup){CRScheduleArray[i] = CRScheduleArray[i+1];}} //Removes the desired group element and collapses the array to account for the removed element
    CRScheduleArray[CRActiveNumb] = MinCRGroup; //Puts the removed element at the now freed-up location
    CRScheduleOut_file.open(CRScheduleFile);
    CRScheduleOut_file << CRGroupNumb << " " << CRActiveNumb << "  " << CRPartialNumb << "  \n";
    for(i=0;i<CRGroupNumb;i++){CRScheduleOut_file << CRScheduleArray[i]+1 << " \n";}
    CRScheduleOut_file.close();
    //Read in the current withdrawal list
    string wholeline;
    char numb [2];
    int value = 0;
    withdrawIn_file.open(withdrawFile);
    if (withdrawIn_file.fail()) {withdrawNumb = 0;} //If file does not yet exist, set number of withdrawn rods to zero
        else{
            getline(withdrawIn_file, wholeline );
            for (l = 0; l < 2; l++){
                numb[l] = wholeline[l];}
            value = atoi(numb);
            withdrawNumb = value;
            for(i=0;i<withdrawNumb;i++){
                getline(withdrawIn_file, wholeline );
                for (l = 0; l < 2; l++){
                    numb[l] = wholeline[l];}
                value = atoi(numb);
                withdrawArray[i] = value-1;}
            withdrawIn_file.close();}
    //Add the removed location to the withdrawal list
    withdrawOut_file.open(withdrawFile);
    withdrawArray[withdrawNumb] = MinCRGroup;
    withdrawNumb++;
    withdrawOut_file << withdrawNumb << " \n";
    for(i=0;i<withdrawNumb;i++){withdrawOut_file << withdrawArray[i]+1 << " \n";}
    withdrawOut_file.close();
    cout << "\nThe file " << withdrawFile << " has been written with updated Control Blade withdrawal information.\n";
}

void resmRead (void) {
    int i = 0;
    int l = 0;
    int k = 0;
    bool doneReading = false;
    int linnumber = 0; // Implicit keff
    int eofnumber = 0;
    string impk = "IMP_KEFF";
    string wholeline;
    string component;
    char kchar [9];
    char knumb [12];
    char ken [8];
    double kvalue = 0;
    double keval = 0;
    resm_file.open(resmFile);
    if (resm_file.fail()) {
        cout << "\nCould not open input.txt_res.m file to extract results from prior simulation.\nNote: this is expected for the first iteration step.\n";}
        else {
            //Find end of file line
            do {
                getline(resm_file, wholeline);
                if (resm_file.eof()) {doneReading = true;}
                    else {eofnumber++;}
            }while (!doneReading);
            resm_file.close();
            doneReading = false;
            linnumber = eofnumber-1;
            //Find last IMP_KEFF line
            resm_file.open(resmFile);
            do {
                k++;
                getline(resm_file, wholeline);
                if (k == linnumber){
                    component = wholeline.substr (0,8);
                    if (component.compare(impk) == 0){doneReading = true;} //Found last instance of IMP_KEFF
                        else{
                            linnumber--;
                            resm_file.close();
                            k = 0;
                            resm_file.open(resmFile);}}
            }while (!doneReading);
            resm_file.close();
            doneReading = false;
            resm_file.open(resmFile);
            k = 0;
            //Extract keff value
            do {
                k++;
                getline(resm_file, wholeline);
                if (k == linnumber){
                    component = wholeline.substr (47,11);
                    for (l = 0; l < 11; l++){
                        knumb[l] = component[l];}
                    kvalue = atof(knumb);
                    keff = kvalue;
                    component = wholeline.substr (59,7);
                    for (l = 0; l < 7; l++){
                        ken[l] = component[l];}
                    keval = atof(ken);
                    kerr = keval;
                    doneReading = true;}
            }while (!doneReading);
            resm_file.close();
        }
}

void depRead (void){
    int i = 0;
    int l = 0;
    string wholeline;
    string component;
    char dopt [1];
    char dtot [6];
    char dstep [18];
    int ivalue = 0;
    double fvalue = 0;
    dep_file.open(depFile);
    if (dep_file.fail()) {cout << "Could not open dep.txt file.\n";}
        else {
            getline(dep_file, wholeline);
            for (l = 0; l < 1; l++){
                dopt[l] = wholeline[l];}
            ivalue = atoi(dopt);
            depType = ivalue;
            if ( (depType < 1) || (depType > 6)) {cout<< "***ERROR: Unsupported depletion type chosen in dep.txt Line 1.***\n";}
            getline(dep_file, wholeline);
            component = wholeline.substr (0,2);
            for (l = 0; l < 2; l++){
                dtot[l] = component[l];}
            ivalue = atoi(dtot);
            depStepTotal = ivalue;
            if (depStepTotal > 99) {cout<< "***ERROR: Only up to 99 depletion steps are supported. Script will not function properly.***";}
            for(i=0;i<depStepTotal;i++){
                getline(dep_file, wholeline);
                for (l = 0; l < 18; l++){
                    dstep[l] = wholeline[l];}
                fvalue = atof(dstep);
                depArray[i] = fvalue;}
            dep_file.close();}
}

void iterationIn (void) {
    int k = 0;
    int l = 0;
    bool doneReading = false;
    int eofnumber = 0;
    string wholeline;
    string component;
    string componentk;
    char tchar[3];
    char kchar[19];
    int ivalue = 0;
    double fvalue = 0;
    iterationIn_file.open(iterationFile);
    if (iterationIn_file.fail()) {
        firstRun = true;
        cout << "\nCould not open iteration.txt file to extract results from prior simulation.\nNote that this is expected for the first iteration step.\n";}
        else{
            //Find end of file line number
            do {
                getline(iterationIn_file, wholeline);
                if (iterationIn_file.eof()) {doneReading = true;}
                    else {eofnumber++;}
            }while (!doneReading);
            iterationIn_file.close();
            doneReading = false;
            //Read-in results from the last line
            iterationIn_file.open(iterationFile);
            do {
                k++;
                getline(iterationIn_file, wholeline);
                if (k == (eofnumber-2)) {
                    component = wholeline.substr (9,2);
                    for (l = 0; l < 2; l++){
                        tchar[l] = component[l];}
                    ivalue = atoi(tchar);
                    CRActivePrev = ivalue;
                    component = wholeline.substr (12,2);
                    for (l = 0; l < 2; l++){
                        tchar[l] = component[l];}
                    ivalue = atoi(tchar);
                    CRPartialPrev = ivalue;
                    componentk = wholeline.substr (15, 18);
                    for (l = 0; l < 18; l++){
                        kchar[l] = componentk[l];}
                    fvalue = atof(kchar);
                    kprev = fvalue;}
                    else if (k == (eofnumber-1)){
                        component = wholeline.substr (0,2);
                        for (l = 0; l < 2; l++){
                            tchar[l] = component[l];}
                        ivalue = atoi(tchar);
                        depStep = ivalue;
                        component = wholeline.substr (3,2);
                        for (l = 0; l < 2; l++){
                            tchar[l] = component[l];}
                        ivalue = atoi(tchar);
                        CRMoveStep = ivalue;
                        component = wholeline.substr (6,2);
                        for (l = 0; l < 2; l++){
                            tchar[l] = component[l];}
                        ivalue = atoi(tchar);
                        THstep = ivalue;
                        component = wholeline.substr (9,2);
                        for (l = 0; l < 2; l++){
                            tchar[l] = component[l];}
                        ivalue = atoi(tchar);
                        CRActiveNumb = ivalue;
                        component = wholeline.substr (12,2);
                        for (l = 0; l < 2; l++){
                            tchar[l] = component[l];}
                        ivalue = atoi(tchar);
                        CRPartialNumb = ivalue;
                        doneReading = true;}
            }while (!doneReading);
            iterationIn_file.close();}
}

void CRSearchMode (void){
    double deltak = 0;
    double deltaCR = 0;
    double deltaCR2 = 0;
    double truncCR = 0;
    double dpartial = 0;
    int cap = CRActivePrev;
    int cpp = CRPartialPrev;
    CRActivePrev = CRActiveNumb;
    CRPartialPrev = CRPartialNumb;
    if (abs(keff - ktarget) > kepsilon){ //Need to move CBs to reach criticality target
        if (CRMoveStep > 1){ //Use previous and current k to estimate the critical position
            deltak = keff - kprev;
            if ((abs(deltak)) < 0.00015){ // If too small, causes numerical instability. Assume movement effect was so slight that it is has already converged.
                nextDepStep = true;}
                else{
                deltaCR = 1.0*(CRActiveNumb - cap) + 1.0*(CRPartialNumb - cpp)/(nL - 2); //1.0 is necessary to force int to float
                deltaCR2 = (keff-ktarget)*deltaCR/deltak;
                truncCR = trunc(deltaCR2); //Compute CR portion
                dpartial = (deltaCR2 - truncCR)*(nL-2);
                dpartial = round(dpartial);
                if ((truncCR == 0) && (dpartial == 0)) { //Search is stuck. Does not want to move CBs. Just deplete to break cycle, likely too small of k window by user.
                    nextDepStep = true;}
                    else{
                        CRActiveNumb -= truncCR;
                        CRPartialNumb -= dpartial;
                        if (CRPartialNumb < 0) { //If negative partial, add number of layers and decrement groups.
                            CRPartialNumb += (nL-2);
                            CRActiveNumb--;}
                            else if (CRPartialNumb >= (nL-2)) { //If partial greater than full, subtract number of layers and increment groups
                                CRPartialNumb -= (nL-2);
                                CRActiveNumb++;}
                        if (CRActiveNumb < 0) { //If negative insertion, just have all withdrawn. Deplete since best possible, if subcritical.
                            CRActiveNumb = 0;
                            CRPartialNumb = 0;
                            if (keff < ktarget) {nextDepStep = true;}}
                            else if (CRActiveNumb > 83){ //If more than all CBs needed, just insert all. Deplete since best possible, if supercritical.
                                CRActiveNumb = 84;
                                CRPartialNumb = 0;
                                if (keff > ktarget) {nextDepStep = true;}}}}}
            else{ //No usable k value available.
                if (CRActiveNumb == 0){
                    if (CRPartialNumb != 0){ //If on last CB group, just fully withdraw it. No other options left for control, so can just proceed with depletion
                        CRPartialNumb = 0;
                        nextDepStep = true;}
                        else{nextDepStep = true;}}
                    else{ //Withdraw an assembly and gauge impact
                        CRActiveNumb--;}}}
        else{ //Satisfied search criteria, can run depletion with current.
            nextDepStep = true;}
    if (firstRun){ //Only for the first run, use user's guess for insertion. Assume that it needs to be tested.
        CRActiveNumb = initialCRGroups;
        CRPartialNumb = 0;
        nextDepStep = false;}
}

void iterationOut (void){
    if (firstRun){iterationOut_file.open(iterationFile);} //If first run, create iteration file
        else{iterationOut_file.open(iterationFile,std::ofstream::app);} //Otherwise, append to existing file
    if (depStep <  10){iterationOut_file << "0";}
    iterationOut_file << depStep << " ";
    if (CRMoveStep < 10){iterationOut_file << "0";}
    iterationOut_file << CRMoveStep << " ";
    if (THstep < 10) {iterationOut_file << "0";}
    iterationOut_file << THstep << " ";
    if (CRActivePrev < 10){iterationOut_file << "0";}
    iterationOut_file << CRActivePrev << " ";
    if (CRPartialPrev < 10){iterationOut_file << "0";}
    iterationOut_file << CRPartialPrev << " " << keff << " \n";
    if (CRCritSearch){
        if (nextDepStep) {
            if (THstep == THiters){
                depStep++;
                CRMoveStep = 1;
                THstep = 0;}
                else{
                    THstep++;
                    nextDepStep = false;}}
            else{CRMoveStep++;}}
        else{ //THMode without crit search
            if(THstep == THiters){
                nextDepStep = true;
                depStep++;
                THstep = 0;}
                else{THstep++;}}
    if (depStep <  10){iterationOut_file << "0";}
    iterationOut_file << depStep << " ";
    if (CRMoveStep < 10){iterationOut_file << "0";}
    iterationOut_file << CRMoveStep << " ";
    if (THstep < 10) {iterationOut_file << "0";}
    iterationOut_file << THstep << " ";
    if (CRActiveNumb < 10){iterationOut_file << "0";}
    iterationOut_file << CRActiveNumb << " ";
    if (CRPartialNumb < 10){iterationOut_file << "0";}
    iterationOut_file << CRPartialNumb << " \n\n";
    iterationOut_file.close();
    cout << "\nThe file " << iterationFile << " has been written with updated criticality and iteration information.\n";
}

void RodInsertionAmount (void){
    int loopCounter = 0;
    for (assm=0;assm<CRActiveNumb;assm++){ //Complete Insertions
        for (layr=(nL-1);layr>0;layr--){
            controlBladeIn[layr][CRScheduleArray[assm]] = true;
            loopCounter++;
            //cout << loopCounter << " " << layr << " " << assm << " " << CRScheduleArray[assm] << endl;
            }}
    if (CRPartialNumb > 0){
        for(layr=(nL-1);layr>(nL-2-CRPartialNumb);layr--){ //Partial Insertion
            controlBladeIn[layr][CRScheduleArray[CRActiveNumb]] = true;
            loopCounter++;
            //cout << loopCounter << " ";
            }}
}

void dims (void){
    dimActiveCoreHeight=0;
    dimAssemblyPitch[0] = dimColdAssemblyPitch*(1+alphag*(inletT-coldT));
    dimAssemblyPitch[nL-1] = dimColdAssemblyPitch*(1+alphasic*(outletT-coldT));
    double plateAverage=0;
    dimReflHoleRadius = dimColdReflHoleRadius*(1+alphag*(inletT-coldT));
    dimReflHex = dimColdReflHex*(1+alphag*(inletT-coldT));
    dimLayrHeight[0]=dimColdAxReflHeight*(1+alphag*(inletT-coldT));
    dimLayrHeight[nL-1]=dimColdAxReflHeight*(1+alphag*(outletT-coldT));
    dimAxSuppPlateHeight[0] = dimColdAxSuppPlateHeight*(1+alphag*(inletT-coldT));
    dimAxSuppPlateHeight[1] = dimColdAxSuppPlateHeight*(1+alphasic*(outletT-coldT));
    for (layr=1; layr<(nL-1);layr++){
        dimZpitch[layr]=dimColdZPitch*(1+alphag*dTmatave[layr]);
        dimZEupitch[layr]=dimColdZEuPitch*(1+alphag*dTmatave[layr]); // This temperature used to ensure that no BP volume is lost due to particles being cut-off
        dimLayrHeight[layr] = (dimColdActiveCoreHeight/(nL-2))*(1+alphag*dTmatave[layr]);
        dimActiveCoreHeight+=dimLayrHeight[layr];
        if (axialFlower == true){
            dimAssemblyPitch[layr] = dimAssemblyPitch[0] + (dimAssemblyPitch[nL-1]-dimAssemblyPitch[0])*(dimLayrHeight[0]+dimLayrHeight[layr]*(layr+0.5))/(dimActiveCoreHeight);
            dimAssemblyPitch[0] = dimAssemblyPitch[0] + (dimAssemblyPitch[nL-1]-dimAssemblyPitch[0])*(dimLayrHeight[0]*0.5)/(dimActiveCoreHeight);
            dimAssemblyPitch[nL-1] = dimAssemblyPitch[0] + (dimAssemblyPitch[nL-1]-dimAssemblyPitch[0])*((dimActiveCoreHeight-dimLayrHeight[0]*0.5))/(dimActiveCoreHeight);
            }
            else{
                plateAverage = 0.5*(dimAssemblyPitch[0] + dimAssemblyPitch[nL-1]);
                dimAssemblyPitch[0] = plateAverage;
                dimAssemblyPitch[nL-1] = plateAverage;
                dimAssemblyPitch[layr] = plateAverage;}
        dimFlibeHex[layr]=0.5*dimAssemblyPitch[layr];
        for (assm=0; assm<nA;assm++){
            for (sect=0;sect<3;sect++){
                dimFuelRadius[layr][assm][sect]=dimColdFuelR*(1+alphaf*dTf[layr][assm][sect]);
                dimOpycRadius[layr][assm][sect]=dimColdOpycR*(1+alphag*dTmat[layr][assm][sect]);
                //Have Opyc, Sic, and IPyc follow from densities only. Fuel freely expands. Buffer accomodates the two.
                denOpyc[layr][assm][sect]=denColdOpyc/pow((1+alphaopyc*dTopyc[layr][assm][sect]),3);
                // mh=mc  VhDenh=VcDenc  Vh=Vc(Denc/Denh) roh3-rih3=(roc3-ric3)*(Dc/Dh)   rih=cbrt(roh3-(roc3-ric3)*(Dc/Dh))
                dimSicRadius[layr][assm][sect]=cbrt(pow(dimOpycRadius[layr][assm][sect],3)-(denColdOpyc/denOpyc[layr][assm][sect])*(pow(dimColdOpycR,3) - pow(dimColdSicR,3)));
                denSic[layr][assm][sect]=denColdSic/pow((1+alphasic*dTsic[layr][assm][sect]),3);
                dimIpycRadius[layr][assm][sect]=cbrt(pow(dimSicRadius[layr][assm][sect],3)-(denColdSic/denSic[layr][assm][sect])*(pow(dimColdSicR,3) - pow(dimColdIpycR,3)));
                denIpyc[layr][assm][sect]=denColdIPyc/pow((1+alphaipyc*dTipyc[layr][assm][sect]),3);
                dimBufferRadius[layr][assm][sect]=cbrt(pow(dimIpycRadius[layr][assm][sect],3)-(denColdIPyc/denIpyc[layr][assm][sect])*(pow(dimColdIpycR,3) - pow(dimColdBufferR,3)));
                //This B value is a "fudge factor" to get equivalent thermal expansion behavior in different directions.
                if (!cubic){
                    B[layr][assm][sect]=sqrt((pow((1+alphag*dTmat[layr][assm][sect]),3))/(1+alphag*dTmatave[layr]));
                    dimXpitch[layr][assm][sect]=dimColdXPitchAsym*B[layr][assm][sect];
                    dimYpitch[layr][assm][sect]=dimColdYPitchAsym*B[layr][assm][sect];}
                    else{
                        B[layr][assm][sect] = 1;
                        dimXpitch[layr][assm][sect] = dimZpitch[layr];
                        dimYpitch[layr][assm][sect] = dimZpitch[layr];}
                dimXLength[layr][assm][sect]=nTW*dimXpitch[layr][assm][sect];
                dimYWidth[layr][assm][sect]=nTL*dimYpitch[layr][assm][sect];
                dimEuRadius[layr][assm][sect]=dimColdEuR*(1+alphaeu*dTeu[layr][assm][sect]);
    }}}
    //Assume Permanent reflector expands freely; assume B4C and CB have free expansion density and force radius to that corresponding value. (In -> out)
    dimPermReflRadius=dimColdPermRadReflRadius*(1+alphag*(averageT-coldT));
    dimBoronCarbideRadius=dimPermReflRadius+(dimColdBoronCarbideRadius-dimColdPermRadReflRadius)*(1+alphabc*(averageT-coldT));
    dimCoreBarrelRadius=dimBoronCarbideRadius+(dimColdCoreBarrelRadius-dimColdBoronCarbideRadius)*(1+alphag*(averageT-coldT));
    //Assume RPV radius expands freely; assume Liner has free expansion density and force its inner radius (downcomer radius) to corresponding value. (Out -> in)
    dimPressureVesselRadius=dimColdPressureVesselRadius*(1+alpha800H*(inletT-coldT));
    dimAlloyNLinerRadius = dimColdAlloyNLinerRadius*(1+alpha800H*(inletT-coldT));
    dimDowncomerRadius=dimAlloyNLinerRadius - (dimColdAlloyNLinerRadius-dimColdDowncomerRadius)*(1+alphaalloyN*(inletT-coldT));
    dimModelHeight = 0;
    for(layr=0;layr<nL;layr++){dimModelHeight += dimLayrHeight[layr];}
    dimModelHeight += dimAxSuppPlateHeight[0] + dimAxSuppPlateHeight[1];
}

void XYLocs (void) {
    for (layr=0; layr<nL;layr++){
        for (assm=0; assm<nA;assm++){
            for (sect=0;sect<nS;sect++){ //Account for expansions here
                //Structural Planes
                locSP[3][layr][assm][sect] = 2*(1+alphacc*dTcc[layr][assm][sect]);
                locSP[4][layr][assm][sect] = 4*(1+alphacc*dTcc[layr][assm][sect]);
                locSP[5][layr][assm][sect] = 0.5*(1+alphacc*dTcc[layr][assm][sect]);
                locSP[6][layr][assm][sect] = 10*(1+alphacb*dTcc[layr][assm][sect]);
                locSP[7][layr][assm][sect] = 20*(1+alphacb*dTcc[layr][assm][sect]);
                locSP[8][layr][assm][sect] = 1*(1+alphacb*dTcc[layr][assm][sect]);
                locSP[9][layr][assm][sect] = 0.88*(1+alphacc*dTcc[layr][assm][sect]);
                locSP[10][layr][assm][sect] = 10.38*(1+alphacc*dTcc[layr][assm][sect]);
                locSP[11][layr][assm][sect] = 20.76*(1+alphacb*dTcc[layr][assm][sect]);
                locSP[12][layr][assm][sect] = 1.76*(1+alphacc*dTcc[layr][assm][sect]);
                locSP[13][layr][assm][sect] = 22.5*(1+alphacc*dTcc[layr][assm][sect]);
                locSP[14][layr][assm][sect] = 21.5*(1+alphacc*dTcc[layr][assm][sect]);
                locSP[15][layr][assm][sect] = 45*(1+alphacc*dTcc[layr][assm][sect]);
                locSP[16][layr][assm][sect] = 42.3*(1+alphacc*dTcc[layr][assm][sect]);
                lociL[layr][assm][sect] = (43-23.1*sqrt3)*(1+alphacc*dTcc[layr][assm][sect]); //2.98...
                lociR[layr][assm][sect] = 43*(1+alphacc*dTcc[layr][assm][sect]);
                //Central Plank Planes (goes with structural expansion due to shift)
                locCP[1][layr][assm][sect] = 3.625*(1+alphacc*dTcc[layr][assm][sect]);
                locCP[2][layr][assm][sect] = 6.875*(1+alphacc*dTcc[layr][assm][sect]);
                locCP[3][layr][assm][sect] = 10.125*(1+alphacc*dTcc[layr][assm][sect]);
                locCP[4][layr][assm][sect] = 13.375*(1+alphacc*dTcc[layr][assm][sect]);
                locCP[5][layr][assm][sect] = 16.625*(1+alphacc*dTcc[layr][assm][sect]);
                locCP[6][layr][assm][sect] = 19.875*(1+alphacc*dTcc[layr][assm][sect]);
                //Fuel Plank Planes (goes with plank expansion, so slightly hotter)
                locPP[1][layr][assm][sect] = locCP[1][layr][assm][sect] - (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                locPP[2][layr][assm][sect] = locCP[1][layr][assm][sect] + (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                locPP[3][layr][assm][sect] = locCP[2][layr][assm][sect] - (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                locPP[4][layr][assm][sect] = locCP[2][layr][assm][sect] + (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                locPP[5][layr][assm][sect] = locCP[3][layr][assm][sect] - (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                locPP[6][layr][assm][sect] = locCP[3][layr][assm][sect] + (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                locPP[7][layr][assm][sect] = locCP[4][layr][assm][sect] - (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                locPP[8][layr][assm][sect] = locCP[4][layr][assm][sect] + (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                locPP[9][layr][assm][sect] = locCP[5][layr][assm][sect] - (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                locPP[10][layr][assm][sect] = locCP[5][layr][assm][sect] + (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                locPP[11][layr][assm][sect] = locCP[6][layr][assm][sect] - (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                locPP[12][layr][assm][sect] = locCP[6][layr][assm][sect] + (0.5*(dimColdPlankWidth-2*dimColdSleeveWidth)*(1+alphag*dTmeat[layr][assm][sect]) + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]));
                //Spacer Radius. Will half the smaller one in file; just carry one value. Just span the gap as it stands.
                dimSpacer[layr][assm][sect] = locPP[3][layr][assm][sect] - locPP[2][layr][assm][sect];
                //Sleeve Planes (expand at sleeve temp)
                locSL[1][layr][assm][sect] = locPP[1][layr][assm][sect] + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                locSL[2][layr][assm][sect] = locPP[2][layr][assm][sect] - dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                locSL[3][layr][assm][sect] = locPP[3][layr][assm][sect] + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                locSL[4][layr][assm][sect] = locPP[4][layr][assm][sect] - dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                locSL[5][layr][assm][sect] = locPP[5][layr][assm][sect] + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                locSL[6][layr][assm][sect] = locPP[6][layr][assm][sect] - dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                locSL[7][layr][assm][sect] = locPP[7][layr][assm][sect] + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                locSL[8][layr][assm][sect] = locPP[8][layr][assm][sect] - dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                locSL[9][layr][assm][sect] = locPP[9][layr][assm][sect] + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                locSL[10][layr][assm][sect] = locPP[10][layr][assm][sect] - dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                locSL[11][layr][assm][sect] = locPP[11][layr][assm][sect] + dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                locSL[12][layr][assm][sect] = locPP[12][layr][assm][sect] - dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
                //Central particle locations - dependent on the structural
                for (plnk=1;plnk<7;plnk++){
                    locBP[plnk][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locCP[plnk][layr][assm][sect] - 2*(1+alphacc*dTcc[layr][assm][sect])));}
                //BP cylinders (for odd stacks per plank)
                for (plnk = 1; plnk < 7; plnk++){
                    for (colm = 1; colm < (nBPC+1); colm++){
                        locBPC[plnk][colm][layr][assm][sect] = locBP[plnk][layr][assm][sect] - nBPW*(nBPC-colm-(0.5*(nBPC-1)))*dimZEupitch[layr];}}
                // Fuel planes. Three for each fuel stripe since sleeve plane already exists.
                locFy[1][layr][assm][sect] = locSL[1][layr][assm][sect] + nTL*dimYpitch[layr][assm][sect];
                locFy[2][layr][assm][sect] = locSL[2][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
                locFy[3][layr][assm][sect] = locSL[3][layr][assm][sect] + nTL*dimYpitch[layr][assm][sect];
                locFy[4][layr][assm][sect] = locSL[4][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
                locFy[5][layr][assm][sect] = locSL[5][layr][assm][sect] + nTL*dimYpitch[layr][assm][sect];
                locFy[6][layr][assm][sect] = locSL[6][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
                locFy[7][layr][assm][sect] = locSL[7][layr][assm][sect] + nTL*dimYpitch[layr][assm][sect];
                locFy[8][layr][assm][sect] = locSL[8][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
                locFy[9][layr][assm][sect] = locSL[9][layr][assm][sect] + nTL*dimYpitch[layr][assm][sect];
                locFy[10][layr][assm][sect] = locSL[10][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
                locFy[11][layr][assm][sect] = locSL[11][layr][assm][sect] + nTL*dimYpitch[layr][assm][sect];
                locFy[12][layr][assm][sect] = locSL[12][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
                locFxL[1][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[1][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxL[2][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[2][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxL[3][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[3][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxL[4][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[4][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxL[5][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[5][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxL[6][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[6][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxL[7][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[7][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxL[8][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[8][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxL[9][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[9][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxL[10][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[10][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxL[11][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[11][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxL[12][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[12][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) - 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[1][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[1][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[2][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[2][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[3][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[3][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[4][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[4][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[5][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[5][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[6][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[6][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[7][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[7][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[8][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[8][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[9][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[9][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[10][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[10][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[11][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[11][layr][assm][sect] - 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                locFxR[12][layr][assm][sect] = (12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*((locFy[12][layr][assm][sect] + 0.5*nTL*dimYpitch[layr][assm][sect]) - 2*(1+alphacc*dTcc[layr][assm][sect]))) + 0.5*nTW*dimXpitch[layr][assm][sect];
                //Spacer. Vertical determined by hotter plate expansion. X determined by cooler structural expansion.
                locSpL[1][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[1][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) - 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpL[2][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[3][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) - 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpL[3][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[5][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) - 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpL[4][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[7][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) - 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpL[5][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[9][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) - 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpL[6][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[11][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) - 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpL[7][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[12][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) - 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpR[1][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[1][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) + 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpR[2][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[3][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) + 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpR[3][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[5][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) + 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpR[4][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[7][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) + 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpR[5][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[9][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) + 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpR[6][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[11][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) + 7*(1+alphacc*dTcc[layr][assm][sect]));
                locSpR[7][layr][assm][sect] = ((12.2109582*(1+alphacc*dTcc[layr][assm][sect]) - (1/sqrt3)*(locPP[12][layr][assm][sect]-2*(1+alphacc*dTcc[layr][assm][sect]))) + 7*(1+alphacc*dTcc[layr][assm][sect]));
                //Fuel Lattice Shifts
                locFLx[1][layr][assm][sect] = locFxL[1][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLx[2][layr][assm][sect] = locFxL[2][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLx[3][layr][assm][sect] = locFxL[3][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLx[4][layr][assm][sect] = locFxL[4][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLx[5][layr][assm][sect] = locFxL[5][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLx[6][layr][assm][sect] = locFxL[6][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLx[7][layr][assm][sect] = locFxL[7][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLx[8][layr][assm][sect] = locFxL[8][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLx[9][layr][assm][sect] = locFxL[9][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLx[10][layr][assm][sect] = locFxL[10][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLx[11][layr][assm][sect] = locFxL[11][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLx[12][layr][assm][sect] = locFxL[12][layr][assm][sect] - dimXpitch[layr][assm][sect];
                locFLy[1][layr][assm][sect] = locSL[1][layr][assm][sect];
                locFLy[2][layr][assm][sect] = locSL[2][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
                locFLy[3][layr][assm][sect] = locSL[3][layr][assm][sect];
                locFLy[4][layr][assm][sect] = locSL[4][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
                locFLy[5][layr][assm][sect] = locSL[5][layr][assm][sect];
                locFLy[6][layr][assm][sect] = locSL[6][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
                locFLy[7][layr][assm][sect] = locSL[7][layr][assm][sect];
                locFLy[8][layr][assm][sect] = locSL[8][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
                locFLy[9][layr][assm][sect] = locSL[9][layr][assm][sect];
                locFLy[10][layr][assm][sect] = locSL[10][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
                locFLy[11][layr][assm][sect] = locSL[11][layr][assm][sect];
                locFLy[12][layr][assm][sect] = locSL[12][layr][assm][sect] - nTL*dimYpitch[layr][assm][sect];
            }}}
}

void densities(void) {
    if (!coldDim){ //Only update densities if not using cold dimensions.
        denbarrel = denColdGraphite/pow((1+alphag*(daverageT)),3);
        denBC = denColdBC/pow((1+alphabc*(daverageT)),3);
        denAlloyN = denColdAlloyN/pow((1+alphaalloyN*(inletT-coldT)),3);
        denH800 = denColdH800/pow((1+alpha800H*(inletT-coldT)),3);;
        for (layr=0; layr<nL; layr++){
            denaxialavegraph[layr] = 0;
            denaxialaveflibe[layr] = 0;
            denaxialavemat[layr] = 0;
            for (assm=0; assm<nA;assm++){
                for (sect=0;sect<nS;sect++){
                    denFuel[layr][assm][sect]=denColdFuel*pow((dimColdFuelR/dimFuelRadius[layr][assm][sect]),3);
                    denBuffer[layr][assm][sect]=denColdBuffer*((pow(dimColdBufferR,3)-pow(dimColdFuelR,3))/(pow(dimBufferRadius[layr][assm][sect],3)-pow(dimFuelRadius[layr][assm][sect],3)));
                    denIpyc[layr][assm][sect] = denColdIPyc*((pow(dimColdIpycR,3)- pow(dimColdBufferR,3))/(pow(dimIpycRadius[layr][assm][sect],3)-pow(dimBufferRadius[layr][assm][sect],3)));
                    denSic[layr][assm][sect] = denColdSic*((pow(dimColdSicR,3)-pow(dimColdIpycR,3))/(pow(dimSicRadius[layr][assm][sect],3)-pow(dimIpycRadius[layr][assm][sect],3)));
                    denOpyc[layr][assm][sect] = denColdOpyc*((pow(dimColdOpycR,3)-pow(dimColdSicR,3))/(pow(dimOpycRadius[layr][assm][sect],3)-pow(dimSicRadius[layr][assm][sect],3)));
                    denMat[layr][assm][sect]=denColdGraphite*((dimColdZPitch*dimColdXPitchAsym*dimColdYPitchAsym-(4/3)*pi*pow(dimColdOpycR,3))/(dimZpitch[layr]*dimXpitch[layr][assm][sect]*dimYpitch[layr][assm][sect]-(4/3)*pi*pow(dimOpycRadius[layr][assm][sect],3)));
                    denaxialavemat[layr] += denMat[layr][assm][sect]/(3*84);
                    denEu[layr][assm][sect] = denColdEu*pow((dimColdEuR/dimEuRadius[layr][assm][sect]),3);
                    denMeat[layr][assm][sect] = denColdGraphite/pow((1+alphag*dTmeat[layr][assm][sect]),3);
                    denSleeve[layr][assm][sect] = denColdGraphite/pow((1+alphag*dTsleeve[layr][assm][sect]),3);
                    denSpacer[layr][assm][sect] = denColdGraphite/pow((1+alphag*dTspacer[layr][assm][sect]),3);
                    denPlank[layr][assm][sect] = denColdGraphite/pow((1+alphag*dTg[layr][assm][sect]),3);
                    denaxialavegraph[layr] += denPlank[layr][assm][sect]/(84*3);
                    denCB[layr][assm][sect] = denColdCB/pow((1+alphacb*dTcb[layr][assm][sect]),3);
                    denCCC[layr][assm][sect] = denColdCCC/pow((1+alphacc*dTcc[layr][assm][sect]),3);
                    denflibe[layr][assm][sect] = 0.001*(2279.7 - 0.4884*(Tflibe[layr][assm][sect] - 273.15));
                    denaxialaveflibe[layr] += denflibe[layr][assm][sect]/(84*3);
                    }}}}
        else{ //If using cold dimensions, will use cold densities.
            denbarrel = denColdGraphite; //update
            denBC = denColdBC;
            denAlloyN = denColdAlloyN;
            denH800 = denColdH800;
            for (layr=0; layr<nL; layr++){
                denaxialavegraph[layr] = denColdGraphite;
                denaxialaveflibe[layr] = 1.95;
                denaxialavemat[layr] = denColdGraphite;
                for (assm=0; assm<nA;assm++){
                    for (sect=0;sect<nS;sect++){
                        denFuel[layr][assm][sect]=denColdFuel;
                        denBuffer[layr][assm][sect]=denColdBuffer;
                        denIpyc[layr][assm][sect] = denColdIPyc;
                        denSic[layr][assm][sect] = denColdSic;
                        denOpyc[layr][assm][sect] = denColdOpyc;
                        denMat[layr][assm][sect]=denColdGraphite;
                        denaxialavemat[layr] += denMat[layr][assm][sect]/(3*84);
                        denEu[layr][assm][sect] = denColdEu;
                        denSleeve[layr][assm][sect] = denColdGraphite;
                        denSpacer[layr][assm][sect] = denColdGraphite;
                        denPlank[layr][assm][sect] = denColdGraphite;
                        denCB[layr][assm][sect] = denColdCB;
                        denCCC[layr][assm][sect] = denColdCCC;
                        denflibe[layr][assm][sect] = 1.95;}}}}
    denSupportPlate[0] = 0.211*denaxialaveflibe[0] + 0.789*denColdCCC/pow((1+alphacc*(inletT-coldT)),3);
    denSupportPlate[1] = 0.43*denaxialaveflibe[nL-1] + 0.57*denColdSic/pow((1+alphasic*(outletT-coldT)),3);
}

void THTemps (void) {
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    double rt[8];
    double Tmatbound[18][84][3][4];
    for (layr=1;layr<(nL-1);layr++){
        for (assm=0;assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
    //Initializations to run
        Tcoolantboundary[0][assm][sect] = inletT;
        //Conductivities [W/(m*K)]
        double kgrph = coldkgrph;
        double kmatrix = coldkmatrix;
        double kfuel = coldkfuel;
        double kbuff = coldkbuff;
        double kipyro = coldkipyro;
        double kopyro = coldkopyro;
        double ksic = coldksic;
        double kmtrx = 0;
        double kflibe = coldkflibe;
        //Geometric Parameters [m]
        double rfuel = 0.01*dimFuelRadius[layr][assm][sect];
        double rbuff = 0.01*dimBufferRadius[layr][assm][sect];
        double ripyro = 0.01*dimIpycRadius[layr][assm][sect];
        double rsic = 0.01*dimSicRadius[layr][assm][sect];
        double ropyro = 0.01*dimOpycRadius[layr][assm][sect];
        double requiv = 0.01*pow(((3*dimXpitch[layr][assm][sect]*dimYpitch[layr][assm][sect]*dimZpitch[layr])/(4*pi)),(0.33333333));
        double xpitch = 0.01*dimXpitch[layr][assm][sect];
        double ypitch = 0.01*dimYpitch[layr][assm][sect];
        double zpitch = 0.01*dimZpitch[layr];
        double loccg = 0.01*0.5*(dimSpacer[layr][assm][sect]);
        double locgm = loccg + 0.01*dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
        double loccenter = 0.01*0.5*(locCP[2][layr][assm][sect] - locCP[1][layr][assm][sect]);
        double channelWidth = loccg;
        double channelLength = 0.01*22.112515*(1+alphacc*dTcc[layr][assm][sect]);
        double nTZ = dimLayrHeight[layr]/dimZpitch[layr];
        double channelPerimeter = 0;
        double Dh = 0;
        double locmc = 0;
        double dyc = 0;
        double dyg = 0;
        double dym = 0;
        double dycenter = 0;
        //Other geometric
        double vol = 0; //m3
        double channelVol = 0;
        double channelArea = 0; //m2
        double TrisoLayers = nTL;
        //Coolant parameters
        double denFlibe = denflibe[layr][assm][sect]; //g/cc
        double viscosityFlibe = flibeVisc; //P*s
        double PrandtlFlibe = flibePr;
        double cpFlibe = flibecp; //J/(kg*K)
        double ReFlibe = 0;
        double NuFlibe = 0;
        double hFlibe = 0; //W/(m2*K)
        double mfrFlibe = coreMFR/(252*12*3); //kg/s in 0.35 channel (coreTotal/(assm*channelPerSection*sect)
        double sectMFRFlibe = coreMFR/(252);
        double spdFlibe = 0; //m/s
        double ff = 0; //friction factor
        //Power parameters
        double powerperpart=0.95*SectionPower[layr][assm][sect][0]/(12*3*nTW*nTL*nTZ); //W
        double power = 0; //W
        double powerdensity = 0; //W/m3
        double heatflux = 0; //W/m2
        double peaking=0;
        //Temperatures [C]
        double dTc=0;
        double dTgrph=0;
        double dTm=0;
        double Tc = 0;
        double Tcg = 0;
        double Tgm = 0;
        double Tmc = 0;
        double Tin=inletT;
        double Tout=outletT;
    //Find homogenized fuel stripe average thermal conductivity (Maxwell's Method)
        double kd = 0; //discontinuous conductivity (particle)
        double kc = 0; //continuous conductivity (surrounding media)
        double pd = 0; //volume fraction
        double km = 0; //mixture conductivity
        //Homogenize fuel into buffer
        kd = kfuel;
        kc = kbuff;
        pd = pow((rfuel/rbuff),3);
        km = kc*(kd + 2*kc - 2*pd*(kc-kd))/(kd+2*kc+pd*(kc-kd));
        //Homogenize fuel/buffer into inner pyrolytic carbon
        kd = km;
        kc = kipyro;
        pd = pow((rbuff/ripyro),3);
        km = kc*(kd + 2*kc - 2*pd*(kc-kd))/(kd+2*kc+pd*(kc-kd));
        //Homogenize f/b/i into silicon carbide
        kd = km;
        kc = ksic;
        pd = pow((ripyro/rsic),3);
        km = kc*(kd + 2*kc - 2*pd*(kc-kd))/(kd+2*kc+pd*(kc-kd));
        //Homogenize f/b/i/s into outer pyrolytic carbon
        kd = km;
        kc = kopyro;
        pd = pow((rsic/ropyro),3);
        km = kc*(kd + 2*kc - 2*pd*(kc-kd))/(kd+2*kc+pd*(kc-kd));
        //Homogenize f/b/i/s/o into matrix carbon
        kd = km;
        kc = kgrph;
        pd = pow((ropyro/requiv),3);
        km = kc*(kd + 2*kc - 2*pd*(kc-kd))/(kd+2*kc+pd*(kc-kd));
        kmtrx = km;
    //Parameters needed to run
        vol=xpitch*ypitch*zpitch;
        channelVol=xpitch*zpitch*channelWidth;
        locmc = TrisoLayers*ypitch+locgm;
        dyc = loccg;
        dyg = locgm - loccg;
        dym = locmc - locgm;
        dycenter = loccenter - locmc;
        channelPerimeter = channelLength+2*channelWidth*(2/sqrt(3));
        channelArea = channelWidth*channelLength;
        Dh=4*channelArea/channelPerimeter;
        power = nTL*powerperpart;
        powerdensity=powerperpart/vol;
        heatflux = power/(xpitch*zpitch);
        spdFlibe = mfrFlibe/(denFlibe*channelArea);
        ReFlibe = denFlibe*spdFlibe*Dh/viscosityFlibe;
        ff=pow((0.79*log(ReFlibe)-1.64),(-2));
        NuFlibe = (ff/8)*(ReFlibe-1000)*(PrandtlFlibe)/(1+12.7*(pow((ff/8),(0.5)))*((pow(PrandtlFlibe,0.666666667))-1));
        hFlibe = kflibe*NuFlibe/Dh;
        dTc = heatflux/hFlibe;
        dTgrph = heatflux*dyg/kgrph;
        dTm = powerdensity*dym*dym/(2*kmtrx);
        Tcoolantboundary[layr][assm][sect] = Tcoolantboundary[layr-1][assm][sect] + 0.95*SectionPower[layr][assm][sect][0]/(cpFlibe*sectMFRFlibe);
        Tflibe[layr][assm][sect] = 0.5*(Tcoolantboundary[layr][assm][sect] + Tcoolantboundary[layr-1][assm][sect]);
        Tc = Tflibe[layr][assm][sect];
        Tcg=Tc+dTc;
        Tgm=Tcg+dTgrph;
        Tmc=Tgm+dTm;
        Tmatbound[layr][assm][sect][0] = Tc;
        Tmatbound[layr][assm][sect][1] = Tcg;
        Tmatbound[layr][assm][sect][2] = Tgm;
        Tmatbound[layr][assm][sect][3] = Tmc;
    //Resolve heterogeneous TRISO layer temperatures
        double currentLoc = 0;
        double qp = (power/nTL)*((3/(4*pi))/pow((rfuel),3)-pow((1/ypitch),3));
        double qn = -(power/nTL)*pow((1/ypitch),3);
        double partloc = 0;
        double kratio=0;
        double Tid [7];
        double kt[7];
        double At[7];
        double Bt[7];
        double Tdecom[7];
        kt[0]=kgrph;
        kt[1]=kmatrix;
        kt[2]=kopyro;
        kt[3]=ksic;
        kt[4]=kipyro;
        kt[5]=kbuff;
        kt[6]=kfuel;
        rt[1]=ypitch/2;
        rt[2]=ropyro;
        rt[3]=rsic;
        rt[4]=ripyro;
        rt[5]=rbuff;
        rt[6]=rfuel;
        rt[7]=0; //Center
        //Find expect deltaT across the fuel particle
        //These are related to the homogenized power shape
        double Am = dTm/(locmc*locmc-locgm*locgm-2*locmc*dym);
        double Bm = -2*Am*locmc;
        double Cm = Tgm - Am*locgm*locgm - Bm*locgm;
        double homoTemp = 0;
        //Heterogeneous interface parameters
        double Tint = 0;
        double Tinti[7];
        Tint += (-qp*pow(rt[6],5)/(30*kt[6]) + (1/3)*Bt[6]*pow(rt[6],3));
        Tinti[6] = (-qp*pow(rt[6],5)/(30*kt[6]) + (1/3)*Bt[6]*pow(rt[6],3));
        for(i=5;i>0;i--){
            Tint += ( -qn*pow(rt[i],5)/(30*kt[i]) + 0.5*At[i]*pow(rt[i],2) + (1/3)*Bt[i]*pow(rt[i],3) ) - ( -qn*pow(rt[i+1],5)/(30*kt[i]) + 0.5*At[i]*pow(rt[i+1],2) + (1/3)*Bt[i]*pow(rt[i+1],3) );
            Tinti[i] = ( -qn*pow(rt[i],5)/(30*kt[i]) + 0.5*At[i]*pow(rt[i],2) + (1/3)*Bt[i]*pow(rt[i],3) ) - ( -qn*pow(rt[i+1],5)/(30*kt[i]) + 0.5*At[i]*pow(rt[i+1],2) + (1/3)*Bt[i]*pow(rt[i+1],3) );
            }
        Tint += (2/pi-1/3)*pow(rt[1],3)*(-qn*pow(rt[1],2)/(6*kt[1]) + At[1]/rt[1] + Bt[1]);
        Tinti[0]=(2/pi-1/3)*pow(rt[1],3)*(-qn*pow(rt[1],2)/(6*kt[1]) + At[1]/rt[1] + Bt[1]);
        Tint /= 2*pow(rt[1],3)/pi;
        At[6]=0;
        Bt[1]=((qn*pow(rt[1],2))/(6*kt[1]))-((At[1])/(rt[1]));
        kratio=kt[6]/kt[5];
        At[5]=pow(rt[6],3)*(kratio*(qp/(3*kt[6]))-(qn/(3*kt[5])));
        for (i=5;i>1;i--){
            kratio=kt[i]/kt[i-1];
            At[i-1]=pow(rt[i],2)*(kratio*(((qn*rt[i])/(3*kt[i]))+((At[i])/(pow(rt[i],2))))-((qn*rt[i])/(3*kt[i-1])));}
        for (i=2;i<6;i++){
            Bt[i]=Bt[i-1] + (1/rt[i])*(At[i-1] - At[i]) - (qn*pow(rt[i],2)/6)*((1/kt[i-1])-(1/kt[i]));}
        Bt[6]=Bt[5] + (At[5]/rt[6]) - (pow(rt[6],2)/6)*((qn/kt[5])-(qp/kt[6]));
        //Heterogeneous step-function behavior
        Tdecom[0] = -qn*pow(ypitch,2)/(6*kt[1]) + At[i]/ypitch + Bt[1] - Tint;
        for (i=1;i<6;i++){
            partloc = rt[i+1];
            Tdecom[i] = -qn*pow(partloc,2)/(6*kt[i]) + At[i]/partloc + Bt[i];
            Tdecom[i] -= Tint;}
        Tdecom[6] = Bt[6] - Tint;
        for (l=0;l<nTL;l++){
            for (i=0;i<11;i++){
            k=1;
            j = i+1;
            if (j>6){
                j=11-j;
                k=-1;}
            currentLoc = locgm + (0.5+l)*ypitch - k*rt[j];
            homoTemp=Am*pow(currentLoc,2)+Bm*currentLoc+Cm;
            THetero[layr][assm][sect][l][i] = homoTemp + Tdecom[j];}}}}}
    //Use linear interpolation with volume-bias to estimate the average temperature of each material in each TRISO
    double lint = 0;
    double lintvol = 0;
    double Aregion1 = 0; //Need two, one for the coolant- and one for the meat-facing halves
    double Bregion1 = 0;
    double inte1 = 0;
    double Aregion2 = 0;
    double Bregion2 = 0;
    double inte2 = 0;
    double intvolume = 0;
    for (layr=1;layr<(nL-1);layr++){
        for (assm=0;assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
                for (l=0;l<nTL;l++){
                    for (i=0;i<5;i++){
                        Aregion1 = THetero[layr][assm][sect][l][5-i] - rt[7-i]*(THetero[layr][assm][sect][l][5-i-1] - THetero[layr][assm][sect][l][5-i])/(rt[6-i]-rt[7-i]);
                        Aregion2 = THetero[layr][assm][sect][l][5+i] - rt[7-i]*(THetero[layr][assm][sect][l][5+i+1] - THetero[layr][assm][sect][l][5+i])/(rt[6-i]-rt[7-i]);
                        Bregion1 = (THetero[layr][assm][sect][l][5-i-1] - THetero[layr][assm][sect][l][5-i])/(rt[6-i]-rt[7-i]);
                        Bregion2 = (THetero[layr][assm][sect][l][5+i+1] - THetero[layr][assm][sect][l][5+i])/(rt[6-i]-rt[7-i]);
                        inte1 = (0.333333333333)*Aregion1*( pow((rt[6-i]),3) - pow((rt[7-i]),3) ) + (0.25)*Bregion1*( pow((rt[6-i]),4) - pow((rt[7-i]),4) );
                        inte2 = (0.333333333333)*Aregion2*( pow((rt[6-i]),3) - pow((rt[7-i]),3) ) + (0.25)*Bregion2*( pow((rt[6-i]),4) - pow((rt[7-i]),4) );
                        intvolume = (0.666666666667)*( pow(rt[6-i],3) - pow(rt[7-i],3) ); //Twice to account for both regions, which were treated like full spheres instead of hemispheres
                        Tregion[layr][assm][sect][l][i] = (inte1 + inte2)/intvolume;
                        //if ( (layr==9) && (assm == 22) && (sect == 2) ) {cout << Aregion1 << " " << Aregion2 << " " << Bregion1 << " " << Bregion2 << " " << inte1 << " " << inte2 << " " << intvolume << " " << Tregion[layr][assm][sect][l][i] << "\n";}
                    }
                    Tregion[layr][assm][sect][l][5] = 0.5*(THetero[layr][assm][sect][l][0] + THetero[layr][assm][sect][l][10]);
                    //if ( (layr==9) && (assm == 22) && (sect == 2) ) {cout << Tregion[layr][assm][sect][l][5] << endl;}
                    }}}}
    //Average materials for all TRISOS
    double matsum = 0;
    double TTriMat[18][84][3][6];
    for (layr=1;layr<(nL-1);layr++){
        for (assm=0;assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
                for (i=0;i<6;i++){
                    matsum = 0;
                    for (l=0;l<nTL;l++){
                        matsum += Tregion[layr][assm][sect][l][i];}
                    TTriMat[layr][assm][sect][i] = matsum/nTL;
                    //if ( (layr==9) && (assm == 22) && (sect == 2) ) {cout << TTriMat[layr][assm][sect][i] << " ";}
                }
                //if ( (layr==9) && (assm == 22) && (sect == 2) ) {cout << endl;}
    }}}
    //Temperatures and temperature differentials (mostly for thermal expansion)
    //Set bottom layer to inlet temperatures
    layr = 0;
    Taxialave[layr] = inletT;
    dTaxialave[layr] = Taxialave[layr] - coldT;
    for (assm=0;assm<nA;assm++){
        for (sect=0;sect<nS;sect++){
            Tflibe[layr][assm][sect] = inletT;
            Tmeat[layr][assm][sect] = inletT; //Graphite meat
            dTmeat[layr][assm][sect] = Tmeat[layr][assm][sect] - coldT;
            Tsleeve[layr][assm][sect] = inletT; //Graphite Sleeve
            dTsleeve[layr][assm][sect] = Tsleeve[layr][assm][sect] - coldT;
            Tspacer[layr][assm][sect] = inletT; //Spacers
            dTspacer[layr][assm][sect] = Tspacer[layr][assm][sect] - coldT;
            Tg[layr][assm][sect] = inletT; // General graphite
            dTg[layr][assm][sect] = Tg[layr][assm][sect] - coldT;
            Tcc[layr][assm][sect] = inletT; //carbon-carbon composite
            dTcc[layr][assm][sect] = Tcc[layr][assm][sect] - coldT;}}
    //Set active core temperatures
    outletT = 0;
    for (layr=1;layr<(nL-1);layr++){
            Tmatave[layr]=0;
            Taxialave[layr]=0;
        for (assm=0;assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
                if (layr == nL-2){outletT+=Tcoolantboundary[layr][assm][sect]/(84*3);}
                Tf[layr][assm][sect] = TTriMat[layr][assm][sect][0]; //fuel
                dTf[layr][assm][sect] = Tf[layr][assm][sect] - coldT;
                Tb[layr][assm][sect] = TTriMat[layr][assm][sect][1]; //buffer
                dTb[layr][assm][sect] = Tb[layr][assm][sect] - coldT;
                Tipyc[layr][assm][sect] = TTriMat[layr][assm][sect][2]; //inner pyrolytic carbon
                dTipyc[layr][assm][sect] = Tipyc[layr][assm][sect] - coldT;
                Tsic[layr][assm][sect] = TTriMat[layr][assm][sect][3]; //sic
                dTsic[layr][assm][sect] = Tsic[layr][assm][sect] - coldT;
                Topyc[layr][assm][sect] = TTriMat[layr][assm][sect][4]; //outer pyrolytic carbon
                dTopyc[layr][assm][sect] = Topyc[layr][assm][sect] - coldT;
                Tmat[layr][assm][sect] = TTriMat[layr][assm][sect][5]; //matrix
                dTmat[layr][assm][sect] = Tmat[layr][assm][sect] - coldT;
                Tmatave[layr] += Tmat[layr][assm][sect]/(84*3); // matrix average
                Tmeat[layr][assm][sect] = Tmatbound[layr][assm][sect][3]; //Graphite meat
                dTmeat[layr][assm][sect] = Tmeat[layr][assm][sect] - coldT;
                Tsleeve[layr][assm][sect] = 0.5*(Tmatbound[layr][assm][sect][1] + Tmatbound[layr][assm][sect][2]); //Graphite Sleeve
                dTsleeve[layr][assm][sect] = Tsleeve[layr][assm][sect] - coldT;
                Tspacer[layr][assm][sect] = Tflibe[layr][assm][sect]; //Spacers
                dTspacer[layr][assm][sect] = Tspacer[layr][assm][sect] - coldT;
                Tg[layr][assm][sect] = Tflibe[layr][assm][sect]; // General graphite
                dTg[layr][assm][sect] = Tg[layr][assm][sect] - coldT;
                Tcc[layr][assm][sect] = Tflibe[layr][assm][sect]; //carbon-carbon composite
                dTcc[layr][assm][sect] = Tcc[layr][assm][sect] - coldT;
                Taxialave[layr] += Tflibe[layr][assm][sect]/(84*3); //Structural graphite average (for use with reflectors)
                Tcb[layr][assm][sect] = Tflibe[layr][assm][sect]; //Control Blade
                dTcb[layr][assm][sect] = Tcb[layr][assm][sect] - coldT;
                Teu[layr][assm][sect] = Tmeat[layr][assm][sect]; //Europia BP spheres
                dTeu[layr][assm][sect] = Teu[layr][assm][sect] - coldT;}}
        dTmatave[layr] = Tmatave[layr] - coldT;
        dTaxialave[layr] = Taxialave[layr] - coldT;}
    //Set top layer to outlet temperatures
    layr = (nL-1);
    Taxialave[layr] = 0;
    for (assm=0;assm<nA;assm++){
        for (sect=0;sect<nS;sect++){
            Tflibe[layr][assm][sect] = Tcoolantboundary[layr-1][assm][sect];
            Tmeat[layr][assm][sect] = Tcoolantboundary[layr-1][assm][sect]; //Graphite meat
            dTmeat[layr][assm][sect] = Tmeat[layr][assm][sect] - coldT;
            Tsleeve[layr][assm][sect] = Tcoolantboundary[layr-1][assm][sect]; //Graphite Sleeve
            dTsleeve[layr][assm][sect] = Tsleeve[layr][assm][sect] - coldT;
            Tspacer[layr][assm][sect] = Tcoolantboundary[layr-1][assm][sect]; //Spacers
            dTspacer[layr][assm][sect] = Tspacer[layr][assm][sect] - coldT;
            Tg[layr][assm][sect] = Tcoolantboundary[layr-1][assm][sect]; // General graphite
            dTg[layr][assm][sect] = Tg[layr][assm][sect] - coldT;
            Tcc[layr][assm][sect] = Tcoolantboundary[layr-1][assm][sect]; //carbon-carbon composite
            dTcc[layr][assm][sect] = Tcc[layr][assm][sect] - coldT;
            Tcb[layr][assm][sect] = Tcoolantboundary[layr-1][assm][sect]; //Control Blade
            dTcb[layr][assm][sect] = Tcb[layr][assm][sect] - coldT;
            Taxialave[layr] += Tcoolantboundary[layr-1][assm][sect]/(84*3);}}
    dTaxialave[layr] = Taxialave[layr] - coldT;
    dTAxial[0]=inletT-coldT;
    dTAxial[1]=outletT-coldT;
    averageT=0.5*(inletT+outletT);
    daverageT = averageT - coldT;
    if ((coldDim) || (firstTH) ){ //If no thermal expansion, set all dT's to zero (but keep temperatures)
        dTAxial[0]=0;
        dTAxial[1]=0;
        for (layr=0;layr<nL;layr++){
            dTmatave[layr]=0;
            dTaxialave[layr]=0;
            for (assm=0; assm<nA;assm++){
                for (sect=0;sect<nS;sect++){
                    dTf[layr][assm][sect]=0; //fuel
                    dTb[layr][assm][sect]=0; //buffer
                    dTipyc[layr][assm][sect]=0; //inner pyrolytic carbon
                    dTsic[layr][assm][sect]=0; //sic
                    dTopyc[layr][assm][sect]=0; //outer pyrolytic carbon
                    dTmat[layr][assm][sect]=0; //matrix
                    dTmeat[layr][assm][sect]=0; //graphite meat
                    dTg[layr][assm][sect]=0; //graphite
                    dTsleeve[layr][assm][sect]=0; //Graphite Sleeve
                    dTspacer[layr][assm][sect]=0; //948 instead of 1110
                    dTcc[layr][assm][sect]=0; //carbon-carbon composite
                    dTcb[layr][assm][sect]=0; //control blade , 948 instead of 1110
                    dTeu[layr][assm][sect]=0;}}}}
}

void THTempsPlot (void){
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    //Initializations to run
        //Select the region of interest
        layr = hotSection[0];
        assm = hotSection[1];
        sect = hotSection[2];
        //Initialize large arrays for data
        double inLow[10000];//Input file array
        double inHigh[10000];//Associated errors with each result
        double inAverage[10000];
        double tempAve[10000];
        double temp[10000];
        double locs[10000];
        double divisions=10000;
        Tcoolantboundary[0][assm][sect] = inletT;
        //Conductivities [W/(m*K)]
        double kgrph = coldkgrph;
        double kmatrix = coldkmatrix;
        double kfuel = coldkfuel;
        double kbuff = coldkbuff;
        double kipyro = coldkipyro;
        double kopyro = coldkopyro;
        double ksic = coldksic;
        double kmtrx = 0;
        double kflibe = coldkflibe;
        //Geometric Parameters [m]
        double rfuel = 0.01*dimFuelRadius[layr][assm][sect];
        double rbuff = 0.01*dimBufferRadius[layr][assm][sect];
        double ripyro = 0.01*dimIpycRadius[layr][assm][sect];
        double rsic = 0.01*dimSicRadius[layr][assm][sect];
        double ropyro = 0.01*dimOpycRadius[layr][assm][sect];
        double requiv = 0.01*pow(((3*dimXpitch[layr][assm][sect]*dimYpitch[layr][assm][sect]*dimZpitch[layr])/(4*pi)),(0.33333333));
        double xpitch = 0.01*dimXpitch[layr][assm][sect];
        double ypitch = 0.01*dimYpitch[layr][assm][sect];
        double zpitch = 0.01*dimZpitch[layr];
        double loccg = 0.01*0.5*(dimSpacer[layr][assm][sect]);
        double locgm = loccg + 0.01*dimColdSleeveWidth*(1+alphag*dTsleeve[layr][assm][sect]);
        double loccenter = 0.01*0.5*(locCP[2][layr][assm][sect] - locCP[1][layr][assm][sect]);
        double channelWidth = loccg;
        double channelLength = 0.01*22.112515*(1+alphacc*dTcc[layr][assm][sect]);
        double nTZ = dimLayrHeight[layr]/dimZpitch[layr];
        double channelPerimeter = 0;
        double Dh = 0;
        double locmc = 0;
        double dyc = 0;
        double dyg = 0;
        double dym = 0;
        double dycenter = 0;
        //Other geometric
        double vol = 0; //m3
        double channelVol = 0;
        double channelArea = 0; //m2
        double TrisoLayers = nTL;
        //Coolant parameters
        double denFlibe = denflibe[layr][assm][sect]; //g/cc
        double viscosityFlibe = flibeVisc; //P*s
        double PrandtlFlibe = flibePr;
        double cpFlibe = flibecp; //J/(kg*K)
        double ReFlibe = 0;
        double NuFlibe = 0;
        double hFlibe = 0; //W/(m2*K)
        double mfrFlibe = coreMFR/(252*12*3); //kg/s in 0.35 channel (coreTotal/(assm*channelPerSection*sect)
        double sectMFRFlibe = coreMFR/(252);
        double spdFlibe = 0; //m/s
        double ff = 0; //friction factor
        //Power parameters
        double powerperpart=0.95*SectionPower[layr][assm][sect][0]/(12*3*nTW*nTL*nTZ); //W
        double power = 0; //W
        double powerdensity = 0; //W/m3
        double heatflux = 0; //W/m2
        double peaking=0;
        //Temperatures [C]
        double dTc=0;
        double dTgrph=0;
        double dTm=0;
        double Tc = 0;
        double Tcg = 0;
        double Tgm = 0;
        double Tmc = 0;
        double Tin=inletT;
        double Tout=outletT;
        double peakTemp = 0;
    //Find homogenized fuel stripe average thermal conductivity (Maxwell's Method)
        double kd = 0; //discontinuous conductivity (particle)
        double kc = 0; //continuous conductivity (surrounding media)
        double pd = 0; //volume fraction
        double km = 0; //mixture conductivity
        //Homogenize fuel into buffer
        kd = kfuel;
        kc = kbuff;
        pd = pow((rfuel/rbuff),3);
        km = kc*(kd + 2*kc - 2*pd*(kc-kd))/(kd+2*kc+pd*(kc-kd));
        //Homogenize fuel/buffer into inner pyrolytic carbon
        kd = km;
        kc = kipyro;
        pd = pow((rbuff/ripyro),3);
        km = kc*(kd + 2*kc - 2*pd*(kc-kd))/(kd+2*kc+pd*(kc-kd));
        //Homogenize f/b/i into silicon carbide
        kd = km;
        kc = ksic;
        pd = pow((ripyro/rsic),3);
        km = kc*(kd + 2*kc - 2*pd*(kc-kd))/(kd+2*kc+pd*(kc-kd));
        //Homogenize f/b/i/s into outer pyrolytic carbon
        kd = km;
        kc = kopyro;
        pd = pow((rsic/ropyro),3);
        km = kc*(kd + 2*kc - 2*pd*(kc-kd))/(kd+2*kc+pd*(kc-kd));
        //Homogenize f/b/i/s/o into matrix carbon
        kd = km;
        kc = kgrph;
        pd = pow((ropyro/requiv),3);
        km = kc*(kd + 2*kc - 2*pd*(kc-kd))/(kd+2*kc+pd*(kc-kd));
        kmtrx = km;
    //Parameters needed to run
        vol=xpitch*ypitch*zpitch;
        channelVol=xpitch*zpitch*channelWidth;
        locmc = TrisoLayers*ypitch+locgm;
        dyc = loccg;
        dyg = locgm - loccg;
        dym = locmc - locgm;
        dycenter = loccenter - locmc;
        channelPerimeter = channelLength+2*channelWidth*(2/sqrt(3));
        channelArea = channelWidth*channelLength;
        Dh=4*channelArea/channelPerimeter;
        power = nTL*powerperpart;
        powerdensity=powerperpart/vol;
        heatflux = power/(xpitch*zpitch);
        spdFlibe = mfrFlibe/(denFlibe*channelArea);
        ReFlibe = denFlibe*spdFlibe*Dh/viscosityFlibe;
        ff=pow((0.79*log(ReFlibe)-1.64),(-2));
        NuFlibe = (ff/8)*(ReFlibe-1000)*(PrandtlFlibe)/(1+12.7*(pow((ff/8),(0.5)))*((pow(PrandtlFlibe,0.666666667))-1));
        hFlibe = kflibe*NuFlibe/Dh;
        dTc = heatflux/hFlibe;
        dTgrph = heatflux*dyg/kgrph;
        dTm = powerdensity*dym*dym/(2*kmtrx);
        Tc = Tflibe[layr][assm][sect];
        Tcg=Tc+dTc;
        Tgm=Tcg+dTgrph;
        Tmc=Tgm+dTm;
    //Average Fuel Profile
        double currentLoc = 0;
        double nvalue=0; //Coolant polynomial fit order
        nvalue = heatflux*dyc/(kflibe*dTc);
        //Functional fit coolant parameters
        double Qc=heatflux/kflibe;
        double Ac=dTc/pow(dyc,nvalue);
        //Parabolic fit matrix parameters
        double Am = dTm/(locmc*locmc-locgm*locgm-2*locmc*dym);
        double Bm = -2*Am*locmc;
        double Cm = Tgm - Am*locgm*locgm - Bm*locgm;
        for(i=0;i<divisions;i++){
            locs[i]=currentLoc;
            if (currentLoc<loccg){temp[i]=Ac*pow(currentLoc,nvalue)+Tc;}
                else if (currentLoc<locgm){temp[i]=(dTgrph/dyg)*currentLoc + (Tcg - (dTgrph/dyg)*loccg);}
                    else if(currentLoc<locmc){temp[i]=Am*pow(currentLoc,2)+Bm*currentLoc+Cm;}
                        else{temp[i]=Tmc;}
            currentLoc += loccenter/divisions;}
        for(i=0;i<divisions;i++){tempAve[i] = temp[i];
        //cout << temp[i] <<  " ";
        }
    //Decomposed fuel profile
        double Tdecom [10000];
        double Tid [7];
        double kt[7];
        double rt[7];
        double At[7];
        double Bt[7];
        if (heteroFuelStripe == true){
            double qp = (power/TrisoLayers)*((3/(4*pi))/pow((rfuel),3)-pow((1/ypitch),3));
            double qn = -(power/TrisoLayers)*pow((1/ypitch),3);
            double partloc = 0;
            double kratio=0;
            kt[0]=kgrph;
            kt[1]=kmatrix;
            kt[2]=kopyro;
            kt[3]=ksic;
            kt[4]=kipyro;
            kt[5]=kbuff;
            kt[6]=kfuel;
            rt[1]=ypitch/2;
            rt[2]=ropyro;
            rt[3]=rsic;
            rt[4]=ripyro;
            rt[5]=rbuff;
            rt[6]=rfuel;
            rt[7]=0; //Center
            //Find expect deltaT across the fuel particle
            At[6]=0;
            Bt[1]=((qn*pow(rt[1],2))/(6*kt[1]))-((At[1])/(rt[1]));
            kratio=kt[6]/kt[5];
            At[5]=pow(rt[6],3)*(kratio*(qp/(3*kt[6]))-(qn/(3*kt[5])));
            for (i=5;i>1;i--){
                kratio=kt[i]/kt[i-1];
                At[i-1]=pow(rt[i],2)*(kratio*(((qn*rt[i])/(3*kt[i]))+((At[i])/(pow(rt[i],2))))-((qn*rt[i])/(3*kt[i-1])));}
            for (i=2;i<6;i++){
                Bt[i]=Bt[i-1] + (1/rt[i])*(At[i-1] - At[i]) - (qn*pow(rt[i],2)/6)*((1/kt[i-1])-(1/kt[i]));
            }
            Bt[6]=Bt[5] + (At[5]/rt[6]) - (pow(rt[6],2)/6)*((qn/kt[5])-(qp/kt[6]));
    //Figure out where you are next and use the right formula
        for (i=0;i<divisions;i++){
            if ((locs[i] > locgm) && (locs[i] < locmc)){ // Only operate in the fuel region
                //Find where in the particle you are
                partloc = locs[i]-locgm-0.5*ypitch;
                for (j=0;j<TrisoLayers;j++){if (partloc > 0.5*ypitch) {partloc -= ypitch;}} //Moves to reference particle space
                if (partloc < 0) {partloc *= -1;} //Move negative values to the positive space
                if (partloc > rt[2]){Tdecom[i] = -qn*pow(partloc,2)/(6*kt[1]) + At[1]/partloc + Bt[1];} //Matrix
                    else if (partloc > rt[3]){Tdecom[i] = -qn*pow(partloc,2)/(6*kt[2]) + At[2]/partloc + Bt[2];}//OPyC
                        else if (partloc > rt[4]){Tdecom[i] = -qn*pow(partloc,2)/(6*kt[3]) + At[3]/partloc + Bt[3];} //Sic
                            else if (partloc > rt[5]){Tdecom[i] = -qn*pow(partloc,2)/(6*kt[4]) + At[4]/partloc + Bt[4];} //IPyC
                                else if (partloc > rt[6]){Tdecom[i] = -qn*pow(partloc,2)/(6*kt[5]) + At[5]/partloc + Bt[5];} //Buff
                                    else {Tdecom[i] = -qp*pow(partloc,2)/(6*kt[6]) + Bt[6];}}}
        double Tint = 0;
        Tint += (-qp*pow(rt[6],5)/(30*kt[6]) + (1/3)*Bt[6]*pow(rt[6],3));
        for(i=5;i>0;i--){
            Tint += ( -qn*pow(rt[i],5)/(30*kt[i]) + 0.5*At[i]*pow(rt[i],2) + (1/3)*Bt[i]*pow(rt[i],3) ) - ( -qn*pow(rt[i+1],5)/(30*kt[i]) + 0.5*At[i]*pow(rt[i+1],2) + (1/3)*Bt[i]*pow(rt[i+1],3) );
        }
        Tint += (2/pi-1/3)*pow(rt[1],3)*(-qn*pow(rt[1],2)/(6*kt[1]) + At[1]/rt[1] + Bt[1]);
        Tint /= 2*pow(rt[1],3)/pi;
        for (i=0;i<divisions;i++){
            if ((locs[i] > locgm) && (locs[i] < locmc)){
                Tdecom[i] -= Tint;
                temp[i]+=Tdecom[i];
                if (temp[i] > peakTemp){peakTemp=temp[i];}}}}
            else{peakTemp = Tmc;}
    //Write results to output file
        th_file.open(THProfileOut);
        th_file << "Highest Power Zone:\n" <<
        "Axial Layer:\t" << hotSection[0] << " of " << (nLt-2) << "\n" <<
        "Assembly Group:\t" << hotSection[1] << " of 84\n" <<
        "Section Number:\t" << hotSection[2] << " of 3\n" <<
        "PPF:\t" << PeakLASPPF << endl <<
        "Peak Temperature:\t" << peakTemp << " K\n\n" <<
        "Loc[m]\tTemp[K]\n";
        for (i=0;i<divisions;i++){th_file << locs[i] << "\t" << temp[i] << endl;}
        th_file.close();
}

void materials (void) {  //Writes all the calculated values to the user-specified file
    string refTemp; // Flibe melts at 459 C (732 K). Assume cases cannot be below 300 K though (even with solid Flibe, allow for it).
    //Choose highest cross section library reference temperature below the inlet temperature (so that it can be used for all materials)
    if (inletT < 300){cout << "***ERROR***: inlet temperature too low, cross sections will have issues.\n";}
        else if (inletT < 600) {
            string refTemp1 = "03";
            refTemp = refTemp1;}
            else if (inletT < 900) {
                string refTemp1 = "06";
                refTemp = refTemp1;}
                else if (inletT < 1200) {
                    string refTemp1 = "09";
                    refTemp = refTemp1;}
                    else { // Flibe boils at 1430 C (1703 K), but not all isotopes have 1500 K data.
                        string refTemp1 = "12";
                        refTemp = refTemp1;}
    mat_file.open(materialsOut);
    mat_file <<
    "% ***** MATERIALS WHICH ONLY REQUIRE A SINGLE DEFINITION *****\n\n" <<
    "% --- Core barrel:\n" <<
    "mat barrel -" << denbarrel << " tms " << inletT << " moder grph1 6012 rgb 50 50 50\n" <<
        "6012." << refTemp << "c 1\n\n" <<
    "% --- Inlet FLiBe:\n" <<
    "mat inletF -" << denaxialaveflibe[0] << " tms " << inletT << " rgb 0 0 100\n" <<
        "3006." << refTemp << "c 1.38322919E-06\n" <<
        "3007." << refTemp << "c 2.37168571E-02\n" <<
        "4009." << refTemp << "c 1.185912017E-02\n" <<
        "9019." << refTemp << "c 4.743648067E-02\n\n" <<
    "% --- Boron Carbide (Density is slightly lower than Theo. Dens. of 2.52 g/cc):\n" <<
    "mat B4C -" << denBC << " tms " << inletT << " rgb 100 255 100\n" <<
        "5010." << refTemp << "c 0.796\n" <<
        "5011." << refTemp << "c 3.204\n" <<
        "6012." << refTemp << "c 1\n\n" <<
    "% --- Reactor Vessel Liner (Alloy-N):\n" <<
    "% Data obtained from: http://www.haynesintl.com/alloys/alloy-portfolio_/Corrosion-resistant-Alloys/hastelloy-n-alloy/nominal-composition\n" <<
    "mat alloyN -" << denAlloyN << " tms " << inletT <<  " rgb 100 100 100\n" <<
        "28000." << refTemp << "c -69.09\n" <<
        "24000." << refTemp << "c -7\n" <<
        "42000." << refTemp << "c -16\n" <<
        "26000." << refTemp << "c -4\n" <<
        "14000." << refTemp << "c -1\n" <<
        "25055." << refTemp << "c -0.8\n" <<
        "23000." << refTemp << "c -0.5\n" <<
        "6012." << refTemp << "c -0.06\n" <<
        "27059." << refTemp << "c -0.2\n" <<
        "29000." << refTemp << "c -0.35\n" <<
        "74000." << refTemp << "c -0.5\n" <<
        "13027." << refTemp << "c -0.25\n" <<
        "22000." << refTemp << "c -0.25\n\n" <<
    "% --- Reactor Vessel (800-H Alloy):\n" <<
    "% Data obtained from: https://www.corrosionmaterials.com/documents/dataSheet/alloy800DataSheet.pdf\n" <<
    "mat vessel -" << denH800 << " tms " << inletT << " rgb 150 150 150\n" <<
        "28000." << refTemp << "c -32.5\n" <<
        "24000." << refTemp << "c -21\n" <<
        "6012." << refTemp << "c -0.075\n" <<
        "25055." << refTemp << "c -1.5\n" <<
        "16000." << refTemp << "c -0.015\n" <<
        "14000." << refTemp << "c -1\n" <<
        "29000." << refTemp << "c -0.75\n" <<
        "15031." << refTemp << "c -0.045\n" <<
        "13027." << refTemp << "c -0.375\n" <<
        "22000." << refTemp << "c -0.375\n" <<
        "26000." << refTemp << "c -42.365\n\n" <<
    "% --- Axial Reflectors\n" <<
    "mat botSuppPlate -" << denSupportPlate[0] << " tms " << inletT << " rgb 255 0 255\n" << //might need to change density to something else later
        "6012." << refTemp << "c 0.034256954\n" <<
        "3006." << refTemp << "c 1.82366E-06\n" <<
        "3007." << refTemp << "c 3.12685E-02\n" <<
        "4009." << refTemp << "c 1.56352E-02\n" <<
        "9019." << refTemp << "c 6.25407E-02\n\n" <<
    "mat topSuppPlate -" << denSupportPlate[1] << " tms " << outletT << " rgb 225 0 225\n" << //might need to change density to something else later
        "14028." << refTemp << "c 2.62835550e-1\n" <<
        "14029." << refTemp << "c 1.33522500e-2\n" <<
        "14030." << refTemp << "c 8.81220000e-3\n" <<
        "6012." << refTemp << "c 2.85000000e-1\n" <<
        "9019." << refTemp << "c 2.45714286e-1\n" <<
        "3006." << refTemp << "c 7.16493230e-6\n" <<
        "3007." << refTemp << "c 1.22849978e-1\n" <<
        "4009." << refTemp << "c 6.14285714e-2\n\n" <<
    "% ***** MATERIALS USING A RADIAL AVERAGE *****\n\n" <<
    "% --- Radial Reflector\n";
    if (singleGraph == false){
        for (layr=0; layr < nL; layr++){
            mat_file << "mat refl_" << layr << " -" << denaxialavegraph[layr] << " tms " << Taxialave[layr] << " moder grph1 6012 rgb 50 50 50\n" <<
            "6012." << refTemp << "c 1\n\n";}}
        else {
            for (layr=0; layr < nL; layr++){
                mat_file << "mat refl_" << layr << " -" << denColdGraphite  << " tms " << inletT << " moder grph1 6012 rgb 50 50 50\n" << "6012." << refTemp << "c 1\n\n";}}
    mat_file << "% ***** MATERIALS USING 1/3 SECTION RESOLUTION *****\n\n" <<
    "% --- UCO Fuel 9% enrichment\n";
    double volfuel = 0;
    double totvolfuel = 0;
    if (singleFuel == false){
        for(layr=1;layr<nL-1;layr++){
            for(assm=0;assm<nA;assm++){
                for(sect=0;sect<nS;sect++){
                    volfuel = (4.0/3.0)*pi*pow(dimFuelRadius[layr][assm][sect],3)*nTL*nTW*2*6*3*(dimLayrHeight[layr]/(dimZpitch[layr]));
                    volumeFuel[layr][assm][sect]= volfuel;
                    totvolfuel += volfuel;
                    mat_file <<
                    "mat fuel_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denFuel[layr][assm][sect] << " tms " << Tf[layr][assm][sect] << " burn 1 vol " << volfuel << " rgb 255 100 100\n" << //Get volume correct (WRT stripe size, expansion, etc._\)
                    "92235." << refTemp << "c 2.25731222E-03\n" <<
                    "92238." << refTemp << "c 2.28239347E-02\n" <<
                    "8016." << refTemp << "c  3.28882742E-02\n" <<
                    "6012." << refTemp << "c  0.011192585\n\n";}}}}
        else{
            volfuel = (4.0/3.0)*pi*dimColdFuelR*dimColdFuelR*dimColdFuelR*nTL*nTW*2*6*3*(dimActiveCoreHeight/(dimColdZPitch*(nL-2)))*(nL-2)*nA*nS;
            mat_file <<
            "mat fuel -10.9 tms 1110 burn 1 vol " << volfuel << " rgb 255 100 100\n" <<
            "92235." << refTemp << "c 2.25731222E-03\n" <<
            "92238." << refTemp << "c 2.28239347E-02\n" <<
            "8016." << refTemp << "c  3.28882742E-02\n" <<
            "6012." << refTemp << "c  0.011192585\n\n";}
    mat_file << "% --- MATERIALS FOR TRISO PARTICLES LAYERS BEYOND FUEL\n\n" <<
    "% --- Buffer\n";
    if (singleTrisoL == false){
        for(layr=1;layr<nL-1;layr++){
            for(assm=0;assm<nA;assm++){
                for(sect=0;sect<nS;sect++){
                    mat_file <<
                    "mat buffer_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denBuffer[layr][assm][sect] << " rgb 0 255 255 tms " << Tb[layr][assm][sect] << "\n" <<
                    "6012." << refTemp << "c 1\n\n";}}}}
        else{
            mat_file << "mat buffer -1.0 rgb 0 255 255 tms 1110\n" <<
            "6012." << refTemp << "c 1\n\n";}
    mat_file << "% --- Inner Pyrolytic Carbon\n";
    if (singleTrisoL == false){
        for(layr=1;layr<nL-1;layr++){
            for(assm=0;assm<nA;assm++){
                for(sect=0;sect<nS;sect++){
                    mat_file <<
                    "mat IPyC_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denIpyc[layr][assm][sect] << " moder grph1 6012 rgb 0 255 255 tms " << Tipyc[layr][assm][sect] << "\n" <<
                    "6012." << refTemp << "c 1\n\n";}}}}
        else{
            mat_file << "mat IPyC -1.9 moder grph1 6012 rgb 0 255 255 tms 1110\n" <<
            "6012." << refTemp << "c 1\n\n";}
    mat_file << "% --- Silicon Carbide\n";
    if (singleTrisoL == false) {
        for(layr=1;layr<nL-1;layr++){
            for(assm=0;assm<nA;assm++){
                for(sect=0;sect<nS;sect++){
                    mat_file <<
                    "mat SiC_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denSic[layr][assm][sect] << " rgb 0 255 255 tms " << Tsic[layr][assm][sect] << "\n" <<
                    "14028." << refTemp << "c 4.43237548E-02\n" <<
                    "14029." << refTemp << "c 2.25168115E-03\n" <<
                    "14030." << refTemp << "c 1.48606150E-03\n" <<
                    "6012." << refTemp << "c  4.806149742E-02\n\n";}}}}
        else{
            mat_file << "mat SiC -3.1 rgb 0 255 255 tms 1110\n" <<
            "14028." << refTemp << "c 4.43237548E-02\n" <<
            "14029." << refTemp << "c 2.25168115E-03\n" <<
            "14030." << refTemp << "c 1.48606150E-03\n" <<
            "6012." << refTemp << "c  4.806149742E-02\n\n";}
    mat_file << "% --- Outer Pyrolytic Carbon\n";
    if (singleTrisoL == false){
        for(layr=1;layr<nL-1;layr++){
            for(assm=0;assm<nA;assm++){
                for(sect=0;sect<nS;sect++){
                    mat_file <<
                    "mat OPyC_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denOpyc[layr][assm][sect] << " moder grph1 6012 rgb 0 255 255 tms " << Topyc[layr][assm][sect] << "\n" <<
                    "6012." << refTemp << "c 1\n\n";}}}}
        else{
            mat_file << "mat OPyC -1.87 moder grph1 6012 rgb 0 255 255 tms 1110\n" <<
            "6012." << refTemp << "c 1\n\n";}
    mat_file << "% --- Matrix\n";
    if (singleTrisoL == false){
        for(layr=1;layr<nL-1;layr++){
            for(assm=0;assm<nA;assm++){
                for(sect=0;sect<nS;sect++){
                    mat_file <<
                    "mat matrix_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denMat[layr][assm][sect] << " tms " << Tmat[layr][assm][sect] << " moder grph1 6012 rgb 0 255 255\n" <<
                    "6012." << refTemp << "c 1\n\n";}}}}
        else {
            mat_file << "mat matrix -1.75 tms 1110 moder grph1 6012 rgb 0 255 255\n" <<
            "6012." << refTemp << "c 1\n\n";}
    mat_file << "% --- GRAPHITE MATERIALS\n\n";
    if (singleGraph == false){
        "% --- Meat Graphite\n";
        for(layr=0;layr<nL;layr++){
            for(assm=0;assm<nA;assm++){
                for(sect=0;sect<nS;sect++){
                    mat_file << "mat meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denMeat[layr][assm][sect] << " tms " << Tmeat[layr][assm][sect] << " moder grph1 6012 rgb 255 0 0\n" <<
                    "6012." << refTemp << "c 1\n\n";}}}
        mat_file << "% --- Stripe Surrounding Sleeve:\n";
        for(layr=0;layr<nL;layr++){
            for(assm=0;assm<nA;assm++){
                for(sect=0;sect<nS;sect++){
                    mat_file << "mat sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denSleeve[layr][assm][sect] << " tms " << Tsleeve[layr][assm][sect] << " moder grph1 6012 rgb 150 0 0\n" <<
                    "6012." << refTemp << "c 1\n\n";}}}
        mat_file << "% --- Spacers\n";
        for(layr=0;layr<nL;layr++){
                for(assm=0;assm<nA;assm++){
                    for(sect=0;sect<nS;sect++){
                        mat_file << "mat spacer_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denSpacer[layr][assm][sect] << " tms " << Tspacer[layr][assm][sect] << " moder grph1 6012 rgb 235 235 235\n" <<
                        "6012." << refTemp << "c 1\n\n";}}}
        mat_file << "% --- Assembly Structural Materials:\n";
        for(layr=0;layr<nL;layr++){
                for(assm=0;assm<nA;assm++){
                    for(sect=0;sect<nS;sect++){
                        mat_file << "mat structural_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denCCC[layr][assm][sect] << " tms " << Tcc[layr][assm][sect] << " rgb 150 150 150\n" <<
                        "6012." << refTemp << "c 9.82741E-02\n\n";}}}}
        else {
            mat_file <<
            "% --- Meat Graphite\n" <<
            "mat meat -1.75 tms 1110 moder grph1 6012 rgb 255 0 0\n" <<
            "6012." << refTemp << "c 1\n\n" <<
            "% --- Stripe Surrounding Sleeve:\n" <<
            "mat sleeve -1.75 tms 948  moder grph1 6012 rgb 150 0 0\n" <<
            "6012." << refTemp << "c 1\n\n" <<
            "% --- Spacers\n" <<
            "mat spacer -1.75 tms 1029 moder grph1 6012 rgb 235 235 235\n" <<
            "6012." << refTemp << "c 1\n\n" <<
            "% --- Assembly Structural Materials:\n" <<
            "mat structural -1.95 tms 948 rgb 150 150 150\n" <<
            "6012." << refTemp << "c 9.82741E-02\n\n";}
    mat_file << "% --- Burnable Poison\n";
    double volbp = 0;
    if (singleBP == false){
        for(layr=1;layr<nL-1;layr++){
            for(assm=0;assm<nA;assm++){
                for(sect=0;sect<nS;sect++){
                    volbp = (4.0/3.0)*pi*pow(dimEuRadius[layr][assm][sect],3)*nBPC*(dimLayrHeight[layr]/dimZEupitch[layr])*6*3;
                    volumeBP[layr][assm][sect] = volbp;
                    mat_file << "mat bPoison_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denEu[layr][assm][sect] << " tms " << Teu[layr][assm][sect] << " burn " << BPBurnZones << " vol ‭" << volbp/BPBurnZones << " rgb 255 50 255\n" << //correct volume
                    "63151." << refTemp << "c 8.181202072E-03\n" <<
                    "63153." << refTemp << "c 8.930703538E-03\n" <<
                    "8016." << refTemp << "c  2.566785841E-02\n\n";}}}}
        else{
            volbp = (4.0/3.0)*pi*pow(dimColdEuR,3)*nBPC*(dimActiveCoreHeight/(dimColdZEuPitch*(nL-2)))*6*3*(nL-2)*nA*nS;
            mat_file << "mat bPoison -" << denColdEu << " tms 1110 burn " << BPBurnZones << " vol " << volbp/BPBurnZones << "‭ rgb 255 50 255\n" <<
            "63151." << refTemp << "c 8.181202072E-03\n" <<
            "63153." << refTemp << "c 8.930703538E-03\n" <<
            "8016." << refTemp << "c  2.566785841E-02\n\n";}
    mat_file << "% --- Control Blade:\n";
    if (singleCB == false){
        for(layr=1;layr<nL;layr++){
            for(assm=0;assm<nA;assm++){
                for(sect=0;sect<nS;sect++){
                    mat_file << "mat cBlade_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denCB[layr][assm][sect] << " tms " << Tcb[layr][assm][sect] << " rgb 0 50 25\n" <<
                    "6012." << refTemp << "c  5.158967334e-4\n" <<
                    "72174." << refTemp << "c 6.659530024e-7\n" <<
                    "72176." << refTemp << "c 2.189320495e-5\n" <<
                    "72177." << refTemp << "c 7.741703653e-5\n" <<
                    "72178." << refTemp << "c 1.135449869e-4\n" <<
                    "72179." << refTemp << "c 5.668924933e-5\n" <<
                    "72180." << refTemp << "c 1.460101958e-4\n" <<
                    "42092." << refTemp << "c 9.252040149e-3\n" <<
                    "42094." << refTemp << "c 5.826301925e-3\n" <<
                    "42095." << refTemp << "c 1.008618830e-2\n" <<
                    "42096." << refTemp << "c 1.061469438e-2\n" <<
                    "42097." << refTemp << "c 6.112841393e-3\n" <<
                    "42098." << refTemp << "c 1.553043766e-2\n" <<
                    "42100." << refTemp << "c 6.252927341e-3\n\n";}}}}
        else{
            mat_file << "mat cBlade -" << denColdCB << " tms 1110 rgb 0 50 25\n" <<
            "6012." << refTemp << "c  5.158967334e-4\n" <<
            "72174." << refTemp << "c 6.659530024e-7\n" <<
            "72176." << refTemp << "c 2.189320495e-5\n" <<
            "72177." << refTemp << "c 7.741703653e-5\n" <<
            "72178." << refTemp << "c 1.135449869e-4\n" <<
            "72179." << refTemp << "c 5.668924933e-5\n" <<
            "72180." << refTemp << "c 1.460101958e-4\n" <<
            "42092." << refTemp << "c 9.252040149e-3\n" <<
            "42094." << refTemp << "c 5.826301925e-3\n" <<
            "42095." << refTemp << "c 1.008618830e-2\n" <<
            "42096." << refTemp << "c 1.061469438e-2\n" <<
            "42097." << refTemp << "c 6.112841393e-3\n" <<
            "42098." << refTemp << "c 1.553043766e-2\n" <<
            "42100." << refTemp << "c 6.252927341e-3\n\n";}
    mat_file << "% --- FLiBe Coolant:\n";
    if (singleFlibe == false){
        for(layr=0;layr<nL;layr++){
            for(assm=0;assm<nA;assm++){
                for(sect=0;sect<nS;sect++){
                    mat_file << "mat flibe_" << layr << "_" << assm+1 << "_" << sect+1 << " -" << denflibe[layr][assm][sect] << " tms " << Tflibe[layr][assm][sect] << " rgb 0 0 125\n" <<
                    "3006." << refTemp << "c 1.38322919E-06\n" <<
                    "3007." << refTemp << "c 2.37168571E-02\n" <<
                    "4009." << refTemp << "c 1.185912017E-02\n" <<
                    "9019." << refTemp << "c 4.743648067E-02\n\n";}}}}
        else{
            mat_file << "mat flibe -1.9503 tms 948 rgb 0 0 125\n" <<
            "3006." << refTemp << "c 1.38322919E-06\n" <<
            "3007." << refTemp << "c 2.37168571E-02\n" <<
            "4009." << refTemp << "c 1.185912017E-02\n" <<
            "9019." << refTemp << "c 4.743648067E-02\n\n";}
    mat_file << "% --- Thermal Scattering Data\n" <<
    "therm grph1 0 gre7.00t gre7.04t  gre7.08t  gre7.12t gre7.16t gre7.18t gre7.20t gre7.22t gre7.24t gre7.26t\n\n";
    //Declare the specified materials for equilibrium xenon treatment in physics file
    if (equilXenon == true){mat_file << "set xenon 1\n\n";}
    //Write out neutron history parameters
    mat_file << "set powdens " << corePowDens << endl;
    if (!depMode){mat_file << "\nset pop " << kCycleSize << " " << kActiveCycles << " " << kInactiveCycles << endl;}
        else{
            if (depStep <= depStepTotal){ //If step > total, then withhold set pop so that SERPENT fails to launch.
                if (!nextDepStep){ //If still doing CB movement search
                    mat_file << "\nset pop " << kCycleSize << " " << kActiveCycles << " " << kInactiveCycles << endl;
                    if (depStep > 0){mat_file << "\nset rfr idx 1 " << "Restart_" << depStep;}}
                    else{
                        mat_file << "\nset pop " << depCycleSize << " " << depActiveCycles << " " << depInactiveCycles << endl;
                        mat_file << "\nset rfw 1 Restart_" << depStep;
                        if (depStep > 1){mat_file << "\nset rfr idx 1 " << "Restart_" << depStep-1;}
                        mat_file << "\ndep ";
                        if (depType == 1) {mat_file << "bustep ";}
                            else if (depType == 2) {mat_file << "butot ";}
                                else if (depType == 3) {mat_file << "daystep ";}
                                    else if (depType == 4) {mat_file << "daytot ";}
                                        else if (depType == 5) {mat_file << "decstep ";}
                                            else if (depType == 6) {mat_file << "dectot ";}
                        mat_file << depArray[depStep-1];}}}
    mat_file << endl;
    mat_file.close();
    cout << "\nThe file " << materialsOut << " has been written with the materials for the Serpent simulation.\n";
}

void geoOut (void) {  //Writes all the calculated values to the user-specified file
    int i = 0;
    int j = 0;
    bool singleTriso = false;
    if (singleTrisoL && singleFuel){singleTriso = true;}
    out_file.open(geometryOut);
    out_file <<
    "% ---------------------------------------\n" <<
    "% ----Geometry---------------------------\n" <<
    "% ---------------------------------------\n\n" <<
    "% --- Core Slice and Whole Core Surfaces\n" << // Likely will use the same dim for each layer (since the cooler downcomer on the far side makes it difficult. Likely just expand to average temperature. Extra tall since will just be cut by the layers later.
    "surf 60 cylz 0 0 " << dimPermReflRadius << " 0 1000 %Permanent Reflector Outer (including radial reflector)\n" <<
    "surf 61 cylz 0 0 " << dimBoronCarbideRadius << " 0 1000 %Boron Carbide\n" <<
    "surf 62 cylz 0 0 " << dimCoreBarrelRadius << " 0 1000 %Core barrel outer diameter\n" <<
    "surf 63 cylz 0 0 " << dimDowncomerRadius << " 0 1000 %Downcomer\n" <<
    "surf 64 cylz 0 0 " << dimAlloyNLinerRadius << " 0 1000 %Alloy N Liner\n" <<
    "surf 65 cylz 0 0 " << dimPressureVesselRadius << " 0 " << dimModelHeight << " %Whole Problem Boundary\n" <<
    "surf bigS cylz 0 0 100 0 1000 %Arbitrary Big Surface\n\n" <<
    "% --- Axial Reflectors\n" <<
    "pin botSupportPlate\n" <<
    "botSuppPlate\n\n" <<
    "pin topSupportPlate\n" <<
    "topSuppPlate\n\n"<<
    "% --- Radial Reflectors\n";
    for (layr=0;layr<nL;layr++){out_file << "pin 50_" << layr << "\nrefl_" << layr << "\n\n";} //Permanent reflector
    out_file << "surf 21 cyl 0 0 " << dimReflHoleRadius << " %Removable Reflector Hole\n" <<
    "surf 22 hexyprism 0 0 " << dimReflHex << " 0 1000 %Removable Reflector Hex (Since along periphery, assume at inlet)\n\n";
    for (layr=0;layr<nL;layr++){
        if (reflHole == true){ //Assume axial average graphite expansion for removable reflector assemblies
            out_file <<
            "cell RR0_" << layr << " 51_" << layr << " inletF -21\n" <<
            "cell RR1_" << layr << " 51_" << layr << " refl_" << layr << " 21 -22\n" <<
            "cell RR2_" << layr << " 51_" << layr << " inletF 22\n\n";}
            else{
                out_file <<
                "cell RR1_" << layr << " 51_" << layr << " refl_" << layr << " -22\n" <<
                "cell RR2_" << layr << " 51_" << layr << " inletF 22\n\n";}}
    out_file <<
    "% --- Fixes planes for dividing assemblies into thirds\n" <<
    "surf div1 py 0\n" <<
    "surf div2 plane -1.73205081 1 0 0\n" <<
    "surf div3 plane 1.73205081 1 0 0\n\n";
    if (cubic == false){
        out_file << "% --- Buffer Particles for lattices (prevents small undefined regions due to rounding)\n";
        for (layr=1; layr<(nL-1);layr++){
            for (assm=0; assm<nA;assm++){
                for (sect=0;sect<nS;sect++){
                    if (singleTrisoL) {
                        out_file <<
                        "particle 98_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "matrix\n\n";}
                        else {
                            out_file <<
                            "particle 98_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                            "matrix_" << layr << "_" << assm+1 << "_" << sect+1 << "\n\n";}}}}}
    out_file <<
    "% --- TRISO Particles\n";
    layr = 1;
    assm = 1;
    sect = 1;
    if (singleTriso){
        out_file <<
        "particle triso\n" <<
        "fuel " << dimFuelRadius[layr][assm][sect] << endl <<
        "buffer " << dimBufferRadius[layr][assm][sect] << endl <<
        "IPyC " << dimIpycRadius[layr][assm][sect] << endl <<
        "SiC " << dimSicRadius[layr][assm][sect] << endl <<
        "OPyC " << dimOpycRadius[layr][assm][sect] << endl <<
        "matrix\n\n";}
        else{
            for (layr=1; layr<(nL-1);layr++){
                for (assm=0; assm<nA;assm++){
                    for (sect=0;sect<nS;sect++){
                        out_file <<
                        "particle triso_" << layr << "_" << assm+1 << "_" << sect+1 << endl;
                        if (singleFuel) {out_file << "fuel " << dimFuelRadius[layr][assm][sect] << endl;}
                            else {out_file << "fuel_" << layr << "_" << assm+1 << "_" << sect+1 << " " << dimFuelRadius[layr][assm][sect] << endl;}
                        if (singleTrisoL) {
                            out_file <<
                            "buffer " << dimBufferRadius[layr][assm][sect] << endl <<
                            "IPyC " << dimIpycRadius[layr][assm][sect] << endl <<
                            "SiC " << dimSicRadius[layr][assm][sect] << endl <<
                            "OPyC " << dimOpycRadius[layr][assm][sect] << endl <<
                            "matrix\n\n";}
                            else{
                                out_file <<
                                "buffer_" << layr << "_" << assm+1 << "_" << sect+1 << " " << dimBufferRadius[layr][assm][sect] << endl <<
                                "IPyC_" << layr << "_" << assm+1 << "_" << sect+1 << " " << dimIpycRadius[layr][assm][sect] << endl <<
                                "SiC_" << layr << "_" << assm+1 << "_" << sect+1 << " " << dimSicRadius[layr][assm][sect] << endl <<
                                "OPyC_" << layr << "_" << assm+1 << "_" << sect+1 << " " << dimOpycRadius[layr][assm][sect] << endl <<
                                "matrix_" << layr << "_" << assm+1 << "_" << sect+1 << "\n\n";}}}}}
    bool singleEu = false;
    if (singleBP && singleGraph) {singleEu = true;}
    out_file << "% --- Europium Spheres\n";
    if (BPUsage){
        if (singleEu){
            out_file <<
            "particle bp\n" <<
            "bPoison " << dimEuRadius[layr][assm][sect] << endl <<
            "meat\n\n";}
            else{
                for (layr=1; layr<(nL-1);layr++){
                    for (assm=0; assm<nA;assm++){
                        for (sect=0;sect<nS;sect++){
                            out_file <<
                            "particle bp_" << layr << "_" << assm+1 << "_" << sect+1 << endl;
                            if (singleBP){out_file << "bPoison " << dimEuRadius[layr][assm][sect] << endl;}
                                else {out_file << "bPoison_" << layr << "_" << assm+1 << "_" << sect+1 << " " << dimEuRadius[layr][assm][sect] << endl;}
                            if (singleGraph){out_file << "meat\n\n";}
                                else {out_file << "meat_" << layr << "_" << assm+1 << "_" << sect+1 << "\n\n";}}}}}}
    out_file << "% --- Fuel Kernel Lattices\n"; //Make fuel lattices large then just cut them to the desired size (max 225 wide in x)
    if (cubic){
        int strp = 0; //Stripe index
        for (layr=1; layr<(nL-1);layr++){
            for (assm=0; assm<nA;assm++){
                for (sect=0;sect<nS;sect++){
                    for (strp=1;strp<13;strp++){
                    out_file << "lat stripe" << strp << "_" << layr << "_" << assm+1 << "_" << sect+1 << " 6 " << locFLx[strp][layr][assm][sect]-0.5*dimXpitch[layr][assm][sect] << " " << locFLy[strp][layr][assm][sect]-0.5*dimYpitch[layr][assm][sect] << " " << dimZpitch[layr] << " ";
                    if (singleTriso){out_file << "triso";}
                        else {out_file << "triso_" << layr << "_" << assm+1 << "_" << sect+1;}
                    out_file << endl;
                }
                out_file << endl;}}}}
        else {
            for (layr=1; layr<(nL-1);layr++){
                for (assm=0; assm<nA;assm++){
                    for (sect=0;sect<nS;sect++){
                        out_file <<
                        "lat 101_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 0 0 0 15 " << nTL+2 << " 1 " << dimXpitch[layr][assm][sect] << " " << dimYpitch[layr][assm][sect] << " " << dimZpitch[layr] << endl;
                        for (i=0;i<15;i++){out_file << "98_" << layr << "_" << assm+1 << "_" << sect+1 << " ";}
                        out_file << endl;
                        for (j=0;j<nTL;j++){
                            for (i=0;i<15;i++){
                                if (singleTriso){out_file << "triso ";}
                                    else {out_file << "triso_" << layr << "_" << assm+1 << "_" << sect+1 << " ";}}
                            out_file << endl;}
                        for (i=0;i<15;i++){out_file << "98_" << layr << "_" << assm+1 << "_" << sect+1 << " ";}
                        out_file << "\n\n";
                        }}}
            for (layr=1; layr<(nL-1);layr++){
                for (assm=0; assm<nA;assm++){
                    for (sect=0;sect<nS;sect++){
                        out_file <<
                        "lat 102_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 0 0 0 15 1 1 " << 15*dimXpitch[layr][assm][sect] << " " << (nTL+2)*dimYpitch[layr][assm][sect] << " " << dimZpitch[layr] << endl;
                        for (i=0;i<15;i++){out_file << "101_" << layr << "_" << assm+1 << "_" << sect+1 << " ";}
                        out_file << "\n\n";
                        }}}
            for (layr=1; layr<(nL-1);layr++){
                for (assm=0; assm<nA;assm++){
                    for (sect=0;sect<nS;sect++){
                        out_file <<
                        "lat 103_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << 0.5*nTW*dimXpitch[layr][assm][sect] << " " << 0.5*nTL*dimYpitch[layr][assm][sect] << " 0 1 1 20 " << nTW*dimXpitch[layr][assm][sect] << " " << nTL*dimYpitch[layr][assm][sect] << " " << dimZpitch[layr] << endl;
                        for (i=0;i<20;i++){out_file << "102_" << layr << "_" << assm+1 << "_" << sect+1 << " ";}
                        out_file << "\n\n";
                        }}}
            for (layr=1; layr<(nL-1);layr++){
                for (assm=0; assm<nA;assm++){
                    for (sect=0;sect<nS;sect++){
                        out_file <<
                        "lat 104_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 0 0 " << 185*dimZpitch[layr] << " 1 1 19 " << nTW*dimXpitch[layr][assm][sect] << " " << nTL*dimYpitch[layr][assm][sect] << " " << 20*dimZpitch[layr] << endl;
                        for (i=0;i<19;i++){out_file << "103_" << layr << "_" << assm+1 << "_" << sect+1 << " ";}
                        out_file << "\n\n";
                        }}}
            out_file << "% --- Shift lattices left and down one unit so that no voids are present when stripes are cut \n";
            for (layr=1; layr<(nL-1);layr++){
                for (assm=0; assm<nA;assm++){
                    for (sect=0;sect<nS;sect++){
                        out_file <<
                        "lat stripe1_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[1][layr][assm][sect] << " " << locFLy[1][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "lat stripe2_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[2][layr][assm][sect] << " " << locFLy[2][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "lat stripe3_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[3][layr][assm][sect] << " " << locFLy[3][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "lat stripe4_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[4][layr][assm][sect] << " " << locFLy[4][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "lat stripe5_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[5][layr][assm][sect] << " " << locFLy[5][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "lat stripe6_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[6][layr][assm][sect] << " " << locFLy[6][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "lat stripe7_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[7][layr][assm][sect] << " " << locFLy[7][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "lat stripe8_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[8][layr][assm][sect] << " " << locFLy[8][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "lat stripe9_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[9][layr][assm][sect] << " " << locFLy[9][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "lat stripe10_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[10][layr][assm][sect] << " " << locFLy[10][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "lat stripe11_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[11][layr][assm][sect] << " " << locFLy[11][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "lat stripe12_" << layr << "_" << assm+1 << "_" << sect+1 << " 11 " << locFLx[12][layr][assm][sect] << " " << locFLy[12][layr][assm][sect] << " 0 1 1 1 20 2 36 104_" << layr << "_" << assm+1 << "_" << sect+1 << "\n\n";
                        }}}}
    if (BPUsage == true){
    out_file << "% --- Poison Lattices\n";
    for (layr=1; layr<(nL-1);layr++){
        for (assm=0; assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
                out_file << "lat BP1_" << layr << "_" << assm+1 << "_" << sect+1 << " 6 " << locBP[1][layr][assm][sect] << " " << locCP[1][layr][assm][sect] << " " << dimZEupitch[layr];
                if (singleEu) {out_file << " bp\n";}
                    else{out_file << " bp_" << layr << "_" << assm+1 << "_" << sect+1 << endl;}
                out_file << "lat BP2_" << layr << "_" << assm+1 << "_" << sect+1 << " 6 " << locBP[2][layr][assm][sect] << " " << locCP[2][layr][assm][sect] << " " << dimZEupitch[layr];
                if (singleEu) {out_file << " bp\n";}
                    else{out_file << " bp_" << layr << "_" << assm+1 << "_" << sect+1 << endl;}
                out_file << "lat BP3_" << layr << "_" << assm+1 << "_" << sect+1 << " 6 " << locBP[3][layr][assm][sect] << " " << locCP[3][layr][assm][sect] << " " << dimZEupitch[layr];
                if (singleEu) {out_file << " bp\n";}
                    else{out_file << " bp_" << layr << "_" << assm+1 << "_" << sect+1 << endl;}
                out_file << "lat BP4_" << layr << "_" << assm+1 << "_" << sect+1 << " 6 " << locBP[4][layr][assm][sect] << " " << locCP[4][layr][assm][sect] << " " << dimZEupitch[layr];
                if (singleEu) {out_file << " bp\n";}
                    else{out_file << " bp_" << layr << "_" << assm+1 << "_" << sect+1 << endl;}
                out_file << "lat BP5_" << layr << "_" << assm+1 << "_" << sect+1 << " 6 " << locBP[5][layr][assm][sect] << " " << locCP[5][layr][assm][sect] << " " << dimZEupitch[layr];
                if (singleEu) {out_file << " bp\n";}
                    else{out_file << " bp_" << layr << "_" << assm+1 << "_" << sect+1 << endl;}
                out_file << "lat BP6_" << layr << "_" << assm+1 << "_" << sect+1 << " 6 " << locBP[6][layr][assm][sect] << " " << locCP[6][layr][assm][sect] << " " << dimZEupitch[layr];
                if (singleEu) {out_file << " bp\n";}
                    else{out_file << " bp_" << layr << "_" << assm+1 << "_" << sect+1 << endl;}
                out_file <<
                "utrans BP1_" << layr << "_" << assm+1 << "_" << sect+1 << " 0 0 -" << 0.5*dimZEupitch[layr] << " 0 0 0\n" <<
                "utrans BP2_" << layr << "_" << assm+1 << "_" << sect+1 << " 0 0 -" << 0.5*dimZEupitch[layr] << " 0 0 0\n" <<
                "utrans BP3_" << layr << "_" << assm+1 << "_" << sect+1 << " 0 0 -" << 0.5*dimZEupitch[layr] << " 0 0 0\n" <<
                "utrans BP4_" << layr << "_" << assm+1 << "_" << sect+1 << " 0 0 -" << 0.5*dimZEupitch[layr] << " 0 0 0\n" <<
                "utrans BP5_" << layr << "_" << assm+1 << "_" << sect+1 << " 0 0 -" << 0.5*dimZEupitch[layr] << " 0 0 0\n" <<
                "utrans BP6_" << layr << "_" << assm+1 << "_" << sect+1 << " 0 0 -" << 0.5*dimZEupitch[layr] << " 0 0 0\n\n";
                }}}
    out_file << "% --- Poison Cylinders\n"; //Since cold size is 0.035, 0.05 should allow for realistic expansions and changes in user input.
    for (layr=1; layr<(nL-1);layr++){
        for (assm=0; assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
                for (plnk=1;plnk<7;plnk++){
                    for(colm=1;colm<(nBPC+1);colm++){
                        out_file << "surf BPC_"<< plnk << "_" << colm << "_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locBPC[plnk][colm][layr][assm][sect] << " " << locCP[plnk][layr][assm][sect] << " 0.05\n";}}}}}
    }
    out_file << "% --- Structural Planes\n";
    for (layr=0; layr<nL;layr++){
        for (assm=0; assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
                out_file <<
                "surf SP3_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSP[3][layr][assm][sect] << " %Top of Bot Y\n" <<
                "surf SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " plane " << sqrt3 << " 1 0 " << locSP[4][layr][assm][sect] << " %Right of Left Y\n" <<
                "surf SP5_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSP[5][layr][assm][sect] << " %Top of Bot CB\n" <<
                "surf SP6_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locSP[6][layr][assm][sect] << " %Right of Bot CB\n" <<
                "surf SP7_" << layr << "_" << assm+1 << "_" << sect+1 << " plane -1 " << sqrt3 << " 0 " << locSP[7][layr][assm][sect] << " %Top of Left CB \n" <<
                "surf SP8_" << layr << "_" << assm+1 << "_" << sect+1 << " plane " << sqrt3 << " 1 0 " << locSP[8][layr][assm][sect] << " %Right of Left CB\n" <<
                "surf SP9_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSP[9][layr][assm][sect] << " %Top of Bot CB Slot\n" <<
                "surf SP10_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locSP[10][layr][assm][sect] << " %Right of Bot CB Slot\n" <<
                "surf SP11_" << layr << "_" << assm+1 << "_" << sect+1 << " plane -1 " << sqrt3 << " 0 " << locSP[11][layr][assm][sect] << " %Top of Left CB Slot\n" <<
                "surf SP12_" << layr << "_" << assm+1 << "_" << sect+1 << " plane " << sqrt3 << " 1 0 " << locSP[12][layr][assm][sect] << " %Right of Left CB Slot\n" <<
                "surf SP13_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSP[13][layr][assm][sect] << " %Top Wrapper Out\n" <<
                "surf SP14_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSP[14][layr][assm][sect] << " %Top Wrapper In\n" <<
                "surf SP15_" << layr << "_" << assm+1 << "_" << sect+1 << " plane " << sqrt3 << " 1 0 " << locSP[15][layr][assm][sect] << " %Right Wrapper Out\n" <<
                "surf SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " plane " << sqrt3 << " 1 0 " << locSP[16][layr][assm][sect] << " %Right Wrapper In\n" <<
                //Likely need some logic here for only writing these if needed (indent/nonindent - if given option)
                "surf iL_" << layr << "_" << assm+1 << "_" << sect+1 << " plane " << sqrt3 << " 1 0 " << lociL[layr][assm][sect] << "\n" <<
                "surf iR_" << layr << "_" << assm+1 << "_" << sect+1 << " plane " << sqrt3 << " 1 0 " << lociR[layr][assm][sect] << "\n\n";
                }}}
    out_file << "% --- Plank Planes\n";
    for (layr=0; layr<nL;layr++){
        for (assm=0; assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
                out_file <<
                "surf PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[1][layr][assm][sect] << "\n" <<
                "surf PP2_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[2][layr][assm][sect] << "\n" <<
                "surf PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[3][layr][assm][sect] << "\n" <<
                "surf PP4_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[4][layr][assm][sect] << "\n" <<
                "surf PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[5][layr][assm][sect] << "\n" <<
                "surf PP6_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[6][layr][assm][sect] << "\n" <<
                "surf PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[7][layr][assm][sect] << "\n" <<
                "surf PP8_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[8][layr][assm][sect] << "\n" <<
                "surf PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[9][layr][assm][sect] << "\n" <<
                "surf PP10_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[10][layr][assm][sect] << "\n" <<
                "surf PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[11][layr][assm][sect] << "\n" <<
                "surf PP12_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locPP[12][layr][assm][sect] << "\n\n";
                }}}
    out_file << "% --- Sleeve Planes\n";
    for (layr=0; layr<nL;layr++){
        for (assm=0; assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
                out_file <<
                "surf SL1_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[1][layr][assm][sect] << "\n" <<
                "surf SL2_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[2][layr][assm][sect] << "\n" <<
                "surf SL3_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[3][layr][assm][sect] << "\n" <<
                "surf SL4_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[4][layr][assm][sect] << "\n" <<
                "surf SL5_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[5][layr][assm][sect] << "\n" <<
                "surf SL6_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[6][layr][assm][sect] << "\n" <<
                "surf SL7_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[7][layr][assm][sect] << "\n" <<
                "surf SL8_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[8][layr][assm][sect] << "\n" <<
                "surf SL9_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[9][layr][assm][sect] << "\n" <<
                "surf SL10_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[10][layr][assm][sect] << "\n" <<
                "surf SL11_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[11][layr][assm][sect] << "\n" <<
                "surf SL12_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locSL[12][layr][assm][sect] << "\n\n";
                }}}
    out_file << "% --- Fuel Planes\n";
    for (layr=1; layr<(nL-1);layr++){
        for (assm=0; assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
            out_file <<
            "surf F1y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[1][layr][assm][sect] << "\n" <<
            "surf F1xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[1][layr][assm][sect] << "\n" <<
            "surf F1xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[1][layr][assm][sect] << "\n" <<
            "surf F2y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[2][layr][assm][sect] << "\n" <<
            "surf F2xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[2][layr][assm][sect] << "\n" <<
            "surf F2xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[2][layr][assm][sect] << "\n" <<
            "surf F3y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[3][layr][assm][sect] << "\n" <<
            "surf F3xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[3][layr][assm][sect] << "\n" <<
            "surf F3xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[3][layr][assm][sect] << "\n" <<
            "surf F4y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[4][layr][assm][sect] << "\n" <<
            "surf F4xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[4][layr][assm][sect] << "\n" <<
            "surf F4xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[4][layr][assm][sect] << "\n" <<
            "surf F5y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[5][layr][assm][sect] << "\n" <<
            "surf F5xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[5][layr][assm][sect] << "\n" <<
            "surf F5xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[5][layr][assm][sect] << "\n" <<
            "surf F6y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[6][layr][assm][sect] << "\n" <<
            "surf F6xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[6][layr][assm][sect] << "\n" <<
            "surf F6xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[6][layr][assm][sect] << "\n" <<
            "surf F7y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[7][layr][assm][sect] << "\n" <<
            "surf F7xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[7][layr][assm][sect] << "\n" <<
            "surf F7xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[7][layr][assm][sect] << "\n" <<
            "surf F8y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[8][layr][assm][sect] << "\n" <<
            "surf F8xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[8][layr][assm][sect] << "\n" <<
            "surf F8xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[8][layr][assm][sect] << "\n" <<
            "surf F9y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[9][layr][assm][sect] << "\n" <<
            "surf F9xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[9][layr][assm][sect] << "\n" <<
            "surf F9xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[9][layr][assm][sect] << "\n" <<
            "surf F10y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[10][layr][assm][sect] << "\n" <<
            "surf F10xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[10][layr][assm][sect] << "\n" <<
            "surf F10xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[10][layr][assm][sect] << "\n" <<
            "surf F11y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[11][layr][assm][sect] << "\n" <<
            "surf F11xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[11][layr][assm][sect] << "\n" <<
            "surf F11xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[11][layr][assm][sect] << "\n" <<
            "surf F12y_" << layr << "_" << assm+1 << "_" << sect+1 << " py " << locFy[12][layr][assm][sect] << "\n" <<
            "surf F12xL_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxL[12][layr][assm][sect] << "\n" <<
            "surf F12xR_" << layr << "_" << assm+1 << "_" << sect+1 << " px " << locFxR[12][layr][assm][sect] << "\n\n";
            }}}
    out_file << "% --- Spacers\n";
    for (layr=0; layr<nL;layr++){
        for (assm=0; assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
                out_file <<
                "surf S1L_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpL[1][layr][assm][sect] << " " << locPP[1][layr][assm][sect] << " " << 0.5*dimSpacer[layr][assm][sect] << "\n" <<
                "surf S1R_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpR[1][layr][assm][sect] << " " << locPP[1][layr][assm][sect] << " " << 0.5*dimSpacer[layr][assm][sect] << "\n" <<
                "surf S2L_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpL[2][layr][assm][sect] << " " << locPP[3][layr][assm][sect] << " " << dimSpacer[layr][assm][sect] << "\n" <<
                "surf S2R_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpR[2][layr][assm][sect] << " " << locPP[3][layr][assm][sect] << " " << dimSpacer[layr][assm][sect] << "\n" <<
                "surf S3L_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpL[3][layr][assm][sect] << " " << locPP[5][layr][assm][sect] << " " << dimSpacer[layr][assm][sect] << "\n" <<
                "surf S3R_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpR[3][layr][assm][sect] << " " << locPP[5][layr][assm][sect] << " " << dimSpacer[layr][assm][sect] << "\n" <<
                "surf S4L_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpL[4][layr][assm][sect] << " " << locPP[7][layr][assm][sect] << " " << dimSpacer[layr][assm][sect] << "\n" <<
                "surf S4R_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpR[4][layr][assm][sect] << " " << locPP[7][layr][assm][sect] << " " << dimSpacer[layr][assm][sect] << "\n" <<
                "surf S5L_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpL[5][layr][assm][sect] << " " << locPP[9][layr][assm][sect] << " " << dimSpacer[layr][assm][sect] << "\n" <<
                "surf S5R_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpR[5][layr][assm][sect] << " " << locPP[9][layr][assm][sect] << " " << dimSpacer[layr][assm][sect] << "\n" <<
                "surf S6L_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpL[6][layr][assm][sect] << " " << locPP[11][layr][assm][sect] << " " << dimSpacer[layr][assm][sect] << "\n" <<
                "surf S6R_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpR[6][layr][assm][sect] << " " << locPP[11][layr][assm][sect] << " " << dimSpacer[layr][assm][sect] << "\n" <<
                "surf S7L_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpL[7][layr][assm][sect] << " " << locPP[12][layr][assm][sect] << " " << 0.5*dimSpacer[layr][assm][sect] << "\n" <<
                "surf S7R_" << layr << "_" << assm+1 << "_" << sect+1 << " cyl " << locSpR[7][layr][assm][sect] << " " << locPP[12][layr][assm][sect] << " " << 0.5*dimSpacer[layr][assm][sect] << "\n\n";
                }}}
    "% --- Create Cells for Each Axial Slice\n"; //Need to update the plank displacements (x,y)
    //int cbCounter = 0;
    for (layr=0; layr<nL;layr++){
        for (assm=0; assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
                if (controlBladeIn[layr][assm]){
                    //cbCounter++;
                    //cout << cbCounter << " ";
                    //cout << layr << " " << assm << " " << sect <<  " " << controlBladeIn[layr][assm];
                    if (singleCB){
                        out_file <<"cell CB_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " cBlade div1 div3 -SP6_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP7_" << layr << "_" << assm+1 << "_" << sect+1 << " (-SP5_" << layr << "_" << assm+1 << "_" << sect+1 << ":-SP8_" << layr << "_" << assm+1 << "_" << sect+1 << ")\n";}
                        else{out_file <<"cell CB_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " cBlade_" << layr << "_" << assm+1 << "_" << sect+1 << " div1 div3 -SP6_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP7_" << layr << "_" << assm+1 << "_" << sect+1 << " (-SP5_" << layr << "_" << assm+1 << "_" << sect+1 << ":-SP8_" << layr << "_" << assm+1 << "_" << sect+1 << ")\n";}}
                    else{
                        if (singleFlibe){out_file <<"cell CB_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe div1 div3 -SP6_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP7_" << layr << "_" << assm+1 << "_" << sect+1 << " (-SP5_" << layr << "_" << assm+1 << "_" << sect+1 << ":-SP8_" << layr << "_" << assm+1 << "_" << sect+1 << ")\n";}
                            else {out_file <<"cell CB_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe_" << layr << "_" << assm+1 << "_" << sect+1 << " div1 div3 -SP6_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP7_" << layr << "_" << assm+1 << "_" << sect+1 << " (-SP5_" << layr << "_" << assm+1 << "_" << sect+1 << ":-SP8_" << layr << "_" << assm+1 << "_" << sect+1 << ")\n";}}
                if (singleFlibe){
                    out_file << "cell YCool_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe div1 div3 -SP10_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP11_" << layr << "_" << assm+1 << "_" << sect+1 << " (-SP9_" << layr << "_" << assm+1 << "_" << sect+1 << ":-SP12_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(div1 div3 -SP6_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP7_" << layr << "_" << assm+1 << "_" << sect+1 << " (-SP5_" << layr << "_" << assm+1 << "_" << sect+1 << ":-SP8_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                    "cell OutFlibe_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe (SP13_" << layr << "_" << assm+1 << "_" << sect+1 << ":SP15_" << layr << "_" << assm+1 << "_" << sect+1 << ")\n" <<
                    "cell Channel1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe SP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S1L_" << layr << "_" << assm+1 << "_" << sect+1 << " S1R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Channel2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe PP2_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S2L_" << layr << "_" << assm+1 << "_" << sect+1 << " S2R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Channel3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe PP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S3L_" << layr << "_" << assm+1 << "_" << sect+1 << " S3R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Channel4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe PP6_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S4L_" << layr << "_" << assm+1 << "_" << sect+1 << " S4R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Channel5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe PP8_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S5L_" << layr << "_" << assm+1 << "_" << sect+1 << " S5R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Channel6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe PP10_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S6L_" << layr << "_" << assm+1 << "_" << sect+1 << " S6R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Channel7_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe PP12_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP14_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S7L_" << layr << "_" << assm+1 << "_" << sect+1 << " S7R_" << layr << "_" << assm+1 << "_" << sect+1 << endl;}
                    else {
                        out_file << "cell YCool_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe_" << layr << "_" << assm+1 << "_" << sect+1 << " div1 div3 -SP10_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP11_" << layr << "_" << assm+1 << "_" << sect+1 << " (-SP9_" << layr << "_" << assm+1 << "_" << sect+1 << ":-SP12_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(div1 div3 -SP6_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP7_" << layr << "_" << assm+1 << "_" << sect+1 << " (-SP5_" << layr << "_" << assm+1 << "_" << sect+1 << ":-SP8_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                        "cell OutFlibe_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe_" << layr << "_" << assm+1 << "_" << sect+1 << " (SP13_" << layr << "_" << assm+1 << "_" << sect+1 << ":SP15_" << layr << "_" << assm+1 << "_" << sect+1 << ")\n" <<
                        "cell Channel1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe_" << layr << "_" << assm+1 << "_" << sect+1 << " SP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S1L_" << layr << "_" << assm+1 << "_" << sect+1 << " S1R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Channel2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe_" << layr << "_" << assm+1 << "_" << sect+1 << " PP2_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S2L_" << layr << "_" << assm+1 << "_" << sect+1 << " S2R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Channel3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe_" << layr << "_" << assm+1 << "_" << sect+1 << " PP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S3L_" << layr << "_" << assm+1 << "_" << sect+1 << " S3R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Channel4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe_" << layr << "_" << assm+1 << "_" << sect+1 << " PP6_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S4L_" << layr << "_" << assm+1 << "_" << sect+1 << " S4R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Channel5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe_" << layr << "_" << assm+1 << "_" << sect+1 << " PP8_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S5L_" << layr << "_" << assm+1 << "_" << sect+1 << " S5R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Channel6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe_" << layr << "_" << assm+1 << "_" << sect+1 << " PP10_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S6L_" << layr << "_" << assm+1 << "_" << sect+1 << " S6R_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Channel7_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " flibe_" << layr << "_" << assm+1 << "_" << sect+1 << " PP12_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP14_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " S7L_" << layr << "_" << assm+1 << "_" << sect+1 << " S7R_" << layr << "_" << assm+1 << "_" << sect+1 << endl;}
                if (singleGraph){
                    out_file <<
                    "cell iL1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP2_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell iL2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP4_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell iL3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP6_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell iL4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP8_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell iL5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP10_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell iL6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP12_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell iR1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP2_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell iR2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP4_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell iR3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP6_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell iR4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP8_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell iR5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP10_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell iR6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP12_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Structural_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " structural div1 div3 -SP13_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP15_" << layr << "_" << assm+1 << "_" << sect+1 << " -(div1 div3 -SP10_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP11_" << layr << "_" << assm+1 << "_" << sect+1 << " (-SP9_" << layr << "_" << assm+1 << "_" << sect+1 << ":-SP12_" << layr << "_" << assm+1 << "_" << sect+1 << ")) -(SP3_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP14_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP2_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP4_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP6_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP8_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP10_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP12_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP2_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP4_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP6_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP8_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP10_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP12_" << layr << "_" << assm+1 << "_" << sect+1 << ")\n" <<
                    "cell Sleeve1B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL1_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Sleeve1T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP2_" << layr << "_" << assm+1 << "_" << sect+1 << " SL2_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Sleeve2B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL3_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Sleeve2T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP4_" << layr << "_" << assm+1 << "_" << sect+1 << " SL4_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Sleeve3B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL5_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Sleeve3T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP6_" << layr << "_" << assm+1 << "_" << sect+1 << " SL6_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Sleeve4B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL7_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Sleeve4T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP8_" << layr << "_" << assm+1 << "_" << sect+1 << " SL8_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Sleeve5B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL9_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Sleeve5T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP10_" << layr << "_" << assm+1 << "_" << sect+1 << " SL10_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Sleeve6B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL11_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Sleeve6T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP12_" << layr << "_" << assm+1 << "_" << sect+1 << " SL12_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell S1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer (-S1L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S1R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " SP3_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell S2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer (-S2L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S2R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " PP2_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell S3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer (-S3L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S3R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " PP4_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell S4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer (-S4L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S4R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " PP6_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell S5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer (-S5L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S5R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " PP8_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell S6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer (-S6L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S6R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " PP10_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell S7_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer (-S7L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S7R_" << layr << "_" << assm+1 << "_" << sect+1 << ") PP12_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP14_" << layr << "_" << assm+1 << "_" << sect+1 << endl ;}
                    else {
                        out_file <<
                        "cell iL1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP2_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell iL2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP4_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell iL3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP6_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell iL4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP8_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell iL5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP10_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell iL6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP12_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell iR1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP2_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell iR2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP4_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell iR3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP6_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell iR4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP8_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell iR5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP10_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell iR6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " -iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP12_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Structural_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " structural_" << layr << "_" << assm+1 << "_" << sect+1 << " div1 div3 -SP13_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP15_" << layr << "_" << assm+1 << "_" << sect+1 << " -(div1 div3 -SP10_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP11_" << layr << "_" << assm+1 << "_" << sect+1 << " (-SP9_" << layr << "_" << assm+1 << "_" << sect+1 << ":-SP12_" << layr << "_" << assm+1 << "_" << sect+1 << ")) -(SP3_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP14_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP2_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP4_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP6_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP8_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP10_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " iL_" << layr << "_" << assm+1 << "_" << sect+1 << " PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP12_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP2_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP4_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP6_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP8_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP10_" << layr << "_" << assm+1 << "_" << sect+1 << ") -(-iR_" << layr << "_" << assm+1 << "_" << sect+1 << " SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP12_" << layr << "_" << assm+1 << "_" << sect+1 << ")\n" <<
                        "cell Sleeve1B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL1_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Sleeve1T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP2_" << layr << "_" << assm+1 << "_" << sect+1 << " SL2_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Sleeve2B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL3_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Sleeve2T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP4_" << layr << "_" << assm+1 << "_" << sect+1 << " SL4_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Sleeve3B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL5_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Sleeve3T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP6_" << layr << "_" << assm+1 << "_" << sect+1 << " SL6_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Sleeve4B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL7_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Sleeve4T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP8_" << layr << "_" << assm+1 << "_" << sect+1 << " SL8_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Sleeve5B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL9_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Sleeve5T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP10_" << layr << "_" << assm+1 << "_" << sect+1 << " SL10_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Sleeve6B_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL11_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell Sleeve6T_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " sleeve_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " -PP12_" << layr << "_" << assm+1 << "_" << sect+1 << " SL12_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell S1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer_" << layr << "_" << assm+1 << "_" << sect+1 << " (-S1L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S1R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP1_" << layr << "_" << assm+1 << "_" << sect+1 << " SP3_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell S2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer_" << layr << "_" << assm+1 << "_" << sect+1 << " (-S2L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S2R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP3_" << layr << "_" << assm+1 << "_" << sect+1 << " PP2_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell S3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer_" << layr << "_" << assm+1 << "_" << sect+1 << " (-S3L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S3R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP5_" << layr << "_" << assm+1 << "_" << sect+1 << " PP4_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell S4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer_" << layr << "_" << assm+1 << "_" << sect+1 << " (-S4L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S4R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP7_" << layr << "_" << assm+1 << "_" << sect+1 << " PP6_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell S5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer_" << layr << "_" << assm+1 << "_" << sect+1 << " (-S5L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S5R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP9_" << layr << "_" << assm+1 << "_" << sect+1 << " PP8_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell S6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer_" << layr << "_" << assm+1 << "_" << sect+1 << " (-S6L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S6R_" << layr << "_" << assm+1 << "_" << sect+1 << ") -PP11_" << layr << "_" << assm+1 << "_" << sect+1 << " PP10_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                        "cell S7_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " spacer_" << layr << "_" << assm+1 << "_" << sect+1 << " (-S7L_" << layr << "_" << assm+1 << "_" << sect+1 << ":-S7R_" << layr << "_" << assm+1 << "_" << sect+1 << ") PP12_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP14_" << layr << "_" << assm+1 << "_" << sect+1 << endl ;}
                if ( (layr>0) && (layr<(nL-1)) ){
                    out_file <<
                    "cell Fuel1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe1_" << layr << "_" << assm+1 << "_" << sect+1 << " SL1_" << layr << "_" << assm+1 << "_" << sect+1 << " -F1y_" << layr << "_" << assm+1 << "_" << sect+1 << " F1xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F1xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Fuel2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe2_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL2_" << layr << "_" << assm+1 << "_" << sect+1 << " F2y_" << layr << "_" << assm+1 << "_" << sect+1 << " F2xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F2xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Fuel3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe3_" << layr << "_" << assm+1 << "_" << sect+1 << " SL3_" << layr << "_" << assm+1 << "_" << sect+1 << " -F3y_" << layr << "_" << assm+1 << "_" << sect+1 << " F3xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F3xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Fuel4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL4_" << layr << "_" << assm+1 << "_" << sect+1 << " F4y_" << layr << "_" << assm+1 << "_" << sect+1 << " F4xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F4xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Fuel5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe5_" << layr << "_" << assm+1 << "_" << sect+1 << " SL5_" << layr << "_" << assm+1 << "_" << sect+1 << " -F5y_" << layr << "_" << assm+1 << "_" << sect+1 << " F5xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F5xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Fuel6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe6_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL6_" << layr << "_" << assm+1 << "_" << sect+1 << " F6y_" << layr << "_" << assm+1 << "_" << sect+1 << " F6xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F6xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Fuel7_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe7_" << layr << "_" << assm+1 << "_" << sect+1 << " SL7_" << layr << "_" << assm+1 << "_" << sect+1 << " -F7y_" << layr << "_" << assm+1 << "_" << sect+1 << " F7xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F7xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Fuel8_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe8_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL8_" << layr << "_" << assm+1 << "_" << sect+1 << " F8y_" << layr << "_" << assm+1 << "_" << sect+1 << " F8xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F8xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Fuel9_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe9_" << layr << "_" << assm+1 << "_" << sect+1 << " SL9_" << layr << "_" << assm+1 << "_" << sect+1 << " -F9y_" << layr << "_" << assm+1 << "_" << sect+1 << " F9xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F9xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Fuel10_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe10_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL10_" << layr << "_" << assm+1 << "_" << sect+1 << " F10y_" << layr << "_" << assm+1 << "_" << sect+1 << " F10xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F10xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Fuel11_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe11_" << layr << "_" << assm+1 << "_" << sect+1 << " SL11_" << layr << "_" << assm+1 << "_" << sect+1 << " -F11y_" << layr << "_" << assm+1 << "_" << sect+1 << " F11xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F11xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                    "cell Fuel12_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill stripe12_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL12_" << layr << "_" << assm+1 << "_" << sect+1 << " F12y_" << layr << "_" << assm+1 << "_" << sect+1 << " F12xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F12xR_" << layr << "_" << assm+1 << "_" << sect+1 << endl;
                    if (BPUsage == true){
                        for (plnk=1;plnk<7;plnk++){
                        out_file <<
                        "cell BPA" << plnk << "_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " fill BP" << plnk << "_" << layr << "_" << assm+1 << "_" << sect+1 << " (";
                        for (colm=1;colm<(nBPC+1);colm++){
                            out_file << "-BPC_" << plnk << "_" << colm << "_" << layr << "_" << assm+1 << "_" << sect+1;
                            if (colm < nBPC){out_file << ":";}}
                        out_file << ")\n";}
                        if (singleGraph){
                            for (plnk=1;plnk<7;plnk++){
                                out_file << "cell Plank" << plnk << "_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL" << (2*plnk-1) << "_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL" << (2*plnk) << "_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL" << (2*plnk-1) << "_" << layr << "_" << assm+1 << "_" << sect+1 << " -F" << (2*plnk-1) << "y_" << layr << "_" << assm+1 << "_" << sect+1 << " F" << (2*plnk-1) << "xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F" << (2*plnk-1) << "xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL" << (2*plnk) << "_" << layr << "_" << assm+1 << "_" << sect+1 << " F" << (2*plnk) << "y_" << layr << "_" << assm+1 << "_" << sect+1 << " F" << (2*plnk) << "xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F" << (2*plnk) << "xR_" << layr << "_" << assm+1 << "_" << sect+1 << ")) -(";
                                for (colm=1;colm<(nBPC+1);colm++){
                                    out_file << "-BPC_" << plnk << "_" << colm << "_" << layr << "_" << assm+1 << "_" << sect+1;
                                    if (colm < nBPC){out_file << ":";}}
                                out_file << ")\n";}}
                            else {
                                for (plnk=1;plnk<7;plnk++){
                                    out_file <<
                                    "cell Plank" << plnk << "_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL" << (2*plnk-1) << "_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL" << (2*plnk) << "_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL" << (2*plnk-1) << "_" << layr << "_" << assm+1 << "_" << sect+1 << " -F" << (2*plnk-1) << "y_" << layr << "_" << assm+1 << "_" << sect+1 << " F" << (2*plnk-1) << "xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F" << (2*plnk-1) << "xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL" << (2*plnk) << "_" << layr << "_" << assm+1 << "_" << sect+1 << " F" << (2*plnk) << "y_" << layr << "_" << assm+1 << "_" << sect+1 << " F" << (2*plnk) << "xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F" << (2*plnk) << "xR_" << layr << "_" << assm+1 << "_" << sect+1 << ")) -(";
                                    for (colm=1;colm<(nBPC+1);colm++){
                                        out_file << "-BPC_" << plnk << "_" << colm << "_" << layr << "_" << assm+1 << "_" << sect+1;
                                        if (colm < nBPC){out_file << ":";}}
                                    out_file << ")\n";}}}
                        else{
                            if (singleGraph){
                                out_file <<
                                "cell Plank1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL1_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL2_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL1_" << layr << "_" << assm+1 << "_" << sect+1 << " -F1y_" << layr << "_" << assm+1 << "_" << sect+1 << " F1xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F1xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL2_" << layr << "_" << assm+1 << "_" << sect+1 << " F2y_" << layr << "_" << assm+1 << "_" << sect+1 << " F2xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F2xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                                "cell Plank2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL3_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL4_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL3_" << layr << "_" << assm+1 << "_" << sect+1 << " -F3y_" << layr << "_" << assm+1 << "_" << sect+1 << " F3xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F3xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL4_" << layr << "_" << assm+1 << "_" << sect+1 << " F4y_" << layr << "_" << assm+1 << "_" << sect+1 << " F4xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F4xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                                "cell Plank3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL5_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL6_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL5_" << layr << "_" << assm+1 << "_" << sect+1 << " -F5y_" << layr << "_" << assm+1 << "_" << sect+1 << " F5xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F5xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL6_" << layr << "_" << assm+1 << "_" << sect+1 << " F6y_" << layr << "_" << assm+1 << "_" << sect+1 << " F6xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F6xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                                "cell Plank4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL7_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL8_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL7_" << layr << "_" << assm+1 << "_" << sect+1 << " -F7y_" << layr << "_" << assm+1 << "_" << sect+1 << " F7xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F7xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL8_" << layr << "_" << assm+1 << "_" << sect+1 << " F8y_" << layr << "_" << assm+1 << "_" << sect+1 << " F8xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F8xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                                "cell Plank5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL9_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL10_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL9_" << layr << "_" << assm+1 << "_" << sect+1 << " -F9y_" << layr << "_" << assm+1 << "_" << sect+1 << " F9xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F9xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL10_" << layr << "_" << assm+1 << "_" << sect+1 << " F10y_" << layr << "_" << assm+1 << "_" << sect+1 << " F10xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F10xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                                "cell Plank6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL11_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL12_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL11_" << layr << "_" << assm+1 << "_" << sect+1 << " -F11y_" << layr << "_" << assm+1 << "_" << sect+1 << " F11xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F11xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL12_" << layr << "_" << assm+1 << "_" << sect+1 << " F12y_" << layr << "_" << assm+1 << "_" << sect+1 << " F12xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F12xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n";}
                                else {
                                    out_file <<
                                    "cell Plank1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL1_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL2_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL1_" << layr << "_" << assm+1 << "_" << sect+1 << " -F1y_" << layr << "_" << assm+1 << "_" << sect+1 << " F1xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F1xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL2_" << layr << "_" << assm+1 << "_" << sect+1 << " F2y_" << layr << "_" << assm+1 << "_" << sect+1 << " F2xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F2xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                                    "cell Plank2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL3_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL4_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL3_" << layr << "_" << assm+1 << "_" << sect+1 << " -F3y_" << layr << "_" << assm+1 << "_" << sect+1 << " F3xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F3xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL4_" << layr << "_" << assm+1 << "_" << sect+1 << " F4y_" << layr << "_" << assm+1 << "_" << sect+1 << " F4xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F4xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                                    "cell Plank3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL5_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL6_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL5_" << layr << "_" << assm+1 << "_" << sect+1 << " -F5y_" << layr << "_" << assm+1 << "_" << sect+1 << " F5xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F5xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL6_" << layr << "_" << assm+1 << "_" << sect+1 << " F6y_" << layr << "_" << assm+1 << "_" << sect+1 << " F6xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F6xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                                    "cell Plank4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL7_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL8_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL7_" << layr << "_" << assm+1 << "_" << sect+1 << " -F7y_" << layr << "_" << assm+1 << "_" << sect+1 << " F7xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F7xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL8_" << layr << "_" << assm+1 << "_" << sect+1 << " F8y_" << layr << "_" << assm+1 << "_" << sect+1 << " F8xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F8xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                                    "cell Plank5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL9_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL10_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL9_" << layr << "_" << assm+1 << "_" << sect+1 << " -F9y_" << layr << "_" << assm+1 << "_" << sect+1 << " F9xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F9xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL10_" << layr << "_" << assm+1 << "_" << sect+1 << " F10y_" << layr << "_" << assm+1 << "_" << sect+1 << " F10xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F10xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n" <<
                                    "cell Plank6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL11_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL12_" << layr << "_" << assm+1 << "_" << sect+1 << " -((SL11_" << layr << "_" << assm+1 << "_" << sect+1 << " -F11y_" << layr << "_" << assm+1 << "_" << sect+1 << " F11xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F11xR_" << layr << "_" << assm+1 << "_" << sect+1 << "):(-SL12_" << layr << "_" << assm+1 << "_" << sect+1 << " F12y_" << layr << "_" << assm+1 << "_" << sect+1 << " F12xL_" << layr << "_" << assm+1 << "_" << sect+1 << " -F12xR_" << layr << "_" << assm+1 << "_" << sect+1 << "))\n";}}}
                    else{
                        if (singleGraph) {
                            out_file <<
                            "cell Plank1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL1_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL2_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                            "cell Plank2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL3_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL4_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                            "cell Plank3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL5_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL6_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                            "cell Plank4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL7_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL8_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                            "cell Plank5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL9_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL10_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                            "cell Plank6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL11_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL12_" << layr << "_" << assm+1 << "_" << sect+1 << endl;}
                            else{
                                out_file <<
                                "cell Plank1_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL1_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL2_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                                "cell Plank2_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL3_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL4_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                                "cell Plank3_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL5_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL6_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                                "cell Plank4_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL7_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL8_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                                "cell Plank5_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL9_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL10_" << layr << "_" << assm+1 << "_" << sect+1 << endl <<
                                "cell Plank6_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " meat_" << layr << "_" << assm+1 << "_" << sect+1 << " SP4_" << layr << "_" << assm+1 << "_" << sect+1 << " -SP16_" << layr << "_" << assm+1 << "_" << sect+1 << " SL11_" << layr << "_" << assm+1 << "_" << sect+1 << " -SL12_" << layr << "_" << assm+1 << "_" << sect+1 << endl;}}
               out_file << "cell below_" << layr << "_" << assm+1 << "_" << sect+1 << " third_" << layr << "_" << assm+1 << "_" << sect+1 << " void (-div1:-div3)\n\n";}}}
    out_file << "% --- Make Assembly Sections\n";
    sect = 0;
    for (layr=0;layr<nL;layr++){
        for (assm=0;assm<nA;assm++){
            out_file <<
            "cell sect1_" << layr << "_" << assm+1 << " third1_" << layr << "_" << assm+1 << " fill third_" << layr << "_" << assm+1 << "_" << sect+1 << " -bigS\n" <<
            "cell sect2_" << layr << "_" << assm+1 << " third2_" << layr << "_" << assm+1 << " fill third_" << layr << "_" << assm+1 << "_" << sect+2 << " -bigS\n" <<
            "cell sect3_" << layr << "_" << assm+1 << " third3_" << layr << "_" << assm+1 << " fill third_" << layr << "_" << assm+1 << "_" << sect+3 << " -bigS\n\n";
        }}
    out_file << "% --- Section Rotations\n";
    for (layr=0;layr<nL;layr++){
        for (assm=0;assm<nA;assm++){
            out_file <<
            "utrans third1_" << layr << "_" << assm+1 << " 0 0 0 0 0 -60\n" <<
            "utrans third2_" << layr << "_" << assm+1 << " 0 0 0 0 0 180\n" <<
            "utrans third3_" << layr << "_" << assm+1 << " 0 0 0 0 0 60\n";}
        out_file << endl;}
    sect = 0;
    out_file << "% --- Make Assemblies\n";
    for (layr=0;layr<nL;layr++){
        for (assm=0;assm<nA;assm++){
                out_file <<
                "cell section1_" << layr << "_" << assm+1 << " assemS_" << layr << "_" << assm+1 << " fill third1_" << layr << "_" << assm+1 << " div1 div2 -bigS" << endl <<
                "cell section2_" << layr << "_" << assm+1 << " assemS_" << layr << "_" << assm+1 << " fill third2_" << layr << "_" << assm+1 << " -div1 -div3 -bigS" << endl <<
                "cell section3_" << layr << "_" << assm+1 << " assemS_" << layr << "_" << assm+1 << " fill third3_" << layr << "_" << assm+1 << " -div2 div3 -bigS\n\n"; }}
    out_file << "% --- Duplicate Assemblies Within Groups\n";
    for (layr=0;layr<nL;layr++){
        for (assm=0;assm<nA;assm++){
            for (sect=0;sect<nS;sect++){
                out_file <<
                "cell Assem_" << layr << "_" << assm+1 << "_" << sect+1 << "_1 A_" << layr << "_" << assm+1 << "_" << sect+1 << " fill assemS_" << layr << "_" << assm+1 << " -bigS\n";}
            out_file << endl;}}
    out_file << "% --- Rotate Group Assemblies\n";
    for (layr=0;layr<nL;layr++){
        for (assm=0;assm<nA;assm++){
            out_file <<
            "utrans A_" << layr << "_" << assm+1 << "_2 0 0 0 0 0 -120\n" <<
            "utrans A_" << layr << "_" << assm+1 << "_3 0 0 0 0 0 120\n\n";}}
    out_file << "% ---  Make Each Axial Core Layout\n"; //May want to include upper axial in with fuel had have it depend on 1/3 assembly outlet conditions
    "%Note: Y-direction definition is reversed since SERPENT builds lattices up, not down\n";
    for (layr=0;layr<nL;layr++){
        out_file <<
        "lat core_" << layr << " 3 0 0 25 25 " << dimAssemblyPitch[layr] << "\n" <<
        "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %25\n" <<
         "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %24\n" <<
          "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " 51_" << layr << " 51_" << layr << " 51_" << layr << " 51_" << layr << " 51_" << layr << " 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %23\n" <<
           "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " 51_" << layr << " A_" << layr << "_75_2 A_" << layr << "_79_2 A_" << layr << "_83_2 A_" << layr << "_82_2 A_" << layr << "_78_2 A_" << layr << "_74_2 51_" << layr << " 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %22\n" <<
            "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " 51_" << layr << " A_" << layr << "_58_2 A_" << layr << "_61_2 A_" << layr << "_65_2 A_" << layr << "_69_2 A_" << layr << "_72_2 A_" << layr << "_68_2 A_" << layr << "_64_2 A_" << layr << "_60_2 A_" << layr << "_57_2 51_" << layr << " 51_" << layr << " 50_" << layr << " 50_" << layr << " %21\n" <<
             "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_76_2 A_" << layr << "_62_2 A_" << layr << "_44_2 A_" << layr << "_47_2 A_" << layr << "_51_2 A_" << layr << "_55_2 A_" << layr << "_54_2 A_" << layr << "_50_2 A_" << layr << "_46_2 A_" << layr << "_43_2 A_" << layr << "_59_2 A_" << layr << "_73_2 51_" << layr << " 50_" << layr << " 50_" << layr << " %20\n" <<
              "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_80_2 A_" << layr << "_66_2 A_" << layr << "_48_2 A_" << layr << "_32_2 A_" << layr << "_37_2 A_" << layr << "_41_2 A_" << layr << "_34_2 A_" << layr << "_40_2 A_" << layr << "_36_2 A_" << layr << "_31_2 A_" << layr << "_45_2 A_" << layr << "_63_2 A_" << layr << "_77_2 51_" << layr << " 50_" << layr << " 50_" << layr << " %19\n" <<
               "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_84_2 A_" << layr << "_70_2 A_" << layr << "_52_2 A_" << layr << "_38_2 A_" << layr << "_22_2 A_" << layr << "_25_2 A_" << layr << "_29_2 A_" << layr << "_28_2 A_" << layr << "_24_2 A_" << layr << "_21_2 A_" << layr << "_35_2 A_" << layr << "_49_2 A_" << layr << "_67_2 A_" << layr << "_81_2 51_" << layr << " 50_" << layr << " 50_" << layr << " %18\n" <<
                "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_81_3 A_" << layr << "_71_3 A_" << layr << "_56_2 A_" << layr << "_42_2 A_" << layr << "_26_2 A_" << layr << "_14_2 A_" << layr << "_17_2 A_" << layr << "_20_2 A_" << layr << "_16_2 A_" << layr << "_13_2 A_" << layr << "_23_2 A_" << layr << "_39_2 A_" << layr << "_53_2 A_" << layr << "_71_2 A_" << layr << "_84_1 51_" << layr << " 50_" << layr << " 50_" << layr << " %17\n" <<
                 "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_77_3 A_" << layr << "_67_3 A_" << layr << "_53_3 A_" << layr << "_33_3 A_" << layr << "_30_2 A_" << layr << "_18_2 A_" << layr << "_8_2 A_" << layr << "_11_2 A_" << layr << "_10_2 A_" << layr << "_7_2 A_" << layr << "_15_2 A_" << layr << "_27_2 A_" << layr << "_33_2 A_" << layr << "_56_1 A_" << layr << "_70_1 A_" << layr << "_80_1 51_" << layr << " 50_" << layr << " 50_" << layr << " %16\n" <<
                  "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_73_3 A_" << layr << "_63_3 A_" << layr << "_49_3 A_" << layr << "_39_3 A_" << layr << "_27_3 A_" << layr << "_19_3 A_" << layr << "_12_2 A_" << layr << "_4_2 A_" << layr << "_5_2 A_" << layr << "_3_2 A_" << layr << "_9_2 A_" << layr << "_19_2 A_" << layr << "_30_1 A_" << layr << "_42_1 A_" << layr << "_52_1 A_" << layr << "_66_1 A_" << layr << "_76_1 51_" << layr << " 50_" << layr << " 50_" << layr << " %15\n" <<
                   "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_59_3 A_" << layr << "_45_3 A_" << layr << "_35_3 A_" << layr << "_23_3 A_" << layr << "_15_3 A_" << layr << "_9_3 A_" << layr << "_6_3 A_" << layr << "_2_2 A_" << layr << "_1_2 A_" << layr << "_6_2 A_" << layr << "_12_1 A_" << layr << "_18_1 A_" << layr << "_26_1 A_" << layr << "_38_1 A_" << layr << "_48_1 A_" << layr << "_62_1 51_" << layr << " 50_" << layr << " 50_" << layr << " %14\n" <<
                    "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_57_3 A_" << layr << "_43_3 A_" << layr << "_31_3 A_" << layr << "_21_3 A_" << layr << "_13_3 A_" << layr << "_7_3 A_" << layr << "_3_3 A_" << layr << "_1_3 51_" << layr << " A_" << layr << "_2_1 A_" << layr << "_4_1 A_" << layr << "_8_1 A_" << layr << "_14_1 A_" << layr << "_22_1 A_" << layr << "_32_1 A_" << layr << "_44_1 A_" << layr << "_58_1 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %13\n" <<
                     "50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_60_3 A_" << layr << "_46_3 A_" << layr << "_36_3 A_" << layr << "_24_3 A_" << layr << "_16_3 A_" << layr << "_10_3 A_" << layr << "_5_3 A_" << layr << "_2_3 A_" << layr << "_1_1 A_" << layr << "_5_1 A_" << layr << "_11_1 A_" << layr << "_17_1 A_" << layr << "_25_1 A_" << layr << "_37_1 A_" << layr << "_47_1 A_" << layr << "_61_1 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %12\n" <<
                      "50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_74_3 A_" << layr << "_64_3 A_" << layr << "_50_3 A_" << layr << "_40_3 A_" << layr << "_28_3 A_" << layr << "_20_3 A_" << layr << "_11_3 A_" << layr << "_4_3 A_" << layr << "_6_1 A_" << layr << "_3_1 A_" << layr << "_10_1 A_" << layr << "_20_1 A_" << layr << "_29_1 A_" << layr << "_41_1 A_" << layr << "_51_1 A_" << layr << "_65_1 A_" << layr << "_75_1 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %11\n" <<
                       "50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_78_3 A_" << layr << "_68_3 A_" << layr << "_54_3 A_" << layr << "_34_3 A_" << layr << "_29_3 A_" << layr << "_17_3 A_" << layr << "_8_3 A_" << layr << "_12_3 A_" << layr << "_9_1 A_" << layr << "_7_1 A_" << layr << "_16_1 A_" << layr << "_28_1 A_" << layr << "_34_1 A_" << layr << "_55_1 A_" << layr << "_69_1 A_" << layr << "_79_1 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %10\n" <<
                        "50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_82_3 A_" << layr << "_72_3 A_" << layr << "_55_3 A_" << layr << "_41_3 A_" << layr << "_25_3 A_" << layr << "_14_3 A_" << layr << "_18_3 A_" << layr << "_19_1 A_" << layr << "_15_1 A_" << layr << "_13_1 A_" << layr << "_24_1 A_" << layr << "_40_1 A_" << layr << "_54_1 A_" << layr << "_72_1 A_" << layr << "_83_1 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %9\n" <<
                         "50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_83_3 A_" << layr << "_69_3 A_" << layr << "_51_3 A_" << layr << "_37_3 A_" << layr << "_22_3 A_" << layr << "_26_3 A_" << layr << "_30_3 A_" << layr << "_27_1 A_" << layr << "_23_1 A_" << layr << "_21_1 A_" << layr << "_36_1 A_" << layr << "_50_1 A_" << layr << "_68_1 A_" << layr << "_82_1 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %8\n" <<
                          "50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_79_3 A_" << layr << "_65_3 A_" << layr << "_47_3 A_" << layr << "_32_3 A_" << layr << "_38_3 A_" << layr << "_42_3 A_" << layr << "_33_1 A_" << layr << "_39_1 A_" << layr << "_35_1 A_" << layr << "_31_1 A_" << layr << "_46_1 A_" << layr << "_64_1 A_" << layr << "_78_1 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << "  %7\n" <<
                           "50_" << layr << " 50_" << layr << " 51_" << layr << " A_" << layr << "_75_3 A_" << layr << "_61_3 A_" << layr << "_44_3 A_" << layr << "_48_3 A_" << layr << "_52_3 A_" << layr << "_56_3 A_" << layr << "_53_1 A_" << layr << "_49_1 A_" << layr << "_45_1 A_" << layr << "_43_1 A_" << layr << "_60_1 A_" << layr << "_74_1 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %6\n" <<
                            "50_" << layr << " 50_" << layr << " 51_" << layr << " 51_" << layr << " A_" << layr << "_58_3 A_" << layr << "_62_3 A_" << layr << "_66_3 A_" << layr << "_70_3 A_" << layr << "_71_1 A_" << layr << "_67_1 A_" << layr << "_63_1 A_" << layr << "_59_1 A_" << layr << "_57_1 51_" << layr << " 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %5\n" <<
                             "50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " 51_" << layr << " A_" << layr << "_76_3 A_" << layr << "_80_3 A_" << layr << "_84_3 A_" << layr << "_81_1 A_" << layr << "_77_1 A_" << layr << "_73_1 51_" << layr << " 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %4\n" <<
                              "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 51_" << layr << " 51_" << layr << " 51_" << layr << " 51_" << layr << " 51_" << layr << " 51_" << layr << " 51_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %3\n" <<
                               "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %2\n" <<
                                "50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " 50_" << layr << " %1\n\n";}
    out_file << "% --- Make Each Axial Slice\n";
    for (layr=0;layr<nL;layr++){
        out_file <<
        "cell S_" << layr << "_1 slice_" << layr << " fill core_" << layr << " -60\n" <<
        "cell S_" << layr << "_2 slice_" << layr << " B4C 60 -61\n" <<
        "cell S_" << layr << "_3 slice_" << layr << " barrel 61 -62\n" <<
        "cell S_" << layr << "_4 slice_" << layr << " inletF 62 -63\n" <<
        "cell S_" << layr << "_5 slice_" << layr << " alloyN 63 -64\n" <<
        "cell S_" << layr << "_6 slice_" << layr << " vessel 64\n\n";}
    double cumuHeight = dimAxSuppPlateHeight[0];
    //Need to adjust the heights of each slice to account for expansion. dimLayrHeight[layr]
    out_file << //If particles get lost between slices, may be due to rounding in code. Can make the sections be higher (e.g. 35cm) and just cut back to the appropriate height later.
    "% --- Make 3D Core by Stacking Slices\n" <<
    "lat stack 9 0 0 " << nL+2 << "\n" << //Change back to 20 after debug
    "0 botSupportPlate\n";
    for(layr=0;layr<nL;layr++){
        out_file << cumuHeight << " slice_" << layr << " %Axial Slice " << layr << "\n";
        cumuHeight += dimLayrHeight[layr];}    //out_file << "\tAxSlice_2\n";
    //dimActiveCoreHeight += 85;
    out_file << cumuHeight << "\ttopSupportPlate\n\n" <<
    "cell whole_1 0 fill stack -65\n" <<
    "cell whole_2 0 outside 65\n\n";
    out_file << "det tri dht 3 0 0 " << 0.5*(dimAssemblyPitch[0]+dimAssemblyPitch[nL-1]) << " 19 19 " << dimAxSuppPlateHeight[0]+dimLayrHeight[0] << " " << dimModelHeight-dimLayrHeight[nL-1]-dimAxSuppPlateHeight[1] << " " << nLt-2 << " dr -8 void\n";
    out_file << "det hex dh 3 0 0 " << 0.5*(dimAssemblyPitch[0]+dimAssemblyPitch[nL-1]) << " 19 19 " << dimAxSuppPlateHeight[0]+dimLayrHeight[0] << " " << dimModelHeight-dimLayrHeight[nL-1]-dimAxSuppPlateHeight[1] << " 1 dr -8 void\n";
    out_file << "det axialPower dz " << dimAxSuppPlateHeight[0]+dimLayrHeight[0] << " " << dimModelHeight-dimLayrHeight[nL-1]-dimAxSuppPlateHeight[1] << " 112 dr -8 void\n";
    out_file << "det hexSlice dh 3 0 0 " << 0.5*(dimAssemblyPitch[0]+dimAssemblyPitch[nL-1]) << " 19 19 " << dimAxSuppPlateHeight[0]+dimLayrHeight[0] << " " << dimModelHeight-dimLayrHeight[nL-1]-dimAxSuppPlateHeight[1] << " " << nLt-2 << " dr -8 void\n";
    out_file.close();
    cout << "\nThe file " << geometryOut << " has been written with the geometry for the Serpent simulation.\n";
}
