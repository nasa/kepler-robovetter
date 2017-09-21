/* The DR25 Kepler Robovetter

  The code is written in C++11.
  It has been tested to work with:
    - The g++ complier (Minimum version 4.7.2 tested - earlier versions unlikely to work.)
    - The clang++ compiler (Minimum version 3.4 tested. Version 3.3 may work, but untested. Earlier than 3.3 will not work.)
    - The intel icpc compiler (Latest version 17 tested. Earlier versions as far back as 11 very likely to work, but untested.)

  Compile with your favorite compiler, using one of the following options:
    - g++     -std=c++11 -O2 -o robovet Kepler-RoboVetter-DR25.cpp
    - clang++ -std=c++11 -O2 -o robovet Kepler-RoboVetter-DR25.cpp
    - icpc    -std=c++11 -O2 -o robovet Kepler-RoboVetter-DR25.cpp

  Run as "./robovet INPUTFILE  OUTFILE  NMC"

  INPUTFILE is the name of the input file for a particular run (OBS, INJ, INV, SCR1).
  OUTFILE is the name of the output file, as desired by the user.
  NMC is the number of Monte Carlo runs desired for computing the scores. A value of at least 100 is reccomemend for useable scores, and preferably at least 1,000. A value of 10,000 was used for calculating the scores archived at NExSCI.

  Example:
  ./robovet RoboVetter-Input-OBS.txt RoboVetter-Output-OBS.txt 10000


  Notices:
  Copyright © 2017 United States Government as represented by the Administrator of the National Aeronautics and Space Administration. All Rights Reserved.
  NASA acknowledges the SETI Institute’s primary role in authoring and producing the Kepler Robovetter under Cooperative Agreement Number NNX13AD01A

  Disclaimers:
  No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE. FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."
  Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT. IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW. RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.

*/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>

using namespace std;

// Initialze the random number generator. Using mt19937 should make this consistent across all platforms and compilers.
mt19937 randgen(1527);  // Initialize seed. Change the number if you want a different set of random numbers, and therefore slightly different scores. The seed 1527 is the number of days between launch and the 2nd reaction wheel failure - seemed like a fitting seed number as any, but it can be any integer as far as the code is concerned.
uniform_real_distribution<double> dist;  // Generates a random number between 0 and 1.

// Declare Various Functions
int READDATA(string);  // Function to read in the input data file. Returns the number of TCEs.

void TRANSITLIKE();     // Function to deterimine if the TCE is transit-like or not
void STELLARECLIPSE();  // Function to determine if the TCE has a secondary eclipse or other features indicating an EB
void SEASONAL();        // Function to look for seasonal depth differences
void APPCENTROID();     // Function to apply the centroid robovetter results from the input file
void APPEPHEMMATCH();   // Function to apply the ephemeris match results from the input file
int  ISSEC();           // Function to determine if the TCE is itself a secondary eclipse. Returns a 1 if it is designated a secondary eclipse, else returns 0
void PERIODALIAS();     // Function to determine if the Modelshift test data incidates the TCE is detected at a period ratio of the true signal
void ADDNOISE();        // Function to add noise to each Robovetter metric's intial value when in the Monte Carlo loop to compute scores
void SETERRORS();       // Function to determine the 'error' on each metric, i.e., the amount of noise to add for the Monte Carlo loop to compute scores

vector<double> COMPPT(double,double,double,double);  // Function to compute several metrics for comparing two sets of periods and epochs

double INVERFC(double);               // The complementary inverse error function, used in several functions when comparing the closeness of two values
double PSIG(double,double);           // A function used in COMPPT to compare two period values
double ESIG(double,double,double);    // A function used in COMPPT to compare two epoch values
double AGRAND(double,double,double);  // A function to add random asymmetrical Gaussian noise to a value given positive and negative error values


// Declare constants used in period/epoch comparisons - these can be adjusted if desired, but likely do not need to be.
const double PSIG_THRESH = 3.5;    // A threshold used to determine how close two periods have to be to be considered 'matching'. A higher value means two periods have to be closer to each other to 'match'.
const double ESIG_THRESH = 2.0;    // A threshold used to determine how close two epochs  have to be to be considered 'matching'. A higher value means two epochs  have to be closer to each other to 'match'.
const double WIDTHFAC = 2.5;       // A threshold used to determine how close, in units of TCE duration, two events have to be to be considered 'overlapping'. A lower value means two events have to be closer to each other to 'overlap'.
const double MISSIONDUR = 1600.0;  // Approximate mission duration in days. Used when comparing the closeness of two TCE's ephemerides over the length of the mission.


// Declare Robovetter thresholds for most metrics.
// The thresholds below are for the nominal run, used to compute the dispositoins displayed at NExScI.
// These thresholds can be adjusted to tune the Robovetter for a desired combination of completeness and reliability.
// The two commented-out blocks of thresholds below are two examples of other threshold combinations - one will produce a more reliable KOI population, and the other will produce a more complete KOI population.
//
const double           SWEET_THRESH = 50.0;  // Threshold for the SWEET test SNR
const double      HALO_GHOST_THRESH =  4.0;  // Threshold for the Ghostbuster core-to-halo SNR test
const double      SES_TO_MES_THRESH =  0.8;  // Threshold for the ses-to-mes ratio transit consistency test
const double ALL_TRAN_CHASES_THRESH =  0.8;  // Threshold for the all transits CHASES test
const double           SHAPE_THRESH = 1.04;  // Threshold for the shape metric test
const double     MOD_VAL1_DV_THRESH =  1.0;  // Threshold for the Modelshift Test #1 utilizing the DV detrending
const double     MOD_VAL2_DV_THRESH =  2.0;  // Threshold for the Modelshift Test #2 utilizing the DV detrending
const double     MOD_VAL3_DV_THRESH =  4.0;  // Threshold for the Modelshift Test #3 utilizing the DV detrending
const double    MOD_VAL1_ALT_THRESH = -3.0;  // Threshold for the Modelshift Test #1 utilizing the ALT detrending
const double    MOD_VAL2_ALT_THRESH =  1.0;  // Threshold for the Modelshift Test #2 utilizing the ALT detrending
const double    MOD_VAL3_ALT_THRESH =  1.0;  // Threshold for the Modelshift Test #3 utilizing the ALT detrending
const double          LPP_DV_THRESH =  2.2;  // Threshold for the LPP test utilizing the DV detrending
const double         LPP_ALT_THRESH =  3.2;  // Threshold for the LPP test utilizing the DV detrending
const double        RV_OE_DV_THRESH =  1.1;  // Threshold for the simple depth-based odd-even test utilizing the DV detrending
const double       RV_OE_ALT_THRESH =  1.1;  // Threshold for the simple depth-based odd-even test utilizing the ALT detrending
const double       MOD_OE_DV_THRESH = 11.2;  // Threshold for the Modelshift odd-even test utilizing the DV detrending
const double      MOD_OE_ALT_THRESH = 19.8;  // Threshold for the Robovetter odd-even test utilizing the DV detrending
const double     MOD_VAL4_DV_THRESH =  1.0;  // Threshold for the Modelshift Test #4 utilizing the DV detrending
const double     MOD_VAL5_DV_THRESH =  0.0;  // Threshold for the Modelshift Test #5 utilizing the DV detrending
const double     MOD_VAL6_DV_THRESH =  0.0;  // Threshold for the Modelshift Test #6 utilizing the DV detrending
const double    MOD_VAL4_ALT_THRESH =  1.0;  // Threshold for the Modelshift Test #4 utilizing the ALT detrending
const double    MOD_VAL5_ALT_THRESH =  0.0;  // Threshold for the Modelshift Test #5 utilizing the ALT detrending
const double    MOD_VAL6_ALT_THRESH =  0.0;  // Threshold for the Modelshift Test #6 utilizing the ALT detrending


// // Higher Reliability Set of Robovetter Thresholds
//
// const double           SWEET_THRESH = 50.00;
// const double      HALO_GHOST_THRESH =  4.00;
// const double      SES_TO_MES_THRESH =  0.75;
// const double ALL_TRAN_CHASES_THRESH =  0.55;
// const double           SHAPE_THRESH =  1.04;
// const double     MOD_VAL1_DV_THRESH = -1.0;
// const double     MOD_VAL2_DV_THRESH = -0.7;
// const double     MOD_VAL3_DV_THRESH = -1.6;
// const double    MOD_VAL1_ALT_THRESH = -4.3;
// const double    MOD_VAL2_ALT_THRESH =  2.5;
// const double    MOD_VAL3_ALT_THRESH =  0.2;
// const double          LPP_DV_THRESH =  2.7;
// const double         LPP_ALT_THRESH =  3.2;
// const double        RV_OE_DV_THRESH =  1.1;
// const double       RV_OE_ALT_THRESH =  1.1;
// const double       MOD_OE_DV_THRESH = 11.2;
// const double      MOD_OE_ALT_THRESH = 19.8;
// const double     MOD_VAL4_DV_THRESH =  1.0;
// const double     MOD_VAL5_DV_THRESH =  0.0;
// const double     MOD_VAL6_DV_THRESH =  0.0;
// const double    MOD_VAL4_ALT_THRESH =  1.0;
// const double    MOD_VAL5_ALT_THRESH =  0.0;
// const double    MOD_VAL6_ALT_THRESH =  0.0;


// //Higher Completness Set of Robovetter Thresholds
//
// const double           SWEET_THRESH = 50.0;
// const double      HALO_GHOST_THRESH =  4.0;
// const double      SES_TO_MES_THRESH =  0.9;
// const double ALL_TRAN_CHASES_THRESH =  1.0;
// const double           SHAPE_THRESH = 1.04;
// const double     MOD_VAL1_DV_THRESH =  2.4;
// const double     MOD_VAL2_DV_THRESH =  5.0;
// const double     MOD_VAL3_DV_THRESH =  7.5;
// const double    MOD_VAL1_ALT_THRESH = -2.5;
// const double    MOD_VAL2_ALT_THRESH = -0.5;
// const double    MOD_VAL3_ALT_THRESH =  0.5;
// const double          LPP_DV_THRESH =  2.8;
// const double         LPP_ALT_THRESH =  3.2;
// const double        RV_OE_DV_THRESH =  1.1;
// const double       RV_OE_ALT_THRESH =  1.1;
// const double       MOD_OE_DV_THRESH = 11.2;
// const double      MOD_OE_ALT_THRESH = 19.8;
// const double     MOD_VAL4_DV_THRESH =  1.0;
// const double     MOD_VAL5_DV_THRESH =  0.0;
// const double     MOD_VAL6_DV_THRESH =  0.0;
// const double    MOD_VAL4_ALT_THRESH =  1.0;
// const double    MOD_VAL5_ALT_THRESH =  0.0;
// const double    MOD_VAL6_ALT_THRESH =  0.0;


// Set up the structs needed to hold all the input data
//
struct valstruct  // Struct to contain a metric's value and associated positive and negative errors. Values and Errors are defaulted to zero. A value of zero defaults to a 'pass' for all metrics.
  {
  double             val=0.0;  // The value of the metric
  double            perr=0.0;  // The positive error on the metric
  double            merr=0.0;  // The negative error on the metric
  };

struct datastruct  // Struct that holds all the data for each individual TCE
  {
  string                 tce;  // TCE string (KIC-PN)
  string            comments;  // All minor flags/comments for each TCE
  string         disp_string;  // Final disposition string that containts the major disposition (PC/FP), major flag values, and minor flag values (comments)

  int                    kic;  // KIC number
  int                     pn;  // Planet Number
  int            num_planets;  // Highest PN number for the given KIC
  int         robo_cent_disp;  // 0 = centroid PC, 1 = centroid FP
  int       ephem_match_disp;  // 0 = no ephem match, 1 = ephem match identified
  int                nitrans;  // Number of individual transit events identified by TPS. Includes both real transits and SES=0 transits due to gaps.
  int             nrealtrans;  // Number of individual transits that were actually included to compute the MES, i.e., no ses=0 events.
  int             ngoodtrans;  // Number of "good" transits after "bad" transits have been ruled out by the Robovetter's individual transit metrics.
  int  itransfailflagints[6];  // 1/0 flags to record when at least one of a TCE's individual transit events failed due to a given individual transit metric.
  int    centroidflagints[9];  // 1/0 flags to record the minor flags of the Robovetter.

  int       not_transit_like;  // Not Transit-Like major flag.
  int        stellar_eclipse;  // Stellar Eclipse major flag.
  int     planet_occultation;  // Planet occultation flag (If set to 1, and stellar_eclipse set to 1, then disposition will be PC.)
  int        centroid_offset;  // Centroid Offset major flag.
  int       period_is_double;  // Period is double flag (If set to 1, and stellar eclipse set to 1, then disposition will be PC.)
  int        ephemeris_match;  // Ephemeris Match major flag.

  double      cent_score=1.0;  // Centroid RV score. Default to 1.0 (PC) in case of missing data.
  double          disp_score;  // Robovetter Score

  valstruct           period;  // Period of the system in days from DV
  valstruct            epoch;  // Epoch from DV
  valstruct              mes;  // MES from TPS/DV
  valstruct            depth;  // Depth of the transit from DV in ppm
  valstruct         duration;  // Duration of the transit from DV in hours
  valstruct            rstar;  // Size of the star in solar radii from DV
  valstruct              sma;  // Semn-major axis of the system in AU from DV
  valstruct           impact;  // Impact parameter of the system from DV
  valstruct   mod_sig_pri_dv;  // Significance of primary from Modelshift test run on DV detrending
  valstruct   mod_sig_sec_dv;  // Significance of secondary from Modelshift test run on DV detrending
  valstruct   mod_sig_ter_dv;  // Significance of tertiary from Modelshift test run on DV detrending
  valstruct   mod_sig_pos_dv;  // Significance of positive feature from Modelshift test run on DV detrending
  valstruct       mod_fa1_dv;  // False Alarm threshold 1 from Modelshift test run on DV detrending
  valstruct       mod_fa2_dv;  // False Alarm threshold 2 from Modelshift test run on DV detrending
  valstruct      mod_fred_dv;  // Ratio of Red noise / Gaussian noise from Modelshift test run on DV detrending
  valstruct    mod_ph_pri_dv;  // Phase of primary from Modelshift test run on DV detrending
  valstruct    mod_ph_sec_dv;  // Phase of secondary from Modelshift test run on DV detrending
  valstruct    mod_ph_ter_dv;  // Phase of tertiary from Modelshift test run on DV detrending
  valstruct    mod_ph_pos_dv;  // Phase of positive feature from Modelshift test run on DV detrending
  valstruct mod_depth_pri_dv;  // Depth of primary from Modelshift test run on DV detrending
  valstruct mod_depth_sec_dv;  // Depth of secondary from Modelshift test run on DV detrending
  valstruct        mod_oe_dv;  // Odd-Even metric from Modelshift test run on DV detrending
  valstruct  mod_sig_pri_alt;  // Significance of primary from Modelshift test run on ALT detrending
  valstruct  mod_sig_sec_alt;  // Significance of secondary from Modelshift test run on ALT detrending
  valstruct  mod_sig_ter_alt;  // Significance of tertiary from Modelshift test run on ALT detrending
  valstruct  mod_sig_pos_alt;  // Significance of positive feature from Modelshift test run on ALT detrending
  valstruct      mod_fa1_alt;  // False Alarm threshold 1 from Modelshift test run on ALT detrending
  valstruct      mod_fa2_alt;  // False Alarm threshold 2 from Modelshift test run on ALT detrending
  valstruct     mod_fred_alt;  // Ratio of Red noise / Gaussian noise from Modelshift test run on ALT detrending
  valstruct   mod_ph_pri_alt;  // Phase of primary from Modelshift test run on ALT detrending
  valstruct   mod_ph_sec_alt;  // Phase of secondary from Modelshift test run on ALT detrending
  valstruct   mod_ph_ter_alt;  // Phase of tertiary from Modelshift test run on ALT detrending
  valstruct   mod_ph_pos_alt;  // Phase of positive feature from Modelshift test run on ALT detrending
  valstruct mod_pridepth_alt;  // Depth of primary from Modelshift test run on ALT detrending
  valstruct mod_secdepth_alt;  // Depth of secondary from Modelshift test run on ALT detrending
  valstruct       mod_oe_alt;  // Odd-Even metric from Modelshift on ALT detrending
  valstruct         oesig_dv;  // Simple depth-based Odd-Even test done on the DV detrending.
  valstruct     sdepthsig_dv;  // Seasonal Depth Difference Test done on the DV detrending
  valstruct        oesig_alt;  // Simple depth-based Odd-Even test done on the ALT detrending.
  valstruct    sdepthsig_alt;  // Seasonal Depth Difference Test done on the ALT detrending
  valstruct           lpp_dv;  // Susan's LPP value based on DV detrending
  valstruct          lpp_alt;  // Susan's LPP value based on ALT detrending
  valstruct           alb_dv;  // Albedo based on secondary eclipse measurements from the model-shift test run on DV detrending
  valstruct            rp_dv;  // Planet Radius based on secondary eclipse measurements from the model-shift test run on DV detrending
  valstruct          alb_alt;  // Albedo based on secondary eclipse measurements from the model-shift test run on ALT detrending
  valstruct           rp_alt;  // Planet Radius based on secondary eclipse measurements from the model-shift test run on ALT detrending
  valstruct  modshiftval1_dv;  // modshiftval1_dv.val  = mod_sig_pri_dv.val/mod_fred_dv.val  -  mod_fa1_dv.val
  valstruct  modshiftval2_dv;  // modshiftval2_dv.val  = mod_sig_pri_dv.val  - mod_sig_ter_dv.val - mod_fa2_dv.val;
  valstruct  modshiftval3_dv;  // modshiftval3_dv.val  = mod_sig_pri_dv.val  - mod_sig_pos_dv.val - mod_fa2_dv.val;
  valstruct  modshiftval4_dv;  // modshiftval1_dv.val  = mod_sig_sec_dv.val/mod_fred_dv.val  -  mod_fa1_dv.val
  valstruct  modshiftval5_dv;  // modshiftval2_dv.val  = mod_sig_sec_dv.val  - mod_sig_ter_dv.val - mod_fa2_dv.val;
  valstruct  modshiftval6_dv;  // modshiftval3_dv.val  = mod_sig_sec_dv.val  - mod_sig_pos_dv.val - mod_fa2_dv.val;
  valstruct modshiftval1_alt;  // modshiftval1_alt.val = mod_sig_pri_alt.val/mod_fred_alt.val - mod_fa1_alt.val;
  valstruct modshiftval2_alt;  // modshiftval2_alt.val = mod_sig_pri_alt.val - mod_sig_ter_alt.val - mod_fa2_alt.val;
  valstruct modshiftval3_alt;  // modshiftval3_alt.val = mod_sig_pri_alt.val - mod_sig_pos_alt.val - mod_fa2_alt.val;
  valstruct modshiftval4_alt;  // modshiftval1_alt.val = mod_sig_sec_alt.val/mod_fred_alt.val - mod_fa1_alt.val;
  valstruct modshiftval5_alt;  // modshiftval2_alt.val = mod_sig_sec_alt.val - mod_sig_ter_alt.val - mod_fa2_alt.val;
  valstruct modshiftval6_alt;  // modshiftval3_alt.val = mod_sig_sec_alt.val - mod_sig_pos_alt.val - mod_fa2_alt.val;
  valstruct       ses_to_mes;  // Ratio of the maximum ses value to the mes value
  valstruct     shape_metric;  // Simple shape metric. shape_metric = impact + Rp/Rs from DV fit
  valstruct  all_tran_chases;  // CHASES test applied to all transits
  valstruct       halo_ghost;  // Ghostbutser test halo-to-core ratio
  valstruct        sweet_snr;  // SNR of SWEET Test
  valstruct        sweet_amp;  // Amplitude of SWEET Test
  valstruct          new_mes;  // New mes value re-calculated from the individual 'good' transits after applying individual transit tests
  valstruct        depth_trap;  // Depth of the trapezoidal fit to the ALT detrending in ppm
  };


// Declare the main "data" vector struct to hold all our data, and another, "dataorig", to backup all the original data for use in the Monte Carlo iterations
vector <datastruct> data, dataorig;  // Use of vector allows for dynamic memory allocation
int n;  // Global counting integer only used to loop through TCEs in the data vector


// The main function where inputs are taken, the Robovetter is run, the Robovetter is re-run through the Monte Carlo to compute scores, and the result file is outputted
int main (int argc, char* argv[])
  {
  int m;               // Integer used to count through the Monte Carlo iterations when calculating the Robovetter scores
  int NMC;             // The number of Monte Carlo iterations to perform when calculating the Robovetter scores
  int ntces;           // Number of TCEs to robo-vet
  int secfound;        // Integer to keep track if we have designated a TCE a secondary eclipse in the system
  string disp_pcfp;    // The PC or FP disposition of the TCE
  string  infilename;  // Name of the input data file
  string outfilename;  // Name of the output results file
  ofstream outfile;    // Stream to print out the results file

  // Get inputs - take from command line if given, else prompt the user
  if(argc>1)
    infilename = argv[1];
  else
    {
    cout << "Name of input file? ";
    cin >> infilename;
    }

  if(argc>2)
    outfilename = argv[2];
  else
    {
    cout << "Name of output file? ";
    cin >> outfilename;
    }

  if(argc>3)
    NMC = atoi(argv[3]);
  else
    {
    cout << "Number of Monte Carlo runs? ";
    cin >> NMC;
    }

  cout << "Reading in data..." << endl;
  ntces = READDATA(infilename);  // Read all Input Data. Returns the total number of TCEs.

  dataorig = data;  // Copy the big vector struct of all the data so that we can do monte carlo witout modifying original values.

  // Okay, let the dispositioning begin!
  cout << "Dispositioning..." << endl;
  outfile.open(outfilename.c_str());
  outfile << "#1:TCE  2:Robovetter_Score  3:Disposition  4:Not_Transit-Like_Flag  5:Stellar_Eclipse_Flag  6:Centroid Offset_Flag  7:Ephemeris_Match_Flag  8:Minor_Descriptive_Flags" << endl;  // Print header line for output file

  // Do the loop for the number of Monte Carlo runs.
  for(m=0;m<=NMC;m++)  // It is <=NMC because the 0th iteration is the 'nominal' run, and the 1->NMC runs are the runs used to compute scores
    {

    // Print out the current Monte Carlo iteration number to the command line so we know how fast the code is running. Step up by powers of 10 so command line is not flooded.
    if(m==0)
      cout << "NMC = " << flush;
    if(m>0 && m<=10)
      cout << m << "..." << flush;
    if(m>10 && m<=100 && m%10==0)
      cout << m << "..." << flush;
    if(m>100 && m<=1000 && m%100==0)
      cout << m << "..." << flush;
    if(m>1000 && m<=10000 && m%1000==0)
      cout << m << "..." << flush;
    if(m>10000 && m<=100000 && m%10000==0)
      cout << m << "..."  << flush;
    if(m==NMC)
      cout << "...DONE" << endl;


    // Do a pass through all TCEs first
    for(n=0;n<ntces;n++)
      {
      if(m==0)  // If it is the first iteration, assign errors to each metric
        SETERRORS();
      else      // If not the first iteration, and we are in Monte Carlo portion to compute scores, add random noise to each variable to compute scores
        ADDNOISE();

      // Set all flags and coments to default values
      data[n].not_transit_like=data[n].stellar_eclipse=data[n].planet_occultation=data[n].period_is_double=data[n].centroid_offset=data[n].ephemeris_match=0;  // Make sure all flags start at 0. PC until proven guilty.
      data[n].comments="";  // Default to no comments / minor flags
      }

    for(n=0;n<ntces;n++)
      {
      // Let's keep track if we already found a seconary eclipse in the system or not
      if(n==0 || data[n].kic!=data[n-1].kic || data[n-1].not_transit_like==0)  // If it's the first TCE we're looking at, or if it's a new system/KIC, or if a transit-like TCE was found in the system since we last found a secondary, start checking if the TCE is a secondary again.
        secfound=0;

      // Check to see if TCE is a secondary eclipse
      if(secfound==0)
        secfound = ISSEC();  // ISSEC() will return a 1 if the TCE was found to be a secondary, else 0. This makes sure we don't designate more than one secondary for every transit-like primary we find.

      // If the TCE is not a secondary, check to see if TCE is Transit-Like. (ISSEC will set both flags to 1 if it found the TCE is a secondary).
      if(data[n].not_transit_like==0 && data[n].stellar_eclipse==0)
        TRANSITLIKE();  // Will set not_transit_like to 1 if found to be not transit-like

      // Now if it is Transit-Like, check to see if there is a significant secondary eclipse
      if(data[n].not_transit_like==0 && data[n].stellar_eclipse==0)
        STELLARECLIPSE();  // Will set stellar_eclipse to 1 if found to be due to a stellar eclipse

      // Check for seasonal depth differences if it is transit-like
      if(data[n].not_transit_like==0)
        SEASONAL();

      // Check for any period aliasing if it transit-like
      if(data[n].not_transit_like==0)
        PERIODALIAS();

      // Apply centroid robovetter module dispositions and flags
      APPCENTROID();

      // Apply ephem match dispositions and flags
      APPEPHEMMATCH();


      // Make final PC/FP determination
      disp_pcfp="PC";  // Default to PC. Innocent until declared FP.
      if(data[n].not_transit_like==1)
        disp_pcfp="FP";
      if(data[n].stellar_eclipse==1 && data[n].planet_occultation==0)
        disp_pcfp="FP";
      if(data[n].centroid_offset==1)
        disp_pcfp="FP";
      if(data[n].ephemeris_match==1)
        disp_pcfp="FP";

      // If m>0, then record the result if it is a PC so we can compute the score (the fraction of Monte Carlo runs the disposition is PC).
      if(disp_pcfp=="PC" && m>0)
        data[n].disp_score += 1.0;

      // If it is the first iteration, save the major and minor flags.
      if(m==0)
        data[n].disp_string = disp_pcfp + " " + to_string(data[n].not_transit_like) + " " + to_string(data[n].stellar_eclipse) + " " + to_string(data[n].centroid_offset) + " " + to_string(data[n].ephemeris_match) + " " + data[n].comments;

      // If it is the last iteration, write to the output file
      if(m==NMC)
        {
        data[n].disp_score /= NMC;  // Compute the score - total number of PC dispositoins divided by the total number of Monte Carlo runs
        outfile << data[n].tce << " " << fixed << setprecision(3) << data[n].disp_score << " " << data[n].disp_string << endl;  // Write out the TCE name, Score, and then Major and Minor flags (comments)
        }
      }
    }
  outfile.close();
  }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to calculate whether the current TCE is the secondary of a previous TCE in the system
//
int ISSEC() {

int i;  // Counting integer
double epochthresh;  // Minimum threshold to consider the TCE distinct enough in time to be considered a secondary, and not a residual systematic left over after the primary was removed

// If a previous TCE in the system had a sec eclipse detected, and this TCE has the same period / diff epoch as it, then this is the sec eclipse
for(i=1;i<data[n].pn;i++)
  {
  vector<double> tmpdobs = COMPPT(data[n].period.val,data[n-data[n].pn+i].period.val,data[n].epoch.val,data[n-data[n].pn+i].epoch.val);  // Computes 5 metrics for comparing the period and epoch
  tmpdobs[2] = data[n-data[n].pn+i].period.val/data[n].period.val;  // Modify the period ratio so that it's always the other TCE divided by the current one for this test. That way when tmpdob>=2, it always means the current TCE is half the period or less of the previous one.
  epochthresh = WIDTHFAC*data[n-data[n].pn+i].duration.val/24.0;  // Calculate the threshold required to say whether or not the epochs overlap or not

  // If the criteria are met, disposition this TCE as a FP with both the N and S flags set to indicate it is a secondary
  if(data[n-data[n].pn+i].not_transit_like==0 && (data[n-data[n].pn+i].stellar_eclipse==1 || data[n-data[n].pn+i].period_is_double==1) && tmpdobs[0] > PSIG_THRESH && ((fabs(tmpdobs[3]) > epochthresh) || (fabs(tmpdobs[3]) < epochthresh && rint(tmpdobs[2])>=2 )) && ((fabs(tmpdobs[4]) > epochthresh) || (fabs(tmpdobs[4]) < epochthresh && rint(tmpdobs[2])>=2)) && ((tmpdobs[3]<0 && tmpdobs[4]<0) || (tmpdobs[3]>0 && tmpdobs[4]>0) || rint(tmpdobs[2])>=2) )  // Either same period and differnt epoch as a previous TCE, or half the period and same epoch. Both end up corresponding to the secondary eclipse.
    {
    data[n].not_transit_like=1;
    data[n].stellar_eclipse=1;
    if(data[n].comments!="")  data[n].comments+="---";
    data[n].comments+="IS_SEC_TCE";
    i=99;  // Only need to trigger once, so set j to high number to end loop
    }
  }

if(i==100)  // If the TCE was declared a secondary eclipse (j gets set to 99, then +1 at end of loop is 100) return 1, else return 0
  return(1);
else
  return(0);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to calculate whether the TCE is Transit-Like
//
void TRANSITLIKE() {

int i;  // Counting integer

// Check if, after rejected transits that fail individual transit metrics, the new (re-computed) MES is either less than 7.1, or there are less than 3 good transits left. If so, fail it and add flags.
if(data[n].new_mes.val < 7.1 || data[n].ngoodtrans < 3)
  {
  data[n].not_transit_like=1;

  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="INDIV_TRANS";

  if(data[n].itransfailflagints[0]==1)
    data[n].comments+= "_RUBBLE";

  if(data[n].itransfailflagints[1]==1)
    data[n].comments+= "_CHASES";

  if(data[n].itransfailflagints[2]==1)
    data[n].comments+= "_MARSHALL";

  if(data[n].itransfailflagints[3]==1)
    data[n].comments+= "_SKYE";

  if(data[n].itransfailflagints[4]==1)
    data[n].comments+= "_ZUMA";

  if(data[n].itransfailflagints[5]==1)
    data[n].comments+= "_TRACKER";
  }


// Check if more than half of the original transit events have failed individual transit events. If so, mark as not transit-like and add flags.
if(1.0*data[n].nrealtrans/data[n].nitrans <= 0.5)
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="TRANS_GAPPED";
  }


// SWEET Test - Look to see if the PDC signal has a sinusoidal signal at the TCE period that is greater than the depth of the TCE.  If so, mark as not transit-like and add flags.
if(data[n].sweet_snr.val > SWEET_THRESH && data[n].sweet_amp.val > data[n].depth.val && data[n].sweet_amp.val > data[n].depth_trap.val && data[n].period.val < 5.0)  // Sweet SNR is high
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="SWEET_NTL";
  }


// DV LPP Test
if(data[n].lpp_dv.val > LPP_DV_THRESH)
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="LPP_DV";
  }


// Alt LPP Test
if(data[n].lpp_alt.val > LPP_ALT_THRESH)
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="LPP_ALT";
  }


// Check the CHASES metric as run on all transits combined
if(data[n].all_tran_chases.val > ALL_TRAN_CHASES_THRESH)
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="ALL_TRANS_CHASES";
  }


// Check is primary is significant in DV Modelshift
if(data[n].modshiftval1_dv.val > MOD_VAL1_DV_THRESH && data[n].mod_sig_pri_dv.val > 0)
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="MOD_NONUNIQ_DV";
  }

// Check is primary is significantly greater than tertiary in DV Modelshift
if(data[n].modshiftval2_dv.val > MOD_VAL2_DV_THRESH && data[n].mod_sig_pri_dv.val > 0 && data[n].mod_sig_ter_dv.val > 0)  // 0 indicates NULL result
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="MOD_TER_DV";
  }

// Check is primary is significantly greater than positive in DV Modelshift
if(data[n].modshiftval3_dv.val > MOD_VAL3_DV_THRESH && data[n].mod_sig_pri_dv.val > 0 && data[n].mod_sig_pos_dv.val > 0)  // 0 indicates NULL result
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="MOD_POS_DV";
  }


// Check if primary is significant in ALT Modelshift
if(data[n].modshiftval1_alt.val > MOD_VAL1_ALT_THRESH && data[n].mod_sig_pri_alt.val > 0)
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="MOD_NONUNIQ_ALT";
  }

// Check is primary is significantly greater than tertiary in ALT Modelshift
if(data[n].modshiftval2_alt.val > MOD_VAL2_ALT_THRESH && data[n].mod_sig_pri_alt.val > 0 && data[n].mod_sig_ter_alt.val > 0)  // 0 indicates NULL result
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="MOD_TER_ALT";
  }

// Check is primary is significantly greater than positive in ALT Modelshift
if(data[n].modshiftval3_alt.val > MOD_VAL3_ALT_THRESH && data[n].mod_sig_pri_alt.val > 0 && data[n].mod_sig_pos_alt.val > 0)  // 0 indicates NULL result
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="MOD_POS_ALT";
  }


// Fail if DV and ALT data[n].mod_sig_pri_dv.val is 0 - means no fit exists for either, so it is very likely a not transit-like TCE
if(data[n].mod_sig_pri_dv.val==0 && data[n].mod_sig_pri_alt.val==0)
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="NO_FITS";
  }


// Check consistency of transits via SES to MES ratio
if(data[n].ses_to_mes.val > SES_TO_MES_THRESH && data[n].period.val > 90)  // Maybe 0.95?  0.99?   1.0 could work.
  {
  data[n].not_transit_like=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="INCONSISTENT_TRANS";
  }


// If a previous TCE in the system has the same period as this one, and it was not transit-like, then this one should be too. Also check if previous TCE was transit-like, then this TCE is triggering off the residual of the earlier TCE.
for(i=1;i<data[n].pn;i++)
  {
  vector<double> tmpdobs = COMPPT(data[n].period.val,data[n-data[n].pn+i].period.val,data[n].epoch.val,data[n-data[n].pn+i].epoch.val);  // Computes 5 metrics for comparing the period and epoch
  if(tmpdobs[0] > PSIG_THRESH)  // This TCE matches the period of a previous TCE in the system
    {
    if(data[n-data[n].pn+i].not_transit_like==1)  // If the previous TCE was deemed not transit-like, then this TCE should be not transit-like as well
      {
      data[n].not_transit_like=1;
      if(data[n].comments!="")  data[n].comments+="---";
      data[n].comments+="SAME_NTL_PERIOD";
      i=99;  // Only need to trigger once
      }
    else
      if(fabs(tmpdobs[3]) < WIDTHFAC*data[n-data[n].pn+i].duration.val/24.0 || fabs(tmpdobs[4]) < WIDTHFAC*data[n-data[n].pn+i].duration.val/24.0 || (tmpdobs[3]<0 && tmpdobs[4]>0) || (tmpdobs[3]>0 && tmpdobs[4]<0))  // If previous TCE was transit-like, check to see if this TCE is triggering on its residuals, i.e., it's epoch is within two transit durations.
        {
        data[n].not_transit_like=1;
        if(data[n].comments!="")  data[n].comments+="---";
        data[n].comments+="RESIDUAL_TCE";
        i=99;  // Only need to trigger once
        }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to calculate whether the TCE is shows signs of being a stellar eclipse
//
void STELLARECLIPSE() {

int i;  // Counting integer

// SWEET Test - Look to see if there is significant, sinusoidal out-of-transit variability that is less than the depth of the TCE.
if(data[n].sweet_snr.val > SWEET_THRESH && (data[n].sweet_amp.val < data[n].depth.val || data[n].sweet_amp.val < data[n].depth_trap.val) && data[n].sweet_amp.val > 5000 && data[n].period.val < 10.0) // Sweet SNR is high but amp is less than depth, so it's an EB
  {
  data[n].stellar_eclipse=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="SWEET_EB";
  }


// Check if planet is inside star via DV fit. Just flag it, don't fail though, as it is possible the stellar parameters are not reliable.
if(data[n].sma.val*214.939469384 < data[n].rstar.val)
  {
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="PLANET_IN_STAR";
  }


// Look for a secondary as indicated by the Modelshift test run on the DV detrending
if(data[n].modshiftval4_dv.val > MOD_VAL4_DV_THRESH && data[n].mod_sig_sec_dv.val > 0)  // See if secondary is significant
  if(data[n].modshiftval5_dv.val > MOD_VAL5_DV_THRESH || data[n].mod_sig_ter_dv.val <= 0)  // If ter measurement exists, check if sec is more significant
    if(data[n].modshiftval6_dv.val > MOD_VAL6_DV_THRESH || data[n].mod_sig_pos_dv.val <= 0)  // If pos measurement exists, check if sec is more significant
      {
      data[n].stellar_eclipse=1;
      if(data[n].comments!="")  data[n].comments+="---";
      data[n].comments+="MOD_SEC_DV";

      // Check to see if sec. eclipse could be due to a planet. Require a radius less than 30 R_earth, an albedo less than 1, a secondary 10 times or more shallow than the primary, and a impact parameter less than 0.95.
      if(data[n].alb_dv.val > 0.0 && data[n].alb_dv.val < 1.0 && data[n].rp_dv.val < 30.0 && data[n].mod_depth_sec_dv.val < 0.10*data[n].mod_depth_pri_dv.val && data[n].impact.val < 0.95)  // Check to see if occultation could be due to planet. Only apply to things with less than 30 earth radii, and if the secodary is less than 10% the depth of the primary, and if impact parameter is less than 0.9.
        {
        data[n].planet_occultation=1;
        if(data[n].comments!="")  data[n].comments+="---";
        data[n].comments+="PLANET_OCCULT_DV";
        }

      // Check is sec eclipse is identical to primary, indicating period is wrong by factor of 2
      if(fabs(0.5 - data[n].mod_ph_sec_dv.val)*data[n].period.val < 0.25*data[n].duration.val/24.0 && fabs(data[n].mod_sig_pri_dv.val - data[n].mod_sig_sec_dv.val) < data[n].mod_fa2_dv.val)  // Check to see if secondary could be identical to the transit so that it's really a PC phased at twice the period
        {
        data[n].period_is_double=1;
        if(data[n].comments!="")  data[n].comments+="---";
        data[n].comments+="PLANET_PERIOD_IS_HALF_DV";
        }
      }


// Look for a secondary as indicated by the Modelshift test run on the ALT detrending
if(data[n].modshiftval4_alt.val > MOD_VAL4_ALT_THRESH && data[n].mod_sig_sec_alt.val > 0)  // See if secondary is significant
  if(data[n].modshiftval5_alt.val > MOD_VAL5_ALT_THRESH || data[n].mod_sig_ter_alt.val <= 0)  // If ter measurement exists, check if sec is more significant
    if(data[n].modshiftval6_alt.val > MOD_VAL6_ALT_THRESH || data[n].mod_sig_pos_alt.val <= 0)  // If pos measurement exists, check if sec is more significant
      {
      data[n].stellar_eclipse=1;
      if(data[n].comments!="")  data[n].comments+="---";
      data[n].comments+="MOD_SEC_ALT";

      // Check to see if sec. eclipse could be due to a planet. Require a radius less than 30 R_earth, an albedo less than 1, a secondary 10 times or more shallow than the primary, and a impact parameter less than 0.95.
      if(data[n].alb_alt.val > 0.0 && data[n].alb_alt.val < 1.0 && data[n].rp_alt.val < 30.0 && data[n].mod_secdepth_alt.val < 0.1*data[n].mod_pridepth_alt.val && data[n].impact.val < 0.95)  // Check to see if occultation could be due to planet
        {
        data[n].planet_occultation=1;
        if(data[n].comments!="")  data[n].comments+="---";
        data[n].comments+="PLANET_OCCULT_ALT";
        }

      // Check is sec eclipse is identical to primary, indicating period is wrong by factor of 2
      if(fabs(0.5 - data[n].mod_ph_sec_alt.val)*data[n].period.val < 0.25*data[n].duration.val/24.0 && fabs(data[n].mod_sig_pri_alt.val - data[n].mod_sig_sec_alt.val) < data[n].mod_fa2_alt.val)  // Check to see if secondary could be identical to the transit so that it's really a PC phased at twice the period
        {
        data[n].period_is_double=1;
        if(data[n].comments!="")  data[n].comments+="---";
        data[n].comments+="PLANET_PERIOD_IS_HALF_ALT";
        }
      }


// Simple Depth-based Odd-Even Test using DV detrending
if(data[n].period.val < 90.0 && data[n].oesig_dv.val > RV_OE_DV_THRESH)
  {
  data[n].stellar_eclipse=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="DEPTH_ODDEVEN_DV";
  }

// Simple Depth-based Odd-Even Test using ALT detrending
if(data[n].period.val < 90.0 && data[n].oesig_alt.val > RV_OE_ALT_THRESH)
  {
  data[n].stellar_eclipse=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="DEPTH_ODDEVEN_ALT";
  }


// Modelshift Odd-Even Test using DV detrending
if(data[n].period.val < 90.0 && data[n].mod_oe_dv.val > MOD_OE_DV_THRESH)
  {
  data[n].stellar_eclipse=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="MOD_ODDEVEN_DV";
  }

// Modelshift Odd-Even Test using ALT detrending
if(data[n].period.val < 90.0 && data[n].mod_oe_alt.val > MOD_OE_ALT_THRESH)
  {
  data[n].stellar_eclipse=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="MOD_ODDEVEN_ALT";
  }


// Shape metric - combination of the impact parameter and radius ratio. I.e., if the transit is too deep and too V-shaped, then it is deemed a stellar eclipse
if(data[n].shape_metric.val > SHAPE_THRESH)
  {
  data[n].stellar_eclipse=1;
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="DEEP_V_SHAPED";
  }


// Look for a secondary detected as a subsequent TCE
for(i=1;i<=data[n].num_planets-data[n].pn;i++)  // Start looking at the next TCE in the system, and every subsequent TCE, through the last TCE in the system
  if(data[n].kic == data[n+i].kic)  // Check that indeed the next TCE (n+i) is the same KIC - we don't want to assume all TCEs from each KIC are in the input file
    {
    vector<double> tmpdobs = COMPPT(data[n].period.val,data[n+i].period.val,data[n].epoch.val,data[n+i].epoch.val);  // Computes 5 metrics for comparing the period and epoch
    if(tmpdobs[0] > PSIG_THRESH && (fabs(tmpdobs[3]) > WIDTHFAC*data[n].duration.val/24.0 || (fabs(tmpdobs[3]) < WIDTHFAC*data[n].duration.val/24.0 && rint(tmpdobs[2])>=2 )) && (fabs(tmpdobs[4]) > WIDTHFAC*data[n].duration.val/24.0 || (fabs(tmpdobs[4]) < WIDTHFAC*data[n].duration.val/24.0 && rint(tmpdobs[2])>=2)) && ((tmpdobs[3]<0 && tmpdobs[4]<0) || (tmpdobs[3]>0 && tmpdobs[4]>0) || rint(tmpdobs[2])>=2))  // Significant period match, at least 2 transit durations away, and either an insigniicant epoch match or more than a day apart. Put in condition so that if close periods, but some drift, will look for drift across the primary transit
      {
      data[n].stellar_eclipse=1;
      if(data[n].comments!="")  data[n].comments+="---";
      data[n].comments+="HAS_SEC_TCE";
      i=99;  // Only need to trigger this once
      }
    }


// Now if a secondary was detected, but the period could be double the true period, mark it as not actually having a secondary
if(data[n].period_is_double==1)
  data[n].stellar_eclipse=0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to apply the centrod disposition and minor flags as provided by the input file
//
void APPCENTROID() {

data[n].centroid_offset = data[n].robo_cent_disp;

if(data[n].centroidflagints[0]==1)
  {
  if(data[n].comments!="")
    data[n].comments+="---";
  data[n].comments+= "CENT_KIC_POS";
  }

if(data[n].centroidflagints[1]==1)
  {
  if(data[n].comments!="")
    data[n].comments+="---";
  data[n].comments+= "CENT_UNCERTAIN";
  }

if(data[n].centroidflagints[2]==1)
  {
  if(data[n].comments!="")
    data[n].comments+="---";
  data[n].comments+= "CENT_RESOLVED_OFFSET";
  }

if(data[n].centroidflagints[3]==1)
  {
  if(data[n].comments!="")
    data[n].comments+="---";
  data[n].comments+= "CENT_CROWDED";
  }

if(data[n].centroidflagints[4]==1)
  {
  if(data[n].comments!="")
    data[n].comments+="---";
  data[n].comments+= "CENT_NOFITS";
  }

if(data[n].centroidflagints[5]==1)
  {
  if(data[n].comments!="")
    data[n].comments+="---";
  data[n].comments+= "CENT_UNRESOLVED_OFFSET";
  }

if(data[n].centroidflagints[6]==1)
  {
  if(data[n].comments!="")
    data[n].comments+="---";
  data[n].comments+= "CENT_SATURATED";
  }

if(data[n].centroidflagints[7]==1)
  {
  if(data[n].comments!="")
    data[n].comments+="---";
  data[n].comments+= "CENT_FEW_MEAS";
  }

if(data[n].centroidflagints[8]==1)
  {
  if(data[n].comments!="")
    data[n].comments+="---";
  data[n].comments+= "CENT_FEW_DIFFS";
  }

// Perform ghostbuster test
if(data[n].halo_ghost.val > HALO_GHOST_THRESH)
  {
  data[n].centroid_offset = 1;
  if(data[n].comments!="") data[n].comments+="---";
    data[n].comments+="HALO_GHOST";
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to apply the ephemeris matching results as provided by the input file
//
void APPEPHEMMATCH() {

data[n].ephemeris_match = data[n].ephem_match_disp;
if(data[n].ephem_match_disp==1)
  {
  if(data[n].comments!="")
    data[n].comments+="---";
  data[n].comments += "EPHEM_MATCH";
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to calculate whether there is significant seasonal depth variations, as provided by the input file.
// Will only flag via the minor flags. Will not fail as seasonal depth variations can arise in PCs if there is a bright star on the edge of the pixel stamp image
//
void SEASONAL() {

if(data[n].sdepthsig_dv.val > 3.6)  // DV
  {
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="SEASONAL_DEPTH_DV";
  }

if(data[n].sdepthsig_alt.val > 3.6) // ALT
  {
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="SEASONAL_DEPTH_ALT";
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to calculate whether the Modelshift results indicate that this TCE is detected at many times the true period of the signal
//
void PERIODALIAS() {

double tmpdob0 = 0;  // Temporary double used in calculations
double tmpdob1 = 0;  // Temporary double used in calculations
double tmpdob2 = 0;  // Temporary double used in calculations
double tmpdob3 = 0;  // Temporary double used in calculations
double tmpdob4 = 0;  // Temporary double used in calculations

// If the primary, secondary, and tertiary events from DV Modelshift are statistically signifant, check to see if they are equidistant in phase. If so, calculate the true period and flag the TCE.
if(data[n].mod_ph_sec_dv.val > 0 && data[n].mod_ph_ter_dv.val > 0 && data[n].mod_sig_sec_dv.val/data[n].mod_fred_dv.val > data[n].mod_fa1_dv.val && data[n].mod_sig_ter_dv.val/data[n].mod_fred_dv.val > data[n].mod_fa1_dv.val)
  {
  if(data[n].mod_ph_sec_dv.val <=0.5)
    tmpdob0 = 1.0/fabs(data[n].mod_ph_sec_dv.val);
  else
    tmpdob0 = 1.0/fabs(data[n].mod_ph_sec_dv.val-1.0);

  if(data[n].mod_ph_ter_dv.val <=0.5)
    tmpdob1 = 1.0/fabs(data[n].mod_ph_ter_dv.val);
  else
    tmpdob1 = 1.0/fabs(data[n].mod_ph_ter_dv.val-1.0);

  tmpdob2 = sqrt(2)*INVERFC(fabs(tmpdob0 - rint(tmpdob0)));
  tmpdob3 = sqrt(2)*INVERFC(fabs(tmpdob1 - rint(tmpdob1)));

  if(tmpdob0 > tmpdob1)
    tmpdob4 = tmpdob0;
  else
    tmpdob4 = tmpdob1;
  }

if(tmpdob2 > ESIG_THRESH && tmpdob3 > ESIG_THRESH)
  {
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="PERIOD_ALIAS_DV";
  }


// If the primary, secondary, and tertiary events from ALT Modelshift are statistically signifant, check to see if they are equidistant in phase.  If so, calculate the true period and flag the TCE.
tmpdob0 = tmpdob1 = tmpdob2 = tmpdob3 = tmpdob4 = 0;
if(data[n].mod_ph_sec_alt.val > 0.0 && data[n].mod_ph_ter_alt.val > 0.0 && data[n].mod_sig_sec_alt.val/data[n].mod_fred_alt.val > data[n].mod_fa1_alt.val && data[n].mod_sig_ter_alt.val/data[n].mod_fred_alt.val > data[n].mod_fa1_alt.val)
  {
  if(data[n].mod_ph_sec_alt.val <=0.5)
    tmpdob0 = 1.0/fabs(data[n].mod_ph_sec_alt.val);
  else
    tmpdob0 = 1.0/fabs(data[n].mod_ph_sec_alt.val-1.0);

  if(data[n].mod_ph_ter_alt.val <=0.5)
    tmpdob1 = 1.0/fabs(data[n].mod_ph_ter_alt.val);
  else
    tmpdob1 = 1.0/fabs(data[n].mod_ph_ter_alt.val-1.0);

  tmpdob2 = sqrt(2)*INVERFC(fabs(tmpdob0 - rint(tmpdob0)));
  tmpdob3 = sqrt(2)*INVERFC(fabs(tmpdob1 - rint(tmpdob1)));

  if(tmpdob0 > tmpdob1)
    tmpdob4 = tmpdob0;
  else
    tmpdob4 = tmpdob1;
  }

if(tmpdob2 > ESIG_THRESH && tmpdob3 > ESIG_THRESH)
  {
  if(data[n].comments!="")  data[n].comments+="---";
  data[n].comments+="PERIOD_ALIAS_ALT";
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to compute how close two sets of period and epoch are, taking into account the mission duration.
// Returns a vector of five different metrics, the components of which are:
// 0 - Period match significance
// 1 - Epoch match significance
// 2 - Period Ratio
// 3 - Difference in epoch
// 4 - Difference in epoch by end of mission
//
vector<double> COMPPT(double P1, double P2, double T1, double T2) {

vector<double> tmpdobs;  // The vector that will be returned
double tmpdob0, tmpdob1, tmpdob2, tmpdob3, tmpdob4;  // Temporary doubles for calculations

if(P1 < P2)
  {
  tmpdob0 = PSIG(P1,P2);  // Period match significance
  tmpdob1 = ESIG(T1,T2,P1);  // Epoch match significance
  tmpdob2 = P2/P1;  // Period Ratio
  tmpdob3 = T1 - T2;  // Difference in epoch
  while(tmpdob3 > 0.5*P1)
    tmpdob3 -= P1;
  while(tmpdob3 < -0.5*P1)
    tmpdob3 += P1;
  tmpdob4 = tmpdob3 + int(MISSIONDUR/P2)*(rint(tmpdob2)*P1-P2);  // Difference in epoch by end of mission

  }
else
  {
  tmpdob0 = PSIG(P2,P1);  // Period match significance
  tmpdob1 = ESIG(T2,T1,P2);  // Epoch match significance
  tmpdob2 = P1/P2;  // Period Ratio
  tmpdob3 = T2 - T1;  // Difference in epoch
  while(tmpdob3 > 0.5*P2)
    tmpdob3 -= P2;
  while(tmpdob3 < -0.5*P2)
    tmpdob3 += P2;
  tmpdob4 = tmpdob3 + int(MISSIONDUR/P1)*(rint(tmpdob2)*P2-P1);  // Difference in epoch by end of mission
  }

tmpdobs.push_back(tmpdob0);
tmpdobs.push_back(tmpdob1);
tmpdobs.push_back(tmpdob2);
tmpdobs.push_back(tmpdob3);
tmpdobs.push_back(tmpdob4);

return tmpdobs;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to compute the closeness of a match between two periods, accounting for any period ratio
//
double PSIG(double P1, double P2) {
  return sqrt(2)*INVERFC(fabs((P1-P2)/P1 - rint((P1-P2)/P1)));
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to compute the closeness of a match between two epochs, accounting for any period ratio
//
double ESIG(double E1, double E2, double P) {
  return sqrt(2)*INVERFC(fabs((E1-E2)/P - rint((E1-E2)/P)));
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The complementary inverse error function
//
double INVERFC(double p) {
  int i;  // Counting integer
  double x, err, t, pp;

  if (p >= 2.) return -100.;
  if (p <= 0.0) return 100.;

  pp=(p < 1.0)? p:2.-p;
  t=sqrt(-2.*log(pp/2));

  x= -0.70711*((2.30753+t*0.27061)/(1+t*(0.99229+t*0.04481))-t);

  for (i=0;i<2;i++) {
    err=erfc(x)-pp;
    x += err/(1.12837916709551257*exp(-x*x)-x*err);
  }
  return (p<1.0 ? x:-x);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to generate a new value, applying asymmetric Gaussian error to a given a value, with given positive and negative errors
//
double AGRAND(double center, double psigma, double msigma)  // Create asymmetric random gaussian number with center and separate postitive and negative sigma values
  {
  int i;  // Counting integer
  double u1,u2,out;  // Doubles used in the computation

  if(psigma < 0)  // Make sure psigma is positive
    psigma = fabs(psigma);
  if(msigma < 0)  // Make sure msigma is positive
    msigma = fabs(msigma);

  i=0;
  while(i==0) // Create Gaussian Noise Array of sigma S
    {
    u1 = 2.0*dist(randgen);
    if(u1>1.0)
      u1 = u1 - 2.0;
    u2 = 2.0*dist(randgen);
    if(u2>1.0)
      u2 = u2 - 2.0;
    if(u1*u1+u2*u2>0 && u1*u1+u2*u2<1)
      {
      if(dist(randgen) < 0.5)
        out = center + psigma*fabs(u1*sqrt(-2*log(u1*u1+u2*u2)/(u1*u1+u2*u2)));
      else
        out = center - msigma*fabs(u1*sqrt(-2*log(u1*u1+u2*u2)/(u1*u1+u2*u2)));
      i=1;
      }
    }
  return out;
  }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function to read in the data input file
//
int READDATA(string infilename)
  {
  int i;  // Counting integer
  string tmpstr;  // String to dump data to
  ifstream infile;  // C++ ifsteram to open the file
  datastruct tmpdata;  // Temp struct for the data

  // Open the file
  infile.open(infilename.c_str());
  if(infile.fail()==1) // If file doesn't exist, exit with warning
    {
    cout << infilename << " doesn't exist or cannot open file..." << endl;
    exit(0);
    }

  for(n=0;n<223;n++)
    infile.ignore(9E9,'\n');  // Ignore header, 223 lines

  n=0;
  infile >> tmpstr;  // First column is data_set so ignore it
  while(!infile.eof())
    {
    infile >> tmpdata.tce;
    infile >> tmpdata.kic;
    infile >> tmpdata.pn;
    infile >> tmpdata.num_planets;
    infile >> tmpdata.nitrans;
    infile >> tmpdata.nrealtrans;
    infile >> tmpdata.ngoodtrans;
    for(i=0;i<6;i++)
      infile >> tmpdata.itransfailflagints[i];
    for(i=0;i<9;i++)
      infile >> tmpdata.centroidflagints[i];
    infile >> tmpdata.ephem_match_disp;
    infile >> tmpdata.robo_cent_disp;
    infile >> tmpdata.cent_score;
    infile >> tmpdata.period.val;
    infile >> tmpdata.period.perr;
    tmpdata.period.merr = tmpdata.period.perr;
    infile >> tmpdata.epoch.val;
    infile >> tmpdata.epoch.perr;
    tmpdata.epoch.merr = tmpdata.epoch.perr;
    infile >> tmpdata.duration.val;
    infile >> tmpdata.duration.perr;
    tmpdata.duration.merr = tmpdata.duration.perr;
    infile >> tmpdata.impact.val;
    infile >> tmpdata.impact.perr;
    tmpdata.impact.merr = tmpdata.impact.perr;
    infile >> tmpdata.depth.val;
    infile >> tmpdata.depth.perr;
    tmpdata.depth.merr = tmpdata.depth.perr;
    infile >> tmpdata.depth_trap.val;
    infile >> tmpdata.sma.val;
    infile >> tmpdata.rstar.val;
    infile >> tmpdata.mes.val;
    infile >> tmpdata.ses_to_mes.val;
    infile >> tmpdata.new_mes.val;
    infile >> tmpdata.lpp_dv.val;
    infile >> tmpdata.lpp_alt.val;
    infile >> tmpdata.all_tran_chases.val;
    infile >> tmpdata.sweet_snr.val;
    infile >> tmpdata.sweet_amp.val;
    infile >> tmpdata.shape_metric.val;
    infile >> tmpdata.halo_ghost.val;
    infile >> tmpdata.mod_sig_pri_dv.val;
    infile >> tmpdata.mod_sig_sec_dv.val;
    infile >> tmpdata.mod_sig_ter_dv.val;
    infile >> tmpdata.mod_sig_pos_dv.val;
    infile >> tmpdata.mod_fred_dv.val;
    infile >> tmpdata.mod_fa1_dv.val;
    infile >> tmpdata.mod_fa2_dv.val;
    infile >> tmpdata.mod_sig_pri_alt.val;
    infile >> tmpdata.mod_sig_sec_alt.val;
    infile >> tmpdata.mod_sig_ter_alt.val;
    infile >> tmpdata.mod_sig_pos_alt.val;
    infile >> tmpdata.mod_fred_alt.val;
    infile >> tmpdata.mod_fa1_alt.val;
    infile >> tmpdata.mod_fa2_alt.val;
    infile >> tmpdata.modshiftval1_dv.val;
    infile >> tmpdata.modshiftval2_dv.val;
    infile >> tmpdata.modshiftval3_dv.val;
    infile >> tmpdata.modshiftval4_dv.val;
    infile >> tmpdata.modshiftval5_dv.val;
    infile >> tmpdata.modshiftval6_dv.val;
    infile >> tmpdata.modshiftval1_alt.val;
    infile >> tmpdata.modshiftval2_alt.val;
    infile >> tmpdata.modshiftval3_alt.val;
    infile >> tmpdata.modshiftval4_alt.val;
    infile >> tmpdata.modshiftval5_alt.val;
    infile >> tmpdata.modshiftval6_alt.val;
    infile >> tmpdata.oesig_dv.val;
    infile >> tmpdata.oesig_alt.val;
    infile >> tmpdata.mod_oe_dv.val;
    infile >> tmpdata.mod_oe_alt.val;
    infile >> tmpdata.rp_dv.val;
    infile >> tmpdata.rp_dv.perr;
    infile >> tmpdata.rp_dv.merr;
    infile >> tmpdata.alb_dv.val;
    infile >> tmpdata.alb_dv.perr;
    infile >> tmpdata.alb_dv.merr;
    infile >> tmpdata.mod_depth_pri_dv.val;
    infile >> tmpdata.mod_depth_sec_dv.val;
    infile >> tmpdata.mod_ph_sec_dv.val;
    infile >> tmpdata.mod_ph_ter_dv.val;
    infile >> tmpdata.rp_alt.val;
    infile >> tmpdata.rp_alt.perr;
    infile >> tmpdata.rp_alt.merr;
    infile >> tmpdata.alb_alt.val;
    infile >> tmpdata.alb_alt.perr;
    infile >> tmpdata.alb_alt.merr;
    infile >> tmpdata.mod_pridepth_alt.val;
    infile >> tmpdata.mod_secdepth_alt.val;
    infile >> tmpdata.mod_ph_sec_alt.val;
    infile >> tmpdata.mod_ph_ter_alt.val;
    infile >> tmpdata.sdepthsig_dv.val;
    infile >> tmpdata.sdepthsig_alt.val;
    data.push_back(tmpdata);  // Add to the vector - dynamically allocates memory as needed
    n++;
    infile >> tmpstr;
    }
  infile.close();

  return(n);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Function to calculate the positive and negative errors on each metric, for use in computing the scores.
  // The hard-coded values from from fitting an asymmetrical Median Absolute Deviation to the on-target planet injection population
  // Note that not all metrics will have calculatable errors - if they are absent below, they default to errors of 0.0
  //
  void SETERRORS() {

  dataorig[n].lpp_dv.perr           = 0.5397772069560*pow(dataorig[n].period.val, 0.02247106554120)*pow(dataorig[n].mes.val,-0.337093838422000);
  dataorig[n].lpp_dv.merr           = 0.3160579035580*pow(dataorig[n].period.val, 0.04550069988450)*pow(dataorig[n].mes.val,-0.293554956782000);
  dataorig[n].lpp_alt.perr          = 0.7490695137030*pow(dataorig[n].period.val,-0.00218980170285)*pow(dataorig[n].mes.val,-0.379157259846000);
  dataorig[n].lpp_alt.merr          = 0.4044933570010*pow(dataorig[n].period.val,-0.00538770927072)*pow(dataorig[n].mes.val,-0.198142651931000);
  dataorig[n].all_tran_chases.perr  = 0.0719123465792*pow(dataorig[n].period.val, 0.41175064199500)*pow(dataorig[n].mes.val,-0.659862269823000);
  dataorig[n].all_tran_chases.merr  = 0.6439450737480*pow(dataorig[n].period.val, 0.20025392291000)*pow(dataorig[n].mes.val,-1.172196362760000);
  dataorig[n].halo_ghost.perr       = 2.9482209369900*pow(dataorig[n].period.val,-0.02190951413390)*pow(dataorig[n].mes.val,-0.857819594071000);
  dataorig[n].halo_ghost.merr       = 0.8937698849180*pow(dataorig[n].period.val,-0.00390653359316)*pow(dataorig[n].mes.val,-0.571733988953000);
  dataorig[n].modshiftval1_dv.perr  = 0.0359851586684*pow(dataorig[n].period.val, 0.34526817970300)*pow(dataorig[n].mes.val, 0.927869889567000);
  dataorig[n].modshiftval1_dv.merr  = 0.1172076101240*pow(dataorig[n].period.val, 0.00849023837934)*pow(dataorig[n].mes.val, 1.032036176590000);
  dataorig[n].modshiftval2_dv.perr  = 0.1236488661910*pow(dataorig[n].period.val, 0.00891441491596)*pow(dataorig[n].mes.val, 1.102582134940000);
  dataorig[n].modshiftval2_dv.merr  = 0.0422097490826*pow(dataorig[n].period.val, 0.26289907325300)*pow(dataorig[n].mes.val, 1.114357224830000);
  dataorig[n].modshiftval3_dv.perr  = 0.1311128208920*pow(dataorig[n].period.val, 0.01381750915310)*pow(dataorig[n].mes.val, 1.081788756960000);
  dataorig[n].modshiftval3_dv.merr  = 0.0398861608861*pow(dataorig[n].period.val, 0.25594508971500)*pow(dataorig[n].mes.val, 1.134482812590000);
  dataorig[n].modshiftval4_dv.perr  = 0.0664371610204*pow(dataorig[n].period.val, 0.42952841123900)*pow(dataorig[n].mes.val,-0.024285611070500);
  dataorig[n].modshiftval4_dv.merr  = 0.2907981857700*pow(dataorig[n].period.val, 0.12742094229400)*pow(dataorig[n].mes.val,-0.091871085127600);
  dataorig[n].modshiftval5_dv.perr  = 0.0789837740045*pow(dataorig[n].period.val, 0.41201629331000)*pow(dataorig[n].mes.val,-0.000324054693769);
  dataorig[n].modshiftval5_dv.merr  = 0.0844302407587*pow(dataorig[n].period.val, 0.22419810187900)*pow(dataorig[n].mes.val, 0.114666081746000);
  dataorig[n].modshiftval6_dv.perr  = 0.1226360702260*pow(dataorig[n].period.val, 0.33528497119100)*pow(dataorig[n].mes.val, 0.113309539340000);
  dataorig[n].modshiftval6_dv.merr  = 0.3613052340200*pow(dataorig[n].period.val, 0.19188256199400)*pow(dataorig[n].mes.val,-0.049283582758000);
  dataorig[n].modshiftval1_alt.perr = 0.0189792764423*pow(dataorig[n].period.val, 0.47829450400900)*pow(dataorig[n].mes.val, 1.050287856470000);
  dataorig[n].modshiftval1_alt.merr = 0.0525292208090*pow(dataorig[n].period.val, 0.43726986699200)*pow(dataorig[n].mes.val, 1.131079878130000);
  dataorig[n].modshiftval2_alt.perr = 0.0631248378436*pow(dataorig[n].period.val, 0.11483007676400)*pow(dataorig[n].mes.val, 1.221681217290000);
  dataorig[n].modshiftval2_alt.merr = 0.1796980205730*pow(dataorig[n].period.val, 0.26046520761600)*pow(dataorig[n].mes.val, 0.860137186619000);
  dataorig[n].modshiftval3_alt.perr = 0.0598379093547*pow(dataorig[n].period.val, 0.12040380816200)*pow(dataorig[n].mes.val, 1.226771599810000);
  dataorig[n].modshiftval3_alt.merr = 0.1899071992290*pow(dataorig[n].period.val, 0.25354584282800)*pow(dataorig[n].mes.val, 0.856120926414000);
  dataorig[n].modshiftval4_alt.perr = 0.1101513129750*pow(dataorig[n].period.val, 0.35132868047500)*pow(dataorig[n].mes.val,-0.048852263483900);
  dataorig[n].modshiftval4_alt.merr = 0.1083669637780*pow(dataorig[n].period.val, 0.27307896478300)*pow(dataorig[n].mes.val,-0.009711959339870);
  dataorig[n].modshiftval5_alt.perr = 0.3935371566150*pow(dataorig[n].period.val, 0.01912753954750)*pow(dataorig[n].mes.val,-0.050275761627100);
  dataorig[n].modshiftval5_alt.merr = 0.2980949739620*pow(dataorig[n].period.val,-0.05578527333360)*pow(dataorig[n].mes.val, 0.002266302822800);
  dataorig[n].modshiftval6_alt.perr = 0.5796868236340*pow(dataorig[n].period.val,-0.02835559900120)*pow(dataorig[n].mes.val, 0.049987624130100);
  dataorig[n].modshiftval6_alt.merr = 0.7411287382880*pow(dataorig[n].period.val,-0.05142749975300)*pow(dataorig[n].mes.val,-0.039861639027500);
  dataorig[n].ses_to_mes.perr       = 0.0369526604942*pow(dataorig[n].period.val, 0.14959656139300)*pow(dataorig[n].mes.val,-0.064601117291500);
  dataorig[n].ses_to_mes.merr       = 0.0226668330445*pow(dataorig[n].period.val, 0.20183714233600)*pow(dataorig[n].mes.val,-0.083435212884100);
  dataorig[n].oesig_dv.perr         = 0.0104399394241*pow(dataorig[n].period.val, 0.73020957218500)*pow(dataorig[n].mes.val, 0.122579886640000);
  dataorig[n].oesig_dv.merr         = 0.0222665648308*pow(dataorig[n].period.val, 0.57648534599100)*pow(dataorig[n].mes.val,-0.032124993308200);
  dataorig[n].oesig_alt.perr        = 0.0033989952765*pow(dataorig[n].period.val, 0.88017571820700)*pow(dataorig[n].mes.val, 0.228459350993000);
  dataorig[n].oesig_alt.merr        = 0.0159254420463*pow(dataorig[n].period.val, 0.62179263457300)*pow(dataorig[n].mes.val, 0.001876976549390);
  dataorig[n].mod_oe_dv.perr        = 0.3049992884340*pow(dataorig[n].period.val, 0.14051796494200)*pow(dataorig[n].mes.val, 0.291823245766000);
  dataorig[n].mod_oe_dv.merr        = 0.2328598160200*pow(dataorig[n].period.val, 0.11908199333700)*pow(dataorig[n].mes.val, 0.269773182986000);
  dataorig[n].mod_oe_alt.perr       = 0.3333666684530*pow(dataorig[n].period.val, 0.18379427556600)*pow(dataorig[n].mes.val, 0.282428580033000);
  dataorig[n].mod_oe_alt.merr       = 0.3999634117150*pow(dataorig[n].period.val, 0.11272637617500)*pow(dataorig[n].mes.val, 0.153009297647000);
  dataorig[n].new_mes.perr          = 0.0750311183399*pow(dataorig[n].period.val,-0.01825752820980)*pow(dataorig[n].mes.val, 0.898854699844000);
  dataorig[n].new_mes.merr          = 0.0434497036450*pow(dataorig[n].period.val, 0.07975608411790)*pow(dataorig[n].mes.val, 0.884972257143000);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Function to add noise to the original value for each metric, given their errors
  //
  void ADDNOISE() {

  // Add noise to parameters with errors calculated based on injection, i.e., those computed in SETERRORS
  data[n].lpp_dv.val           = AGRAND(dataorig[n].lpp_dv.val,          dataorig[n].lpp_dv.perr,          dataorig[n].lpp_dv.merr);
  data[n].lpp_alt.val          = AGRAND(dataorig[n].lpp_alt.val,         dataorig[n].lpp_alt.perr,         dataorig[n].lpp_alt.merr);
  data[n].all_tran_chases.val  = AGRAND(dataorig[n].all_tran_chases.val, dataorig[n].all_tran_chases.perr, dataorig[n].all_tran_chases.merr);
  data[n].halo_ghost.val       = AGRAND(dataorig[n].halo_ghost.val,      dataorig[n].halo_ghost.perr,      dataorig[n].halo_ghost.merr);
  data[n].modshiftval1_dv.val  = AGRAND(dataorig[n].modshiftval1_dv.val, dataorig[n].modshiftval1_dv.perr, dataorig[n].modshiftval1_dv.merr);
  data[n].modshiftval2_dv.val  = AGRAND(dataorig[n].modshiftval2_dv.val, dataorig[n].modshiftval2_dv.perr, dataorig[n].modshiftval2_dv.merr);
  data[n].modshiftval3_dv.val  = AGRAND(dataorig[n].modshiftval3_dv.val, dataorig[n].modshiftval3_dv.perr, dataorig[n].modshiftval3_dv.merr);
  data[n].modshiftval4_dv.val  = AGRAND(dataorig[n].modshiftval4_dv.val, dataorig[n].modshiftval4_dv.perr, dataorig[n].modshiftval4_dv.merr);
  data[n].modshiftval5_dv.val  = AGRAND(dataorig[n].modshiftval5_dv.val, dataorig[n].modshiftval5_dv.perr, dataorig[n].modshiftval5_dv.merr);
  data[n].modshiftval6_dv.val  = AGRAND(dataorig[n].modshiftval6_dv.val, dataorig[n].modshiftval6_dv.perr, dataorig[n].modshiftval6_dv.merr);
  data[n].modshiftval1_alt.val = AGRAND(dataorig[n].modshiftval1_alt.val,dataorig[n].modshiftval1_alt.perr,dataorig[n].modshiftval1_alt.merr);
  data[n].modshiftval2_alt.val = AGRAND(dataorig[n].modshiftval2_alt.val,dataorig[n].modshiftval2_alt.perr,dataorig[n].modshiftval2_alt.merr);
  data[n].modshiftval3_alt.val = AGRAND(dataorig[n].modshiftval3_alt.val,dataorig[n].modshiftval3_alt.perr,dataorig[n].modshiftval3_alt.merr);
  data[n].modshiftval4_alt.val = AGRAND(dataorig[n].modshiftval4_alt.val,dataorig[n].modshiftval4_alt.perr,dataorig[n].modshiftval4_alt.merr);
  data[n].modshiftval5_alt.val = AGRAND(dataorig[n].modshiftval5_alt.val,dataorig[n].modshiftval5_alt.perr,dataorig[n].modshiftval5_alt.merr);
  data[n].modshiftval6_alt.val = AGRAND(dataorig[n].modshiftval6_alt.val,dataorig[n].modshiftval6_alt.perr,dataorig[n].modshiftval6_alt.merr);
  data[n].ses_to_mes.val       = AGRAND(dataorig[n].ses_to_mes.val,      dataorig[n].ses_to_mes.perr,      dataorig[n].ses_to_mes.merr);
  data[n].oesig_dv.val         = AGRAND(dataorig[n].oesig_dv.val,        dataorig[n].oesig_dv.perr,        dataorig[n].oesig_dv.merr);
  data[n].oesig_alt.val        = AGRAND(dataorig[n].oesig_alt.val,       dataorig[n].oesig_alt.perr,       dataorig[n].oesig_alt.merr);
  data[n].mod_oe_dv.val        = AGRAND(dataorig[n].mod_oe_dv.val,       dataorig[n].mod_oe_dv.perr,       dataorig[n].mod_oe_dv.merr);
  data[n].mod_oe_alt.val       = AGRAND(dataorig[n].mod_oe_alt.val,      dataorig[n].mod_oe_alt.perr,      dataorig[n].mod_oe_alt.merr);
  data[n].new_mes.val          = AGRAND(dataorig[n].new_mes.val,         dataorig[n].new_mes.perr,         dataorig[n].new_mes.merr);

  // Add noise to parameters with errors from the input file
  data[n].period.val   = AGRAND(dataorig[n].period.val,  dataorig[n].period.perr,  dataorig[n].period.merr);
  while(data[n].period.val<=0.0)  // If period randomly got perturbed to a negative value, try again
    data[n].period.val = AGRAND(dataorig[n].period.val,  dataorig[n].period.perr,  dataorig[n].period.merr);  // Ensure period is always positive
  data[n].epoch.val    = AGRAND(dataorig[n].epoch.val,   dataorig[n].epoch.perr,   dataorig[n].epoch.merr);
  data[n].duration.val = AGRAND(dataorig[n].duration.val,dataorig[n].duration.perr,dataorig[n].duration.merr);
  data[n].depth.val    = AGRAND(dataorig[n].depth.val,   dataorig[n].depth.perr,   dataorig[n].depth.merr);
  data[n].impact.val   = AGRAND(dataorig[n].impact.val,  dataorig[n].impact.perr,  dataorig[n].impact.merr);
  data[n].alb_dv.val   = AGRAND(dataorig[n].alb_dv.val,  dataorig[n].alb_dv.perr,  dataorig[n].alb_dv.merr);
  data[n].rp_dv.val    = AGRAND(dataorig[n].rp_dv.val,   dataorig[n].rp_dv.perr,   dataorig[n].rp_dv.merr);
  data[n].alb_alt.val  = AGRAND(dataorig[n].alb_alt.val, dataorig[n].alb_alt.perr, dataorig[n].alb_alt.merr);
  data[n].rp_alt.val   = AGRAND(dataorig[n].rp_alt.val,  dataorig[n].rp_alt.perr,  dataorig[n].rp_alt.merr);

  // Randomly perturb the centroid disposition, according to the centroid score value
  if(data[n].cent_score > dist(randgen))  // dist(randgen) is a random number between 0 and 1.0
    data[n].robo_cent_disp = 0;  // If the CRV score is greater than a random number 0->1, then it's a pass. E.g., a score of 0.9 will greater than a random number between 0.0 and 1.0 90% of the time.
  else
    data[n].robo_cent_disp = 1;  // If the CRV score is less than a random number 0->1, it's a fail. E.g., a score of 0.05 will be less than random number between 0.0 and 1.0 95% of the time.
  }
