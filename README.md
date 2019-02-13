# The Kepler DR25 Robovetter

The DR25 Kepler Robovetter is a robotic decision-making code that dispositions each Threshold Crossing Event (TCE) from the Kepler pipeline into Planet Candidates (PCs) and False Positives (FPs). The Robovetter also provides four major flags to designate each FP TCE as Not Transit-Like (NTL), a Stellar Eclipse (SS), a Centroid Offset (CO), and/or an Ephemeris Match (EM). It also produces a score ranging from 0.0 to 1.0 that indicates the Robovetter's disposition confidence, where 1.0 indicates strong confidence in PC, and 0.0 indicates strong confidence in FP. Finally, the Robovetter provides comments in a text string that indicate the specific tests each FP TCE fails and provides supplemental information on all TCEs as necessary.

More information can be found at: https://exoplanetarchive.ipac.caltech.edu/docs/PurposeOfKOITable.html#q1-q17dr25


## Compiling and Running the Code

The DR25 Robovetter code is provided, along with the necessary input files that contain every metric the Robovetter uses to make its decisions, and the resultant output files, so that users may validate their implementation. (Note that input files are compressed using tar/gzip and and should be uncompressed via "tar xvzf FILENAME" before being used.)


### Prerequisites

The code is written in C++11 and only requires the standard C++ library (specifically, the required libraries are iomainip, iostream, fstream, vector, and random). It has been tested on Linux and Mac to work with:
  - The g++ complier (Minimum version 4.7.2 tested - earlier versions unlikely to work.)
  - The clang++ compiler (Versions 3.4 and 3.5 tested. Version 3.3 may work, but untested. Earlier than 3.3 will not work.)
  - The Intel icpc compiler (Version 17 tested. Earlier versions as far back as 11 very likely to work, but untested.)

Note that Danley Hsu was able to compile and run the code on Windows, and he has graciously provided detailed instructions in the file https://github.com/nasa/kepler-robovetter/blob/master/windows-compile-instructions. These instructions have not been verified to work on another Windows machine, but we provide them in case anyone using a Windows machine finds them helpful.


### Compiling

To compile the code, use your available C++ compiler with the -std=c++11 option, and a recommended O2 level of optimization. For example, use one of the commands below based on your available compiler:

```
g++     -std=c++11 -O2 -o robovet Kepler-RoboVetter-DR25.cpp
```
```
clang++ -std=c++11 -O2 -o robovet Kepler-RoboVetter-DR25.cpp
```
```
icpc    -std=c++11 -O2 -o robovet Kepler-RoboVetter-DR25.cpp
```

### Running the Code

Run as "./robovet INPUTFILE  OUTFILE  NMC"

INPUTFILE is the name of the input file for a particular run. OBS is the Observed run. INV is the Inverted run. SCR1 is the Scrambled (ordering #1) run. SCR2 is the Scrambled (ordering #2) run. SCR3 is the Scrambled (ordering #3) run. INJ1 is the Injected group 1 run (on-target injected planets.) INJ2 is the Injected group 2 run (off-target planets.) INJ3 is the Injected group 3 run (Eclipsing Binaries). For more details on each run, see https://exoplanetarchive.ipac.caltech.edu/docs/KSCI-19114-001.pdf

OUTFILE is the name of the output file, as desired by the user.

NMC is the number of Monte Carlo iterations desired for computing the scores. A value of at least 100 is recommended for useable scores, and preferably at least 1,000. A value of 10,000 was used for calculating the scores archived at NExScI and the output files provided here.

For example:

```
./robovet kplr_dr25_obs_robovetter_input.txt kplr_dr25_obs_robovetter_output.txt 10000
```

will run the Robovetter on the OBS data, and perform 10,000 Monte Carlo runs to compute the score. The resulting output file should exactly match that provided in this GitHub repository. To replicate the other data sets (INV, SCR1, SCR2, SCR3, INJ1, INJ2, INJ3) one would use the following commands:

```
./robovet kplr_dr25_inv_robovetter_input.txt kplr_dr25_inv_robovetter_output.txt 10000
```
```
./robovet kplr_dr25_scr1_robovetter_input.txt kplr_dr25_scr1_robovetter_output.txt 10000
```
```
./robovet kplr_dr25_scr2_robovetter_input.txt kplr_dr25_scr2_robovetter_output.txt 10000
```
```
./robovet kplr_dr25_scr3_robovetter_input.txt kplr_dr25_scr3_robovetter_output.txt 10000
```
```
./robovet kplr_dr25_inj1_robovetter_input.txt kplr_dr25_inj1_robovetter_output.txt 10000
```
```
./robovet kplr_dr25_inj2_robovetter_input.txt kplr_dr25_inj2_robovetter_output.txt 10000
```
```
./robovet kplr_dr25_inj3_robovetter_input.txt kplr_dr25_inj3_robovetter_output.txt 10000
```


### Terminal output

In addition to the output file that is created, the code writes to the command line the currently executing task (e.g., reading in the data or vetting the TCEs) and the current Monte-Carlo iteration number so that users can monitor progress.


## Acknowledgments

Please reference Thompson et al. 2018, ApJS, 235, 38 (http://adsabs.harvard.edu/abs/2018ApJS..235...38T) if you make use of this code or the files provided.


## Notices

Copyright © 2017 United States Government as represented by the Administrator of the National Aeronautics and Space Administration. All Rights Reserved.

NASA acknowledges the SETI Institute’s primary role in authoring and producing the Kepler Robovetter under Cooperative Agreement Number NNX13AD01A


## Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE. FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT. IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW. RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT.
