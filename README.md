# PeakCompPoisson

Peak comparison program using Poisson ML ratio chi-square

Maintainers: Aaron Chester, Jonathan Williams

## Description

Uses Poisson maximum likelihood ratio chi-square for the comparison of experimental to Geant4-simulated RDM lineshapes for TIP experiment analysis. Utilizes the ROOT Minuit function minimization libraries for numerical minimization of chi-square statistic w.r.t. model parameters. Particularly useful for low statistics data sets.

### Useful References: 

* A. Chester et al. Nucl. Inst. and Meth. A 882 (2017) 69-83.
* A. Chester et al. Phys. Rev. C. 96 (2017) 011302.
* S. Baker and R. D. Cousins. Nucl. Inst. and Meth. 221 (1984) 437-442.
* G. Cowan. Statistical Data Analysis. Clarendon Press, 1998.
* F. James. Minuit Function Minimization and Error Analysis Reference Manual. CERN 3 1994. Version 94.1. CERN Program Library entry D506.

## Installation

Requires ROOT (tested with v5.x) to be installed with environment variables set up properly.  Environment variable setup can be done by adding to your `.bashrc` (and then reloading the terminal):

```
#ROOT configuration in .bashrc
export ROOTSYS=/path/to/root
export ROOTINC=$ROOTSYS/include
export ROOTLIB=$ROOTSYS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTLIB
```

Once this is done, use `make` to compile.  Tested using g++ and GNU make on Ubuntu 16.04.
