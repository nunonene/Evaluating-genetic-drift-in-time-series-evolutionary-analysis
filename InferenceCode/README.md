Inference code for "Evaluating-genetic-drift-in-time-series-evolutionary-analysis"
=================================================================================

Forward-backward/predict-update HMM algorithm for drift model parameter estimation (N or sigma) through likelihood maximization and posterior calculation. Adapted for time-series analysis from Andrej Fischer, Ignacio Vazquez-Garcia, Christopher J.R. Illingworth and Ville Mustonen. High-definition reconstruction of subclonal composition in cancer. Cell Reports (2014), http://dx.doi.org/10.1016/j.celrep.2014.04.055

Additional features for time-resolved data: neutral exact Wright-Fisher (drift parameter N) and Gaussian (drift paramter sigma) propagation with absorbing boundaries drift model on a frequency grid with matrix exponentiation routine.

Executable: DMS
-----------

Required external C libraries:
------------------------------

--gsl (see makefile) 

Required Arguments:
--------------------

  --data    [file]   The input observed datafile with four columns (extra columns of "read depth" correspond to more replicates which                     are analysed independently from each other), locus time reads depth

  --mode    [int]    1 Binomial 2 Beta-Binomial

  --GWF     [int]    0 Gaussian propagation with absorbing boundaries 1 Wright-Fisher propagation by matrix exponentiation (2N<=1000)                     or pre-computed matrix powers (2N>1000)
  		
Optional Arguments:
-------------------

  --pre        [str]     The prefix to put before all output files
  
  --grid       [int]     The grid size for the distributions (partitions [0,1] into [grid] bins)
  
  --dist                 Prints all the posterior distributions as well
  
  --sigmai     [double]  Initial drift amplitude (default minimum value sigmai=0.01)
  
  --sigmaf     [double]  Final drift amplitude (default minimum value sigmaf=0.05)
  
  --sigmagrid  [int]     Number of initial points for drift parameter optimization
  
  --shape      [double]  To fix the shape parameter in mode 2

Output:
-------

  -- Parameter estimate, corresponding likelihood and goodness-of-fit calculated through the posterior distribution at each locus and sampling instant
  
  -- The posterior mean and standard deviation of each locus frequency at each sampling instant
  
  -- On request (--dist), the posterior distribution at each point (may generate large files and be slow to compute!)
  


Example:
--------

./DMS --data [file] --mode 1 --GWF 1 -pre [str] --sigmai 0.01 --sigmaf 0.05 --sigmagrid 5 --grid 400 

Additional folders:

Examples-examples of Wright-Fisher trajectories with duration T=300, sampling period 10, population size N=1000

MatPow- pre-computed transition matrices


Notes:
------

--data files should have line ending in Unix/Linux format

--Each locus is treated as an independent sample from the same process (linkage is not taken into account). Later versions will maximize likelihood across all replicates taken together. In order to achieve this in the current version just add the additional replicates to the bottom of the 4 column file 

--For the Wright-Fisher model if during the optimization the current population size (N) is above 1000, the routine reads precomputed matrix powers from folder MatPow located in the same directory where the executable is.
