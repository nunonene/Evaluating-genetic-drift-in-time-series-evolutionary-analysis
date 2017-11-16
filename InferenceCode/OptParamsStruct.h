/*
  OptParamsStruct.h
 
  Nuno Rocha Nene   

*/



#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <vector>
#include <list>


// GSL headers...

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_blas.h"




//#include "emission.h"
//#include "drift.h"

// Optimization struct

struct fpar{ // structure of parameters passed into the gsl optimizer
  Drift * myD;
  vector<int> to_opt;
  double sigm_i;
  double sigm_f;
  int sigmgrid;
};
