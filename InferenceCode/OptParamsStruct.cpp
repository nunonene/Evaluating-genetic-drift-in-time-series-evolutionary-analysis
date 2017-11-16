/*
  OptParamsStruct.cpp
 
  Nuno Rocha Nene   

*/



#include "OptParamsStruct.h"



#include "emission.h"
#include "drift.h"

// Optimization struct

struct fpar{ // structure of parameters passed into the gsl optimizer
  Drift * myD;
  vector<int> to_opt;
  double sigm_i;
  double sigm_f;
  int sigmgrid;
};
