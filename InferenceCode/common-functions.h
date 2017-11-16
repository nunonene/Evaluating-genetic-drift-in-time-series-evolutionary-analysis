//common-functions.h

// Adapted from http://dx.doi.org/10.1016/j.celrep.2014.04.055
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


//own headers

class Emission;

using namespace std;

void get_dims( const char * data_fn, int& nRep, vector<int>& Loci, vector<int>& sTimes, int keep);
void get_data( const char * data_fn, Emission * myEmit);

double get_mean( gsl_vector * dist, double xmin, double xmax);
double get_var(  gsl_vector * dist, double xmin, double xmax, double mean);



