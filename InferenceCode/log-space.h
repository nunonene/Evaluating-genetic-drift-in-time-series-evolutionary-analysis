//log_space.h

// Adapted from http://dx.doi.org/10.1016/j.celrep.2014.04.055

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <vector>
#include <algorithm>


// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"

double log_add(double one, double two);
double log_sub(double one, double two);
void log_vector_add(gsl_vector * one, gsl_vector * two);
void log_matrix_add(gsl_matrix * one, gsl_matrix * two);
void log_vector_invert(gsl_vector * vec);
double log_vector_norm(const gsl_vector * x);
double log_matrix_norm(const gsl_matrix * M);
void log_vector_normalize(gsl_vector * x);
