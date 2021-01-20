#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "read_input.h"

// --------------------------------------------------------
// The module "rng.c" is used to initialize, update and
// eventually free the random number generator. The default
// choice for the random number generator is the Marsenne-
// Twister implemented in the GNU Scientific Library and 
// described at the following link (Accessed January 20, 
// 2021): https://www.gnu.org/software/gsl/doc/html/rng.html
// --------------------------------------------------------


// Global variable for random number generator
static gsl_rng *rng_mt;
static long unsigned int r_num_max;

// ------ Initialize random number generator ------
void rng_init(){
  // Set-up random number generator (Marsenne-Twister)
  rng_mt = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng_mt,G_IN.seed);
  r_num_max = gsl_rng_max(rng_mt);
}

// ------ Get one random number in [0,1] ------
double rng_get_double(){
  return (double)gsl_rng_get(rng_mt)/(double)r_num_max;
}

// ------ Get one random integer in [0,MM-1] ------
int rng_get_int(int MM){
  return (int)gsl_rng_uniform_int(rng_mt, MM);
}

// ------ Free random number generator ------
void rng_free(){
  gsl_rng_free(rng_mt);
}

// ------ Write status of random number generator to file ------
void rng_write(FILE *fid){
  gsl_rng_fwrite(fid, rng_mt);
}

// ------ Read status of random number generator to file ------
void rng_read(FILE *fid){
  rng_mt = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_fread(fid, rng_mt);
  r_num_max = gsl_rng_max(rng_mt);
}
