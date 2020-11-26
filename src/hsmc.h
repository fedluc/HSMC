#ifndef HSMC_H
#define HSMC_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_fit.h>
#include "read_input.h"
#include "init.h"

void hs_nvt();

void run_opt();

void run();

void sweep();

bool check_overlap();

void compute_hist(double *hist, int *counter,
                  double *pos, double cutoff);

void average_hist(double *hist, int nn);

void reset_hist(double *hist, int *counter, int nn);

void compute_rdf(double *rdf, double *hist, double *rr, int nn);

void compute_pressure(double *press, double *rdf, double *rr, int nn);


#endif
