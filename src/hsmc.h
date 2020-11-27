#ifndef HSMC_H
#define HSMC_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_fit.h>
#include <time.h>
#include "cell_list.h"

void hs_nvt();

void run_opt();

void run();

void sweep();

double compute_dist(int idx1, int idx2);

bool check_overlap();

void compute_hist(double *hist, int *counter,
                  double *pos, double cutoff);

void average_hist(double *hist, int nn);

void reset_hist(double *hist, int *counter, int nn);

void compute_rdf(double *rdf, double *hist, double *rr, int nn);

void compute_pressure(double *press, double *rdf, double *rr, int nn);


#endif
