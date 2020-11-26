#ifndef HSMC_H
#define HSMC_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_fit.h>
#include "read_input.h"
#include "init.h"

void hs_nvt();

double run_opt(double **part, const gsl_rng *r);

void run(double **part, const gsl_rng *r);


void sweep(double **part, int *acc_moves, int *rej_moves,
           const gsl_rng *r, 
	   long unsigned int max_r);

void metropolis(double **part, int *accept, int *reject,
                int r_idx, double r_x, double r_y, double r_z);

bool check_overlap(double **part, int target_idx, 
		   double r_x, double r_y, double r_z);

void compute_hist(double *hist, int *counter, double **part,
                  double *pos, double cutoff);

void average_hist(double *hist, int nn);

void reset_hist(double *hist, int *counter, int nn);

void compute_rdf(double *rdf, double *hist, double *rr, int nn);

void compute_pressure(double *press, double *rdf, double *rr, int nn);


#endif
