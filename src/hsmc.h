#ifndef HSMC_H
#define HSMC_H

#include <stdbool.h>

void hs_nvt();

void run_opt();

void run(bool prod_flag);

void sweep();

double compute_dist(int idx1, int idx2, 
		    double sf_x, double sf_y, double sf_z);

bool check_overlap(int idx_ref, 
		   double sf_x, double sf_y, double sf_z);

#endif
