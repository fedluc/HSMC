#ifndef COMPUTE_CHEM_POT_H
#define COMPUTE_CHEM_POT_H

#include <stdbool.h>

void compute_mu(bool init);

void widom_insertion();

void widom_rand_pos();

bool widom_check_overlap();

int widom_cell_part_idx();

double widom_compute_dist(int idx);

void mu_output(bool init);

#endif
