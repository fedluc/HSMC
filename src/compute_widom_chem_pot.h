#ifndef COMPUTE_CHEM_POT_H
#define COMPUTE_CHEM_POT_H

#include <stdbool.h>

void compute_mu(bool init,
		int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
                int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void widom_insertion(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		     int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void widom_rand_pos();

bool widom_check_overlap(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
			 int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

int widom_cell_part_idx();

double widom_compute_dist(int idx);

void mu_output(bool init);

#endif
