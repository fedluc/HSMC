#ifndef COMPUTE_PRESS_H
#define COMPUTE_PRESS_H

#include <stdbool.h>

void compute_pressv(bool init,
		    int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		    int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void pressv_hist_init();

void pressv_compute_hist(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
			 int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void pressv_compute_rdf();

void pressv_output(bool init);

void compute_presst(bool init,
		    int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		    int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void presst_hist_init();

void presst_compute_hist(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
			 int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void presst_output(bool init);

#endif
