#ifndef COMPUTE_ORDER_PARAMETER_H
#define COMPUTE_ORDER_PARAMETER_H

#include <stdbool.h>

void compute_op(bool init,
		int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void global_ql_compute(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		       int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void ql_compute(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void qlm2_compute(int ref_idx,
		  int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		  int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void global_ql_output(bool init);

#endif
