#ifndef NPT_H
#define NPT_H

void hs_npt();

void run_npt(bool prod_flag, int sweep_offset,
	     int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
	     int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void sweep_npt(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
               int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

#endif
