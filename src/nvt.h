#ifndef NVT_H
#define NVT_H

void hs_nvt();

void run_nvt(bool prod_flag, int sweep_offset, 
	     int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
	     int cl_neigh_num, int cl_neigh[cl_num_tot][27]);

void sweep_nvt(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
	       int cl_neigh_num, int cl_neigh[cl_num_tot][27]);


#endif
