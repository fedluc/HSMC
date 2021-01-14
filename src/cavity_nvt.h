#ifndef CAVITY_NVT_H
#define CAVITY_NVT_H

void cavity_hs_nvt();

void cavity_run_nvt(bool prod_flag, int sweep_offset,
		    int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
                    int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void cavity_set_distance();

void cavity_sweep_nvt(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		      int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void cavity_psi_output();

void cavity_dist_output(bool init);

#endif
