#ifndef OPTIMIZER_H
#define OPTIMIZER_H

void opt_nvt(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
	     int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void get_sample_nvt(double *dr, double *acc_ratio, 
		    int sample_iter,
		    int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		    int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void opt_npt(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
	     int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void get_sample_npt(double *dr, double *acc_ratio_part,
                    double *dv, double *acc_ratio_vol,
                    int sample_iter,
		    int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		    int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);


void opt_cavity_nvt(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		    int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

void get_sample_cavity_nvt(double *dr, double *acc_ratio, 
			   int sample_iter,
			   int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
			   int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);


#endif
