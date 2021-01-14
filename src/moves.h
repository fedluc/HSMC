#ifndef MOVES_H
#define MOVES_H

#include <stdbool.h>


void part_move();


void vol_move();

bool check_overlap(int idx_ref,
                   double sf_x, double sf_y, double sf_z);

double compute_dist(int idx1, int idx2,
                    double sf_x, double sf_y, double sf_z);

void cavity_part_move(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		      int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

bool cavity_check_move(int idx_ref, int move_type, double en_old,
		       int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		       int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]);

double cavity_interaction(double xx, bool cavity_init);

void get_moves_counters(int *pm, int *apm, int *rpm,
                        int* vm, int *avm, int *rvm);

void reset_moves_counters();

#endif
