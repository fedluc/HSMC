#ifndef MOVES_H
#define MOVES_H

#include <stdbool.h>


void part_move();


void vol_move();

bool check_overlap(int idx_ref,
                   double sf_x, double sf_y, double sf_z);

double compute_dist(int idx1, int idx2,
                    double sf_x, double sf_y, double sf_z);

void apply_pbc(int idx);

void cavity_part_move();

bool cavity_check_move(int idx_ref, int move_type, double en_old);

double cavity_interaction(double xx, bool cavity_init);

void get_moves_counters(int *pm, int *apm, int *rpm,
                        int* vm, int *avm, int *rvm);

void reset_moves_counters();

#endif
