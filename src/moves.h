#ifndef MOVES_H
#define MOVES_H

#include <stdbool.h>
#include <gsl/gsl_rng.h>

extern gsl_rng *rng_mt;
extern long unsigned int r_num_max;
extern int part_moves, vol_moves;
extern int acc_part_moves, rej_part_moves;
extern int acc_vol_moves, rej_vol_moves;

void rng_init();

void part_move();

void vol_move();

double compute_dist(int idx1, int idx2,
                    double sf_x, double sf_y, double sf_z);

bool check_overlap(int idx_ref,
                   double sf_x, double sf_y, double sf_z);

void cavity_rng_init();

void cavity_part_move();

bool cavity_check_move(int idx_ref, int move_type, double en_old);

double cavity_interaction(double xx);

double hspy(double xx, double eta);

#endif
