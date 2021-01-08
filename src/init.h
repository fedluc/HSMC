#ifndef INIT_H
#define INIT_H

#include <gsl/gsl_rng.h>

struct p_info {int NN, Ncell; };

struct box_info {
  double vol;
  double lx, ly, lz;
  double min_size;
  double cell_size;
  int cell_x, cell_y, cell_z;
  int cell_type;
};

extern gsl_rng *rng_mt;
extern long unsigned int r_num_max;
extern struct box_info sim_box_info;
extern struct p_info part_info;
extern double (*part)[4];

void sim_box_init(int cell_type, int nx, int ny,
		  int nz, double rho);

void part_alloc();

void part_init();

void part_init_sc();

void part_init_fcc();

void part_init_err();

void add_particle(int id, double xx, double yy, double zz);


#endif
