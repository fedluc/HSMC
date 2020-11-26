#ifndef INIT_H
#define INIT_H

#include "read_input.h"
#include <math.h>

struct p_info {int NN, Ncell; };

struct box_info {
  double vol;
  double lx, ly, lz;
  double min_size;
  double cell_size;
  int cell_x, cell_y, cell_z;
  int cell_type;
};

extern struct box_info sim_box_info;
extern struct p_info part_info;

void sim_box_init(int cell_type, int nx, int ny,
		  int nz, double rho);

void part_alloc(double ***part);

void part_init(double **part);

void part_init_sc(double **part);

void part_init_fcc(double **part);

void part_init_err();

void add_particle(double **part, int id,
                  double xx, double yy, double zz);


#endif
