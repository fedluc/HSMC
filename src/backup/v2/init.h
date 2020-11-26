#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
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


struct box_info sim_box_init(int cell_type, int nx, int ny,
                             int nz, double rho);

struct p_info part_alloc(double ***part,
                         struct box_info sim_box_info);

void part_init(double **part, struct box_info sim_box_info,
               struct p_info part_info);

void part_init_sc(double **part, struct box_info sim_box_info);

void part_init_fcc(double **part, struct box_info sim_box_info);

void part_init_err();

void add_particle(double **part, int id,
                  double xx, double yy, double zz);
