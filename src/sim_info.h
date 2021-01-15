#ifndef SIM_INFO_H
#define SIM_INFO_H

#include <stdio.h>

struct p_info {int NN, Ncell; };

struct box_info {
  double vol;
  double lx, ly, lz;
  double min_size;
  double cell_size;
  int cell_x, cell_y, cell_z;
  int cell_type;
};

typedef double part_array[][4];

extern double (*part)[4]; // THIS SHOULD BE REMOVED

void sim_box_init(int cell_type, int nx, int ny,
		  int nz, double rho);

void part_alloc();

void part_init();

void part_init_sc();

void part_init_fcc();

void part_init_err();

void add_particle(int id, double xx, double yy, double zz);

struct box_info sim_box_info_get();

struct p_info part_info_get();

//double (*part_config_get())[4];
part_array *part_config_get();

void part_free();

void print_sim_info();

void sim_box_info_write(FILE *fid);

void sim_box_info_read(FILE *fid);

void part_info_write(FILE *fid);

void part_info_read(FILE *fid);


#endif
