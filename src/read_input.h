#ifndef READ_INPUT_H
#define READ_INPUT_H

struct input {

  // Density
  double rho;

  // Number of cells
  int nx, ny, nz, type;

  // Maximum displacement
  double dr_max;

  // Number of sweeps 
  int sweep_eq, sweep_stat;

  // Output on screen
  int output_int;

  // NpT simulations
  double press, dv_max;

  // Optimizer
  int opt_flag, opt_samples, opt_sweeps;
  double opt_part_target, opt_vol_target;

  // Seed for random number generator
  unsigned long int seed;

  // Compute: pressure (virial and thermodynamic)
  double pressv_dr, presst_xi_max, presst_dxi;
  int pressv_sample_int, presst_sample_int;
  
  // Compute: order parameter
  int ql_order, ql_sample_int;
  double ql_rmax;

};

extern struct input in;

void print_example();

void  read_input_file(char *filename);

void read_input_file_err(int err_id, char *last_string);

#endif
