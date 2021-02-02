#ifndef READ_INPUT_H
#define READ_INPUT_H

struct input {

  // Density
  double rho;

  // Number of cells
  int nx, ny, nz, type;

  // Parameters for neighbor list
  double neigh_dr;
  int neigh_max_part;

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

  // Cavity simulations
  double cavity_pcav, cavity_maxdr, cavity_mindr;
  double cavity_out_dr;
  int cavity_sample_int;

  // Cluster moves
  int cluster_flag, cluster_moves_sweep;

  // Read restart
  int restart_read;
  char restart_name[100];

  // Write restart
  int restart_write;

  // Write configuration
  int config_write, config_samples;
  
  // Compute: pressure (virial and thermodynamic)
  double pressv_dr, presst_xi_max, presst_dxi;
  int pressv_sample_int, presst_sample_int;
  
  // Compute: order parameter
  int ql_order, ql_sample_int;
  double ql_rmax;

  // Compute: chemical potential
  int mu_insertions, mu_sample_int;

  // Compute: radial distribution function
  int rdf_sample_int, rdf_samples;
  double rdf_dr, rdf_rmax;

};

extern struct input G_IN;

void print_example();

void  read_input_file(char *filename);

void read_input_file_err(int err_id, char *last_string);

#endif
