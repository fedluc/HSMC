#ifndef READ_INPUT_H
#define READ_INPUT_H

struct input {
  double rho;
  int nx, ny, nz, type;
  double dr_max;
  int sweep_eq, sweep_stat;
  int output_int;
  double pressv_dr, presst_xi_max, presst_dxi;
  int pressv_sample_int, presst_sample_int;
};

extern struct input in;

void print_example();

void  read_input_file(char *filename);

void read_input_file_err(int err_id, char *last_string);

#endif
