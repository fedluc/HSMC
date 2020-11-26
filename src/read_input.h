#ifndef READ_INPUT_H
#define READ_INPUT_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

struct input {
  double rho;
  int nx, ny, nz, type;
  double dr_max;
  int sweep_eq, sweep_stat;
  int dt_output;
  double rdf_dr;
  int rdf_tcompute, rdf_tave;
};

extern struct input in;

void print_example();

void  read_input_file(char *filename);

void read_input_file_err(int err_id, char *last_string);

#endif
