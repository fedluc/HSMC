#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_input.h"

// ----------------------------------------
// ------------- Read Input ---------------
// ----------------------------------------


void print_example(){

  char example[] = 
    
    "# This is an input file for hsmc\n"
    "# Lines starting with # or empty lines are skipped\n\n"
    "# Density\n"
    "# (for npt calculations the density is only used to specify\n" 
    "# the starting configuration)\n"
    "rho 0.8\n\n"
    "# Number of cells along the x, y and z directions\n"
    "# The cells are used a building block of the simulation box\n"
    "cells_x 5\n" 
    "cells_y 5\n"
    "cells_z 5\n\n"
    "# Type of cells\n"
    "# 1: simple cubic, 1 particle per cell\n"
    "# 2: face-centered cubic, 4 particles per cell\n"
    "type 2\n\n"
    "# Neighbor list (minumum size of cell list, default 1.0)\n"
    "neigh_list 1.05\n\n"
    "# Maximum displacement for particles moves\n"
    "dr_max 0.05\n\n"
    "# Parameters for NpT calculations (pressure, maximum volume deformation\n"
    "# the default is npt 0, which corresponds to an NVT simulations)\n"
    "npt 2.5 0.001\n\n"
    "# Optmization\n"
    "# (activation flag (0 or 1), sweeps used during optimization\n"
    "# samples collected (<= sweeps), target acceptance ratio for particles moves\n"
    "# target acceptance ratio for volume deformations, default: 1 1000 10 0.5 0.5)\n"
    "opt 1 1000 10 0.5 0.5\n\n"
    "# Pressure via virial (resolution, saving interval\n"
    "# note: available only for NVT simulations)\n"
    "press_virial 0.002 10 \n\n"
    "# Pressure via thermodynamics (resolution, max relative compression, saving interval)\n"
    "press_thermo 0.0001 0.002 10 \n\n"
    "# Order parameter (order, cutoff distance, saving interval)\n"
    "ql 6 1.5 1 \n\n"
    "# Chemical potential via Widom insertions (insertions, saving interval)\n"
    "widom 100 5 \n\n"
    "# Cavity simulations (probability of moving a cavity, maximum and minimum distance,\n"
    "# saving interval, potential resolution)\n"
    "cavity 0.1 1.2 0.0 4 0.01\n\n"
    "# Seed for random number generator\n"
    "seed 124787\n\n"
    "# Write restart data (saving interval)\n"
    "restart_write 1024\n\n"
    "# Read restart data (activation flag (0 or 1), name with restart file\n"
    "# restart_read 1 restart_04096.bin\n\n"
    "# Write configuration to file (saving interval)\n"
    "config_write 1024\n\n"
    "# Number of sweeps for equilibration (for N particles, 1 sweep = N moves)\n"
    "sweep_eq 10000\n\n"
    "# Number of sweeps for statistics\n"
    "sweep_stat 10000\n\n"
    "# Output interval (how often to print status on screen)\n"
    "out 100\n\n";

   
    printf(example);

}

struct input in;

void read_input_file(char *filename){

  FILE *in_file;
  char *line_buf = NULL;
  char *key, *value;
  size_t line_buf_size = 0;
  ssize_t line_size;

  // Initialize input structure 
  in.rho = 0;
  in.nx = 0;
  in.ny = 0;
  in.nz = 0;
  in.type = 0;
  in.neigh_dr = 1.0;
  in.neigh_max_part = 10;
  in.dr_max = 0;
  in.sweep_eq = 0;
  in.sweep_stat = 0;
  in.output_int = 0;
  in.pressv_dr = 0;
  in.pressv_sample_int = 0;
  in.presst_dxi = 0;
  in.presst_xi_max = 0;
  in.presst_sample_int = 0;
  in.press = 0;
  in.dv_max = 0;
  in.opt_flag = 1;
  in.opt_sweeps = 1000;
  in.opt_samples = 10;
  in.opt_part_target = 0.5;
  in.opt_vol_target = 0.5;
  in.seed = 0;
  in.ql_order = -1;
  in.ql_rmax = 0;
  in.ql_sample_int = 0;
  in.mu_sample_int = 0;
  in.mu_insertions = 0;
  in.cavity_pcav = 0;
  in.cavity_maxdr = 0;
  in.cavity_mindr = 0;
  in.cavity_sample_int = 0;
  in.restart_read = 0;
  in.restart_write = 0;
  in.config_write = 0;

  // Open file
  printf("Reading input data from %s ...\n",filename);
  in_file = fopen(filename,"r");
    
  // Test if the file exists 
  if (!in_file) {  
      printf("Error! Could not open file %s\n", filename); 
      exit(EXIT_FAILURE);
  } 

  // Get first line of file
  line_size = getline(&line_buf, &line_buf_size, in_file);

  // Loop through file until EOF
  while (line_size >= 0){

    if (line_buf[0] != '#' && line_buf[0] != '\n'){
      
      // Key-word
      key = strtok(line_buf, " ");
      
      if (strcmp(key,"rho") == 0 || strcmp(key,"rho\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.rho = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"cells_x") == 0 || strcmp(key,"cells_x\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.nx = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }


      else if (strcmp(key,"cells_y") == 0 || strcmp(key,"cells_y\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.ny = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }


      else if (strcmp(key,"cells_z") == 0 || strcmp(key,"cells_z\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.nz = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }


      else if (strcmp(key,"type") == 0 || strcmp(key,"type\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.type = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"neigh_list") == 0 || strcmp(key,"neigh_list\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.neigh_dr = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"dr_max") == 0 || strcmp(key,"dr_max\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.dr_max = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"sweep_eq") == 0 || strcmp(key,"sweep_eq\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.sweep_eq = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"sweep_stat") == 0 || strcmp(key,"sweep_stat\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.sweep_stat = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"out") == 0 || strcmp(key,"out\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.output_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"press_virial") == 0 || strcmp(key,"press_virial\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.pressv_dr = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.pressv_sample_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"press_thermo") == 0 || strcmp(key,"press_thermo\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.presst_dxi = atof(value);
	}
	else read_input_file_err(1,line_buf);
        value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.presst_xi_max = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.presst_sample_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"npt") == 0 || strcmp(key,"npt\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.press = atof(value);
	}
	else read_input_file_err(1,line_buf);
        value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.dv_max = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"opt") == 0 || strcmp(key,"opt\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.opt_flag = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.opt_sweeps = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.opt_samples = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.opt_part_target = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.opt_vol_target = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"seed") == 0 || strcmp(key,"seed\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  in.seed = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"ql") == 0 || strcmp(key,"ql\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.ql_order = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.ql_rmax = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.ql_sample_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"widom") == 0 || strcmp(key,"widom\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.mu_insertions = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.mu_sample_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"cavity") == 0 || strcmp(key,"cavity\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.cavity_pcav = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.cavity_maxdr = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.cavity_mindr = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.cavity_sample_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.cavity_out_dr = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"restart_read") == 0 || strcmp(key,"restart_read\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.restart_read = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  int len = strlen(value);
	  if (value[len-1] == '\n') value[len-1] = '\0';
	  if (len < 100) strcpy(in.restart_name,value);
	  else read_input_file_err(3,line_buf);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"restart_write") == 0 || strcmp(key,"restart_write\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.restart_write = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"config_write") == 0 || strcmp(key,"config_write\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  in.config_write = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }


      else read_input_file_err(2,line_buf);

    }

    line_size = getline(&line_buf, &line_buf_size, in_file);
  }

  // Free memory
  free(line_buf);
  line_buf = NULL;

  // Close file
  fclose(in_file);

  // Print content of input structure
  /* printf("Density: %.8f\n", in.rho); */
  /* printf("Number of cells along x: %d\n", in.nx); */
  /* printf("Number of cells along y: %d\n", in.ny); */
  /* printf("Number of cells along z: %d\n", in.nz); */
  /* printf("Cell type: %d\n", in.type); */
  /* printf("Maximum MC displacement: %.8f\n", in.dr_max); */
  /* printf("Number of sweeps (equilibration): %d\n", in.sweep_eq); */
  /* printf("Number of sweeps (statistics): %d\n", in.sweep_stat); */
  /* printf("Output interval (sweeps): %d\n", in.output_int); */
  /* printf("Pressure - virial (resolution): %.8f\n", in.pressv_dr); */
  /* printf("Pressure - virial (sweeps/sample): %d\n", in.pressv_sample_int); */
  /* printf("Pressure - thermo (resolution): %.8f\n", in.presst_dxi); */
  /* printf("Pressure - thermo (max compression): %.8f\n", in.presst_xi_max); */
  /* printf("Pressure - thermo (sweeps/sample): %d\n", in.presst_sample_int); */
  fflush(stdout);

}

void read_input_file_err(int err_id, char *last_string){

  switch (err_id)
    {

    case 1: 
      printf("Missing value to key\n");
      break;
    
    case 2:
      printf("Unknown key\n");
      break;

    case 3:
      printf("Name of restart file is too long, maximum 100 characters\n");
      break;

    }

  printf("Last read line in the input file:\n%s\n", last_string);
  exit(EXIT_FAILURE);

}
