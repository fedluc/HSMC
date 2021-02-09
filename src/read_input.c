#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_input.h"

// --------------------------------------------------------
// The module "read_input.c" is used to read the input file
// --------------------------------------------------------

// ------ Global variable ------
// WARNING: This is a truly global variable in the sense 
// that its scope is not limited to this module
struct input G_IN; 


// ------ Print example of input file on screen ------

void print_example(){

  char example[] = 
    
    "# This is an input file to compute the pressure via the virial route in the   \n"
    "# canonical (NVT) ensemble. The simulation is performed at density 0.5/sigma^3,\n"
    "# with 1000 particles 1e6 sweeps for equilibration and 1e6 sweeps for statistics\n"
    "# Configuration and restart files are written every 1e5 sweeps, output on screen\n"
    "# is printed every 1e4 sweeps\n\n"
    "# Density\n"
    "rho 0.5\n\n"
    "# Simulation box\n"
    "cells_x 10\n" 
    "cells_y 10\n"
    "cells_z 10\n"
    "type 1\n\n"
    "# Neighbor list\n"
    "neigh_list 1.0 10\n\n"
    "# Maximum particle displacement\n"
    "dr_max 0.05\n\n"
    "# Optmizer\n"
    "opt 1 1000 10 0.5 0.5\n\n"
    "# Pressure via virial\n"
    "press_virial 0.002 10 \n\n"
    "# Seed for random number generator\n"
    "seed 124787\n\n"
    "# Write restart data\n"
    "restart_write 100000\n\n"
    "# Write configuration to file\n"
    "config_write 100000 100\n\n"
    "# Number of sweeps for equilibration\n"
    "sweep_eq 1000000\n\n"
    "# Number of sweeps for statistics\n"
    "sweep_stat 1000000\n\n"
    "# Output interval\n"
    "out 10000\n\n";

   
    printf(example);

}

// ------ Read input file ------

void read_input_file(char *filename){

  FILE *in_file;
  char *line_buf = NULL;
  char *key, *value;
  size_t line_buf_size = 0;
  ssize_t line_size;

  // Initialize input structure 
  G_IN.rho = 0.5;
  G_IN.nx = 3;
  G_IN.ny = 3;
  G_IN.nz = 3;
  G_IN.type = 1;
  G_IN.neigh_dr = 1.0;
  G_IN.neigh_max_part = 10;
  G_IN.dr_max = 0.05;
  G_IN.sweep_eq = 0;
  G_IN.sweep_stat = 0;
  G_IN.output_int = 0;
  G_IN.cavity_pcav = 0;
  G_IN.cavity_maxdr = 1.2;
  G_IN.cavity_mindr = 0;
  G_IN.cavity_out_dr = 0.01;
  G_IN.cavity_sample_int = 100;
  G_IN.cluster_flag = 0;
  G_IN.cluster_moves_sweep = 1;
  G_IN.cluster_init_step = 10000;
  G_IN.restart_read = 0;
  G_IN.restart_write = 0;
  G_IN.config_write = 0;
  G_IN.config_samples = 128;
  G_IN.pressv_dr = 0.01;
  G_IN.pressv_sample_int = 0;
  G_IN.presst_dxi = 0.0001;
  G_IN.presst_xi_max = 0.002;
  G_IN.presst_sample_int = 0;
  G_IN.press = 0;
  G_IN.dv_max = 0.001;
  G_IN.opt_flag = 1;
  G_IN.opt_sweeps = 1000;
  G_IN.opt_samples = 10;
  G_IN.opt_part_target = 0.5;
  G_IN.opt_vol_target = 0.5;
  G_IN.seed = 0;
  G_IN.ql_order = 6;
  G_IN.ql_rmax = 1.5;
  G_IN.ql_sample_int = 0;
  G_IN.mu_sample_int = 0;
  G_IN.mu_insertions = 100;
  G_IN.rdf_dr = 0.01;
  G_IN.rdf_rmax = 10;
  G_IN.rdf_sample_int = 0;
  G_IN.rdf_samples = 100;

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
	  G_IN.rho = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"cells_x") == 0 || strcmp(key,"cells_x\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.nx = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }


      else if (strcmp(key,"cells_y") == 0 || strcmp(key,"cells_y\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.ny = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }


      else if (strcmp(key,"cells_z") == 0 || strcmp(key,"cells_z\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.nz = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }


      else if (strcmp(key,"type") == 0 || strcmp(key,"type\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.type = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"neigh_list") == 0 || strcmp(key,"neigh_list\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.neigh_dr = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.neigh_max_part = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"dr_max") == 0 || strcmp(key,"dr_max\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.dr_max = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"sweep_eq") == 0 || strcmp(key,"sweep_eq\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.sweep_eq = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"sweep_stat") == 0 || strcmp(key,"sweep_stat\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.sweep_stat = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"out") == 0 || strcmp(key,"out\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.output_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"npt") == 0 || strcmp(key,"npt\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.press = atof(value);
	}
	else read_input_file_err(1,line_buf);
        value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.dv_max = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"opt") == 0 || strcmp(key,"opt\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.opt_flag = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.opt_sweeps = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.opt_samples = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.opt_part_target = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.opt_vol_target = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"seed") == 0 || strcmp(key,"seed\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.seed = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"cavity") == 0 || strcmp(key,"cavity\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.cavity_pcav = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.cavity_maxdr = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.cavity_mindr = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.cavity_sample_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.cavity_out_dr = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"cluster") == 0 || strcmp(key,"cluster\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.cluster_flag = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.cluster_moves_sweep = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.cluster_init_step = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"restart_read") == 0 || strcmp(key,"restart_read\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.restart_read = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  int len = strlen(value);
	  if (value[len-1] == '\n') value[len-1] = '\0';
	  if (len < 100) strcpy(G_IN.restart_name,value);
	  else read_input_file_err(3,line_buf);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"restart_write") == 0 || strcmp(key,"restart_write\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.restart_write = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"config_write") == 0 || strcmp(key,"config_write\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.config_write = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.config_samples = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"press_virial") == 0 || strcmp(key,"press_virial\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.pressv_dr = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.pressv_sample_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"press_thermo") == 0 || strcmp(key,"press_thermo\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.presst_dxi = atof(value);
	}
	else read_input_file_err(1,line_buf);
        value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.presst_xi_max = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.presst_sample_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"ql") == 0 || strcmp(key,"ql\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.ql_order = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.ql_rmax = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.ql_sample_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"widom") == 0 || strcmp(key,"widom\n") == 0){
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.mu_insertions = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
        if(value != NULL ) {
	  G_IN.mu_sample_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"rdf") == 0 || strcmp(key,"rdf\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.rdf_dr = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.rdf_rmax = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.rdf_sample_int = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  G_IN.rdf_samples = atoi(value);
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
  printf("Done\n");

  // Print content of input structure
  /* printf("Density: %.8f\n", G_IN.rho); */
  /* printf("Number of cells along x: %d\n", G_IN.nx); */
  /* printf("Number of cells along y: %d\n", G_IN.ny); */
  /* printf("Number of cells along z: %d\n", G_IN.nz); */
  /* printf("Cell type: %d\n", G_IN.type); */
  /* printf("Maximum MC displacement: %.8f\n", G_IN.dr_max); */
  /* printf("Number of sweeps (equilibration): %d\n", G_IN.sweep_eq); */
  /* printf("Number of sweeps (statistics): %d\n", G_IN.sweep_stat); */
  /* printf("Output interval (sweeps): %d\n", G_IN.output_int); */
  /* printf("Pressure - virial (resolution): %.8f\n", G_IN.pressv_dr); */
  /* printf("Pressure - virial (sweeps/sample): %d\n", G_IN.pressv_sample_int); */
  /* printf("Pressure - thermo (resolution): %.8f\n", G_IN.presst_dxi); */
  /* printf("Pressure - thermo (max compression): %.8f\n", G_IN.presst_xi_max); */
  /* printf("Pressure - thermo (sweeps/sample): %d\n", G_IN.presst_sample_int); */
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

