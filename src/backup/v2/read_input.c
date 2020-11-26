#include "hsmc.h"

// ----------------------------------------
// ------------- Read Input ---------------
// ----------------------------------------


void print_example(){

  char example[] = 
    
    "# This is an input file for hsmc\n"
    "# Lines starting with # or empty lines are skipped\n\n"
    "# Density\n"
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
    "# Radial distribution function (points, how often to compute, how often to average)\n"
    "rdf 1000 1000 10000\n\n"
    "# Maximum displacement for MC moves\n"
    "dr_max 0.05\n\n"
    "# Number of sweeps for equilibration (for N particles, 1 sweep = N moves)\n"
    "sweep_eq 10000\n\n"
    "# Number of sweeps for statistics\n"
    "sweep_stat 10000\n\n"
    "# Output interval (how often output gets printed)\n"
    "out 100\n\n";

    
    printf(example);

}

struct input  read_input_file(char *filename){

  FILE *in_file;
  char *line_buf = NULL;
  char *key, *value;
  size_t line_buf_size = 0;
  ssize_t line_size;
  struct input input_data;

  printf("Reading input data from %s ...\n",filename);
  
  // Open file
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
	  input_data.rho = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"cells_x") == 0 || strcmp(key,"cells_x\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  input_data.nx = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }


      else if (strcmp(key,"cells_y") == 0 || strcmp(key,"cells_y\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  input_data.ny = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }


      else if (strcmp(key,"cells_z") == 0 || strcmp(key,"cells_z\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  input_data.nz = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }


      else if (strcmp(key,"type") == 0 || strcmp(key,"type\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  input_data.type = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"dr_max") == 0 || strcmp(key,"dr_max\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  input_data.dr_max = atof(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"sweep_eq") == 0 || strcmp(key,"sweep_eq\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  input_data.sweep_eq = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"sweep_stat") == 0 || strcmp(key,"sweep_stat\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  input_data.sweep_stat = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"out") == 0 || strcmp(key,"out\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  input_data.dt_output = atoi(value);
	}
	else read_input_file_err(1,line_buf);
      }

      else if (strcmp(key,"rdf") == 0 || strcmp(key,"rdf\n") == 0){
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  input_data.rdf_dr = atof(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  input_data.rdf_tcompute = atoi(value);
	}
	else read_input_file_err(1,line_buf);
	value = strtok(NULL, " ");
	if(value != NULL ) {
	  input_data.rdf_tave = atoi(value);
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

  // Print content of input_data structure
  printf("Density: %.8f\n", input_data.rho);
  printf("Number of cells along x: %d\n", input_data.nx);
  printf("Number of cells along y: %d\n", input_data.ny);
  printf("Number of cells along z: %d\n", input_data.nz);
  printf("Cell type: %d\n", input_data.type);
  printf("Maximum MC displacement: %.8f\n", input_data.dr_max);
  printf("Number of sweeps (equilibration): %d\n", input_data.sweep_eq);
  printf("Number of sweeps (statistics): %d\n", input_data.sweep_stat);
  printf("Output interval (sweeps): %d\n", input_data.dt_output);
  printf("rdf (resolution): %.8f\n", input_data.rdf_dr);
  printf("rdf (sweeps/sample): %d\n", input_data.rdf_tcompute);
  printf("rdf (samples/average): %.d\n", input_data.rdf_tave);
  fflush(stdout);

  // Output
  return input_data;

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
    }

  printf("Last read line in the input file:\n%s\n", last_string);
  exit(EXIT_FAILURE);

}
