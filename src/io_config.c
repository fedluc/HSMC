#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "init.h"
#include "read_input.h"
#include "io_config.h"

void write_restart(int sweep){
  
  // Name of restart file (defined by the sweep number)
  int max_out_length = ceil(log10(in.sweep_stat+in.sweep_eq));
  char restart_file_template[20], restart_file[max_out_length+20];
  sprintf(restart_file_template, "restart_%%0%dd.bin", max_out_length);
  sprintf(restart_file, restart_file_template, sweep);
  
  // Open binary file
  FILE *fid = NULL;
  fid = fopen(restart_file, "wb");
  if (fid == NULL) {
    perror("Error while creating restart file\n");
    exit(EXIT_FAILURE);
  }

  // Write Simulation box information
  fwrite(&sim_box_info, sizeof(struct box_info), 1, fid);

  // Write particles information
  fwrite(&part_info, sizeof(struct p_info), 1, fid);

  // Write configuration
  fwrite(part, sizeof(double), part_info.NN*4, fid);

  // Write random number generator status
  gsl_rng_fwrite(fid, rng_mt);

  // Close binary file
  fclose(fid);

}


void read_restart(char *restart_file){

  // Declare the type of random number generator
  rng_mt = gsl_rng_alloc(gsl_rng_mt19937);
  
  // Open binary file
  FILE *fid = NULL;
  fid = fopen(restart_file, "rb");
  if (fid == NULL) {
    perror("Error while reading restart file");
    exit(EXIT_FAILURE);
  }

  // Read simulation box
  fread(&sim_box_info, sizeof(struct box_info), 1, fid);

  // Read particles information
  fread(&part_info, sizeof(struct p_info), 1, fid);

  // Allocate particles
  part_alloc();
 
  // Read configuration
  fread(part, sizeof(double), part_info.NN*4, fid);

  // Read random number generator status
  gsl_rng_fread(fid, rng_mt);

  // Close binary file
  fclose(fid);

  // Compute the density
  in.rho = part_info.NN/sim_box_info.vol;

  // Initialize maximum number produced by random number generator
  r_num_max = gsl_rng_max(rng_mt);


}
