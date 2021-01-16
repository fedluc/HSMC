#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include <gsl/gsl_rng.h>
#include "sim_info.h"
#include "rng.h"
#include "read_input.h"
#include "io_config.h"

void write_restart(int sweep){
  
  // Name of restart file (defined by the sweep number)
  int max_out_length = ceil(log10(G_IN.sweep_stat+G_IN.sweep_eq));
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

  // Write maximum step information (could be affected by optimization)  
  fwrite(&G_IN.dr_max, sizeof(double), 1, fid); 

  // Write maximum volume step information (could be affected by optimization)  
  fwrite(&G_IN.dv_max, sizeof(double), 1, fid); 

  // Write Simulation box information
  sim_box_info_write(fid);

  // Write particles information
  part_info_write(fid);

  // Write configuration
  part_conf_write(fid);

  // Write random number generator status
  rng_write(fid);

  // Close binary file
  fclose(fid);

}


void read_restart(char *restart_file){

  // Open binary file
  printf("Reading data from restart file %s...\n", restart_file);
  FILE *fid = NULL;
  fid = fopen(restart_file, "rb");
  if (fid == NULL) {
    perror("Error while reading restart file");
    exit(EXIT_FAILURE);
  }

  // Read maximum displacement 
  fread(&G_IN.dr_max, sizeof(double), 1, fid);

  // Read maximum volume displacement 
  fread(&G_IN.dv_max, sizeof(double), 1, fid);
  
  // Read simulation box
  sim_box_info_read(fid);

  // Read particles information
  part_info_read(fid);

  // Allocate particles
  part_alloc();
 
  // Read configuration
  part_conf_read(fid);

  // Read random number generator status
  rng_read(fid);
  
  // Close binary file
  fclose(fid);

  // Compute the density
  box_info sim_box_info = sim_box_info_get();
  p_info part_info = part_info_get();
  G_IN.rho = part_info.NN/sim_box_info.vol;

  // Print message on screen
  printf("The following data was initialized via the restart file:\n"
	 "- Maximum displacement (override with optimization)\n"
	 "- Maximum volume displacement (only for NpT, override with optimization)\n"
	 "- Dimensions of the simulation box\n"
	 "- Cell list information\n"
	 "- Number of particles\n"
	 "- Particle's positions\n"
	 "- Density\n"
	 "- Status of the random number generator\n");

}

void write_config(int sweep){
  
  // Name of restart file (defined by the sweep number)
  int max_out_length = ceil(log10(G_IN.sweep_stat+G_IN.sweep_eq));
  char config_file_template[20], config_file[max_out_length+20];
  sprintf(config_file_template, "config_%%0%dd.dat.gz", max_out_length);
  sprintf(config_file, config_file_template, sweep);
  
  // Open binary file
  gzFile fid = gzopen(config_file, "w");
  if (fid == NULL) {
    perror("Error while creating configuration file\n");
    exit(EXIT_FAILURE);
  }

  // Write configuration
  p_info part_info = part_info_get();
  box_info sim_box_info = sim_box_info_get();
  config part_conf = part_config_get();
  gzprintf(fid, "# Sweep number\n");
  gzprintf(fid, "%d\n", sweep);
  gzprintf(fid, "# Number of particles\n");
  gzprintf(fid, "%d\n", part_info.NN);
  gzprintf(fid, "# Simulation box size\n");
  gzprintf(fid, "%.8f\n", sim_box_info.lx);
  gzprintf(fid, "%.8f\n", sim_box_info.ly);
  gzprintf(fid, "%.8f\n", sim_box_info.lz);
  gzprintf(fid, "# Configuration\n");
  for (int ii=0; ii<part_info.NN; ii++){
    gzprintf(fid, "%d %.8f %.8f %.8f\n", 
	     (int)part_conf[ii][0], part_conf[ii][1],
	     part_conf[ii][2], part_conf[ii][3]);
  }

  // Close binary file
  gzclose(fid);


}
