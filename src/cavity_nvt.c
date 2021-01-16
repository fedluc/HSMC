#include <stdlib.h>
#include <time.h>
#include "sim_info.h"
#include "rng.h"
#include "read_input.h"
#include "cell_list.h"
#include "moves.h"
#include "io_config.h"
#include "optimizer.h"
#include "cavity_nvt.h"
#include "analytic.h"

// ------ Cavity simulation for hard-spheres in the NVT ensemble ------

void cavity_hs_nvt() {

  if (G_IN.restart_read == 0){ // Read from input file

    // Simulation box
    sim_box_init(G_IN.type, G_IN.nx, G_IN.ny, G_IN.nz, G_IN.rho);

    // Particles
    part_alloc();

    // Initialize particle's positions
    part_init();

    // Set cavities to a distance within the allowed interval
    cavity_set_distance();

    // Set-up random number generator (Marsenne-Twister)
    rng_init();

  }
  else { // Read from restart file
    read_restart(G_IN.restart_name);
  }

  // Print simulation info on screen  
  print_sim_info();

  // Initialize cavity interaction potential (0.5 is just a dummy distance)
  cavity_interaction(0.5, true);

  // Write cavity interaction potential to file 
  cavity_psi_output();

  // Set-up the neighbor list
  cell_list_init(true);

  // Optmize maximum displacement
  if (G_IN.opt_flag == 1){
    opt_cavity_nvt();
    part_init();
    cavity_set_distance();
    cell_list_init(false);
  }

  // Start timing
  clock_t start = clock();

  // Initialize move counters
  reset_moves_counters();

  // Run equilibration
  printf("---------------------------------------------------\n");
  printf("Equilibration...\n");
  cavity_run_nvt(false,0);
  printf("Equilibration completed.\n");

  // Run statistics
  printf("---------------------------------------------------\n");
  printf("Production...\n");
  cavity_run_nvt(true,G_IN.sweep_eq);
  printf("Production completed.\n");
  clock_t end = clock();

  // Print acceptance and rejection percentages
  int part_moves, acc_part_moves, rej_part_moves;
  get_moves_counters(&part_moves, &acc_part_moves, &rej_part_moves,
                     NULL, NULL, NULL); 
  printf("---------------------------------------------------\n");
  printf("-- Particle moves: %.8e\n", (double)part_moves);
  printf("   Acceptance percentage: %f\n", (double)acc_part_moves/((double)part_moves));
  printf("   Rejection percentage: %f\n", (double)rej_part_moves/((double)part_moves));

  
  // Stop timing
  printf("Elapsed time: %f seconds\n",
	 (double)(end - start) / CLOCKS_PER_SEC);
  
  // Free memory
  part_free();
  rng_free();
  cell_list_free();

}

void cavity_run_nvt(bool prod_flag, int sweep_offset){

  // Variable declaration
  bool cavity_init = true;
  int n_sweeps;

  // Number of sweeps
  if (prod_flag) n_sweeps = G_IN.sweep_stat;
  else n_sweeps = G_IN.sweep_eq;
  
  // Run MC simulation
  for (int ii=sweep_offset; ii<n_sweeps+sweep_offset; ii++){

    // Output on screen
    if (ii == 0){
      printf("Sweep number\n");
    }
    if (ii % G_IN.output_int == 0) {
      printf("%d\n", ii);
      fflush(stdout);
    }

    // Write restart file
    if (G_IN.restart_write > 0){
      if (ii % G_IN.restart_write == 0) {
        write_restart(ii);
      }
    }

    // Save samples for production runs
    if (prod_flag){

      // Write configuration
      if (G_IN.config_write > 0){
        if (ii % G_IN.config_write == 0) {
          write_config(ii);
        }
      }


      // Output distance between the cavities
      if (G_IN.cavity_sample_int > 0){
        if (ii % G_IN.cavity_sample_int == 0) {
          cavity_dist_output(cavity_init);
          if (cavity_init) cavity_init = false;
        }
      }

    }
    
    // Generate new configuration
    cavity_sweep_nvt();

  }  

}

void cavity_sweep_nvt(){

  // Create N trial moves (N = number of particles)
  p_info part_info = part_info_get();
  for (int ii=0; ii<part_info.NN; ii++){
    cavity_part_move();
  }

}

// ------ Set cavity distance to a value allowed from input ------

void cavity_set_distance(){
  
  double cavity_avedr = (G_IN.cavity_maxdr + G_IN.cavity_mindr)/ 2.0;
  box_info sim_box_info = sim_box_info_get();
  config part_conf = part_config_get();
  part_conf[1][1] = part_conf[0][1];
  part_conf[1][2] = part_conf[0][2];
  part_conf[1][3] = part_conf[0][3] + cavity_avedr;
  if (part_conf[1][3] > sim_box_info.lz) part_conf[1][3] -= sim_box_info.lz;
  else if (part_conf[1][3] < 0.0)        part_conf[1][3] += sim_box_info.lz;

}


// ------ Write cavity interaction potential to file ------

void cavity_psi_output(){

  FILE* fid;
  fid = fopen("cavity_psi.dat", "w");
  if (fid == NULL) {
    perror("Error while creating the file for the cavity potential\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fid, "#########################################################\n");
  fprintf(fid, "# Interaction potential  between the cavities \n");
  fprintf(fid, "# Distance, potential\n");
  fprintf(fid, "#########################################################\n");
  double xx = G_IN.cavity_mindr;
  int n_points = (int)((G_IN.cavity_maxdr - G_IN.cavity_mindr)/G_IN.cavity_out_dr);
  for (int ii=0; ii < n_points; ii++){
    fprintf(fid, "%.8e %.8e\n", xx, cavity_interaction(xx, false));
    xx += G_IN.cavity_out_dr;
  }
  fclose(fid);

}

// ------ Write cavity distance to file ------

void cavity_dist_output(bool init){

  FILE* fid;
  if (init) fid = fopen("cavity_distance.dat", "w");
  else fid = fopen("cavity_distance.dat", "a");
  if (fid == NULL) {
    perror("Error while creating the file for the cavity distance\n");
    exit(EXIT_FAILURE);
  }
  if (init){
    fprintf(fid, "#########################################################\n");
    fprintf(fid, "# Distance between the cavities (each line is one sample)\n");
    fprintf(fid, "#########################################################\n");
  }
  fprintf(fid, "%.8e\n", compute_dist(0,1,1.0,1.0,1.0));
  fclose(fid);

}



