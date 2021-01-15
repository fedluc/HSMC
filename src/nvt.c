#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "sim_info.h"
#include "rng.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_press.h"
#include "compute_order_parameter.h"
#include "compute_widom_chem_pot.h"
#include "moves.h"
#include "io_config.h"
#include "optimizer.h"
#include "nvt.h"

// Hard-sphere simulation in the NVT ensemble
void hs_nvt() {

  if (G_IN.restart_read == 0){ // Read from input file  

    // Simulation box
    sim_box_init(G_IN.type, G_IN.nx, G_IN.ny, G_IN.nz, G_IN.rho);

    // Particles
    part_alloc();

    // Initialize particle's positions
    part_init();
    
    // Set-up random number generator (Marsenne-Twister)
    rng_init();

  }
  else { // Read from restart file
    read_restart(G_IN.restart_name);
  }

  // ---------------- TESTING ----------------------
  double (*part_test)[4] = part_config_get();
  for (int ii=0; ii<100; ii++){
    printf("%f %f %f %f\n", part_test[ii][0], part_test[ii][1],
	   part_test[ii][2], part_test[ii][3]);
  }
  // ---------------- TESTING ----------------------  
  // Print simulation info on screen
  print_sim_info();

  // Set-up the neighbor list
  cell_list_init(true);

  // Optmize maximum displacement
  if (G_IN.opt_flag == 1){
    opt_nvt();
    part_init();
    cell_list_init(false);
  }

  // Start timing
  clock_t start = clock();

  // Initialize move counters
  reset_moves_counters();

  // Run equilibration
  printf("---------------------------------------------------\n");
  printf("Equilibration...\n");
  run_nvt(false,0);
  printf("Equilibration completed.\n");

  // Run statistics
  printf("---------------------------------------------------\n");
  printf("Production...\n");
  run_nvt(true,G_IN.sweep_eq);
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

  /* // Free memory  */
  part_free();
  cell_list_free();
  rng_free();

}

void run_nvt(bool prod_flag, int sweep_offset){

  // Variable declaration
  bool pressv_init = true, presst_init = true;
  bool ql_ave_init = true, mu_ave_init = true;
  int n_sweeps;

  // Number of sweeps
  if (prod_flag) n_sweeps = G_IN.sweep_stat;
  else n_sweeps = G_IN.sweep_eq;
  
  // Run MC simulation
  for (int ii=sweep_offset; ii<n_sweeps+sweep_offset; ii++){

    // Output on screen
    if (ii == sweep_offset){
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

      // Compute pressure via virial route
      if (G_IN.pressv_sample_int > 0){
      	if (ii % G_IN.pressv_sample_int == 0) {
      	  compute_pressv(pressv_init);
      	  if (pressv_init) pressv_init = false;
      	}
      }

      // Compute pressure via thermodynamic route
      if (G_IN.presst_sample_int > 0){
      	if (ii % G_IN.presst_sample_int == 0) {
      	  compute_presst(presst_init);
      	  if (presst_init) presst_init = false;
      	}
      }

      // Compute order parameter
      if (G_IN.ql_sample_int > 0){
      	if (ii % G_IN.ql_sample_int == 0) {
      	  compute_op(ql_ave_init);
      	  if (ql_ave_init) ql_ave_init = false;
      	}
      }

      // Compute chemical potential via Widom insertions
      if (G_IN.mu_sample_int > 0){
      	if (ii % G_IN.mu_sample_int == 0) {
      	  compute_mu(mu_ave_init);
      	  if (mu_ave_init) mu_ave_init = false;
      	}
      }

    }

    // Generate new configuration
    sweep_nvt();

  }

}

void sweep_nvt(){

  // Create N trial moves (N = number of particles)
  struct p_info part_info = part_info_get();
  for (int ii=0; ii<part_info.NN; ii++){
    part_move();
  }

}
