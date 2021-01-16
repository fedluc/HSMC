#include <stdlib.h>
#include <time.h>
#include "sim_info.h"
#include "rng.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_press.h"
#include "compute_order_parameter.h"
#include "moves.h"
#include "io_config.h"
#include "optimizer.h"
#include "npt.h"

// Hard-sphere simulation in the NpT ensemble
void hs_npt() {

  if (G_IN.restart_read == 0){ // Read from input file

    // Initialize simulation box
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

  // Print simulation info on screen
  print_sim_info();

  // Set-up the neighbor list
  cell_list_init(true);


  // Optmize maximum displacement
  if (G_IN.opt_flag == 1){
    double rho_start = G_IN.rho;
    opt_npt();
    sim_box_init(G_IN.type, G_IN.nx, G_IN.ny, G_IN.nz, rho_start);
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
  run_npt(false,0);
  printf("Equilibration completed.\n");

  // Run statistics
  printf("---------------------------------------------------\n");
  printf("Production...\n");
  run_npt(true,G_IN.sweep_eq);
  printf("Production completed.\n");
  clock_t end = clock();

  // Print acceptance and rejection percentages
  int part_moves, acc_part_moves, rej_part_moves;
  int vol_moves, acc_vol_moves, rej_vol_moves;
  get_moves_counters(&part_moves, &acc_part_moves, &rej_part_moves,
                     &vol_moves, &acc_vol_moves, &rej_vol_moves);
  printf("---------------------------------------------------\n");
  printf("-- Particle moves: %.8e\n", (double)part_moves);
  printf("   Acceptance percentage: %f\n", (double)acc_part_moves/((double)part_moves));
  printf("   Rejection percentage: %f\n", (double)rej_part_moves/((double)part_moves));
  printf("-- Volume moves: %.8e\n", (double)vol_moves);
  printf("   Acceptance percentage: %f\n", (double)acc_vol_moves/((double)vol_moves));
  printf("   Rejection percentage: %f\n", (double)rej_vol_moves/((double)vol_moves));
  
  // Stop timing
  printf("Elapsed time: %f seconds\n",
	 (double)(end - start) / CLOCKS_PER_SEC);
  
  // Free memory
  part_free();
  cell_list_free();
  rng_free();

}


void run_npt(bool prod_flag, int sweep_offset){

  // Variable declaration
  bool presst_init = true, ql_ave_init = true;
  int n_sweeps;

  if (prod_flag) n_sweeps = G_IN.sweep_stat;
  else n_sweeps = G_IN.sweep_eq;

  // Run MC simulation
  for (int ii=sweep_offset; ii<n_sweeps+sweep_offset; ii++){

    // Output on screen
    if (ii == 0){
      printf("Sweep number  Density\n");
    }
    if (ii % G_IN.output_int == 0) {
      printf("%d  %.8f\n", ii,G_IN.rho);
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


    }

    // Generate a new configuration
    sweep_npt();

  }  

}

void sweep_npt(){

  int r_move_id;
  p_info part_info = part_info_get();
  // Create N trial moves (N = number of particles)
  for (int ii=0; ii<part_info.NN; ii++){

    // Perform (on average) one volume move every N particles moves
    r_move_id = rng_get_int(part_info.NN+1);

    if (r_move_id < part_info.NN){
      part_move(); // Move one particle
    }
    else{
      vol_move(); // Change the volume
    }
    
  }

}
