#include <time.h>
#include "init.h"
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

  if (in.restart_read == 0){ // Read from input file

    // Initialize simulation box
    sim_box_init(in.type, in.nx, in.ny, in.nz, in.rho);
  
    // Particles
    part_alloc();

    // Initialize particle's positions
    part_init();

    // Set-up random number generator (Marsenne-Twister)
    rng_init();

  }
  else { // Read from restart file
    read_restart(in.restart_name);
  }

  // Print simulation info on screen
  printf("Simulation box size (x, y, z): %.5f %.5f %.5f\n", sim_box_info.lx,
         sim_box_info.ly, sim_box_info.lz);
  printf("Number of particles: %d\n", part_info.NN);
  printf("Pressure: %.8f\n", in.press);

  // Initialize cell lists
  cell_list_init();

  // Optmize maximum displacement
  if (in.opt_flag == 1){
    double rho_start = in.rho;
    opt_npt();
    sim_box_init(in.type, in.nx, in.ny, in.nz, rho_start);
    part_init();
    cell_list_init();

  }

  // Start timing
  clock_t start = clock();

  // Initialize move counters
  part_moves = 0;
  acc_part_moves = 0;
  rej_part_moves = 0;
  vol_moves = 0;
  acc_vol_moves = 0;
  rej_vol_moves = 0;

  // Run equilibration
  printf("---------------------------------------------------\n");
  printf("Equilibration...\n");
  run_npt(false,0);
  printf("Equilibration completed.\n");

  // Run statistics
  printf("---------------------------------------------------\n");
  printf("Production...\n");
  run_npt(true,in.sweep_eq);
  printf("Production completed.\n");
  clock_t end = clock();

  // Print acceptance and rejection percentages 
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
  free(part);
  gsl_rng_free(rng_mt);
  free(cl_neigh);
  free(cl_part_cell);

}


void run_npt(bool prod_flag, int sweep_offset){

  // Variable declaration
  bool presst_init = true, ql_ave_init = true;
  int n_sweeps;

  if (prod_flag) n_sweeps = in.sweep_stat;
  else n_sweeps = in.sweep_eq;

  // Run MC simulation
  for (int ii=sweep_offset; ii<n_sweeps+sweep_offset; ii++){

    // Output on screen
    if (ii == 0){
      printf("Sweep number  Density\n");
    }
    if (ii % in.output_int == 0) {
      printf("%d  %.8f\n", ii,in.rho);
      fflush(stdout);
    }

    // Write restart file
    if (in.restart_write > 0){
      if (ii % in.restart_write == 0) {
        write_restart(ii);
      }
    }

    
    // Save samples for production runs
    if (prod_flag){      

      // Write configuration
      if (in.config_write > 0){
        if (ii % in.config_write == 0) {
          write_config(ii);
        }
      }


      // Compute pressure via thermodynamic route
      if (in.presst_sample_int > 0){
      	if (ii % in.presst_sample_int == 0) {
      	  compute_presst(presst_init);
      	  if (presst_init) presst_init = false;
      	}
      }

      // Compute order parameter
      if (in.ql_sample_int > 0){
        if (ii % in.ql_sample_int == 0) {
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

  // Create N trial moves (N = number of particles)
  for (int ii=0; ii<part_info.NN; ii++){

    // Perform (on average) one volume move every N particles moves
    r_move_id = gsl_rng_uniform_int(rng_mt, part_info.NN+1);

    if (r_move_id < part_info.NN){
      part_move(); // Move one particle
      part_moves+=1;
    }
    else{
      vol_move(); // Change the volume
      vol_moves+=1;
    }
    
  }

}
