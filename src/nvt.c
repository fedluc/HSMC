#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "init.h"
#include "rng.h"
#include "read_input.h"
#include "cell_list.h"
//#include "compute_press.h"
//#include "compute_order_parameter.h"
//#include "compute_widom_chem_pot.h"
#include "moves.h"
#include "io_config.h"
//#include "optimizer.h"
#include "nvt.h"

// Hard-sphere simulation in the NVT ensemble
void hs_nvt() {

  if (in.restart_read == 0){ // Read from input file  

    // Simulation box
    sim_box_init(in.type, in.nx, in.ny, in.nz, in.rho);

    // Particles
    part_alloc();

    // Initialize particle's positions
    part_init();
    
    // Set-up random number generator (Marsenne-Twister)
    rng_init();

  }
  else { // Read from restart file
    //read_restart(in.restart_name);
  }

  // Print simulation info on screen
  printf("Simulation box size (x, y, z): %.5f %.5f %.5f\n", sim_box_info.lx,
	 sim_box_info.ly, sim_box_info.lz);
  printf("Number of particles: %d\n", part_info.NN);

  // Set-up the neighbor list
  cell_list_init();

  /* // Optmize maximum displacement */
  /* if (in.opt_flag == 1){ */
  /*   opt_nvt(cl_num_tot, in.neigh_max_part, cl_part_cell, */
  /* 	    cl_neigh_num, cl_neigh); */
  /*   part_init(); */
  /*   cell_list_init(cl_neigh_num, cl_neigh, in.neigh_max_part, cl_part_cell); */
  /* } */

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
  run_nvt(true,in.sweep_eq);
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
  free(part);
  cell_list_free();
  rng_free();

}

void run_nvt(bool prod_flag, int sweep_offset){

  // Variable declaration
  bool pressv_init = true, presst_init = true;
  bool ql_ave_init = true, mu_ave_init = true;
  int n_sweeps;

  // Number of sweeps
  if (prod_flag) n_sweeps = in.sweep_stat;
  else n_sweeps = in.sweep_eq;
  
  // Run MC simulation
  for (int ii=sweep_offset; ii<n_sweeps+sweep_offset; ii++){

    // Output on screen
    if (ii == sweep_offset){
      printf("Sweep number\n");
    }
    if (ii % in.output_int == 0) {
      printf("%d\n", ii);
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

    /*   // Compute pressure via virial route */
    /*   if (in.pressv_sample_int > 0){ */
    /*   	if (ii % in.pressv_sample_int == 0) { */
    /*   	  compute_pressv(pressv_init, */
    /* 			 cl_num_tot, cl_max_part, cl_part_cell, */
    /* 			 cl_neigh_num, cl_neigh); */
    /*   	  if (pressv_init) pressv_init = false; */
    /*   	} */
    /*   } */

    /*   // Compute pressure via thermodynamic route */
    /*   if (in.presst_sample_int > 0){ */
    /*   	if (ii % in.presst_sample_int == 0) { */
    /*   	  compute_presst(presst_init, */
    /* 			 cl_num_tot, cl_max_part, cl_part_cell, */
    /* 			 cl_neigh_num, cl_neigh); */
    /*   	  if (presst_init) presst_init = false; */
    /*   	} */
    /*   } */

    /*   // Compute order parameter */
    /*   if (in.ql_sample_int > 0){ */
    /*   	if (ii % in.ql_sample_int == 0) { */
    /*   	  compute_op(ql_ave_init, */
    /* 		     cl_num_tot, cl_max_part, cl_part_cell, */
    /* 		     cl_neigh_num, cl_neigh); */
    /*   	  if (ql_ave_init) ql_ave_init = false; */
    /*   	} */
    /*   } */

    /*   // Compute chemical potential via Widom insertions */
    /*   if (in.mu_sample_int > 0){ */
    /*   	if (ii % in.mu_sample_int == 0) { */
    /*   	  compute_mu(mu_ave_init, */
    /* 		     cl_num_tot, cl_max_part, cl_part_cell, */
    /* 		     cl_neigh_num, cl_neigh); */
    /*   	  if (mu_ave_init) mu_ave_init = false; */
    /*   	} */
    /*   } */

    }

    // Generate new configuration
    sweep_nvt();

  }

}

void sweep_nvt(){

  // Create N trial moves (N = number of particles)
    for (int ii=0; ii<part_info.NN; ii++){
      part_move();
    }

}
