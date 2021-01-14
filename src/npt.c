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

  // Set-up the neighbor list
  int cl_neigh_num, cl_num_tot;
  compute_cell_list_info();
  get_cell_list_info(&cl_neigh_num, &cl_num_tot, NULL, NULL, NULL,
                     NULL, NULL, NULL);
  int (*cl_neigh)[cl_neigh_num] = (int (*)[cl_neigh_num])cell_list_alloc(cl_num_tot, cl_neigh_num);
  int (*cl_part_cell)[in.neigh_max_part] = (int (*)[in.neigh_max_part])cell_list_alloc(cl_num_tot, in.neigh_max_part);
  cell_list_init(cl_neigh_num, cl_neigh, in.neigh_max_part, cl_part_cell);

  // Optmize maximum displacement
  if (in.opt_flag == 1){
    double rho_start = in.rho;
    opt_npt(cl_num_tot, in.neigh_max_part, cl_part_cell,
	    cl_neigh_num, cl_neigh);
    sim_box_init(in.type, in.nx, in.ny, in.nz, rho_start);
    part_init();
    cell_list_init(cl_neigh_num, cl_neigh, in.neigh_max_part, cl_part_cell);
  }

  // Start timing
  clock_t start = clock();

  // Initialize move counters
  reset_moves_counters();

  // Run equilibration
  printf("---------------------------------------------------\n");
  printf("Equilibration...\n");
  run_npt(false,0, cl_num_tot, in.neigh_max_part, cl_part_cell,
          cl_neigh_num, cl_neigh);
  printf("Equilibration completed.\n");

  // Run statistics
  printf("---------------------------------------------------\n");
  printf("Production...\n");
  run_npt(true,in.sweep_eq, cl_num_tot, in.neigh_max_part, cl_part_cell,
          cl_neigh_num, cl_neigh);
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
  free(part);
  gsl_rng_free(rng_mt);
  free(cl_neigh);
  free(cl_part_cell);

}


void run_npt(bool prod_flag, int sweep_offset,
	     int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
             int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]){

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
      	  compute_presst(presst_init,
			 cl_num_tot, cl_max_part, cl_part_cell,
			 cl_neigh_num, cl_neigh);
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
    sweep_npt(cl_num_tot, cl_max_part, cl_part_cell,
              cl_neigh_num, cl_neigh);

  }  

}

void sweep_npt(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
               int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]){

  int r_move_id;

  // Create N trial moves (N = number of particles)
  for (int ii=0; ii<part_info.NN; ii++){

    // Perform (on average) one volume move every N particles moves
    r_move_id = gsl_rng_uniform_int(rng_mt, part_info.NN+1);

    if (r_move_id < part_info.NN){
      part_move(cl_num_tot, cl_max_part, cl_part_cell,
		cl_neigh_num, cl_neigh); // Move one particle
    }
    else{
      vol_move(cl_num_tot, cl_max_part, cl_part_cell,
	       cl_neigh_num, cl_neigh); // Change the volume
    }
    
  }

}
