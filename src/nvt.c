#include <time.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_press.h"
#include "moves.h"
#include "nvt.h"

// Hard-sphere simulation in the NVT ensemble
void hs_nvt() {

  // Simulation box
  sim_box_init(in.type, in.nx, in.ny, in.nz, in.rho);
  printf("Simulation box size (x, y, z): %.5f %.5f %.5f\n", sim_box_info.lx,
	 sim_box_info.ly, sim_box_info.lz);

  // Particles
  part_alloc();
  printf("Number of particles: %d\n", part_info.NN);

  // Initialize particle's positions
  part_init();

  // Initialize cell lists
  cell_list_init();

  // Set-up random number generator (Marsenne-Twister)
  rng_init();

  // Run to determine the optimal maximum displacement
  if (in.dr_max < 0){
    in.dr_max *= -1;
    run_opt_nvt();
    part_init();
  }


  // Start timing
  clock_t start = clock();

  // Initialize move counters
  part_moves = 0;
  acc_part_moves = 0;
  rej_part_moves = 0;  

  // Run equilibration
  printf("---------------------------------------------------\n");
  printf("Equilibration...\n");
  run_nvt(false);
  printf("Equilibration completed.\n");

  // Run statistics
  printf("---------------------------------------------------\n");
  printf("Production...\n");
  run_nvt(true);
  printf("Production completed.\n");
  clock_t end = clock();

  // Print acceptance and rejection percentages 
  printf("Acceptance percentage: %f\n", (double)acc_part_moves/((double)part_moves));
  printf("Rejection percentage: %f\n", (double)rej_part_moves/((double)part_moves));
  
  // Stop timing
  printf("Elapsed time: %f seconds\n",
	 (double)(end - start) / CLOCKS_PER_SEC);
  
  // Free memory
  free(part);
  gsl_rng_free(rng_mt);
  free(cl_neigh);
  free(cl_head);
  free(cl_link);

}


void run_opt_nvt(){

  int max_iter=1000, n_samples=10, sample_iter = max_iter/n_samples;
  double acc_ratio_1, acc_ratio_2, dr_1, dr_2;

  printf("---------------------------------------------------\n");
  printf("Maximum displacement optimization started ...\n");
  printf("Sweeps for optimization: %d\n", max_iter);
  printf("Number of samples: %d\n", n_samples);

  // First step
  acc_part_moves=0, rej_part_moves=0;
  for (int ii=0; ii<sample_iter; ii++){
    sweep_nvt();
  }
  dr_1 = in.dr_max;
  acc_ratio_1 = (double)acc_part_moves/((double)sample_iter*part_info.NN);

  // Second step
  if (acc_ratio_1 > 0.5){
    in.dr_max *= 2;
  }
  else {
    in.dr_max /= 2;
  }
  acc_part_moves=0, rej_part_moves=0;
  for (int ii=0; ii<sample_iter; ii++){
    sweep_nvt();
  }
  dr_2 = in.dr_max;
  acc_ratio_2 = (double)acc_part_moves/((double)sample_iter*part_info.NN);

  // Secant-method to find optimum value
  for (int ii=0; ii<n_samples; ii++){

    in.dr_max = dr_2 - (acc_ratio_2 - 0.5) * (dr_2 - dr_1)/(acc_ratio_2 - acc_ratio_1);
    if (in.dr_max > 1.0) {
      in.dr_max = 1.0;
    }
    else if (in.dr_max <= 0.0) {
      printf("Error: maximum displacement is zero!\n");
      exit(EXIT_FAILURE);
    }
    acc_part_moves=0, rej_part_moves=0;
    for (int jj=0; jj<sample_iter; jj++){
      sweep_nvt();
    }
    dr_1 = dr_2;
    dr_2 = in.dr_max;
    acc_ratio_1 = acc_ratio_2;
    acc_ratio_2 = (double)acc_part_moves/((double)sample_iter*part_info.NN);
  }

  printf("Optimal maximum displacement: %.8f. Acceptance ratio: %.8f \n", in.dr_max, acc_ratio_2);
  printf("Maximum displacement optimization completed\n");

}

void run_nvt(bool prod_flag){

  // Variable declaration
  bool pressv_init = true, presst_init = true;
  int n_sweeps;

  // Number of sweeps
  if (prod_flag) n_sweeps = in.sweep_stat;
  else n_sweeps = in.sweep_eq;
  
  // Run MC simulation
  for (int ii=0; ii<n_sweeps; ii++){

    // Sweep
    sweep_nvt();

    if (prod_flag){

      // Compute pressure via virial route
      if (in.pressv_sample_int > 0){
	if (ii % in.pressv_sample_int == 0) {
	  compute_pressv(pressv_init);
	  if (pressv_init) pressv_init = false;
	}
      }


      // Compute pressure via thermodynamic route
      if (in.presst_sample_int > 0){
	if (ii % in.presst_sample_int == 0) {
	  compute_presst(presst_init);
	  if (presst_init) presst_init = false;
	}
      }

    }

    // Output on screen
    if (ii % in.output_int == 0) {
      printf("Sweep number: %d\n", ii);
      fflush(stdout);
    }


  }  

}

void sweep_nvt(){

  // Create N trial moves (N = number of particles)
  for (int ii=0; ii<part_info.NN; ii++){
    part_move();
    part_moves+=1;
  }

}
