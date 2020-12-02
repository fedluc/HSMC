#include <time.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_press.h"
#include "moves.h"
#include "npt.h"

// Hard-sphere simulation in the NpT ensemble
void hs_npt() {

  // Initialize simulation box
  sim_box_init(in.type, in.nx, in.ny, in.nz, in.rho);

  // Particles
  part_alloc();
  printf("Number of particles: %d\n", part_info.NN);

  // Pressure 
  printf("Pressure: %.8f\n", in.npt_press);
  
  // Initialize particle's positions
  part_init();

  // Initialize cell lists
  cell_list_init();

  // Set-up random number generator (Marsenne-Twister)
  rng_init();

  
  /* // Run to determine the optimal maximum displacement */
  /* if (in.dr_max < 0){ */
  /*   in.dr_max *= -1; */
  /*   run_opt_npt(); */
  /*   part_init(); */
  /* } */

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
  run_npt(false);
  printf("Equilibration completed.\n");

  // Run statistics
  printf("---------------------------------------------------\n");
  printf("Production...\n");
  run_npt(true);
  printf("Production completed.\n");
  clock_t end = clock();

  // Print acceptance and rejection percentages 
  printf("---------------------------------------------------\n");
  printf("-- Particle moves: %d\n", part_moves);
  printf("   Acceptance percentage: %f\n", (double)acc_part_moves/((double)part_moves));
  printf("   Rejection percentage: %f\n", (double)rej_part_moves/((double)part_moves));
  printf("-- Volume moves: %d\n", vol_moves);
  printf("   Acceptance percentage: %f\n", (double)acc_vol_moves/((double)vol_moves));
  printf("   Rejection percentage: %f\n", (double)rej_vol_moves/((double)vol_moves));
  
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


/* void run_opt(){ */

/*   int max_iter=1000, n_samples=10, sample_iter = max_iter/n_samples; */
/*   double acc_ratio_1, acc_ratio_2, dr_1, dr_2; */

/*   printf("---------------------------------------------------\n"); */
/*   printf("Maximum displacement optimization started ...\n"); */
/*   printf("Sweeps for optimization: %d\n", max_iter); */
/*   printf("Number of samples: %d\n", n_samples); */

/*   // First step */
/*   acc_part_moves=0, rej_part_moves=0; */
/*   for (int ii=0; ii<sample_iter; ii++){ */
/*     sweep(); */
/*   } */
/*   dr_1 = in.dr_max; */
/*   acc_ratio_1 = (double)acc_part_moves/((double)sample_iter*part_info.NN); */

/*   // Second step */
/*   if (acc_ratio_1 > 0.5){ */
/*     in.dr_max *= 2; */
/*   } */
/*   else { */
/*     in.dr_max /= 2; */
/*   } */
/*   acc_part_moves=0, rej_part_moves=0; */
/*   for (int ii=0; ii<sample_iter; ii++){ */
/*     sweep(); */
/*   } */
/*   dr_2 = in.dr_max; */
/*   acc_ratio_2 = (double)acc_part_moves/((double)sample_iter*part_info.NN); */

/*   // Secant-method to find optimum value */
/*   for (int ii=0; ii<n_samples; ii++){ */

/*     in.dr_max = dr_2 - (acc_ratio_2 - 0.5) * (dr_2 - dr_1)/(acc_ratio_2 - acc_ratio_1); */
/*     if (in.dr_max > 1.0) { */
/*       in.dr_max = 1.0; */
/*     } */
/*     else if (in.dr_max <= 0.0) { */
/*       printf("Error: maximum displacement is zero!\n"); */
/*       exit(EXIT_FAILURE); */
/*     } */
/*     acc_part_moves=0, rej_part_moves=0; */
/*     for (int jj=0; jj<sample_iter; jj++){ */
/*       sweep(); */
/*     } */
/*     dr_1 = dr_2; */
/*     dr_2 = in.dr_max; */
/*     acc_ratio_1 = acc_ratio_2; */
/*     acc_ratio_2 = (double)acc_part_moves/((double)sample_iter*part_info.NN); */
/*   } */

/*   printf("Optimal maximum displacement: %.8f. Acceptance ratio: %.8f \n", in.dr_max, acc_ratio_2); */
/*   printf("Maximum displacement optimization completed\n"); */

/* } */

void run_npt(bool prod_flag){

  // Variable declaration
  bool presst_init = true;
  int n_sweeps;

  if (prod_flag) n_sweeps = in.sweep_stat;
  else n_sweeps = in.sweep_eq;

  // Run MC simulation
  for (int ii=0; ii<n_sweeps; ii++){

    // Sweep
    sweep_npt();

    if (prod_flag){      

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
      printf("Sweep number: %d, Density: %.8f\n", ii,in.rho);
      fflush(stdout);
    }


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
