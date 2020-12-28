#include <time.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_press.h"
#include "compute_order_parameter.h"
#include "compute_widom_chem_pot.h"
#include "moves.h"
#include "optimizer.h"
#include "cavity_nvt.h"

// Cavity simulation for hard-spheres in the NVT ensemble
void cavity_hs_nvt() {

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

  // Optmize maximum displacement
  if (in.opt_flag == 1){
    opt_cavity_nvt();
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
  cavity_run_nvt(false);
  printf("Equilibration completed.\n");

  // Run statistics
  printf("---------------------------------------------------\n");
  printf("Production...\n");
  cavity_run_nvt(true);
  printf("Production completed.\n");
  clock_t end = clock();

  // Print acceptance and rejection percentages 
  printf("---------------------------------------------------\n");
  printf("-- Particle moves: %.8e\n", (double)part_moves);
  printf("   Acceptance percentage: %f\n", (double)acc_part_moves/((double)part_moves));
  printf("   Rejection percentage: %f\n", (double)rej_part_moves/((double)part_moves));

  
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

void cavity_run_nvt(bool prod_flag){

  // Variable declaration
  bool cavity_init = true;
  int n_sweeps;

  // Number of sweeps
  if (prod_flag) n_sweeps = in.sweep_stat;
  else n_sweeps = in.sweep_eq;
  
  // Run MC simulation
  for (int ii=0; ii<n_sweeps; ii++){

    // Sweep
    cavity_sweep_nvt();

    if (prod_flag){
      // Output distance between the cavities
      if (in.cavity_sample_int > 0){
        if (ii % in.cavity_sample_int == 0) {
          cavity_dist_output(cavity_init);
          if (cavity_init) cavity_init = false;
        }
      }
    }

    // Output on screen
    if (ii == 0){
      printf("Sweep number\n");
    }
    if (ii % in.output_int == 0) {
      printf("%d\n", ii);
      fflush(stdout);
    }


  }  

}

void cavity_sweep_nvt(){

  // Create N trial moves (N = number of particles)
  for (int ii=0; ii<part_info.NN; ii++){
    cavity_part_move();
    part_moves+=1;
  }

}


void cavity_dist_output(bool init){

  // Print to file the sample needed to compute the pressure
  FILE* fid;
  if (init) fid = fopen("cavity_distance.dat", "w");
  else fid = fopen("cavity_distance.dat", "a");
  if (fid == NULL) {
    perror("Error while creating the file for the cavity\n");
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
