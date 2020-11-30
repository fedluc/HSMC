#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <time.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_press.h"
#include "hsmc.h"

// ----------------------------------------
// ------------- Monte-Carlo --------------
// ----------------------------------------

// Global variables for random number generator
gsl_rng *rng_mt;
long unsigned int r_num_max;

// Global variables for particles moves
int acc_moves, rej_moves;

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
  rng_mt = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng_mt,0);
  r_num_max = gsl_rng_max(rng_mt);

  // Run to determine the optimal maximum displacement
  if (in.dr_max < 0){
    in.dr_max *= -1;
    run_opt();
    part_init();
  }

  // Start timing
  clock_t start = clock();

  // Run equilibration
  printf("---------------------------------------------------\n");
  printf("Equilibration...\n");
  run(false);
  printf("Equilibration completed.\n");

  // Run statistics
  printf("---------------------------------------------------\n");
  printf("Production...\n");
  run(true);
  printf("Production completed.\n");
  clock_t end = clock();

  // Print acceptance and rejection percentages 
  int moves_tot = (in.sweep_eq+in.sweep_stat)*part_info.NN;
  printf("Acceptance percentage: %f\n", (double)acc_moves/((double)moves_tot));
  printf("Rejection percentage: %f\n", (double)rej_moves/((double)moves_tot));
  
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


void run_opt(){

  int max_iter=1000, n_samples=10, sample_iter = max_iter/n_samples;
  double acc_ratio_1, acc_ratio_2, dr_1, dr_2;

  printf("---------------------------------------------------\n");
  printf("Maximum displacement optimization started ...\n");
  printf("Sweeps for optimization: %d\n", max_iter);
  printf("Number of samples: %d\n", n_samples);

  // First step
  acc_moves=0, rej_moves=0;
  for (int ii=0; ii<sample_iter; ii++){
    sweep();
  }
  dr_1 = in.dr_max;
  acc_ratio_1 = (double)acc_moves/((double)sample_iter*part_info.NN);

  // Second step
  if (acc_ratio_1 > 0.5){
    in.dr_max *= 2;
  }
  else {
    in.dr_max /= 2;
  }
  acc_moves=0, rej_moves=0;
  for (int ii=0; ii<sample_iter; ii++){
    sweep();
  }
  dr_2 = in.dr_max;
  acc_ratio_2 = (double)acc_moves/((double)sample_iter*part_info.NN);

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
    acc_moves=0, rej_moves=0;
    for (int jj=0; jj<sample_iter; jj++){
      sweep();
    }
    dr_1 = dr_2;
    dr_2 = in.dr_max;
    acc_ratio_1 = acc_ratio_2;
    acc_ratio_2 = (double)acc_moves/((double)sample_iter*part_info.NN);
  }

  printf("Optimal maximum displacement: %.8f. Acceptance ratio: %.8f \n", in.dr_max, acc_ratio_2);
  printf("Maximum displacement optimization completed\n");

}

void run(bool prod_flag){

  // Variable declaration
  bool pressv_init = true, presst_init = true;
  int n_sweeps;

  if (prod_flag) n_sweeps = in.sweep_stat;
  else n_sweeps = in.sweep_eq;

  // Run MC simulation
  for (int ii=0; ii<n_sweeps; ii++){

    // Sweep
    sweep();

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

void sweep(){

  int ii, r_idx;
  double r_x, r_y, r_z;
  double x_old, y_old, z_old;

  // Create N trial moves (N = number of particles)
  for (ii=0; ii<part_info.NN; ii++){

    // Random number in [0, NN-1] to select the particle
    r_idx = gsl_rng_uniform_int(rng_mt, part_info.NN);

    // Three random numbers in [0,1] for the displacement
    r_x = (double)gsl_rng_get(rng_mt)/(double)r_num_max;
    r_y = (double)gsl_rng_get(rng_mt)/(double)r_num_max;
    r_z = (double)gsl_rng_get(rng_mt)/(double)r_num_max;

    // Store old coordinates (if the move gets rejected)
    x_old = part[r_idx][1];
    y_old = part[r_idx][2];
    z_old = part[r_idx][3];

    // Proposed move
    part[r_idx][1] += (r_x - 0.5)*in.dr_max;
    part[r_idx][2] += (r_y - 0.5)*in.dr_max;
    part[r_idx][3] += (r_z - 0.5)*in.dr_max;

    // Periodic boundary conditions
    if (part[r_idx][1] > sim_box_info.lx) part[r_idx][1] -= sim_box_info.lx;
    else if (part[r_idx][1] < 0.0)        part[r_idx][1] += sim_box_info.lx;
    if (part[r_idx][2] > sim_box_info.ly) part[r_idx][2] -= sim_box_info.ly;
    else if (part[r_idx][2] < 0.0)        part[r_idx][2] += sim_box_info.ly;
    if (part[r_idx][3] > sim_box_info.lz) part[r_idx][3] -= sim_box_info.lz;
    else if (part[r_idx][3] < 0.0)        part[r_idx][3] += sim_box_info.lz;
      
    // Accept or reject move according to metropolis algorithm
    if (check_overlap(r_idx, 1.0, 1.0, 1.0)){
      // Reject move
      part[r_idx][1] = x_old;
      part[r_idx][2] = y_old;
      part[r_idx][3] = z_old;
      rej_moves+=1;
    }
    else {
      // Accept move
      acc_moves+=1;
      // Update cell list
      cell_list_update();
      
    }
 
  }

}


bool check_overlap(int idx_ref,
                   double sf_x, double sf_y, double sf_z){
  
  // Variable declaration
  int cell_idx, neigh_idx, part_idx;
  double dr;

  // Cell that contains the particle
  cell_idx = cell_part_idx(idx_ref);
    
  // Loop over the neighboring cells
  for (int ii=0; ii<cl_neigh_num; ii++){

    neigh_idx = cl_neigh[cell_idx][ii];
    part_idx = cl_head[neigh_idx];

    // Loop over the particles in the neighboring cells
    while (part_idx > 0){
      
      // Compute inter-particle distance
      dr = compute_dist(idx_ref, part_idx-1, sf_x, sf_y, sf_z);

      // Signal that there is overlap
      if (dr < 1.0 && (part_idx-1) != idx_ref){
	return true;
      }

      // Update index
      part_idx = cl_link[part_idx];

    }

  }

  return false;


}


double compute_dist(int idx1, int idx2, 
		    double sf_x, double sf_y, double sf_z){

  double lx = sim_box_info.lx*sf_x;
  double ly = sim_box_info.ly*sf_y;
  double lz = sim_box_info.lz*sf_z;
  double lx_2 = lx/2.0;
  double ly_2 = ly/2.0;
  double lz_2 = lz/2.0;
  double dx, dy, dz, dr;

  // Cartesian components of the distance
  dx = (part[idx1][1] - part[idx2][1])*sf_x;
  dy = (part[idx1][2] - part[idx2][2])*sf_y;
  dz = (part[idx1][3] - part[idx2][3])*sf_z;
  
  // Periodic boundary conditions
  if (dx > lx_2)       dx -= lx;
  else if (dx < -lx_2) dx += lx;
  if (dy > ly_2)       dy -= ly;
  else if (dy < -ly_2) dy += ly;
  if (dz > lz_2)       dz -= lz;
  else if (dz< -lz_2) dz += lz;
  
  // Radial distance
  dr =  sqrt(dx*dx + dy*dy + dz*dz);
  
  return dr;

}
