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
#include "compute_rdf.h"
#include "moves.h"
#include "io_config.h"
#include "nvt.h"
#include "cluster_nvt.h"

// --------------------------------------------------------
// The module "cluster_nvt.c" is used to ...
// --------------------------------------------------------


// Global variable to track the size of the clusters
static int cluster_size;

// ------ Hard-sphere NVT simulation with cluster moves  ------
void cluster_hs_nvt() {

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

  // Print simulation info on screen
  print_sim_info();

  // Set-up the neighbor list
  cell_list_init(true);

  // Start timing
  clock_t start = clock();

  // Initialize move counters
  reset_moves_counters();

  // Run short NVT simulation to avoid rotating a perfect crystal
  for (int ii=0; ii<G_IN.cluster_init_step; ii++){
    sweep_nvt();
  }

  // Run equilibration
  printf("---------------------------------------------------\n");
  printf("Equilibration...\n");
  cluster_run_nvt(false,0);
  printf("Equilibration completed.\n");

  // Run statistics
  printf("---------------------------------------------------\n");
  printf("Production...\n");
  cluster_run_nvt(true,G_IN.sweep_eq);
  printf("Production completed.\n");
  clock_t end = clock();
  
  // Stop timing
  printf("Elapsed time: %f seconds\n",
  	 (double)(end - start) / CLOCKS_PER_SEC);

  // Free memory
  part_free();
  cell_list_free();
  rng_free();


}


void cluster_run_nvt(bool prod_flag, int sweep_offset){

  // Variable declaration
  bool pressv_init = true, presst_init = true;
  bool ql_init = true, mu_ave_init = true;
  bool rdf_init = true;
  int n_sweeps;

  // Number of sweeps
  if (prod_flag) n_sweeps = G_IN.sweep_stat;
  else n_sweeps = G_IN.sweep_eq;
  
  // Run MC simulation
  for (int ii=sweep_offset; ii<n_sweeps+sweep_offset; ii++){

    // Output on screen
    if (ii == sweep_offset){
      printf("Sweep number, cluster size\n ");
    }
    if (ii % G_IN.output_int == 0) {
      printf("%d %d\n", ii, cluster_size);
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
      	  compute_op(ql_init);
      	  if (ql_init) ql_init = false;
      	}
      }

      // Compute chemical potential via Widom insertions
      if (G_IN.mu_sample_int > 0){
      	if (ii % G_IN.mu_sample_int == 0) {
      	  compute_mu(mu_ave_init);
      	  if (mu_ave_init) mu_ave_init = false;
      	}
      }

      // Compute radial distribution function
      if (G_IN.rdf_sample_int > 0){
      	if (ii % G_IN.rdf_sample_int == 0) {
      	  compute_rdf(rdf_init,ii);
      	  if (rdf_init) rdf_init = false;
      	}
      }

    }

    // Generate new configuration
    cluster_sweep_nvt();

  }

  // Free memory from compute
  if (!pressv_init) pressv_hist_free();
  if (!presst_init) presst_hist_free();
  if (!ql_init) ql_free();
  if (!rdf_init) rdf_hist_free();

}

void cluster_sweep_nvt(){

  // Create N trial moves (N = number of particles)
  for (int ii=0; ii<G_IN.cluster_moves_sweep; ii++){
    cluster_move();
  }

}


// ------ Cluster moves ------

void cluster_move(){

  int r_idx;
  double r_x_p1, r_y_p1, r_z_p1;
  double uu, vv, ww;
  
  // Select random pivot for the rotation
  box_info sim_box_info = sim_box_info_get();
  r_x_p1 = rng_get_double()*sim_box_info.lx;
  r_y_p1 = rng_get_double()*sim_box_info.ly;
  r_z_p1 = rng_get_double()*sim_box_info.lz;

  // Select rotation axis
  r_idx = rng_get_int(3);
  if (r_idx == 0){
    uu = 1.0;
    vv = 0.0;
    ww = 0.0;
  }
  else if (r_idx == 1){
    uu = 0.0;
    vv = 1.0;
    ww = 0.0;
  }
  else {
    uu = 0.0;
    vv = 0.0;
    ww = 1.0;
  }
 
  // Move one particle
  p_info part_info = part_info_get();
  r_idx = rng_get_int(part_info.NN);
  cluster_move_part(r_idx, r_x_p1, r_y_p1, r_z_p1,
		    uu, vv, ww);
  
  // --- Pocket algorithm to resolve overlaps ---

  // Allocate pocket
  int *pocket = malloc( sizeof(int) * part_info.NN);

  // Initialize pocket
  for (int ii=0; ii<part_info.NN; ii++) pocket[ii] = -1;
  int part_pocket = 0;

  // Update pocket
  pocket_add(r_idx, pocket, &part_pocket);

  // Resolve overlaps
  cluster_size = 1;
  while (part_pocket > 0){
    
    r_idx = pocket[rng_get_int(part_pocket)];
    cluster_move_part(r_idx, r_x_p1, r_y_p1, r_z_p1,
  		      uu, vv, ww);
    pocket_del(r_idx, pocket, &part_pocket);
    pocket_add(r_idx, pocket, &part_pocket);
    cluster_size++;

  }

  // Free pocket
  free(pocket);

}

void cluster_move_part(int idx, double aa, double bb, double cc,
		       double uu, double vv, double ww){

  int cell_idx_old, cell_idx_new;
  double x_old, y_old, z_old;  
    
  config part_conf = part_config_get();
  
  // Get information of the particle to move
  x_old = part_conf[idx][1];
  y_old = part_conf[idx][2];
  z_old = part_conf[idx][3];
  cell_idx_old = cell_part_idx(idx);

  // Rotate particle of pi around the rotation axis
  rotate(x_old, y_old, z_old, aa, bb, cc, 
	 uu, vv, ww, M_PI, &part_conf[idx][1],
	 &part_conf[idx][2], &part_conf[idx][3]);

  // Periodic boundary conditions
  apply_pbc(idx);

  // Update cell lists
  cell_idx_new = cell_part_idx(idx);
  if (cell_idx_new != cell_idx_old) {
     cell_list_update(cell_idx_old, cell_idx_new, idx);
  }

}

// ----- Rotate point (x,y,z) around axis passing by (a,b,c) with direction (u,v,w) ------
void rotate(double xx, double yy, double zz, 
	    double aa, double bb, double cc, 
	    double uu, double vv, double ww,
	    double theta, double *xx_rot, 
	    double *yy_rot, double *zz_rot){

  
  double cost = cos(theta), sint = sin(theta);
  double uu2 = uu*uu, vv2 = vv*vv, ww2 = ww*ww;

  *xx_rot = (aa*(vv2 + ww2) - uu*(bb*vv + cc*ww - uu*xx - vv*yy - ww*zz))
    *(1 - cost) + xx*cost + (-cc*vv + bb*ww - ww*yy + vv*zz)*sint;

  *yy_rot = (bb*(uu2 + ww2) - vv*(aa*uu + cc*ww - uu*xx - vv*yy - ww*zz))
    *(1 - cost) + yy*cost + (cc*uu - aa*ww + ww*xx - uu*zz)*sint;

  *zz_rot = (cc*(uu2 + vv2) - ww*(aa*uu + bb*vv - uu*xx - vv*yy - ww*zz))
    *(1 - cost) + zz*cost + (-bb*uu + aa*vv - vv*xx + uu*yy)*sint;


}

// ------ Manage the pocket ------

void pocket_add(int idx_ref, int *pocket, int *part_pocket){

  // Variable declaration
  int cell_idx, neigh_idx, part_idx;
  int n_part_cell;
  double dr;

  // Neighbor list
  cl_info nl = get_cell_list_info();

  // Cell that contains the particle
  cell_idx = cell_part_idx(idx_ref);

  // Loop over the neighboring cells
  for (int ii=0; ii<nl.neigh_num; ii++){

    neigh_idx = nl.neigh_mat[cell_idx][ii];

    // Loop over the particles in the neighboring cell
    n_part_cell = nl.part_cell[neigh_idx][0];
    if (n_part_cell > 0){
      for (int jj=1; jj<=n_part_cell; jj++){
        
	// Particle index
        part_idx = nl.part_cell[neigh_idx][jj];

        //Compute inter-particle distance
        dr = compute_dist(idx_ref, part_idx, 1.0, 1.0, 1.0);

        // Signal that there is overlap
        if (dr < 1.0 && part_idx != idx_ref){
	  if (!in_pocket(part_idx, pocket, *part_pocket)){
	    pocket[*part_pocket] = part_idx;
	    *part_pocket += 1;
	  }
	}

      }

    }

  }
  
}

void pocket_del(int idx_ref, int *pocket, int *part_pocket){
    
  int idx_remove = *part_pocket;
  bool shift_flag = false;

  for (int ii=0; ii<*part_pocket; ii++){
    if (pocket[ii] == idx_ref){
      idx_remove = ii;
      pocket[ii] = -1;
      shift_flag = true;
    }
    if (shift_flag && ii > idx_remove) {
      pocket[ii-1] = pocket[ii];
      pocket[ii] = -1;
    }
  }
  *part_pocket -= 1;

}

bool in_pocket(int idx_ref, int *pocket, int part_pocket){

  for (int ii=0; ii<part_pocket; ii++){
    if (pocket[ii] == idx_ref){
      return true;
    }
  }
  
  return false;

}
