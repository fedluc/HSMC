#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "sim_info.h"
#include "rng.h"
#include "read_input.h"
#include "cell_list.h"
#include "moves.h"
#include "analytic.h"

// --------------------------------------------------------
// The module "moves.c" contains a set of functions used 
// to perform particle moves and volume moves (for NPT 
// simulations). The module contains functions to generate
// new positions (or new volumes) and to check wheter the
// resulting configuration can be accepted or not. 
// --------------------------------------------------------


// Global variables
static int part_moves=0, vol_moves=0;
static int acc_part_moves=0, rej_part_moves=0;
static int acc_vol_moves=0, rej_vol_moves=0;

// ------ Move one particle ------
void part_move(){

  int r_idx, cell_idx_old, cell_idx_new;
  double r_x, r_y, r_z;
  double x_old, y_old, z_old;

  // Get current configuration
  config part_conf = part_config_get();  

  // Random number in [0, NN-1] to select the particle
  p_info part_info = part_info_get();
  r_idx = rng_get_int(part_info.NN);
    
  // Three random numbers in [0,1] for the displacement
  r_x = rng_get_double();
  r_y = rng_get_double();
  r_z = rng_get_double();
 
  // Store old coordinates (if the move gets rejected)
  x_old = part_conf[r_idx][1];
  y_old = part_conf[r_idx][2];
  z_old = part_conf[r_idx][3];
  cell_idx_old = cell_part_idx(r_idx);

  // Proposed move
  part_conf[r_idx][1] += (r_x - 0.5)*G_IN.dr_max;
  part_conf[r_idx][2] += (r_y - 0.5)*G_IN.dr_max;
  part_conf[r_idx][3] += (r_z - 0.5)*G_IN.dr_max;

  // Periodic boundary conditions
  box_info sim_box_info = sim_box_info_get();
  if (part_conf[r_idx][1] > sim_box_info.lx) part_conf[r_idx][1] -= sim_box_info.lx;
  else if (part_conf[r_idx][1] < 0.0)        part_conf[r_idx][1] += sim_box_info.lx;
  if (part_conf[r_idx][2] > sim_box_info.ly) part_conf[r_idx][2] -= sim_box_info.ly;
  else if (part_conf[r_idx][2] < 0.0)        part_conf[r_idx][2] += sim_box_info.ly;
  if (part_conf[r_idx][3] > sim_box_info.lz) part_conf[r_idx][3] -= sim_box_info.lz;
  else if (part_conf[r_idx][3] < 0.0)        part_conf[r_idx][3] += sim_box_info.lz;
      
  // Accept or reject move according to metropolis algorithm
  if (check_overlap(r_idx, 1.0, 1.0, 1.0)){
    // Reject move
    part_conf[r_idx][1] = x_old;
    part_conf[r_idx][2] = y_old;
    part_conf[r_idx][3] = z_old;
    rej_part_moves+=1;
  }
  else {
    // Update cell list if the particle left its original box
    cell_idx_new = cell_part_idx(r_idx);
    if (cell_idx_new != cell_idx_old) {
      cell_list_update(cell_idx_old, cell_idx_new, r_idx);
    }
    // Accept move
    acc_part_moves+=1;
  }

  // Increment counter for particles moves
  part_moves += 1;
  
}

// ------ Change volume (only for NpT) ------
void vol_move(){

  double r_dv, r_acc;
  double log_vol_new, vol_new, vol_ratio, sf;
  double boltz_fact;
  bool overlap = false;

  // Random number for volume perturbation
  r_dv = rng_get_double();

  // Proposed volume move
  box_info sim_box_info = sim_box_info_get();
  log_vol_new = log(sim_box_info.vol) + (r_dv - 0.5)*G_IN.dv_max;
  vol_new = exp(log_vol_new);

  // Volume ratio
  vol_ratio = vol_new/sim_box_info.vol;

  // Scaling factor for coordinates and simulation box
  sf = pow(vol_ratio, 1./3.);
     
  // Check if there is overlap with re-scaled coordinates
  p_info part_info = part_info_get();
  for (int ii=0; ii<part_info.NN; ii++){

    overlap = check_overlap(ii, sf, sf, sf);

    if (overlap) break;

  }

  // If there is no overlap accept move according to npt acceptance rule
  if (!overlap){

    boltz_fact = exp(G_IN.press*(sim_box_info.vol - vol_new) +
		     (part_info.NN + 1)*log(vol_ratio));
    r_acc = rng_get_double();
    
    if (r_acc >= boltz_fact){
      // Reject move
      rej_vol_moves+=1;
    }
    else{
      // Accept move
      acc_vol_moves+=1;
      // Update density
      G_IN.rho = part_info.NN / vol_new;
      // Update the simulation box
      sim_box_init(G_IN.type, G_IN.nx, G_IN.ny, G_IN.nz, G_IN.rho);
      sim_box_info = sim_box_info_get();
      // Re-scale particle's positions
      config part_conf = part_config_get();
      for (int ii=0; ii<part_info.NN; ii++){
	part_conf[ii][1] *= sf;
	part_conf[ii][2] *= sf;
	part_conf[ii][3] *= sf;
	if (part_conf[ii][1] > sim_box_info.lx) part_conf[ii][1] -= sim_box_info.lx;
	else if (part_conf[ii][1] < 0.0)        part_conf[ii][1] += sim_box_info.lx;
	if (part_conf[ii][2] > sim_box_info.ly) part_conf[ii][2] -= sim_box_info.ly;
	else if (part_conf[ii][2] < 0.0)        part_conf[ii][2] += sim_box_info.ly;
	if (part_conf[ii][3] > sim_box_info.lz) part_conf[ii][3] -= sim_box_info.lz;
	else if (part_conf[ii][3] < 0.0)        part_conf[ii][3] += sim_box_info.lz;
      }
      // Construct new  cell-list (the simulation box has changed)
      cell_list_new();
    }
  }
  else{
    // Reject move
    rej_vol_moves+=1;
  }

  // Increment counter for volume moves
  vol_moves += 1;

}


// ------ Check overlap between particles ------
bool check_overlap(int idx_ref,
                   double sf_x, double sf_y, double sf_z){

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
  	dr = compute_dist(idx_ref, part_idx, sf_x, sf_y, sf_z);
	
  	// Signal that there is overlap
  	if (dr < 1.0 && part_idx != idx_ref){
  		return true;
  	}
      }
    }

  }

  /* // Loop over the other particles (no neighbor list are used) */
  /* for (int ii=0; ii<part_info.NN; ii++){ */

  /*   //Compute inter-particle distance */
  /*   dr = compute_dist(idx_ref, ii, sf_x, sf_y, sf_z); */
	
  /*   // Signal that there is overlap */
  /*   if (dr < 1.0 && ii != idx_ref){ */
  /* 		return true; */
  /*   } */

  /* } */

  return false;


}


void get_moves_counters(int *pm, int *apm, int *rpm, 
			int* vm, int *avm, int *rvm){

  if (pm != NULL) *pm = part_moves;
  if (apm != NULL) *apm = acc_part_moves;
  if (rpm != NULL) *rpm = rej_part_moves;
  if (vm != NULL) *vm = vol_moves;
  if (avm != NULL) *avm = acc_vol_moves;
  if (rvm != NULL) *rvm = rej_vol_moves;

}

void reset_moves_counters(){

  part_moves = 0;
  acc_part_moves = 0;
  rej_part_moves = 0;
  vol_moves = 0;
  acc_vol_moves = 0;
  rej_vol_moves = 0;
}

// ------ Move particle in cavity simulations ------
void cavity_part_move(){

  int r_idx, cell_idx_old, cell_idx_new;
  int move_type;
  double r_type, r_x, r_y, r_z;
  double x_old, y_old, z_old, dr_old, en_old=0.0;
  
  // Get current configuration
  config part_conf = part_config_get();

  // Select particle to move (cavities are moved with probability G_IN.cavity_pcav)
  // Note: The cavities have indexes 0 and 1
  p_info part_info = part_info_get();
  r_type = rng_get_double();  
  if (r_type > G_IN.cavity_pcav ) {
    // Move standard particles
    r_idx = rng_get_int(part_info.NN-2) + 2;
    move_type = 0;
  }
  else {
    // Move cavity
    r_idx = rng_get_int(2);
    move_type = 1;
  }

  // Three random numbers in [0,1] for the displacement
  r_x = rng_get_double();
  r_y = rng_get_double();
  r_z = rng_get_double();

  // Store old coordinates (if the move gets rejected)
  x_old = part_conf[r_idx][1];
  y_old = part_conf[r_idx][2];
  z_old = part_conf[r_idx][3];
  cell_idx_old = cell_part_idx(r_idx);

  // If a cavity has to be moved, compute old interaction energy
  if (move_type == 1){
    dr_old = compute_dist(0,1,1.0,1.0,1.0);
    if (dr_old < G_IN.cavity_mindr) {
	printf("Error, cavities are too close!, Distance: %.8e\n", dr_old);
	exit(EXIT_FAILURE);
    }
    if (dr_old > G_IN.cavity_maxdr) {
      printf("Error, cavities are too far apart!, Distance: %.8e\n", dr_old);
      exit(EXIT_FAILURE);
    }
    en_old = cavity_interaction(dr_old,false);
  }

  // Proposed move
  part_conf[r_idx][1] += (r_x - 0.5)*G_IN.dr_max;
  part_conf[r_idx][2] += (r_y - 0.5)*G_IN.dr_max;
  part_conf[r_idx][3] += (r_z - 0.5)*G_IN.dr_max;

  // Periodic boundary conditions
  box_info sim_box_info = sim_box_info_get();
  if (part_conf[r_idx][1] > sim_box_info.lx) part_conf[r_idx][1] -= sim_box_info.lx;
  else if (part_conf[r_idx][1] < 0.0)        part_conf[r_idx][1] += sim_box_info.lx;
  if (part_conf[r_idx][2] > sim_box_info.ly) part_conf[r_idx][2] -= sim_box_info.ly;
  else if (part_conf[r_idx][2] < 0.0)        part_conf[r_idx][2] += sim_box_info.ly;
  if (part_conf[r_idx][3] > sim_box_info.lz) part_conf[r_idx][3] -= sim_box_info.lz;
  else if (part_conf[r_idx][3] < 0.0)        part_conf[r_idx][3] += sim_box_info.lz;

  // Accept or reject move according to metropolis algorithm
  if (cavity_check_move(r_idx,move_type,en_old)){
    // Update cell list if the particle left its original box
    cell_idx_new = cell_part_idx(r_idx);
    if (cell_idx_new != cell_idx_old) {
      cell_list_update(cell_idx_old, cell_idx_new, r_idx);
    }
    // Accept move
    acc_part_moves+=1;
  }
  else {
    // Reject move
    part_conf[r_idx][1] = x_old;
    part_conf[r_idx][2] = y_old;
    part_conf[r_idx][3] = z_old;
    rej_part_moves+=1;
  }

  // Increment counter for particles moves
  part_moves += 1;

}


// ------ Accept or reject move in cavity simulations ------
bool cavity_check_move(int idx_ref, int move_type, double en_old){

  // Variable declaration
  int cell_idx, neigh_idx, part_idx;
  int n_part_cell;
  double dr, dr_cavity;
  double boltz_fact, r_acc;

  // If a cavity is moved, check if the move should be rejected
  if (move_type == 1) {
    dr_cavity = compute_dist(0,1,1.0,1.0,1.0);
    // Reject if the cavities are too far or too close
    if (dr_cavity < G_IN.cavity_mindr || dr_cavity > G_IN.cavity_maxdr){
      return false;
    }
    // Reject according to metropolis algorithm
    boltz_fact =  exp(-(cavity_interaction(dr_cavity,false) - en_old));
    r_acc = rng_get_double();
    if (r_acc > boltz_fact) return false;
  }


  // Check for overlaps

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

	// Reject move if an overlap between standard particles occurs
	if (move_type == 0 && part_idx != idx_ref && dr < 1.0) return false;

	// Reject move if an overlap between a cavity and the standard particles occurs
	if (move_type == 1 && part_idx >= 2 && dr < 1.0) return false;

      }
    }

  }

  // Accept move if it was not rejected
  return true;

}

// ------ Compute distance between two particles ------
double compute_dist(int idx1, int idx2,
		    double sf_x, double sf_y, double sf_z){
  
  box_info sim_box_info = sim_box_info_get();
  config part_conf = part_config_get();
  double lx = sim_box_info.lx*sf_x;
  double ly = sim_box_info.ly*sf_y;
  double lz = sim_box_info.lz*sf_z;
  double lx_2 = lx/2.0;
  double ly_2 = ly/2.0;
  double lz_2 = lz/2.0;
  double dx, dy, dz, dr;

  // Cartesian components of the distance
  dx = (part_conf[idx1][1] - part_conf[idx2][1])*sf_x;
  dy = (part_conf[idx1][2] - part_conf[idx2][2])*sf_y;
  dz = (part_conf[idx1][3] - part_conf[idx2][3])*sf_z;
  
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

// ------ Interaction potential between cavities ------

double cavity_interaction(double xx, bool cavity_init){
 
  return lny_gh(xx, M_PI*G_IN.rho/6.0, 1.0, cavity_init);
  //return 0.0;

}
