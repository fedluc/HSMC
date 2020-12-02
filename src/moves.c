#include <math.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "moves.h"

// Global variables for random number generator
gsl_rng *rng_mt;
long unsigned int r_num_max;

// Global variables for particles moves
int part_moves, vol_moves;
int acc_part_moves, rej_part_moves;
int acc_vol_moves, rej_vol_moves;


void rng_init(){

  // Set-up random number generator (Marsenne-Twister)
  rng_mt = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng_mt,0);
  r_num_max = gsl_rng_max(rng_mt);

}

void part_move(){

  int r_idx;
  double r_x, r_y, r_z;
  double x_old, y_old, z_old;

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
    rej_part_moves+=1;
  }
  else {
    // Accept move
    acc_part_moves+=1;
    // Update cell list
    cell_list_update();      
  }

}

void vol_move(){

  double r_dv, r_acc;
  double log_vol_new, vol_new, vol_ratio, sf;
  double boltz_fact;
  bool overlap = false;

  // Random number for volume perturbation
  r_dv = (double)gsl_rng_get(rng_mt)/(double)r_num_max;

  // Proposed volume move
  log_vol_new = log(sim_box_info.vol) + (r_dv - 0.5)*in.npt_dv_max;
  vol_new = exp(log_vol_new);

  // Volume ratio
  vol_ratio = vol_new/sim_box_info.vol;

  // Scaling factor for coordinates and simulation box
  sf = pow(vol_ratio, 1./3.);
     
  // Check if there is overlap with re-scaled coordinates
  for (int ii=0; ii<part_info.NN; ii++){

    overlap = check_overlap(ii, sf, sf, sf);

    if (overlap) break;

  }

  // If there is no overlap accept move according to npt acceptance rule
  if (!overlap){

    boltz_fact = exp(in.npt_press*(sim_box_info.vol - vol_new) + 
		     (part_info.NN + 1)*log(vol_ratio));
    r_acc = (double)gsl_rng_get(rng_mt)/(double)r_num_max;
    
    if (r_acc >= boltz_fact){
      // Reject move
      rej_vol_moves+=1;
    }
    else{
      // Accept move
      acc_vol_moves+=1;
      // Update density
      in.rho = part_info.NN / vol_new;
      // Update the simulation box
      sim_box_init(in.type, in.nx, in.ny, in.nz, in.rho);
      // Re-scale particle's positions
      for (int ii=0; ii<part_info.NN; ii++){
	part[ii][1] *= sf;
	part[ii][2] *= sf;
	part[ii][3] *= sf;
	if (part[ii][1] > sim_box_info.lx) part[ii][1] -= sim_box_info.lx;
	else if (part[ii][1] < 0.0)        part[ii][1] += sim_box_info.lx;
	if (part[ii][2] > sim_box_info.ly) part[ii][2] -= sim_box_info.ly;
	else if (part[ii][2] < 0.0)        part[ii][2] += sim_box_info.ly;
	if (part[ii][3] > sim_box_info.lz) part[ii][3] -= sim_box_info.lz;
	else if (part[ii][3] < 0.0)        part[ii][3] += sim_box_info.lz;
      }
      // Update cell-list (not sure if necessary)
      cell_list_update();
    }
  }
  else{
    // Reject move
    rej_vol_moves+=1;
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
