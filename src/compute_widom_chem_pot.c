#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_widom_chem_pot.h"
#include "moves.h"

// Variables for the random positions generated during the insertions
static double r_x, r_y, r_z;
static double mu;
static int wtest;

void compute_mu(bool init){

  // Compute the chemical potential via widom insertion method
  widom_insertion();

  // Export output
  mu_output(init);
 
}

void widom_insertion(){

  // Initialize insertion counter
  wtest = 0;

  // Perform the insertions prescribed in input
  for (int ii=0; ii<in.mu_insertions; ii++){

    // Generate one random position inside the simulation box
    widom_rand_pos();

    // Update average chemical potential if there is no overlap
    if (!widom_check_overlap()) wtest++;

  }

  // Chemical potential
  if (wtest > 0){
    mu = -log((double)wtest/in.mu_insertions);
  }
  else {
    mu = 0.0; // No insertion was accepted
  }

}


void widom_rand_pos(){

    r_x = (double)gsl_rng_get(rng_mt)/(double)r_num_max;
    r_y = (double)gsl_rng_get(rng_mt)/(double)r_num_max;
    r_z = (double)gsl_rng_get(rng_mt)/(double)r_num_max;
    r_x *= sim_box_info.lx;
    r_y *= sim_box_info.ly;
    r_z *= sim_box_info.lz;
  
}


bool widom_check_overlap(){
  
  // Variable declaration
  int cell_idx, neigh_idx, part_idx;
  double dr;

  // Cell that contains the virtual particle
  cell_idx = widom_cell_part_idx();
    
  // Loop over the neighboring cells
  for (int ii=0; ii<cl_neigh_num; ii++){

    neigh_idx = cl_neigh[cell_idx][ii];
    part_idx = cl_head[neigh_idx];

    // Loop over the particles in the neighboring cells
    while (part_idx > 0){
      
      // Compute inter-particle distance
      dr = widom_compute_dist(part_idx-1);

      // Signal that there is overlap
      if (dr < 1.0){
	return true;
      }

      // Update index
      part_idx = cl_link[part_idx];

    }

  }

  return false;


}

int widom_cell_part_idx(){

  return  (int)(r_x/sim_box_info.cell_size)*sim_box_info.cell_x*sim_box_info.cell_x
    + (int)(r_y/sim_box_info.cell_size)*sim_box_info.cell_y
    + (int)(r_z/sim_box_info.cell_size);

}

double widom_compute_dist(int idx){

  double lx = sim_box_info.lx;
  double ly = sim_box_info.ly;
  double lz = sim_box_info.lz;
  double lx_2 = lx/2.0;
  double ly_2 = ly/2.0;
  double lz_2 = lz/2.0;
  double dx, dy, dz, dr;

  // Cartesian components of the distance
  dx = r_x - part[idx][1];
  dy = r_y - part[idx][2];
  dz = r_z - part[idx][3];
  
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

void mu_output(bool init){

  // Print to file the global order parameter 
  FILE *fid;
  if (init) fid = fopen("chem_pot.dat", "w");
  else fid = fopen("chem_pot.dat", "a");
  if (fid == NULL) {
    perror("Error while creating the file for the chemical potential\n");
    exit(EXIT_FAILURE);
  }
  if (init){
    fprintf(fid, "##################################################################################\n");
    fprintf(fid, "# Chemical potenital (average over %d insertions, Fraction of accepted insertions)\n", in.mu_insertions);
    fprintf(fid, "##################################################################################\n");
  }
  fprintf(fid, "%.8e %.8e\n", mu, (double)wtest/in.mu_insertions);
  fclose(fid);

}
