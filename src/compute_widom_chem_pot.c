#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sim_info.h"
#include "rng.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_widom_chem_pot.h"
#include "moves.h"

// --------------------------------------------------------
// The module "compute_widom_chem_pot.c" is used to
// compute the chemical potential via the Widom insertion 
// method. For a discussion of the method see the book
// "Understanding Molecular Simulation: From Algorithms to 
// Applications" by B. Smit and D. Frenkel
// --------------------------------------------------------


// ------ Global variables ------

// Variables for the random positions generated during the insertions
static double r_x, r_y, r_z;
static double mu;
static int wtest;

// Variables to assign the virtual particles to the cell lists
static cl_info nl;

// ------ Caller to compute the chemical potential ------

void compute_mu(bool init){

  // Compute the chemical potential via widom insertion method
  widom_insertion();

  // Export output
  mu_output(init);
 
}

// ----- Functions to perform the Widom insertion ------

void widom_insertion(){

  // Initialize insertion counter
  wtest = 0;

  // Neighbor list information
  nl = get_cell_list_info();

  // Perform the insertions prescribed in input
  for (int ii=0; ii<G_IN.mu_insertions; ii++){

    // Generate one random position inside the simulation box
    widom_rand_pos();

    // Update average chemical potential if there is no overlap
    if (!widom_check_overlap()) wtest++;

  }

  // Chemical potential
  if (wtest > 0){
    mu = -log((double)wtest/G_IN.mu_insertions);
  }
  else {
    mu = 0.0; // No insertion was accepted
  }

}

void widom_rand_pos(){
  
  box_info sim_box_info = sim_box_info_get();
  r_x = rng_get_double() * sim_box_info.lx;
  r_y = rng_get_double() * sim_box_info.ly;
  r_z = rng_get_double() * sim_box_info.lz;
  
}

bool widom_check_overlap(){
  
  // Variable declaration
  int cell_idx, neigh_idx, part_idx;
  int n_part_cell;
  double dr;

  // Cell that contains the virtual particle
  cell_idx = widom_cell_part_idx();

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
        dr = widom_compute_dist(part_idx);

        // Signal that there is overlap
        if (dr < 1.0){
	  return true;
        }
      }
    }
  }


  return false;


}

int widom_cell_part_idx(){

  return  (int)(r_x/nl.size_x)*nl.num_x*nl.num_x
    + (int)(r_y/nl.size_y)*nl.num_y
    + (int)(r_z/nl.size_z);

}

double widom_compute_dist(int idx){

  box_info sim_box_info = sim_box_info_get();
  config part_conf = part_config_get();
  double lx = sim_box_info.lx;
  double ly = sim_box_info.ly;
  double lz = sim_box_info.lz;
  double lx_2 = lx/2.0;
  double ly_2 = ly/2.0;
  double lz_2 = lz/2.0;
  double dx, dy, dz, dr;

  // Cartesian components of the distance
  dx = r_x - part_conf[idx][1];
  dy = r_y - part_conf[idx][2];
  dz = r_z - part_conf[idx][3];
  
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

// ------- Write chemical potential to file ------

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
    fprintf(fid, "# Chemical potenital (average over %d insertions, Fraction of accepted insertions)\n", G_IN.mu_insertions);
    fprintf(fid, "##################################################################################\n");
  }
  fprintf(fid, "%.8e %.8e\n", mu, (double)wtest/G_IN.mu_insertions);
  fclose(fid);

}
