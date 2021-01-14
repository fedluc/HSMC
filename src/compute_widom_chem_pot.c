#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "init.h"
#include "rng.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_widom_chem_pot.h"
#include "moves.h"

// Variables for the random positions generated during the insertions
static double r_x, r_y, r_z;
static double mu;
static int wtest;

// Variables to assign the virtual particles to the cell lists
static int cell_num_x;
static int cell_num_y;
static int cell_num_z;
static double cell_size_x;
static double cell_size_y;
static double cell_size_z;

void compute_mu(bool init,
		int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
                int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]){

  // Compute the chemical potential via widom insertion method
  widom_insertion(cl_num_tot, cl_max_part, cl_part_cell,
		  cl_neigh_num, cl_neigh);

  // Export output
  mu_output(init);
 
}

void widom_insertion(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		     int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]){

  // Initialize insertion counter
  wtest = 0;

  // Size of the cell list
  get_cell_list_info(NULL, NULL, &cell_num_x, &cell_num_y, &cell_num_z,
		     &cell_size_x, &cell_size_y, &cell_size_z);

  // Perform the insertions prescribed in input
  for (int ii=0; ii<in.mu_insertions; ii++){

    // Generate one random position inside the simulation box
    widom_rand_pos();

    // Update average chemical potential if there is no overlap
    if (!widom_check_overlap(cl_num_tot, cl_max_part, cl_part_cell,
			     cl_neigh_num, cl_neigh)) wtest++;

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

    r_x = rng_get_double() * sim_box_info.lx;
    r_y = rng_get_double() * sim_box_info.ly;
    r_z = rng_get_double() * sim_box_info.lz;
  
}


bool widom_check_overlap(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
			 int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]){
  
  // Variable declaration
  int cell_idx, neigh_idx, part_idx;
  int n_part_cell;
  double dr;

  // Cell that contains the virtual particle
  cell_idx = widom_cell_part_idx();

  // Loop over the neighboring cells
  for (int ii=0; ii<cl_neigh_num; ii++){

    neigh_idx = cl_neigh[cell_idx][ii];

    // Loop over the particles in the neighboring cell
    n_part_cell = cl_part_cell[neigh_idx][0];
    if (n_part_cell > 0){
      for (int jj=1; jj<=n_part_cell; jj++){

        // Particle index
        part_idx = cl_part_cell[neigh_idx][jj];

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

  return  (int)(r_x/cell_size_x)*cell_num_x*cell_num_x
    + (int)(r_y/cell_size_y)*cell_num_y
    + (int)(r_z/cell_size_z);

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
