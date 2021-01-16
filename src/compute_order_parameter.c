#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_legendre.h>
#include "sim_info.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_order_parameter.h"

// ------ Global variables ------
static double ql_ave;
static double *ql, *qlm2;
static double lx_2, ly_2, lz_2;

// ------ Caller for the other functions in the file ------

void compute_op(bool init){


  // Check if the neighbor list allows for a correct calculation of the order parameter
  if (init){
    
    // Get info on the neighbor list
    cl_info nl = get_cell_list_info();
    double nl_size_min = nl.size_x;
    if (nl_size_min < nl.size_y) nl_size_min = nl.size_y;
    if (nl_size_min < nl.size_z) nl_size_min = nl.size_z;
    if (nl_size_min < G_IN.ql_rmax){
      printf("WARNING: The cutoff for the order parameter was reduced to %f in order to be consistent with the neighbor list size\n", nl_size_min);
      G_IN.ql_rmax = nl_size_min;
    }

  }

  // Compute the global order parameter
  global_ql_compute();

  // Export output
  global_ql_output(init);

}

// ------ Compute the average (over all the particles) order paramter ------

void global_ql_compute(){
    

  // Allocate array for order parameter per particle
  p_info part_info = part_info_get();
  ql = (double*)malloc(sizeof(double) * part_info.NN);
  
  // Compute the order parameter per particle
  ql_compute();

  // Average in order to obtain the global order parameter
  ql_ave = 0.0;
  for (int ii=0; ii<part_info.NN; ii++){
    ql_ave += ql[ii]/part_info.NN;
  }

  // Free memory
  free(ql);

}


// ------ Compute the ql order parameter for each particle ------

void ql_compute(){

   // Variable declaration
  int tlp1 = 2*G_IN.ql_order + 1;

  // Allocate array to store the various m-components of ql
  qlm2 = (double*)malloc(sizeof(double) * tlp1);

  // Half-box size
  box_info sim_box_info = sim_box_info_get();
  lx_2 = sim_box_info.lx/2.0;
  ly_2 = sim_box_info.ly/2.0;
  lz_2 = sim_box_info.lz/2.0;  

  // Compute ql
  p_info part_info = part_info_get();
  for (int ii=0; ii<part_info.NN; ii++){

    // Obtain all the m-components of ql
    qlm2_compute(ii);

    ql[ii] = 0.0;
    for (int jj=0; jj<tlp1; jj++){
      ql[ii] += qlm2[jj];
    }
    ql[ii] *= 4 * M_PI / tlp1;
    ql[ii] = sqrt(ql[ii]);

  }

  // Free memory
  free(qlm2);

}


// ------ Compute the qlm order parameter for each particle ------

void qlm2_compute(int ref_idx){
  
  int num_bonds = 0;
  int cell_idx, neigh_idx, part_idx;
  int n_part_cell;
  double qlm_real_tmp, qlmm_real_tmp;
  double qlm_imag_tmp, qlmm_imag_tmp;
  double plm;
  double dr, dx, dy, dz;
  double phi;

  // Neighbor list
  cl_info nl = get_cell_list_info();

  // Simulation box
  box_info sim_box_info = sim_box_info_get();

  // Configuration
  config part_conf = part_config_get();
  
  // Cell with the reference particle
  cell_idx = cell_part_idx(ref_idx);
 
  // Loop over all possible values of m
  for (int mm=0; mm<G_IN.ql_order+1; mm++){

    num_bonds = 0;
    qlm_real_tmp = 0.;
    qlm_imag_tmp = 0.;
    qlmm_real_tmp = 0.;
    qlmm_imag_tmp = 0.;

    // Loop over the neighboring cells
    for (int ii=0; ii<nl.neigh_num; ii++){

      neigh_idx = nl.neigh_mat[cell_idx][ii];

      // Loop over the particles in the neighboring cell
      n_part_cell = nl.part_cell[neigh_idx][0];
      if (n_part_cell > 0){
    	for (int jj=1; jj<=n_part_cell; jj++){

    	  // Particle index
    	  part_idx = nl.part_cell[neigh_idx][jj];

    	  // Cartesian components of the distance
    	  dx = (part_conf[ref_idx][1] - part_conf[part_idx][1]);
    	  dy = (part_conf[ref_idx][2] - part_conf[part_idx][2]);
    	  dz = (part_conf[ref_idx][3] - part_conf[part_idx][3]);
      
    	  // Periodic boundary conditions
    	  if (dx > lx_2)       dx -= sim_box_info.lx;
    	  else if (dx < -lx_2) dx += sim_box_info.lx;
    	  if (dy > ly_2)       dy -= sim_box_info.ly;
    	  else if (dy < -ly_2) dy += sim_box_info.ly;
    	  if (dz > lz_2)       dz -= sim_box_info.lz;
    	  else if (dz < -lz_2) dz += sim_box_info.lz;
	  
    	  // Radial distance
    	  dr = sqrt(dx*dx + dy*dy + dz*dz);
	  
    	  // Update count if the radial distance falls within the cutoff
    	  if (part_idx != ref_idx && dr <= G_IN.ql_rmax){
	    
    	    // Spherical coordinates
    	    phi = atan2(dy,dx);
    	    if(phi<0) phi += 2.*M_PI;
    	    // Update number of bonds
    	    num_bonds++;
    	    // Legendre polynomial
    	    plm =  gsl_sf_legendre_sphPlm(G_IN.ql_order, mm, dz/dr);
    	    // Spherical harmonics
    	    qlm_real_tmp += plm * cos(mm*phi);
    	    qlm_imag_tmp += plm * sin(mm*phi);
    	    if (mm % 2 != 0) plm *= -1.0;
    	    qlmm_real_tmp += plm * cos(-mm*phi);
    	    qlmm_imag_tmp += plm * sin(-mm*phi);
	    
    	  }
    	}
      }
    }
     
    if (num_bonds != 0) {
      qlm_real_tmp /= num_bonds;
      qlm_imag_tmp /= num_bonds;
      qlmm_real_tmp /= num_bonds;
      qlmm_imag_tmp /= num_bonds;
    }
    qlm2[G_IN.ql_order+mm] = qlm_real_tmp*qlm_real_tmp +
                           qlm_imag_tmp*qlm_imag_tmp;
    qlm2[G_IN.ql_order-mm] = qlmm_real_tmp*qlmm_real_tmp +
                           qlmm_imag_tmp*qlmm_imag_tmp;
                           
    
  }
  
}
 

// ------ Print output to file ------

void global_ql_output(bool init){

  // Print to file the global order parameter 
  FILE *fid;
  if (init) fid = fopen("order_param.dat", "w");
  else fid = fopen("order_param.dat", "a");
  if (fid == NULL) {
    perror("Error while creating the file for the order parameter\n");
    exit(EXIT_FAILURE);
  }
  if (init){
    fprintf(fid, "###############################################################\n");
    fprintf(fid, "# Average order parameter of order %d (each line is one sample)\n", G_IN.ql_order);
    fprintf(fid, "###############################################################\n");
  }
  fprintf(fid, "%.8e\n", ql_ave);
  fclose(fid);

}

