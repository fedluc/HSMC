#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_legendre.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_order_parameter.h"

double ql_ave = 0.0;
double *ql, *qlm2;
double lx_2, ly_2, lz_2;

void compute_op(bool init){
  
  // Compute the global order parameter
  global_ql_compute();

  // Export output
  global_ql_output(init);

}

void global_ql_compute(){

  // Allocate array for order parameter per particle
  ql = (double*)malloc(sizeof(double) * part_info.NN);
  
  // Compute the order parameter per particle
  ql_compute();

  // Average in order to obtain the global order parameter
  for (int ii=0; ii<part_info.NN; ii++){
    ql_ave += ql[ii]/part_info.NN;
  }

  // Free memory
  free(ql);

}

void ql_compute(){

   // Variable declaration
  int tlp1 = 2*in.ql_order + 1;

  // Allocate array to store the various m-components of ql
  qlm2 = (double*)malloc(sizeof(double) * tlp1);

  // Check for consistency 
  // (The cell lists limit the number of neighbor interactions 
  // that can be compute)  
  if (in.ql_rmax > 1.5*sim_box_info.cell_size) {
    printf("WARNING: ql_rmax reduced to %f in order to be consistent with cell list", 1.5*sim_box_info.cell_size); 
    in.ql_rmax = 1.5*sim_box_info.cell_size;
  }

  // Half-box size
  lx_2 = sim_box_info.lx/2.0;
  ly_2 = sim_box_info.ly/2.0;
  lz_2 = sim_box_info.lz/2.0;  

  // Compute ql
  for (int ii=0; ii<part_info.NN; ii++){

    // Obtain all the m-components of ql
    qlm_compute(ii);

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


void qlm_compute(int ref_idx){
  
  int cell_idx, neigh_idx, part_idx;
  int num_bonds = 0;
  double qlm_real_tmp;
  double qlm_imag_tmp;
  double plm;
  double dr, dx, dy, dz;
  double phi;
  
  // Cell index for the reference particle
  cell_idx = cell_part_idx(ref_idx);
    
  // Loop over all possible values of m
  for (int mm=-in.ql_order; mm<in.ql_order+1; mm++){

    num_bonds = 0;
    qlm_real_tmp = 0;
    qlm_imag_tmp = 0;
    
    // Loop over the neighbors
    for (int jj=0; jj<cl_neigh_num; jj++){

      neigh_idx = cl_neigh[cell_idx][jj];
      part_idx = cl_head[neigh_idx];

      while (part_idx > 0){
	
	// Cartesian components of the distance 
	dx = (part[ref_idx][1] - part[part_idx-1][1]);
	dy = (part[ref_idx][2] - part[part_idx-1][2]);
	dz = (part[ref_idx][3] - part[part_idx-1][3]);
	
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
	if (part_idx-1 != ref_idx && dr <= in.ql_rmax){
	  
	  // Spherical coordinates
	  phi = atan2(dy,dx);
	  if(phi<0) phi += 2.*M_PI;
	  // Update number of bonds
	  num_bonds++;
	  // Legendre polynomial
	  plm =  gsl_sf_legendre_sphPlm(in.ql_order, abs(mm), dz/dr);
	  if (mm < 0 && mm % 2 != 0) plm *= -1.0;
	  // Spherical harmonics
	  qlm_real_tmp += plm * cos(mm*phi);
	  qlm_imag_tmp += plm * sin(mm*phi);
	}
	
	part_idx = cl_link[part_idx];
	
      }
      
    }
    
    qlm_real_tmp /= num_bonds;
    qlm_imag_tmp /= num_bonds;
    qlm2[mm+in.ql_order] = qlm_real_tmp*qlm_real_tmp + 
                           qlm_imag_tmp*qlm_imag_tmp;
    
  }
  

}
 


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
    fprintf(fid, "# Average order parameter of order %d (each line is one sample)\n", in.ql_order);
    fprintf(fid, "###############################################################\n");
  }
  fprintf(fid, "%.8e\n", ql_ave);
  fclose(fid);

}

