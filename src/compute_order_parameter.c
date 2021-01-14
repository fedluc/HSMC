#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_legendre.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_order_parameter.h"

static double ql_ave;
static double *ql, *qlm2;
static double lx_2, ly_2, lz_2;

void compute_op(bool init,
		int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]){

  // Check if the neighbor list allows for a correct calculation of the pressure
  if (init){
    if (in.neigh_dr < in.ql_rmax){
      printf("WARNING: The cutoff for the order parameter was reduced to %f in order to be consistent with the neighbor list size\n", in.neigh_dr);
      in.ql_rmax = in.neigh_dr;
    }
  }

  // Compute the global order parameter
  global_ql_compute(cl_num_tot, cl_max_part, cl_part_cell,
		    cl_neigh_num, cl_neigh);

  // Export output
  global_ql_output(init);

}

void global_ql_compute(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		       int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]){
    

  // Allocate array for order parameter per particle
  ql = (double*)malloc(sizeof(double) * part_info.NN);
  
  // Compute the order parameter per particle
  ql_compute(cl_num_tot, cl_max_part, cl_part_cell,
		    cl_neigh_num, cl_neigh);

  // Average in order to obtain the global order parameter
  ql_ave = 0.0;
  for (int ii=0; ii<part_info.NN; ii++){
    ql_ave += ql[ii]/part_info.NN;
  }

  // Free memory
  free(ql);

}

void ql_compute(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]){

   // Variable declaration
  int tlp1 = 2*in.ql_order + 1;

  // Allocate array to store the various m-components of ql
  qlm2 = (double*)malloc(sizeof(double) * tlp1);

  // Half-box size
  lx_2 = sim_box_info.lx/2.0;
  ly_2 = sim_box_info.ly/2.0;
  lz_2 = sim_box_info.lz/2.0;  

  // Compute ql
  for (int ii=0; ii<part_info.NN; ii++){

    // Obtain all the m-components of ql
    qlm2_compute(ii,
		 cl_num_tot, cl_max_part, cl_part_cell,
		 cl_neigh_num, cl_neigh);

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


void qlm2_compute(int ref_idx,
		  int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		  int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]){
  
  int cell_idx, neigh_idx, n_part_cell, part_idx;
  int num_bonds = 0;
  double qlm_real_tmp, qlmm_real_tmp;
  double qlm_imag_tmp, qlmm_imag_tmp;
  double plm;
  double dr, dx, dy, dz;
  double phi;
  
  // Cell index for the reference particle
  cell_idx = cell_part_idx(ref_idx);
    
  // Loop over all possible values of m
  for (int mm=0; mm<in.ql_order+1; mm++){

    num_bonds = 0;
    qlm_real_tmp = 0.;
    qlm_imag_tmp = 0.;
    qlmm_real_tmp = 0.;
    qlmm_imag_tmp = 0.;
    
    // Loop over the neighbors
    for (int jj=0; jj<cl_neigh_num; jj++){

      neigh_idx = cl_neigh[cell_idx][jj];
      n_part_cell = cl_part_cell[neigh_idx][0];

      if (n_part_cell > 0){
    	for (int kk=1; kk<=n_part_cell; kk++){

    	  // Particle index
    	  part_idx = cl_part_cell[neigh_idx][kk];

    	  // Cartesian components of the distance
    	  dx = (part[ref_idx][1] - part[part_idx][1]);
    	  dy = (part[ref_idx][2] - part[part_idx][2]);
    	  dz = (part[ref_idx][3] - part[part_idx][3]);

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
    	  if (part_idx != ref_idx && dr <= in.ql_rmax){

    	    // Spherical coordinates
    	    phi = atan2(dy,dx);
    	    if(phi<0) phi += 2.*M_PI;
    	    // Update number of bonds
    	    num_bonds++;
    	    // Legendre polynomial
    	    plm =  gsl_sf_legendre_sphPlm(in.ql_order, mm, dz/dr);
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
    qlm2[in.ql_order+mm] = qlm_real_tmp*qlm_real_tmp +
                           qlm_imag_tmp*qlm_imag_tmp;
    qlm2[in.ql_order-mm] = qlmm_real_tmp*qlmm_real_tmp +
                           qlmm_imag_tmp*qlmm_imag_tmp;
    qlm2[in.ql_order+mm] = num_bonds; 
    qlm2[in.ql_order-mm] = num_bonds;
                           
    
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

