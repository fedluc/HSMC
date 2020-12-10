#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_legendre.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_order_parameter.h"


void compute_ql_ave(bool init){

  // Variable declaration
  int cell_idx, neigh_idx, part_idx;
  int ll = 6, tlp1 = 2*ll + 1;
  double neigh_cutoff = 1.5; // Rather arbitrary
  double *ql = (double*)malloc(sizeof(double) * part_info.NN);
  double *qlm2 = (double*)malloc(sizeof(double) * tlp1);
  double ql_ave;
  int num_bonds = 0;
  double qlm_real_tmp;
  double qlm_imag_tmp;
  double plm;
  double dr, dx, dy, dz;
  double phi;
  double lx_2 = sim_box_info.lx/2.0;
  double ly_2 = sim_box_info.ly/2.0;
  double lz_2 = sim_box_info.lz/2.0;  
  // Compute ql
  for (int ii=0; ii<part_info.NN; ii++){

    cell_idx = cell_part_idx(ii);

    // Compute qlm
    for (int mm=-ll; mm<ll+1; mm++){
       
      num_bonds = 0;
      qlm_real_tmp = 0;
      qlm_imag_tmp = 0;

      for (int jj=0; jj<cl_neigh_num; jj++){

	neigh_idx = cl_neigh[cell_idx][jj];
	part_idx = cl_head[neigh_idx];

	while (part_idx > 0){
	
	  // Cartesian components of the distance 
	  dx = (part[ii][1] - part[part_idx-1][1]);
	  dy = (part[ii][2] - part[part_idx-1][2]);
	  dz = (part[ii][3] - part[part_idx-1][3]);
	
	  // Periodic boundary conditions
	  if (dx > lx_2)       dx -= sim_box_info.lx;
	  else if (dx < -lx_2) dx += sim_box_info.lx;
	  if (dy > ly_2)       dy -= sim_box_info.ly;
	  else if (dy < -ly_2) dy += sim_box_info.ly;
	  if (dz > lz_2)       dz -= sim_box_info.lz;
	  else if (dz < -lz_2) dz += sim_box_info.lz;
	  
	  // Radial distance
	  dr = sqrt(dx*dx + dy*dy + dz*dz);
	  
	  // ...
	  if (part_idx-1 != ii && dr <= neigh_cutoff){
	    
	    phi = atan2(dy,dx);
	    if(phi<0) phi += 2.*M_PI;
	    num_bonds++;
	    plm =  gsl_sf_legendre_sphPlm(ll, abs(mm), dz/dr);
	    if (mm < 0 && mm % 2 != 0) plm *= -1.0;
	    qlm_real_tmp += plm * cos(mm*phi);
	    qlm_imag_tmp += plm * sin(mm*phi);
	  }
	  
	  part_idx = cl_link[part_idx];

	}
	
      }
      
      qlm_real_tmp /= num_bonds;
      qlm_imag_tmp /= num_bonds;
      qlm2[mm+ll] = qlm_real_tmp*qlm_real_tmp + 
                    qlm_imag_tmp*qlm_imag_tmp;

    }

    ql[ii] = 0.0;
    for (int jj=0; jj<tlp1; jj++){
      ql[ii] += qlm2[jj]; 
    }
    ql[ii] *= 4 * M_PI / tlp1;
    ql[ii] = sqrt(ql[ii]);

  }

  // Average ql
  ql_ave = 0.0;
  for (int ii=0; ii<part_info.NN; ii++){
    ql_ave += ql[ii]/part_info.NN;
  }

  // For NpT calculations print to file also the density
  FILE *fid;
  if (init) fid = fopen("ql6.dat", "w");
  else fid = fopen("ql6.dat", "a");
  if (fid == NULL) {
    perror("Error while creating the file for the order parameter\n");
    exit(EXIT_FAILURE);
  }
  if (init){
    fprintf(fid, "###################################################\n");
    fprintf(fid, "# Average order parameter (each line is one sample)\n");
    fprintf(fid, "###################################################\n");
  }
  fprintf(fid, "%.8e\n", ql_ave);
  fclose(fid);

  
  // Free memory
  free(ql);
  free(qlm2);

}

