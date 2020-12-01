#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_press.h"
#include "hsmc.h"

// Variables for the calculation of the pressure via the virial
int pressv_hist_nn;
double pressv_rmax = 1.05;
double *pressv_rr, *pressv_hist;

// Variables for the calculation of the pressure via virtual
// volume perturbations
int presst_hist_nn;
double *presst_xi, *presst_hist;


void compute_pressv(bool init){

  // Initialize histogram
  pressv_hist_init();

  // Fill histogram
  pressv_compute_hist();

  // compute Radial distribution function close to contact
  pressv_compute_rdf();
  
  // Write output
  FILE* fid;
  if (init) fid = fopen("press_virial.dat", "w");
  else fid = fopen("press_virial.dat", "a");
  if (fid == NULL) {
    perror("Error while creating the file for the virial pressure");
    exit(EXIT_FAILURE);
  }
  fprintf(fid, "######################################\n");
  fprintf(fid, "# Bins, volume, number of particles\n");
  fprintf(fid, "######################################\n");
  fprintf(fid, "%d %.8e %d\n", pressv_hist_nn, 
	  sim_box_info.vol, part_info.NN);
  fprintf(fid, "###############################\n");
  fprintf(fid, "# rr, rdf\n");
  fprintf(fid, "###############################\n");
  for (int ii = 0; ii < pressv_hist_nn; ii++)
    {
      fprintf(fid, "%.8e %.8e\n", pressv_rr[ii], pressv_hist[ii]);
    }
  fclose(fid);


  // Free memory
  free(pressv_rr);
  free(pressv_hist);
  
}

void pressv_hist_init(){

  pressv_hist_nn = (int)((pressv_rmax - 1.0)/in.pressv_dr);
  pressv_rr = (double*)malloc(sizeof(double) * pressv_hist_nn);
  pressv_hist = (double*)malloc(sizeof(double) * pressv_hist_nn);
  if (pressv_rr == NULL ||  pressv_hist == NULL){
    printf("ERROR: Failed histogram allocation\n");
    exit(EXIT_FAILURE);
  }
  for (int ii=0; ii<pressv_hist_nn; ii++){
    pressv_rr[ii] = (ii+1./2.)*in.pressv_dr + 1.0;
    pressv_hist[ii] = 0.0;
  }

  
}

void pressv_compute_hist(){
  
  // Variable declaration
  int cell_idx, neigh_idx, part_idx;
  double dr;
  int bin;

  // Fill histograms
  for (int ii=0; ii<part_info.NN; ii++){

    cell_idx = cell_part_idx(ii);

    for (int jj=0; jj<cl_neigh_num; jj++){

      neigh_idx = cl_neigh[cell_idx][jj];
      part_idx = cl_head[neigh_idx];

      while (part_idx > 0){
	
	dr = compute_dist(ii, part_idx-1, 
			  1.0, 1.0, 1.0);

	if (dr <= pressv_rmax && (part_idx-1) > ii) {
	  bin = (int)((dr-1.0)/in.pressv_dr);
	  pressv_hist[bin] += 2.0;
	}

	part_idx = cl_link[part_idx];

      }
      
    }

  }

}

void pressv_compute_rdf(){

  // Variable declaration
  double r1,r2, dr, bin_vol;
  
  // Normalize histograms
  dr = pressv_rr[1] - pressv_rr[0];
  for (int ii=0; ii<pressv_hist_nn; ii++){
    r1 = pressv_rr[ii] - dr/2.;
    r2 = pressv_rr[ii] + dr/2.;
    bin_vol = (4.*M_PI/3.) * (pow(r2,3.) - pow(r1,3.));
    pressv_hist[ii] = pressv_hist[ii]/(bin_vol * in.rho * part_info.NN);
  }

}


void compute_presst(bool init){

  //Initialize histogram
  presst_hist_init();

  // Fill histogram
  presst_compute_hist();
  
  // Write output
  FILE* fid;
  if (init) fid = fopen("press_thermo.dat", "w");
  else fid = fopen("press_thermo.dat", "a");
  if (fid == NULL) {
    perror("Error while creating the file for the thermo pressure");
    exit(EXIT_FAILURE);
  }
  fprintf(fid, "######################################\n");
  fprintf(fid, "# Bins, volume, number of particles\n");
  fprintf(fid, "######################################\n");
  fprintf(fid, "%d %.8e %d\n", presst_hist_nn, 
	  sim_box_info.vol, part_info.NN);
  fprintf(fid, "######################################\n");
  fprintf(fid, "# abs(xi), exp(-beta*U)\n");
  fprintf(fid, "######################################\n");
  for (int ii = 0; ii < presst_hist_nn; ii++)
    {
      fprintf(fid, "%.8e %.8e\n", presst_xi[ii], presst_hist[ii]);
    }
  fclose(fid);


  // Free memory
  free(presst_xi);
  free(presst_hist);
  
}


void presst_hist_init(){

  presst_hist_nn = (int)(in.presst_xi_max/in.presst_dxi);
  presst_xi = (double*)malloc(sizeof(double) * presst_hist_nn);
  presst_hist = (double*)malloc(sizeof(double) * presst_hist_nn);
  if (presst_xi == NULL ||  presst_hist == NULL){
    printf("ERROR: Failed histogram allocation\n");
    exit(EXIT_FAILURE);
  }
  for (int ii=0; ii<presst_hist_nn; ii++){
    presst_xi[ii] = (ii+1)*in.presst_dxi;
    presst_hist[ii] = 0.0;
  }

  
}

void presst_compute_hist(){
  
  // Variable declaration
  int cell_idx;
  double vol_ratio, sf;
  bool overlap = false;
  
  // Fill histograms
  for (int ii=0; ii<presst_hist_nn; ii++){

    // Box length after re-scaling
    vol_ratio = 1 - presst_xi[ii];
    
    sf = pow(vol_ratio, 1./3.);
    
    // Check if there is overlap
    for (int jj=0; jj<part_info.NN; jj++){

      cell_idx = cell_part_idx(jj);
      
      overlap = check_overlap(cell_idx, sf, sf, sf);

      if (overlap) break;

    }

    // Update histogram if there is no overlap
    if (!overlap) {
      presst_hist[ii] += 1.0;
    }

  }
    
}

/* void average_hist(double *hist, int nn){ */
  
/*   for (int ii=0; ii<nn; ii++){ */
/*      hist[ii] /= in.rdf_tave; */
/*   } */

/* } */


/* void reset_hist(double *hist, int *counter, int nn){ */
  

/*   // Reset histograms */
/*   for (int ii=0; ii<nn; ii++){ */
/*     hist[ii] = 0.0; */
/*   } */

/*   // Reset counter */
/*   *counter = 0; */

/* } */




/* } */

/* void compute_pressure(double *press, double *rdf, double *rr, int nn){ */

/*   // Check how many points are available for linear regression */
/*   int ii=0, fit_pts=0; */
/*   while (rr[ii] < 1.05) { */
/*     fit_pts++; */
/*     ii++; */
/*   } */

/*   // Arrays for linear regression */
/*   double *rr_fit, *rdf_fit; */
/*   rr_fit = (double*)malloc(sizeof(double) * fit_pts); */
/*   rdf_fit = (double*)malloc(sizeof(double) * fit_pts); */
/*   if (rr_fit == NULL || rdf_fit == NULL){ */
/*     printf("ERROR: Failed  allocation for linear regression arrays\n"); */
/*     exit(EXIT_FAILURE); */
/*   } */
/*   for (ii=0; ii<fit_pts; ii++){ */
/*     rr_fit[ii] = rr[ii]; */
/*     rdf_fit[ii] = rdf[ii]; */
/*   } */

/*   // Linear regression */
/*   double c0, c1, cov00, cov01, cov11, chisq; */
/*   gsl_fit_linear(rr_fit, 1, rdf_fit, 1, fit_pts, */
/* 		 &c0, &c1, &cov00, */
/* 		 &cov01,  &cov11, &chisq); */

/*   //Pressure (ideal + excess) */
/*   *press = 1. + (2.*M_PI/3.) * in.rho * (c0 + c1); */

/*   // Free memory */
/*   free(rr_fit); */
/*   free(rdf_fit); */
  
/* } */
