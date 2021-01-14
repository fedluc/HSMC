#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_press.h"
#include "moves.h"

// Variables for the calculation of the pressure via the virial
static int pressv_hist_nn;
static double pressv_rmax = 1.05;
static double *pressv_rr, *pressv_hist;

// Variables for the calculation of the pressure via virtual
// volume perturbations
static int presst_hist_nn;
static double *presst_xi, *presst_hist;


void compute_pressv(bool init){

  // Check if the neighbor list allows for a correct calculation of the pressure
  if (init){

    // Get info on the neighbor list
    double cell_size_x, cell_size_y, cell_size_z, cell_size_min;
    get_cell_list_info(NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, &cell_size_x,
                       &cell_size_y, &cell_size_z);
    cell_size_min = cell_size_x;
    if (cell_size_min < cell_size_y) cell_size_min = cell_size_y;
    if (cell_size_min < cell_size_z) cell_size_min = cell_size_z;
    if (cell_size_min < pressv_rmax){
      printf("ERROR: The size of the cells in the neighbor list does not allow a correct calculation of the pressure, increase to %f\n", 
	     pressv_rmax);
      exit(EXIT_FAILURE);
    }

  }



  if (in.neigh_dr < pressv_rmax){

  }

  // Initialize histogram
  pressv_hist_init();

  // Fill histogram
  pressv_compute_hist();

  // compute Radial distribution function close to contact
  pressv_compute_rdf();
  
  // Write output
  pressv_output(init);

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
  int **cl_part_cell, **cl_neigh, cl_neigh_num;
  int cell_idx, neigh_idx, n_part_cell, part_idx;
  double dr;
  int bin;

  // Neighbor list
  get_cell_list_info(&cl_part_cell, &cl_neigh,
                     &cl_neigh_num, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL);

  // Fill histograms
  for (int ii=0; ii<part_info.NN; ii++){

    cell_idx = cell_part_idx(ii);

    for (int jj=0; jj<cl_neigh_num; jj++){

      neigh_idx = cl_neigh[cell_idx][jj];
      n_part_cell = cl_part_cell[neigh_idx][0];

      if (n_part_cell > 0){
	for (int kk=1; kk<=n_part_cell; kk++){

	  part_idx = cl_part_cell[neigh_idx][kk];

	  dr = compute_dist(ii, part_idx, 1.0, 1.0, 1.0);

	  if (dr <= pressv_rmax && part_idx > ii) {
	    bin = (int)((dr-1.0)/in.pressv_dr);
	    pressv_hist[bin] += 2.0;
	  }

	}

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
  presst_output(init);

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
  double vol_ratio, sf;
  bool overlap = false;
  
  // Fill histograms
  for (int ii=0; ii<presst_hist_nn; ii++){

    // Box length after re-scaling
    vol_ratio = 1 - presst_xi[ii];
    
    sf = pow(vol_ratio, 1./3.);
    
    // Check if there is overlap
    for (int jj=0; jj<part_info.NN; jj++){
      
      overlap = check_overlap(jj, sf, sf, sf);

      if (overlap) break;

    }

    // Update histogram if there is no overlap
    if (!overlap) presst_hist[ii] += 1.0;


  }
    
}

void pressv_output(bool init){

  // Print to file the sample needed to compute the pressure
  FILE* fid;
  if (init) fid = fopen("press_virial.dat", "w");
  else fid = fopen("press_virial.dat", "a");
  if (fid == NULL) {
    perror("Error while creating the file for the virial pressure\n");
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

}

void presst_output(bool init){

  // Print to file the sample needed to compute the pressure
  FILE* fid;
  if (init) fid = fopen("press_thermo.dat", "w");
  else fid = fopen("press_thermo.dat", "a");
  if (fid == NULL) {
    perror("Error while creating the file for the thermo pressure\n");
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

  // For NpT calculations print to file also the density 
  if (in.press > 0) {
    
    if (init) fid = fopen("density.dat", "w");
    else fid = fopen("density.dat", "a");
    if (fid == NULL) {
      perror("Error while creating the file for the density\n");
      exit(EXIT_FAILURE);
    }
    if (init){
      fprintf(fid, "######################################\n");
      fprintf(fid, "# Density (each line is one sample)\n");
      fprintf(fid, "######################################\n");
    }
    fprintf(fid, "%.8e\n", in.rho);
    fclose(fid);
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
