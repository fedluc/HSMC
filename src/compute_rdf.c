#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <zlib.h>
#include "sim_info.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_rdf.h"
#include "moves.h"

// --------------------------------------------------------
// The module "compute_rdf.c" is used to compute the 
// the radial distribution function as defined in any 
// classical textbook.
// --------------------------------------------------------

// ------------------------------------------------------
// NOTE: This compute scales like N^2, it can have a 
// significant impact on performance if used with large
// number of particles 
// ------------------------------------------------------

// ------ Global variables ------

// Variables for the calculation of the radial distribution function
static int rdf_hist_nn;
static double *rdf_rr, *rdf_hist;

// Variables for the output
static int rdf_samples = 0;
static int rdf_file_id = 0;


// ------ Caller to compute the radial distribution function ------
 
void compute_rdf(bool init, int sweep){
  
  // Limit the cutoff to half the simulation box
  if (init){

    // Get info on the simulation box
    box_info sim_box_info = sim_box_info_get();
    double sim_box_size_min = sim_box_info.lx;
    if (sim_box_size_min < sim_box_info.ly) sim_box_size_min = sim_box_info.ly;
    if (sim_box_size_min < sim_box_info.lz) sim_box_size_min = sim_box_info.lz;
    if (sim_box_size_min < 2.0*G_IN.rdf_rmax){
      G_IN.rdf_rmax = sim_box_size_min/2.0;
      printf("WARNING: Cutoff for the rdf extraction reduced to %f in order to be consistent with minimum image convention\n",
             G_IN.rdf_rmax);
      
    }

    // Allocate vectors for rdf calculation
    rdf_hist_alloc();

  }

  // Initialize histogram
  rdf_hist_init();

  // Fill histogram with interparticle distances
  rdf_hist_compute();

  // Normalize histogram
  rdf_hist_norm();
  
  // Write output
  rdf_output(init,sweep);
  
}

// ------ Allocate and free the histograms for the radial distribution function calculation  ------

void rdf_hist_alloc(){

  rdf_hist_nn = (int)((G_IN.rdf_rmax - 1.0)/G_IN.rdf_dr);
  rdf_rr = (double*)malloc(sizeof(double) * rdf_hist_nn);
  rdf_hist = (double*)malloc(sizeof(double) * rdf_hist_nn);
  if (rdf_rr == NULL ||  rdf_hist == NULL){
    printf("ERROR: Failed histogram allocation\n");
    exit(EXIT_FAILURE);
  }

}

void rdf_hist_free(){

  free(rdf_hist);
  free(rdf_rr);

}

// ------ Initialize histogram for the calculation of the pressure via the virial route ------

void rdf_hist_init(){
  
  // Fill the histograms
  for (int ii=0; ii<rdf_hist_nn; ii++){
    rdf_rr[ii] = (ii+1./2.)*G_IN.rdf_dr + 1.0;
    rdf_hist[ii] = 0.0;
  }

  // Adjust the cutoff to the histograms
  G_IN.rdf_rmax = G_IN.rdf_dr*rdf_hist_nn + 1.0;
  
}

// ------ Fill histogrmas for the calculation of radial distribution function ------

void rdf_hist_compute(){
  
  // Variable declaration
  double dr;
  int bin;

  // Fill histograms
  p_info part_info = part_info_get();
  for (int ii=0; ii<part_info.NN; ii++){
    for (int jj=ii+1; jj<part_info.NN; jj++){
      dr = compute_dist(ii, jj, 1.0, 1.0, 1.0);
      if (dr < G_IN.rdf_rmax) {
	bin = (int)((dr-1.0)/G_IN.rdf_dr);
	rdf_hist[bin] += 2.0;
      }
    }
  }

}


// ------ Normalize the histograms for the calculation of the pressure from the virial route ------

void rdf_hist_norm(){

  // Variable declaration
  double r1,r2, dr, bin_vol;
  
  // Number of particles
  p_info part_info = part_info_get();
  
  // Normalize histograms
  dr = rdf_rr[1] - rdf_rr[0];
  for (int ii=0; ii<rdf_hist_nn; ii++){
    r1 = rdf_rr[ii] - dr/2.;
    r2 = rdf_rr[ii] + dr/2.;
    bin_vol = (4.*M_PI/3.) * (pow(r2,3.) - pow(r1,3.));
    rdf_hist[ii] = rdf_hist[ii]/(bin_vol * G_IN.rho * part_info.NN);
  }

}


// ------ Print to file the pressure obtained from the virial route ------

void rdf_output(bool init, int sweep){

  // Check that not too many output files are produced
  if (init){
    if ((double)(G_IN.sweep_stat+G_IN.sweep_eq)/(G_IN.rdf_sample_int * G_IN.rdf_samples) > 100000){
      printf("ERROR: Too many (> 100000) rdf files will be produced. Consider increasing number of samples per file\n");
    }
  }

  // Name of rdf file
  char rdf_file[30];
  sprintf(rdf_file, "rdf_%06d.dat.gz", rdf_file_id);

  // Open file
  gzFile fid = Z_NULL;
  if (rdf_samples == 0){
    fid = gzopen(rdf_file, "w");
  }
  else if (rdf_samples > 0 && rdf_samples < G_IN.rdf_samples){
    fid = gzopen(rdf_file, "a");
  }
  if (fid == Z_NULL) {
    perror("Error while creating rdf file");
    exit(EXIT_FAILURE);
  }

    // Write rdf
  p_info part_info = part_info_get();
  box_info sim_box_info = sim_box_info_get();
  gzprintf(fid, "######################################\n");
  gzprintf(fid, "# Sweep, Bins, volume, number of particles\n");
  gzprintf(fid, "######################################\n");
  gzprintf(fid, "%d %d %.8e %d\n", sweep, rdf_hist_nn,
  	   sim_box_info.vol, part_info.NN);
  gzprintf(fid, "###############################\n");
  gzprintf(fid, "# rr, rdf\n");
  gzprintf(fid, "###############################\n");
  for (int ii = 0; ii < rdf_hist_nn; ii++){
      gzprintf(fid, "%.8e %.8e\n", rdf_rr[ii], rdf_hist[ii]);
  }
  // Close binary file
  gzclose(fid);

  // Update counters
  rdf_samples += 1;
  if (rdf_samples == G_IN.rdf_samples){
    rdf_samples = 0;
    rdf_file_id += 1;
  }

}
