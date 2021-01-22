#include <stdlib.h>
#include <stdio.h>
#include "sim_info.h"
#include "read_input.h"
#include "compute_press.h"
#include "moves.h"
#include "npt.h"
#include "nvt.h"
#include "cavity_nvt.h"
#include "optimizer.h"

// --------------------------------------------------------
// The module "optimizer.c" is used to find the optimal 
// maximum displacements for particles motions and volume
// displacements. Here optimal is intended as the maximum
// displacement that allows approximately 50% of the moves
// to be accepted, but the user is free to specify a 
// different percentage as optimal
// --------------------------------------------------------

// ------ Optimizer for NVT simulations ------

void opt_nvt(){

  int sample_iter = G_IN.opt_sweeps/G_IN.opt_samples;
  double acc_ratio_1, acc_ratio_2, dr_1, dr_2;

  printf("---------------------------------------------------\n");
  printf("Maximum displacement optimization started ...\n");
  printf("Sweeps for optimization: %d\n", G_IN.opt_sweeps);
  printf("Number of samples: %d\n", G_IN.opt_samples);

  // Collect two samples to get started
  get_sample_nvt(&dr_1, &acc_ratio_1, sample_iter);
  if (acc_ratio_1 > G_IN.opt_part_target){
    G_IN.dr_max *= 2;
  }
  else {
    G_IN.dr_max /= 2;
  }
  get_sample_nvt(&dr_2, &acc_ratio_2, sample_iter);


  // Secant-method to find optimum value
  for (int ii=0; ii<G_IN.opt_samples; ii++){

    G_IN.dr_max = dr_2 - (acc_ratio_2 - G_IN.opt_part_target) * (dr_2 - dr_1)/(acc_ratio_2 - acc_ratio_1);
    if (G_IN.dr_max > 1.0) {
      G_IN.dr_max = 1.0;
    }
    else if (G_IN.dr_max <= 0.0) {
      G_IN.dr_max = -G_IN.dr_max;
      if (G_IN.dr_max > 1.0) G_IN.dr_max = dr_2/2;
    }
    dr_1 = dr_2;
    acc_ratio_1 = acc_ratio_2;
    get_sample_nvt(&dr_2, &acc_ratio_2, sample_iter);
  }

  printf("Optimal maximum displacement: %.8f\n", G_IN.dr_max);
  printf("Acceptance ratio: %.8f\n", acc_ratio_2);
  printf("Maximum displacement optimization completed\n");

}

void get_sample_nvt(double *dr, double *acc_ratio, int sample_iter){

  reset_moves_counters();
  for (int ii=0; ii<sample_iter; ii++){
    sweep_nvt();
  }
  int part_moves, acc_part_moves;
  get_moves_counters(&part_moves, &acc_part_moves, NULL,
		     NULL, NULL, NULL);
  *dr = G_IN.dr_max;
  *acc_ratio = (double)acc_part_moves/((double)part_moves);

}


// ------ Optimizer for NpT simulations ------

void opt_npt(){

  int sample_iter = G_IN.opt_sweeps/G_IN.opt_samples;
  double acc_ratio_part_1, acc_ratio_part_2, dr_1, dr_2;
  double acc_ratio_vol_1, acc_ratio_vol_2, dv_1, dv_2;

  printf("---------------------------------------------------\n");
  printf("Maximum displacement optimization started ...\n");
  printf("Sweeps for optimization: %d\n", G_IN.opt_sweeps);
  printf("Number of samples: %d\n", G_IN.opt_samples);

  // Collect two samples to get started
  get_sample_npt(&dr_1, &acc_ratio_part_1, 
		 &dv_1, &acc_ratio_vol_1, 
		 sample_iter);
  if (acc_ratio_part_1 > G_IN.opt_part_target){
    G_IN.dr_max *= 2;
  }
  else {
    G_IN.dr_max /= 2;
  }
  if (acc_ratio_vol_1 > G_IN.opt_vol_target){
    G_IN.dv_max *= 2;
  }
  else {
    G_IN.dv_max /= 2;
  }
  get_sample_npt(&dr_2, &acc_ratio_part_2, 
		 &dv_2, &acc_ratio_vol_2, 
		 sample_iter);

  // Secant-method to find optimum value
  for (int ii=0; ii<G_IN.opt_samples; ii++){

    G_IN.dr_max = dr_2 - (acc_ratio_part_2 - G_IN.opt_part_target) * (dr_2 - dr_1)/(acc_ratio_part_2 - acc_ratio_part_1);
    G_IN.dv_max = dv_2 - (acc_ratio_vol_2 - G_IN.opt_vol_target) * (dv_2 - dv_1)/(acc_ratio_vol_2 - acc_ratio_vol_1);
    if (G_IN.dr_max > 1.0) {
      G_IN.dr_max = 1.0;
    }
    // Stop if the optimal values are negative or zero
    if (G_IN.dr_max <= 0.0) {
      G_IN.dr_max = -G_IN.dr_max;
      if (G_IN.dr_max > 1.0) G_IN.dr_max = dr_2/2;
    }
    if (G_IN.dv_max <= 0.0) {
      G_IN.dv_max = -G_IN.dv_max;
      if (G_IN.dv_max > 0.1) G_IN.dv_max = dv_2/2;
    }
    dr_1 = dr_2;
    dv_1 = dv_2;
    acc_ratio_part_1 = acc_ratio_part_2;
    acc_ratio_vol_1 = acc_ratio_vol_2;
    get_sample_npt(&dr_2, &acc_ratio_part_2,
    		   &dv_2, &acc_ratio_vol_2,
    		   sample_iter);

  }

  // Print result of optimization on screen
  printf("Optimal maximum particle displacement: %.8f\n", G_IN.dr_max);
  printf("Acceptance ratio: %.8f \n", acc_ratio_part_2);
  printf("Optimal maximum volume deformation: %.8f\n", G_IN.dv_max);
  printf("Acceptance ratio: %.8f \n", acc_ratio_vol_2);
  printf("Maximum displacement optimization completed\n");

}

void get_sample_npt(double *dr, double *acc_ratio_part, 
		    double *dv, double *acc_ratio_vol,
		    int sample_iter){

  reset_moves_counters();
  for (int ii=0; ii<sample_iter; ii++){
    sweep_npt();
  }
  int part_moves, acc_part_moves, vol_moves, acc_vol_moves;
  get_moves_counters(&part_moves, &acc_part_moves, NULL,
		     &vol_moves, &acc_vol_moves, NULL);
  *dr = G_IN.dr_max;
  *dv = G_IN.dv_max;
  *acc_ratio_part = (double)acc_part_moves/((double)part_moves);
  *acc_ratio_vol = (double)acc_vol_moves/((double)vol_moves);

}



// ------ Optimizer for cavity simulations ------

void opt_cavity_nvt(){

  int sample_iter = G_IN.opt_sweeps/G_IN.opt_samples;
  double acc_ratio_1, acc_ratio_2, dr_1, dr_2;

  printf("---------------------------------------------------\n");
  printf("Maximum displacement optimization started ...\n");
  printf("Sweeps for optimization: %d\n", G_IN.opt_sweeps);
  printf("Number of samples: %d\n", G_IN.opt_samples);

  // Collect two samples to get started
  get_sample_cavity_nvt(&dr_1, &acc_ratio_1, sample_iter);
  if (acc_ratio_1 > G_IN.opt_part_target){
    G_IN.dr_max *= 2;
  }
  else {
    G_IN.dr_max /= 2;
  }
  get_sample_cavity_nvt(&dr_2, &acc_ratio_2, sample_iter);


  // Secant-method to find optimum value
  for (int ii=0; ii<G_IN.opt_samples; ii++){

    G_IN.dr_max = dr_2 - (acc_ratio_2 - G_IN.opt_part_target) * (dr_2 - dr_1)/(acc_ratio_2 - acc_ratio_1);
    if (G_IN.dr_max > 1.0) {
      G_IN.dr_max = 1.0;
    }
    else if (G_IN.dr_max <= 0.0) {
      G_IN.dr_max = -G_IN.dr_max;
      if (G_IN.dr_max > 1.0) G_IN.dr_max = dr_2/2;
    }
    dr_1 = dr_2;
    acc_ratio_1 = acc_ratio_2;
    get_sample_cavity_nvt(&dr_2, &acc_ratio_2, sample_iter);
  }

  printf("Optimal maximum displacement: %.8f\n", G_IN.dr_max);
  printf("Acceptance ratio: %.8f \n", acc_ratio_2);
  printf("Maximum displacement optimization completed\n");

}

void get_sample_cavity_nvt(double *dr, double *acc_ratio, int sample_iter){

  reset_moves_counters();
  for (int ii=0; ii<sample_iter; ii++){
    cavity_sweep_nvt();
  }
  int part_moves, acc_part_moves;
  get_moves_counters(&part_moves, &acc_part_moves, NULL,
                     NULL, NULL, NULL);
  *dr = G_IN.dr_max;
  *acc_ratio = (double)acc_part_moves/((double)part_moves);

}
