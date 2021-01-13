#include <stdlib.h>
#include "init.h"
#include "read_input.h"
#include "compute_press.h"
#include "moves.h"
#include "npt.h"
#include "nvt.h"
#include "cavity_nvt.h"
#include "optimizer.h"


void opt_nvt(){

  int sample_iter = in.opt_sweeps/in.opt_samples;
  double acc_ratio_1, acc_ratio_2, dr_1, dr_2;

  printf("---------------------------------------------------\n");
  printf("Maximum displacement optimization started ...\n");
  printf("Sweeps for optimization: %d\n", in.opt_sweeps);
  printf("Number of samples: %d\n", in.opt_samples);

  // Collect two samples to get started
  get_sample_nvt(&dr_1, &acc_ratio_1, sample_iter);
  if (acc_ratio_1 > in.opt_part_target){
    in.dr_max *= 2;
  }
  else {
    in.dr_max /= 2;
  }
  get_sample_nvt(&dr_2, &acc_ratio_2, sample_iter);


  // Secant-method to find optimum value
  for (int ii=0; ii<in.opt_samples; ii++){

    in.dr_max = dr_2 - (acc_ratio_2 - in.opt_part_target) * (dr_2 - dr_1)/(acc_ratio_2 - acc_ratio_1);
    if (in.dr_max > 1.0) {
      in.dr_max = 1.0;
    }
    else if (in.dr_max <= 0.0) {
      printf("Error: maximum displacement is zero!\n");
      exit(EXIT_FAILURE);
    }
    dr_1 = dr_2;
    acc_ratio_1 = acc_ratio_2;
    get_sample_nvt(&dr_2, &acc_ratio_2, sample_iter);
  }

  printf("Optimal maximum displacement: %.8f\n", in.dr_max);
  printf("Acceptance ratio: %.8f \n", acc_ratio_2);
  printf("Maximum displacement optimization completed\n");

}

void get_sample_nvt(double *dr, double *acc_ratio, int sample_iter){

  part_moves = 0;
  acc_part_moves=0;
  rej_part_moves=0;
  for (int ii=0; ii<sample_iter; ii++){
    sweep_nvt();
  }
  *dr = in.dr_max;
  *acc_ratio = (double)acc_part_moves/((double)part_moves);

}


void opt_npt(){

  int sample_iter = in.opt_sweeps/in.opt_samples;
  double acc_ratio_part_1, acc_ratio_part_2, dr_1, dr_2;
  double acc_ratio_vol_1, acc_ratio_vol_2, dv_1, dv_2;

  printf("---------------------------------------------------\n");
  printf("Maximum displacement optimization started ...\n");
  printf("Sweeps for optimization: %d\n", in.opt_sweeps);
  printf("Number of samples: %d\n", in.opt_samples);

  // Collect two samples to get started
  get_sample_npt(&dr_1, &acc_ratio_part_1, 
		  &dv_1, &acc_ratio_vol_1, 
		  sample_iter);
  if (acc_ratio_part_1 > in.opt_part_target){
    in.dr_max *= 2;
  }
  else {
    in.dr_max /= 2;
  }
  if (acc_ratio_vol_1 > in.opt_vol_target){
    in.dv_max *= 2;
  }
  else {
    in.dv_max /= 2;
  }
  get_sample_npt(&dr_2, &acc_ratio_part_2, 
		  &dv_2, &acc_ratio_vol_2, 
		  sample_iter);

  // Secant-method to find optimum value
  for (int ii=0; ii<in.opt_samples; ii++){

    in.dr_max = dr_2 - (acc_ratio_part_2 - in.opt_part_target) * (dr_2 - dr_1)/(acc_ratio_part_2 - acc_ratio_part_1);
    in.dv_max = dv_2 - (acc_ratio_vol_2 - in.opt_vol_target) * (dv_2 - dv_1)/(acc_ratio_vol_2 - acc_ratio_vol_1);
    if (in.dr_max > 1.0) {
      in.dr_max = 1.0;
    }
    dr_1 = dr_2;
    dv_1 = dv_2;
    acc_ratio_part_1 = acc_ratio_part_2;
    acc_ratio_vol_1 = acc_ratio_vol_2;
    get_sample_npt(&dr_2, &acc_ratio_part_2, 
		  &dv_2, &acc_ratio_vol_2, 
		  sample_iter);

  }

  // Stop if the optimal values are negative
  if (in.dr_max <= 0.0) {
      printf("Error: maximum particle displacement is zero!\n");
      exit(EXIT_FAILURE);
  }
  if (in.dv_max <= 0.0) {
      printf("Error: maximum volume displacement is zero!\n");
      exit(EXIT_FAILURE);
  }

  // Print result of optimization on screen
  printf("Optimal maximum particle displacement: %.8f\n", in.dr_max);
  printf("Acceptance ratio: %.8f \n", acc_ratio_part_2);
  printf("Optimal maximum volume deformation: %.8f\n", in.dv_max);
  printf("Acceptance ratio: %.8f \n", acc_ratio_vol_2);
  printf("Maximum displacement optimization completed\n");

}

void get_sample_npt(double *dr, double *acc_ratio_part, 
		    double *dv, double *acc_ratio_vol,
		    int sample_iter){

  part_moves = 0;
  vol_moves = 0;
  acc_part_moves= 0;
  rej_part_moves= 0;
  acc_vol_moves= 0;
  rej_vol_moves= 0;
  for (int ii=0; ii<sample_iter; ii++){
    sweep_npt();
  }
  *dr = in.dr_max;
  *dv = in.dv_max;
  *acc_ratio_part = (double)acc_part_moves/((double)part_moves);
  *acc_ratio_vol = (double)acc_vol_moves/((double)vol_moves);

}


void opt_cavity_nvt(){

  int sample_iter = in.opt_sweeps/in.opt_samples;
  double acc_ratio_1, acc_ratio_2, dr_1, dr_2;

  printf("---------------------------------------------------\n");
  printf("Maximum displacement optimization started ...\n");
  printf("Sweeps for optimization: %d\n", in.opt_sweeps);
  printf("Number of samples: %d\n", in.opt_samples);

  // Collect two samples to get started
  get_sample_cavity_nvt(&dr_1, &acc_ratio_1, sample_iter);
  if (acc_ratio_1 > in.opt_part_target){
    in.dr_max *= 2;
  }
  else {
    in.dr_max /= 2;
  }
  get_sample_cavity_nvt(&dr_2, &acc_ratio_2, sample_iter);


  // Secant-method to find optimum value
  for (int ii=0; ii<in.opt_samples; ii++){

    in.dr_max = dr_2 - (acc_ratio_2 - in.opt_part_target) * (dr_2 - dr_1)/(acc_ratio_2 - acc_ratio_1);
    if (in.dr_max > 1.0) {
      in.dr_max = 1.0;
    }
    else if (in.dr_max <= 0.0) {
      printf("Error: maximum displacement is zero!\n");
      exit(EXIT_FAILURE);
    }
    dr_1 = dr_2;
    acc_ratio_1 = acc_ratio_2;
    get_sample_cavity_nvt(&dr_2, &acc_ratio_2, sample_iter);
  }

  printf("Optimal maximum displacement: %.8f\n", in.dr_max);
  printf("Acceptance ratio: %.8f \n", acc_ratio_2);
  printf("Maximum displacement optimization completed\n");

}

void get_sample_cavity_nvt(double *dr, double *acc_ratio, int sample_iter){

  part_moves = 0;
  acc_part_moves=0;
  rej_part_moves=0;
  for (int ii=0; ii<sample_iter; ii++){
    cavity_sweep_nvt();
  }
  *dr = in.dr_max;
  *acc_ratio = (double)acc_part_moves/((double)part_moves);

}
