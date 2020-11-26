#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_fit.h>

// Declare structures
struct p_info {int NN, Ncell; };

struct box_info {
  double vol;
  double lx, ly, lz;
  double min_size;
  double cell_size;
  int cell_x, cell_y, cell_z;
  int cell_type;
};


// Declare functions
void hs_nvt(struct input in);

struct box_info sim_box_init(int cell_type, int nx, int ny,
                             int nz, double rho);

struct p_info part_alloc(double ***part,
                         struct box_info sim_box_info);

void part_init(double **part, struct box_info sim_box_info,
               struct p_info part_info);

void part_init_sc(double **part, struct box_info sim_box_info);

void part_init_fcc(double **part, struct box_info sim_box_info);

void part_init_err();

void add_particle(double **part, int id,
                  double xx, double yy, double zz);

double run_opt(double **part, struct input in, 
	     struct box_info sim_box_info, 
	     struct p_info part_info, const gsl_rng *r);

void run(double **part, struct input in, 
	     struct box_info sim_box_info, 
	     struct p_info part_info, const gsl_rng *r);


void sweep(double **part, int *acc_moves, int *rej_moves,
           struct input in, struct box_info sim_box_info,
           struct p_info part_info, const gsl_rng *r, 
	   long unsigned int max_r);

void metropolis(double **part, int *accept, int *reject,
                int r_idx, double r_x, double r_y, double r_z,
                struct input in, struct box_info sim_box_info,
                struct p_info part_info);

bool check_overlap(double **part, int target_idx, 
		   double r_x, double r_y, double r_z,
                   struct box_info sim_box_info, struct p_info part_info);

void compute_hist(double *hist, int *counter, double **part,
                  double *pos, double cutoff,
                  struct input in, struct box_info sim_box_info,
                  struct p_info part_info);

void average_hist(double *hist, int nn, struct input in);

void reset_hist(double *hist, int *counter, int nn, struct input in);

void compute_rdf(double *rdf, double *hist, double *rr, int nn, 
		 struct input in, struct p_info part_info);

void compute_pressure(double *press, double *rdf, double *rr, int nn,
                        struct input in, struct p_info part_info);
