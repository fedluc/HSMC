#include "hsmc.h"

// ----------------------------------------
// ------------- Monte-Carlo --------------
// ----------------------------------------

// Global variables for random number generator
gsl_rng *rng_mt;
long unsigned int r_num_max;

// Global variables for particles moves
int acc_moves, rej_moves;
int r_idx;

// Hard-sphere simulation in the NVT ensemble
void hs_nvt() {

  // Simulation box
  sim_box_init(in.type, in.nx, in.ny, in.nz, in.rho);
  printf("Simulation box size (x, y, z): %.5f %.5f %.5f\n", sim_box_info.lx,
	 sim_box_info.ly, sim_box_info.lz);

  // Particles
  part_alloc();
  printf("Number of particles: %d\n", part_info.NN);

  // Initialize particle's positions
  part_init();

  // Initialize cell lists
  cell_list_init();

  // Set-up random number generator (Marsenne-Twister)
  rng_mt = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng_mt,0);
  r_num_max = gsl_rng_max(rng_mt);

  // Run to determine the optimal maximum displacement
  if (in.dr_max < 0){
    in.dr_max *= -1;
    run_opt();
    part_init();
  }

  // Run equilibration
  clock_t start = clock();
  run();
  clock_t end = clock();
  printf("Elapsed time: %f seconds\n",
	 (double)(end - start) / CLOCKS_PER_SEC);
  
  // Free memory
  free(part);
  gsl_rng_free(rng_mt);
  free(cl_neigh);
  free(cl_head);
  free(cl_link);

}


void run_opt(){

  int ii,jj;
  int max_iter=1000, n_samples=10, sample_iter = max_iter/n_samples;
  double acc_ratio_1, acc_ratio_2, dr_1, dr_2;

  printf("---------------------------------------------------\n");
  printf("Maximum displacement optimization started ...\n");
  printf("Sweeps for optimization: %d\n", max_iter);
  printf("Number of samples: %d\n", n_samples);

  // First step
  acc_moves=0, rej_moves=0;
  for (ii=0; ii<sample_iter; ii++){
    sweep();
  }
  dr_1 = in.dr_max;
  acc_ratio_1 = (double)acc_moves/((double)sample_iter*part_info.NN);

  // Second step
  if (acc_ratio_1 > 0.5){
    in.dr_max *= 2;
  }
  else {
    in.dr_max /= 2;
  }
  acc_moves=0, rej_moves=0;
  for (ii=0; ii<sample_iter; ii++){
    sweep();
  }
  dr_2 = in.dr_max;
  acc_ratio_2 = (double)acc_moves/((double)sample_iter*part_info.NN);

  // Secant-method to find optimum value
  for (ii=0; ii<n_samples; ii++){

    in.dr_max = dr_2 - (acc_ratio_2 - 0.5) * (dr_2 - dr_1)/(acc_ratio_2 - acc_ratio_1);
    if (in.dr_max > 1.0) {
      in.dr_max = 1.0;
    }
    else if (in.dr_max <= 0.0) {
      printf("Error: maximum displacement is zero!\n");
      exit(EXIT_FAILURE);
    }
    acc_moves=0, rej_moves=0;
    for (jj=0; jj<sample_iter; jj++){
      sweep();
    }
    dr_1 = dr_2;
    dr_2 = in.dr_max;
    acc_ratio_1 = acc_ratio_2;
    acc_ratio_2 = (double)acc_moves/((double)sample_iter*part_info.NN);
  }

  printf("Optimal maximum displacement: %.8f. Acceptance ratio: %.8f \n", in.dr_max, acc_ratio_2);
  printf("Maximum displacement optimization completed\n");

}

void run(){

  // Variable declaration
  int ii, sample_counter=0;
  double * pos_hist;
  int rdf_nn;
  double *rr, *rdf;
  double press=0.0;

  // Arrays for histograms of particle's positions
  rdf_nn = (int)((sim_box_info.min_size - 2.)/(2.*in.rdf_dr)) - 1;
  rr = (double*)malloc(sizeof(double) * rdf_nn);
  rdf = (double*)malloc(sizeof(double) * rdf_nn);
  pos_hist = (double*)malloc(sizeof(double) * rdf_nn);
  if (rr == NULL || rdf == NULL || pos_hist == NULL){
    printf("ERROR: Failed histogram allocation\n");
    exit(EXIT_FAILURE);
  }
  for (ii=0; ii<rdf_nn; ii++){
    rr[ii] = (ii+1./2.)*in.rdf_dr + 1.0;
    rdf[ii] = 0.0;
    pos_hist[ii] = 0.0;
  }

  // Run MC simulation
  printf("---------------------------------------------------\n");
  printf("Equilibration started ...\n");
  for (ii=0; ii<in.sweep_eq; ii++){

    // Sweep
    sweep();

    /* // Histograms of particles positions */
    /* if (ii % in.rdf_tcompute == 0) { */
    /*   compute_hist(pos_hist, &sample_counter, rr,  */
    /* 		   rr[rdf_nn-1]+in.rdf_dr/2.);  */
    /* } */

    /* // Radial distribution function and pressure */
    /* if (sample_counter % in.rdf_tave == 0 && sample_counter > 0) {       */
    /*   average_hist(pos_hist, rdf_nn); */
    /*   compute_rdf(rdf, pos_hist, rr, rdf_nn); */
    /*   compute_pressure(&press, rdf, rr, rdf_nn); */
    /*   reset_hist(pos_hist, &sample_counter, rdf_nn); */
    /* } */

    
    // Output on screen
    if (ii % in.dt_output == 0) {
      //printf("Sweep number: %d, pressure: %.8f, g[sigma]: %.8f\n",
      //ii, press, rdf[0]);
      printf("Sweep number: %d\n", ii, press, rdf[0]);
      fflush(stdout);
    }


  }

  printf("Equilibration completed\n");
  printf("Acceptance percentage: %f\n", (double)acc_moves/((double)in.sweep_eq*part_info.NN));
  printf("Rejection percentage: %f\n", (double)rej_moves/((double)in.sweep_eq*part_info.NN));

  // Output
  FILE* fid;
  fid = fopen("rdf.dat", "w");
  if (fid == NULL) {
    perror("Error while creating the rdf file");
    exit(EXIT_FAILURE);
  }
  for (int ii = 0; ii < rdf_nn; ii++)
    {
      fprintf(fid, "%.8e %.8e\n", rr[ii], rdf[ii]);
    }
  fclose(fid);
  

  /* // Free memory */
  free(rr);
  free(rdf);
  free(pos_hist);

}

void sweep(){

  int ii;
  double r_x, r_y, r_z;
  double x_old, y_old, z_old;

  // Create N trial moves (N = number of particles)
  for (ii=0; ii<part_info.NN; ii++){

    // Random number in [0, NN-1] to select the particle
    r_idx = gsl_rng_uniform_int(rng_mt, part_info.NN);

    // Three random numbers in [0,1] for the displacement
    r_x = (double)gsl_rng_get(rng_mt)/(double)r_num_max;
    r_y = (double)gsl_rng_get(rng_mt)/(double)r_num_max;
    r_z = (double)gsl_rng_get(rng_mt)/(double)r_num_max;

    // Store old coordinates (if the move gets rejected)
    x_old = part[r_idx][1];
    y_old = part[r_idx][2];
    z_old = part[r_idx][3];

    // Proposed move
    part[r_idx][1] += (2.0 * r_x - 1.0)*in.dr_max;
    part[r_idx][2] += (2.0 * r_y - 1.0)*in.dr_max;
    part[r_idx][3] += (2.0 * r_z - 1.0)*in.dr_max;

    // Periodic boundary conditions
    if (part[r_idx][1] > sim_box_info.lx) part[r_idx][1] -= sim_box_info.lx;
    else if (part[r_idx][1] < 0.0)        part[r_idx][1] += sim_box_info.lx;
    if (part[r_idx][2] > sim_box_info.ly) part[r_idx][2] -= sim_box_info.ly;
    else if (part[r_idx][2] < 0.0)        part[r_idx][2] += sim_box_info.ly;
    if (part[r_idx][3] > sim_box_info.lz) part[r_idx][3] -= sim_box_info.lz;
    else if (part[r_idx][3] < 0.0)        part[r_idx][3] += sim_box_info.lz;
      
    // Accept or reject move according to metropolis algorithm
    if (check_overlap()){
      // Reject move
      part[r_idx][1] = x_old;
      part[r_idx][2] = y_old;
      part[r_idx][3] = z_old;
      rej_moves+=1;
    }
    else {
      // Accept move
      acc_moves+=1;
      // Update cell list
      cell_list_update();
      
    }
 
  }

}


bool check_overlap(){
  
  // Variable declaration
  int cell_idx, neigh_idx, part_idx;
  double dr;

  // Cell that contains the particle
  cell_idx = cell_part_idx(r_idx);
    
  // Loop over the neighboring cells
  for (int ii=0; ii<26; ii++){

    neigh_idx = cl_neigh[cell_idx][ii];
    part_idx = cl_head[neigh_idx];

    // Loop over the particles in the neighboring cells
    while (part_idx > 0){
      
      // Compute inter-particle distance
      dr = compute_dist(r_idx, part_idx);

      // Signal that there is overlap
      if (dr < 1.0 && part_idx != r_idx)
	return true;

      // Update index
      part_idx = cl_link[part_idx];
	
    }

  }

  return false;

}


double compute_dist(int idx1, int idx2){

  double lx_2 = sim_box_info.lx/2.0;
  double ly_2 = sim_box_info.ly/2.0;
  double lz_2 = sim_box_info.lz/2.0;
  double dx, dy, dz, dr;

  // Cartesian components of the distance
  dx = part[idx1][1] - part[idx2][1];
  dy = part[idx1][2] - part[idx2][2];
  dz = part[idx1][3] - part[idx2][3];
  
  // Periodic boundary conditions
  if (dx > lx_2)       dx -= sim_box_info.lx;
  else if (dx < -lx_2) dx += sim_box_info.lx;
  if (dy > ly_2)       dy -= sim_box_info.ly;
  else if (dy < -ly_2) dy += sim_box_info.ly;
  if (dz > lz_2)       dz -= sim_box_info.lz;
  else if (dz< -lz_2) dz += sim_box_info.lz;
  
  // Radial distance
  dr =  sqrt(dx*dx + dy*dy + dz*dz);
  
  return dr;

}



void compute_hist(double *hist, int *counter,
		  double *pos, double cutoff){
  
  // Variable declaration
  double dr;
  int bin;

  // Fill histograms
  for (int ii=0; ii<part_info.NN; ii++){
    for (int jj=ii+1; jj<part_info.NN; jj++){

      // Inter-particle distance
      dr = compute_dist(ii, jj);
      
      // Update histogram count
      if (dr <= cutoff) {
      	bin = (int)((dr-1.0)/in.rdf_dr);
	hist[bin] += 2.0;
      }
 
    }
  }

  // Update counter
  *counter += 1;

}


void average_hist(double *hist, int nn){
  
  for (int ii=0; ii<nn; ii++){
     hist[ii] /= in.rdf_tave;
  }

}


void reset_hist(double *hist, int *counter, int nn){
  

  // Reset histograms
  for (int ii=0; ii<nn; ii++){
    hist[ii] = 0.0;
  }

  // Reset counter
  *counter = 0;

}


void compute_rdf(double *rdf, double *hist, double *rr, int nn){


  // Variable declaration
  double r1,r2, dr, bin_vol;
  
  // Normalize histograms
  dr = rr[1] - rr[0];
  for (int ii=0; ii<nn; ii++){
    r1 = rr[ii] - dr/2.;
    r2 = rr[ii] + dr/2.;
    bin_vol = (4.*M_PI/3.) * (pow(r2,3.) - pow(r1,3.));
    rdf[ii] = hist[ii]/(bin_vol * in.rho * part_info.NN);
  }


}

void compute_pressure(double *press, double *rdf, double *rr, int nn){

  // Check how many points are available for linear regression
  int ii=0, fit_pts=0;
  while (rr[ii] < 1.05) {
    fit_pts++;
    ii++;
  }

  // Arrays for linear regression
  double *rr_fit, *rdf_fit;
  rr_fit = (double*)malloc(sizeof(double) * fit_pts);
  rdf_fit = (double*)malloc(sizeof(double) * fit_pts);
  if (rr_fit == NULL || rdf_fit == NULL){
    printf("ERROR: Failed  allocation for linear regression arrays\n");
    exit(EXIT_FAILURE);
  }
  for (ii=0; ii<fit_pts; ii++){
    rr_fit[ii] = rr[ii];
    rdf_fit[ii] = rdf[ii];
  }

  // Linear regression
  double c0, c1, cov00, cov01, cov11, chisq;
  gsl_fit_linear(rr_fit, 1, rdf_fit, 1, fit_pts,
		 &c0, &c1, &cov00,
		 &cov01,  &cov11, &chisq);

  //Pressure (ideal + excess)
  *press = 1. + (2.*M_PI/3.) * in.rho * (c0 + c1);

  // Free memory
  free(rr_fit);
  free(rdf_fit);
  
}
