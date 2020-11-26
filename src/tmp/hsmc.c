#include "hsmc.h"

// ----------------------------------------
// ------------- Monte-Carlo --------------
// ----------------------------------------

// Hard-sphere simulation in the NVT ensemble
void hs_nvt(struct input in) {

  // Simulation box
  struct box_info sim_box_info = sim_box_init(in.type, in.nx, in.ny, in.nz, in.rho);
  printf("Simulation box size (x, y, z): %.5f %.5f %.5f\n", sim_box_info.lx,
	 sim_box_info.ly, sim_box_info.lz);

  // Particles
  double **part;
  struct p_info part_info = part_alloc(&part, sim_box_info);
  printf("Number of particles: %d\n", part_info.NN);

  // Initialize particle's positions
  part_init(part, sim_box_info, part_info);

  // Set-up random number generator (Marsenne-Twister)
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,0);

  // Run to determine the optimal maximum displacement
  if (in.dr_max < 0){
    in.dr_max *= -1;
    in.dr_max = run_opt(part, in, sim_box_info, part_info, r);
    part_init(part, sim_box_info, part_info);
  }
  

  // Run equilibration
  run(part, in, sim_box_info, part_info, r);


  /* // Free memory */
  free(part);
  gsl_rng_free(r);
  

}


struct box_info sim_box_init(int cell_type, int nx, int ny, int nz, double rho){

  // Variable declaration
  int part_cell;
  double cell_vol, cell_size;
  int min_n;
  struct box_info out;

  // Particles per cell
  if (cell_type == 1){
    part_cell = 1; // Simple cubic lattice
  }
  else if (cell_type == 2){
    part_cell = 4; // Face-centered cubic lattice
  }
  else {
    printf("Unknown lattice type, default to fcc\n");
    part_cell = 4;
  }

  cell_vol = part_cell/ rho;
  cell_size = pow(cell_vol,1./3.);

  // Smallest size of the simulation box
  min_n =  nx;
  if (ny < min_n) min_n = ny;
  if (nz < min_n) min_n = nz;

  // Write output
  out.cell_x = nx;
  out.cell_y = ny;
  out.cell_z = nz;
  out.lx = nx * cell_size;
  out.ly = ny * cell_size;
  out.lz = nz * cell_size;
  out.min_size = min_n * cell_size;
  out.vol = nx * ny *nz * cell_vol;
  out.cell_size = cell_size;
  out.cell_type = cell_type;

  return out;

}


struct p_info part_alloc(double ***part, struct box_info sim_box_info){

  // Variable declaration
  int part_cell=1, part_tot=1;
  struct p_info out;

  // Particles per cell
  if (sim_box_info.cell_type == 1){
    part_cell = 1; // Simple cubic lattice
  }
  else if (sim_box_info.cell_type == 2){
    part_cell = 4; // Face-centered cubic lattice
  }
  else {
    printf("Unknown lattice type, default to fcc\n");
    part_cell = 4;
  }

  // Total number of particles
  part_tot = sim_box_info.cell_x *
             sim_box_info.cell_y *
             sim_box_info.cell_z *  part_cell;
  
  // Allocate matrix to store particle information
  *part = malloc( sizeof(double *) * part_tot);

  if (*part == NULL){
    printf("ERROR: Failed particle allocation\n");
    exit(EXIT_FAILURE);
  }

  for (int ii=0; ii<part_tot; ii++){
    (*part)[ii] = malloc( sizeof(double) * 4);
    if ((*part)[ii] == NULL){
      printf("ERROR: Failed particle allocation\n");
      exit(EXIT_FAILURE);
    }
  }

  // Output
  out.Ncell = part_cell;
  out.NN = part_tot;
  return out;

}


void part_init(double **part, struct box_info sim_box_info, struct p_info part_info){

  if (sim_box_info.cell_type == 1){
    part_init_sc(part, sim_box_info);
  }
  else {
    part_init_fcc(part, sim_box_info);
  }

}

void part_init_sc(double **part, struct box_info sim_box_info){

  double aa = sim_box_info.cell_size;
  if (aa < 1) part_init_err(); // Check nearest-neighbor distance
  int partid = 0;
  for (int ii=0; ii<sim_box_info.cell_x; ii++){
    for (int jj=0; jj<sim_box_info.cell_y; jj++){
      for (int kk=0; kk<sim_box_info.cell_z; kk++){
	add_particle(part, partid, ii*aa, jj*aa, kk*aa);
	partid ++;
      }
    }
  }
}

void part_init_fcc(double **part, struct box_info sim_box_info){

  double aa = sim_box_info.cell_size;
  if (aa/sqrt(2.0) < 1) part_init_err(); // Check nearest neighbor distance
  int partid = 0;
  for (int ii=0; ii<sim_box_info.cell_x; ii++){
    for (int jj=0; jj<sim_box_info.cell_y; jj++){
      for (int kk=0; kk<sim_box_info.cell_z; kk++){
        add_particle(part, partid, ii*aa, jj*aa, kk*aa);
	add_particle(part, partid+1, (ii+0.5)*aa, (jj+0.5)*aa, kk*aa);
	add_particle(part, partid+2, (ii+0.5)*aa, jj*aa, (kk+0.5)*aa);
	add_particle(part, partid+3, ii*aa, (jj+0.5)*aa, (kk+0.5)*aa);
	partid += 4;
      }
    }
  }
}


void add_particle(double **part, int id, double xx, double yy, double zz){

  part[id][0] = id;
  part[id][1] = xx;
  part[id][2] = yy;
  part[id][3] = zz;

}

void part_init_err(){

  printf("Overlap in the initial configuration. Possible solutions:\n");
  printf("-- If SC lattice was selected, try to change to FCC\n");
  printf("-- If FCC lattice was selected, the selected value of density is unphysical\n");
  exit(EXIT_FAILURE);

}

double run_opt(double **part, struct input in, struct box_info sim_box_info, 
	     struct p_info part_info, const gsl_rng *r){

  int ii,jj;
  int max_iter=1000, n_samples=10, sample_iter = max_iter/n_samples;
  int acc_moves, rej_moves;
  double acc_ratio_1, acc_ratio_2, dr_1, dr_2;
  long unsigned int max_r = gsl_rng_max(r);

  printf("---------------------------------------------------\n");
  printf("Maximum displacement optimization started ...\n");
  printf("Sweeps for optimization: %d\n", max_iter);
  printf("Number of samples: %d\n", n_samples);

  
  acc_moves=0, rej_moves=0;
  for (ii=0; ii<sample_iter; ii++){
    sweep(part, &acc_moves, &rej_moves, in, sim_box_info, part_info, 
	  r, max_r);
  }
  dr_1 = in.dr_max;
  acc_ratio_1 = (double)acc_moves/((double)sample_iter*part_info.NN);
 
  if (acc_ratio_1 > 0.5){
    in.dr_max *= 2;
  }
  else {
    in.dr_max /= 2;
  }
  acc_moves=0, rej_moves=0;
  for (ii=0; ii<sample_iter; ii++){
    sweep(part, &acc_moves, &rej_moves, in, sim_box_info, part_info, 
	  r, max_r);
  }
  dr_2 = in.dr_max;
  acc_ratio_2 = (double)acc_moves/((double)sample_iter*part_info.NN);

  for (ii=0; ii<n_samples; ii++){

    in.dr_max = dr_2 - (acc_ratio_2 - 0.5) * (dr_2 - dr_1)/(acc_ratio_2 - acc_ratio_1);
    if (in.dr_max > 1.0) {
      in.dr_max = 1.0;
    }
    else if (in.dr_max <= 0.0) {
      exit(EXIT_FAILURE);
    }
    acc_moves=0, rej_moves=0;
    for (jj=0; jj<sample_iter; jj++){
      sweep(part, &acc_moves, &rej_moves, in, sim_box_info, part_info, 
	    r, max_r);
    }
    dr_1 = dr_2;
    dr_2 = in.dr_max;
    acc_ratio_1 = acc_ratio_2;
    acc_ratio_2 = (double)acc_moves/((double)sample_iter*part_info.NN);
  }

  printf("Optimal maximum displacement: %.8f. Acceptance ratio: %.8f \n", in.dr_max, acc_ratio_2);
  printf("Maximum displacement optimization completed\n");


  return in.dr_max;

}

void run(double **part, struct input in, struct box_info sim_box_info,
	     struct p_info part_info, const gsl_rng *r){

  // Variable declaration
  int ii, sample_counter=0;
  int acc_moves=0, rej_moves=0;
  long unsigned int max_r = gsl_rng_max(r);
  double * pos_hist;
  int rdf_nn;
  double *rr, *rdf;
  double press;

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
    sweep(part, &acc_moves, &rej_moves, in, sim_box_info, part_info,
	  r, max_r);

    // Histograms of particles positions
    if (ii % in.rdf_tcompute == 0) {
      compute_hist(pos_hist, &sample_counter, part, rr, 
		   rr[rdf_nn-1]+in.rdf_dr/2., in, sim_box_info, part_info); 
    }

    // Radial distribution function and pressure
    if (sample_counter % in.rdf_tave == 0 && sample_counter > 0) {      
      average_hist(pos_hist, rdf_nn,in);
      compute_rdf(rdf, pos_hist, rr, rdf_nn, in, part_info);
      compute_pressure(&press, rdf, rr, rdf_nn, in, part_info);
      reset_hist(pos_hist, &sample_counter, rdf_nn, in);
    }

    
    // Output on screen
    if (ii % in.dt_output == 0) {
      printf("Sweep number: %d, pressure: %.8f, g[sigma]: %.8f\n",
	     ii, press, rdf[0]);
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

void sweep(double **part, int *acc_moves, int *rej_moves, 
	   struct input in, struct box_info sim_box_info, 
	   struct p_info part_info, const gsl_rng *r, long unsigned int max_r){

  int ii;
  long unsigned int r_idx;
  double r_x, r_y, r_z;

  // Loop over the particles
  for (ii=0; ii<part_info.NN; ii++){

    // Random number in [0, NN-1] to select the particle
    r_idx = gsl_rng_uniform_int(r,part_info.NN);

    // Three random numbers in [0,1] for the displacement
    r_x = (double)gsl_rng_get(r)/(double)max_r;
    r_y = (double)gsl_rng_get(r)/(double)max_r;
    r_z = (double)gsl_rng_get(r)/(double)max_r;
    
    // Metropolis alogorithm to check if the move can be accepted
    metropolis(part, acc_moves, rej_moves, r_idx, r_x, r_y, r_z,
	       in, sim_box_info, part_info);
     
  }

}

void metropolis(double **part, int *accept, int *reject, 
		int r_idx, double r_x, double r_y, double r_z,
                struct input in, struct box_info sim_box_info,
                struct p_info part_info){

  // Proposed move
  double x_new = part[r_idx][1] + (2.0 * r_x - 1.0)*in.dr_max;
  double y_new = part[r_idx][2] + (2.0 * r_y - 1.0)*in.dr_max;
  double z_new = part[r_idx][3] + (2.0 * r_z - 1.0)*in.dr_max;

  // Periodic boundary conditions
  if (x_new > sim_box_info.lx) x_new -= sim_box_info.lx;
  else if (x_new < 0.0)        x_new += sim_box_info.lx;
  if (y_new > sim_box_info.ly) y_new -= sim_box_info.ly;
  else if (y_new < 0.0)        y_new += sim_box_info.ly;
  if (z_new > sim_box_info.lz) z_new -= sim_box_info.lz;
  else if (z_new < 0.0)        z_new += sim_box_info.lz;
      
  // Accept or reject move
  if (check_overlap(part, r_idx, x_new, y_new, z_new, sim_box_info, part_info)){
    // Reject move
    *reject+=1;
  }
  else {
    // Accept move
    part[r_idx][1] = x_new;
    part[r_idx][2] = y_new;
    part[r_idx][3] = z_new;
    *accept+=1;
  }

}

bool check_overlap(double **part, int target_idx, double r_x, double r_y, double r_z,
		   struct box_info sim_box_info, struct p_info part_info){
  
  double lx_2 = sim_box_info.lx/2.0;
  double ly_2 = sim_box_info.ly/2.0;
  double lz_2 = sim_box_info.lz/2.0;
  double dx, dy, dz, rr;

  for (int ii=0; ii<part_info.NN; ii++){

    if (ii != target_idx) {

      // Cartesian components of the distance
      dx = r_x - part[ii][1];
      dy = r_y - part[ii][2];
      dz = r_z - part[ii][3];
      
      // Consider periodic boundary conditions
      if (dx > lx_2)       dx -= sim_box_info.lx;
      else if (dx < -lx_2) dx += sim_box_info.lx;
      if (dy > ly_2)       dy -= sim_box_info.ly;
      else if (dy < -ly_2) dy += sim_box_info.ly;
      if (dz > lz_2)       dz -= sim_box_info.lz;
      else if ( dz< -lz_2) dz += sim_box_info.lz;

      // Radial distance
      rr = sqrt(dx*dx + dy*dy + dz*dz);
      
      // Signal that there is overlap
      if (rr < 1.0) return true;

    }

  }

  return false;

}

void compute_hist(double *hist, int *counter, double **part, 
		  double *pos, double cutoff,
		  struct input in, struct box_info sim_box_info, 
		  struct p_info part_info){
  
  // Variable declaration
  double lx_2 = sim_box_info.lx/2.0;
  double ly_2 = sim_box_info.ly/2.0;
  double lz_2 = sim_box_info.lz/2.0;
  double dx, dy, dz, rr;
  int bin;

  // Fill histograms
  for (int ii=0; ii<part_info.NN; ii++){
    for (int jj=ii+1; jj<part_info.NN; jj++){
 
      // Cartesian components of the distance
      dx = part[jj][1] - part[ii][1];
      dy = part[jj][2] - part[ii][2];
      dz = part[jj][3] - part[ii][3];
      
      // Consider periodic boundary conditions
      if (dx > lx_2)       dx -= sim_box_info.lx;
      else if (dx < -lx_2) dx += sim_box_info.lx;
      if (dy > ly_2)       dy -= sim_box_info.ly;
      else if (dy < -ly_2) dy += sim_box_info.ly;
      if (dz > lz_2)       dz -= sim_box_info.lz;
      else if ( dz< -lz_2) dz += sim_box_info.lz;

      // Radial distance
      rr = sqrt(dx*dx + dy*dy + dz*dz);
      
      // Update histogram count
      if (rr <= cutoff) {
      	bin = (int)((rr-1.0)/in.rdf_dr);
	hist[bin] += 2.0;
      }
 
    }
  }

  // Update counter
  *counter += 1;

}


void average_hist(double *hist, int nn, struct input in){
  
  for (int ii=0; ii<nn; ii++){
     hist[ii] /= in.rdf_tave;
  } 

}


void reset_hist(double *hist, int *counter, int nn, struct input in){
  

  // Reset histograms
  for (int ii=0; ii<nn; ii++){
    hist[ii] = 0.0;
  } 

  // Reset counter
  *counter = 0;

}


void compute_rdf(double *rdf, double *hist, double *rr, int nn, struct input in,
		 struct p_info part_info){


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

void compute_pressure(double *press, double *rdf, double *rr, int nn, 
			struct input in, struct p_info part_info){

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
  
}
