#include <time.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"
#include "compute_press.h"
#include "compute_order_parameter.h"
#include "compute_widom_chem_pot.h"
#include "moves.h"
#include "io_config.h"
#include "optimizer.h"
#include "cavity_nvt.h"
#include "analytic.h"

// Cavity simulation for hard-spheres in the NVT ensemble
void cavity_hs_nvt() {

  if (in.restart_read == 0){ // Read from input file

    // Simulation box
    sim_box_init(in.type, in.nx, in.ny, in.nz, in.rho);

    // Particles
    part_alloc();

    // Initialize particle's positions
    part_init();

    // Set cavities to a distance within the allowed interval
    cavity_set_distance();

    // Set-up random number generator (Marsenne-Twister)
    rng_init();

  }
  else { // Read from restart file
    read_restart(in.restart_name);
  }

  // Print simulation info on screen
  printf("Simulation box size (x, y, z): %.5f %.5f %.5f\n", sim_box_info.lx,
         sim_box_info.ly, sim_box_info.lz);
  printf("Number of particles: %d\n", part_info.NN);

  // Initialize cavity interaction potential (0.5 is just a dummy distance)
  cavity_interaction(0.5, true);

  // Write cavity interaction potential to file 
  cavity_psi_output();

  // Set-up the neighbor list
  int cl_neigh_num, cl_num_tot;
  compute_cell_list_info();
  get_cell_list_info(&cl_neigh_num, &cl_num_tot, NULL, NULL, NULL,
                     NULL, NULL, NULL);
  int (*cl_neigh)[cl_neigh_num] = (int (*)[cl_neigh_num])cell_list_alloc(cl_num_tot, cl_neigh_num);
  int (*cl_part_cell)[in.neigh_max_part] = (int (*)[in.neigh_max_part])cell_list_alloc(cl_num_tot, in.neigh_max_part);
  cell_list_init(cl_neigh_num, cl_neigh, in.neigh_max_part, cl_part_cell);


  // Optmize maximum displacement
  if (in.opt_flag == 1){
    opt_cavity_nvt(cl_num_tot, in.neigh_max_part, cl_part_cell,
		   cl_neigh_num, cl_neigh);
    part_init();
    cavity_set_distance();
    cell_list_init(cl_neigh_num, cl_neigh, in.neigh_max_part, cl_part_cell);
  }

  // Start timing
  clock_t start = clock();

  // Initialize move counters
  reset_moves_counters();

  // Run equilibration
  printf("---------------------------------------------------\n");
  printf("Equilibration...\n");
  cavity_run_nvt(false,0, cl_num_tot, in.neigh_max_part, cl_part_cell,
		 cl_neigh_num, cl_neigh);
  printf("Equilibration completed.\n");

  // Run statistics
  printf("---------------------------------------------------\n");
  printf("Production...\n");
  cavity_run_nvt(true,in.sweep_eq, cl_num_tot, in.neigh_max_part, cl_part_cell,
		 cl_neigh_num, cl_neigh);
  printf("Production completed.\n");
  clock_t end = clock();

  // Print acceptance and rejection percentages
  int part_moves, acc_part_moves, rej_part_moves;
  get_moves_counters(&part_moves, &acc_part_moves, &rej_part_moves,
                     NULL, NULL, NULL); 
  printf("---------------------------------------------------\n");
  printf("-- Particle moves: %.8e\n", (double)part_moves);
  printf("   Acceptance percentage: %f\n", (double)acc_part_moves/((double)part_moves));
  printf("   Rejection percentage: %f\n", (double)rej_part_moves/((double)part_moves));

  
  // Stop timing
  printf("Elapsed time: %f seconds\n",
	 (double)(end - start) / CLOCKS_PER_SEC);
  
  // Free memory
  free(part);
  gsl_rng_free(rng_mt);
  free(cl_neigh);
  free(cl_part_cell);

}

void cavity_run_nvt(bool prod_flag, int sweep_offset,
		    int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		    int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]){

  // Variable declaration
  bool cavity_init = true;
  int n_sweeps;

  // Number of sweeps
  if (prod_flag) n_sweeps = in.sweep_stat;
  else n_sweeps = in.sweep_eq;
  
  // Run MC simulation
  for (int ii=sweep_offset; ii<n_sweeps+sweep_offset; ii++){

    // Output on screen
    if (ii == 0){
      printf("Sweep number\n");
    }
    if (ii % in.output_int == 0) {
      printf("%d\n", ii);
      fflush(stdout);
    }

    // Write restart file
    if (in.restart_write > 0){
      if (ii % in.restart_write == 0) {
        write_restart(ii);
      }
    }

    // Save samples for production runs
    if (prod_flag){

      // Write configuration
      if (in.config_write > 0){
        if (ii % in.config_write == 0) {
          write_config(ii);
        }
      }


      // Output distance between the cavities
      if (in.cavity_sample_int > 0){
        if (ii % in.cavity_sample_int == 0) {
          cavity_dist_output(cavity_init);
          if (cavity_init) cavity_init = false;
        }
      }

    }
    
    // Generate new configuration
    cavity_sweep_nvt(cl_num_tot, in.neigh_max_part, cl_part_cell,
		     cl_neigh_num, cl_neigh);

  }  

}

void cavity_sweep_nvt(int cl_num_tot, int cl_max_part, int cl_part_cell[cl_num_tot][cl_max_part],
		      int cl_neigh_num, int cl_neigh[cl_num_tot][cl_neigh_num]){

  // Create N trial moves (N = number of particles)
  for (int ii=0; ii<part_info.NN; ii++){
    cavity_part_move(cl_num_tot, in.neigh_max_part, cl_part_cell,
		     cl_neigh_num, cl_neigh);
  }

}

void cavity_set_distance(){
  
  double cavity_avedr = (in.cavity_maxdr + in.cavity_mindr)/ 2.0;
  part[1][1] = part[0][1];
  part[1][2] = part[0][2];
  part[1][3] = part[0][3] + cavity_avedr;
  if (part[1][3] > sim_box_info.lz) part[1][3] -= sim_box_info.lz;
  else if (part[1][3] < 0.0)        part[1][3] += sim_box_info.lz;

}

void cavity_psi_output(){

  FILE* fid;
  fid = fopen("cavity_psi.dat", "w");
  if (fid == NULL) {
    perror("Error while creating the file for the cavity potential\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fid, "#########################################################\n");
  fprintf(fid, "# Interaction potential  between the cavities \n");
  fprintf(fid, "# Distance, potential\n");
  fprintf(fid, "#########################################################\n");
  double xx = in.cavity_mindr;
  int n_points = (int)((in.cavity_maxdr - in.cavity_mindr)/in.cavity_out_dr);
  for (int ii=0; ii < n_points; ii++){
    fprintf(fid, "%.8e %.8e\n", xx, cavity_interaction(xx, false));
    xx += in.cavity_out_dr;
  }
  fclose(fid);

}

void cavity_dist_output(bool init){

  FILE* fid;
  if (init) fid = fopen("cavity_distance.dat", "w");
  else fid = fopen("cavity_distance.dat", "a");
  if (fid == NULL) {
    perror("Error while creating the file for the cavity distance\n");
    exit(EXIT_FAILURE);
  }
  if (init){
    fprintf(fid, "#########################################################\n");
    fprintf(fid, "# Distance between the cavities (each line is one sample)\n");
    fprintf(fid, "#########################################################\n");
  }
  fprintf(fid, "%.8e\n", compute_dist(0,1,1.0,1.0,1.0));
  fclose(fid);

}



