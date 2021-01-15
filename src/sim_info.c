#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "read_input.h"
#include "sim_info.h"

// ------ Global variables ------

// Global variable for simulation box information
static struct box_info sim_box_info;

// Global variables for particles information
static struct p_info part_info;
double (*part)[4];
static part_array *part_test;

// ------ Initialize simulation box ------

void sim_box_init(int cell_type, int nx, int ny, int nz, double rho){

  // Variable declaration
  int part_cell;
  double cell_vol, cell_size;
  int min_n;

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
  sim_box_info.cell_x = nx;
  sim_box_info.cell_y = ny;
  sim_box_info.cell_z = nz;
  sim_box_info.lx = nx * cell_size;
  sim_box_info.ly = ny * cell_size;
  sim_box_info.lz = nz * cell_size;
  sim_box_info.min_size = min_n * cell_size;
  sim_box_info.vol = nx * ny *nz * cell_vol;
  sim_box_info.cell_size = cell_size;
  sim_box_info.cell_type = cell_type;

}


// ------ Allocate multidimensional array for particles positions  ------

void part_alloc(){

  // Variable declaration
  int part_cell=1, part_tot=1;

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
  part = malloc(part_tot * sizeof(*part));
  if (part == NULL){
    printf("ERROR: Failed particle allocation\n");
    exit(EXIT_FAILURE);
  }

  // Allocate matrix to store particle information
  part_test = malloc(part_tot * sizeof(*part_test));
  if (part_test == NULL){
    printf("ERROR: Failed particle allocation\n");
    exit(EXIT_FAILURE);
  }

  // Output
  part_info.Ncell = part_cell;
  part_info.NN = part_tot;

}


// ------ Initialize particles positions  ------

void part_init(){

  if (sim_box_info.cell_type == 1){
    part_init_sc();
  }
  else {
    part_init_fcc();
  }

}

void part_init_sc(){

  double aa = sim_box_info.cell_size;
  if (aa < 1) part_init_err(); // Check nearest-neighbor distance
  int partid = 0;
  for (int ii=0; ii<sim_box_info.cell_x; ii++){
    for (int jj=0; jj<sim_box_info.cell_y; jj++){
      for (int kk=0; kk<sim_box_info.cell_z; kk++){
	add_particle(partid, ii*aa, jj*aa, kk*aa);
	partid ++;
      }
    }
  }
}

void part_init_fcc(){

  double aa = sim_box_info.cell_size;
  if (aa/sqrt(2.0) < 1) part_init_err(); // Check nearest neighbor distance
  int partid = 0;
  for (int ii=0; ii<sim_box_info.cell_x; ii++){
    for (int jj=0; jj<sim_box_info.cell_y; jj++){
      for (int kk=0; kk<sim_box_info.cell_z; kk++){
        add_particle(partid, ii*aa, jj*aa, kk*aa);
	add_particle(partid+1, (ii+0.5)*aa, (jj+0.5)*aa, kk*aa);
	add_particle(partid+2, (ii+0.5)*aa, jj*aa, (kk+0.5)*aa);
	add_particle(partid+3, ii*aa, (jj+0.5)*aa, (kk+0.5)*aa);
	partid += 4;
      }
    }
  }
}


void add_particle(int id, double xx, double yy, double zz){

  part[id][0] = id;
  part[id][1] = xx;
  part[id][2] = yy;
  part[id][3] = zz;

  part_test[id][0] = id;
  part_test[id][1] = xx;
  part_test[id][2] = yy;
  part_test[id][3] = zz;

}

void part_init_err(){

  printf("Overlap in the initial configuration. Possible solutions:\n");
  printf("-- If SC lattice was selected, try to change to FCC\n");
  printf("-- If FCC lattice was selected, the selected value of density is unphysical\n");
  exit(EXIT_FAILURE);

}

// ------ Functions to access the simulation box and the particle's information ------

struct box_info sim_box_info_get(){
  return sim_box_info;
}

struct p_info part_info_get(){
  return part_info;
}

/* double (*part_config_get())[4]{ */
/*   return part_test;  */
/* } */
part_array *part_config_get(){
  return part_test;
}


void print_sim_info(){

  if (G_IN.press > 0){
    // Output for NpT simulations
    printf("Simulation box size (x, y, z): %.5f %.5f %.5f\n", sim_box_info.lx,
	   sim_box_info.ly, sim_box_info.lz);
    printf("Number of particles: %d\n", part_info.NN);
    printf("Pressure: %.8f\n", G_IN.press);   
  }
  else if (G_IN.cavity_pcav > 0){
    // Output for NVT cavity simulations
      printf("Simulation box size (x, y, z): %.5f %.5f %.5f\n", sim_box_info.lx,
	     sim_box_info.ly, sim_box_info.lz);
      printf("Number of particles: %d\n", part_info.NN);  
  }
  else {
    // Output for NVT simulations
      printf("Simulation box size (x, y, z): %.5f %.5f %.5f\n", sim_box_info.lx,
         sim_box_info.ly, sim_box_info.lz);
      printf("Number of particles: %d\n", part_info.NN);   
  }


}

// ------ Functions to write and read simulation box information to file ------

void sim_box_info_write(FILE *fid){
  fwrite(&sim_box_info, sizeof(struct box_info), 1, fid);
}

void sim_box_info_read(FILE *fid){
  fread(&sim_box_info, sizeof(struct p_info), 1, fid);
}


// ------ Functions to write and read particles information to file ------

void part_info_write(FILE *fid){
  fwrite(&part_info, sizeof(struct p_info), 1, fid);
}

void part_info_read(FILE *fid){
  fread(&part_info, sizeof(struct p_info), 1, fid);
}


// ------ Function to free the matrix that stores the configuration -------

void part_free(){
  free(part);
  free(part_test);
}





