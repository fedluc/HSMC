#include "init.h"

struct box_info sim_box_info;
struct p_info part_info;

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


void part_alloc(double ***part){

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
  part_info.Ncell = part_cell;
  part_info.NN = part_tot;

}


void part_init(double **part){

  if (sim_box_info.cell_type == 1){
    part_init_sc(part);
  }
  else {
    part_init_fcc(part);
  }

}

void part_init_sc(double **part){

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

void part_init_fcc(double **part){

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
