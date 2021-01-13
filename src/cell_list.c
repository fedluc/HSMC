#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"

int cl_neigh_num;
int (*cl_neigh)[27];
int (*cl_part_cell)[11];

static int cell_num_tot;
static int cell_num_x;
static int cell_num_y;
static int cell_num_z;
static double cell_size_x;
static double cell_size_y;
static double cell_size_z;

void cell_list_init(){

  // Number of cells
  cell_num_x = compute_cell_num(sim_box_info.lx, in.neigh_dr);
  cell_num_y = compute_cell_num(sim_box_info.ly, in.neigh_dr);
  cell_num_z = compute_cell_num(sim_box_info.lz, in.neigh_dr);
  cell_num_tot = cell_num_x * cell_num_y * cell_num_z;

  // Number of neighbor per cell
  cl_neigh_num = 27;
  if (cell_num_tot < 27) cl_neigh_num = cell_num_tot;

  // Allocate linked list and neighbor matrix
  cell_list_alloc();

  // Initialize linked list
  cell_list_new();
  
  // Initialize neighbor matrix
  cell_neigh_init();

}

int compute_cell_num(double sim_box_len, double dr){
  
  double rem = fmod(sim_box_len,dr); 
  return (int)(sim_box_len - rem); 

}

void cell_list_alloc(){

  // Allocate matrix to store cell neighbors
  cl_neigh = malloc(cell_num_tot * sizeof(*cl_neigh));
  if (cl_neigh == NULL){
    printf("ERROR: Failed cell list neighbor allocation\n");
    exit(EXIT_FAILURE);
  }

  // Allocate matrix to store particles in each cell
  cl_part_cell = malloc(cell_num_tot * sizeof(*cl_part_cell));
  if (cl_part_cell == NULL){
    printf("ERROR: Failed cell list particles allocation\n");
    exit(EXIT_FAILURE);
  }

}

void cell_list_new(){

  // Cell size
  cell_size_x = sim_box_info.lx/cell_num_x;
  cell_size_y = sim_box_info.ly/cell_num_y;
  cell_size_z = sim_box_info.lz/cell_num_z;

  // Check for consistency
  if (cell_size_x < 1.0 || cell_size_y < 1.0 || cell_size_z < 1.0) {
    printf("ERROR: size of cells in cell list is smaller than the size of the particles\n");
    exit(EXIT_FAILURE);
  }

  // Initialize the cell lists
  for (int ii=0; ii<cell_num_tot; ii++){
    cl_part_cell[ii][0] = 0;
    for (int jj=1; jj<=10; jj++){
      cl_part_cell[ii][jj] = -1;
    }
  }

  // Fill the cell list with the particles indexes
  int idx_row;
  int *idx_col  = (int*)malloc(sizeof(int) * cell_num_tot);
  for (int ii=0; ii<cell_num_tot; ii++){
    idx_col[ii] = 1;
  }
  for (int ii=0; ii<part_info.NN; ii++){
    // Index of the cell that contains the particle
    idx_row = cell_part_idx(ii);
    // Update number of particles in cell
    cl_part_cell[idx_row][0] += 1;
    cell_list_check(idx_row);
    // Update id of particles in cell
    cl_part_cell[idx_row][idx_col[idx_row]] = ii;
    idx_col[idx_row] += 1;
  }
  free(idx_col);
  
}

void cell_list_update(int cell_idx_del, int cell_idx_add, int part_idx){
  cell_list_del(cell_idx_del, part_idx);
  cell_list_add(cell_idx_add, part_idx);
}

void cell_list_add(int cell_idx, int part_idx){

  int n_part_cell = cl_part_cell[cell_idx][0];
  cl_part_cell[cell_idx][0] += 1;
  cell_list_check(cell_idx);
  cl_part_cell[cell_idx][n_part_cell+1] = part_idx;
  
}

void cell_list_del(int cell_idx, int part_idx){

  int n_part_cell = cl_part_cell[cell_idx][0];
  int idx_remove = n_part_cell;
  bool shift_flag = false;


  cl_part_cell[cell_idx][0] -= 1;
  cell_list_check(cell_idx);
  for (int ii=1; ii<=n_part_cell; ii++){
    if (cl_part_cell[cell_idx][ii] == part_idx){ 
      idx_remove = ii;
      cl_part_cell[cell_idx][ii] = -1;
      shift_flag = true;
    }
    if (shift_flag && ii > idx_remove) { 
      cl_part_cell[cell_idx][ii-1] = cl_part_cell[cell_idx][ii];
      cl_part_cell[cell_idx][ii] = -1;
    }
  }
  if (!shift_flag){
    printf("ERROR: Attempt to delete particle %d from cell %d, but particle is not in cell\n", 
	   part_idx, cell_idx);
    exit(EXIT_FAILURE);
  }
  
}


void cell_list_check(int cell_idx){
  if (cl_part_cell[cell_idx][0] < 0){ 
    printf("ERROR: Trying to remove one particle from an empty cell\n");
    exit(EXIT_FAILURE);
  }
  if (cl_part_cell[cell_idx][0] > 10){ 
    printf("ERROR: More than %d particles in one cell", 10);
    exit(EXIT_FAILURE);
  }
}

int cell_part_idx(int id){

  return  (int)(part[id][1]/cell_size_x)*cell_num_x*cell_num_x 
          + (int)(part[id][2]/cell_size_y)*cell_num_y 
          + (int)(part[id][3]/cell_size_z);

}

void cell_neigh_init(){

  for (int ii=0; ii<cell_num_x; ii++){
    for (int jj=0; jj<cell_num_y; jj++){
      for (int kk=0; kk<cell_num_z; kk++){
  	neigh_id(ii,jj,kk);
      }
    }
  }

}

void neigh_id(int ref_idx_x, int ref_idx_y, int ref_idx_z){

  int idx_x, idx_y, idx_z;
  int counter=0;
  int ref_idx = ref_idx_x*cell_num_x*cell_num_x + 
                ref_idx_y*cell_num_y + 
                ref_idx_z;
  
  // Loop over neighboring cells (including the reference cell)
  for (int ii=-1; ii<2; ii++){
    for (int jj=-1; jj<2; jj++){
      for (int kk=-1; kk<2; kk++){

	// x,y,z components of the neighbor index
	idx_x = ref_idx_x + ii;
	idx_y = ref_idx_y + jj;
	idx_z = ref_idx_z + kk;
	// Apply periodic boundary conditions
	if (idx_x>cell_num_x-1) idx_x -= cell_num_x;
	else if (idx_x<0) idx_x += cell_num_x;
	if (idx_y>cell_num_y-1) idx_y -= cell_num_y;
	else if (idx_y<0) idx_y += cell_num_y;
	if (idx_z>cell_num_z-1) idx_z -= cell_num_z;
	else if (idx_z<0) idx_z += cell_num_z;
	// Update neighbor matrix
	cl_neigh[ref_idx][counter] = 
	  idx_x*cell_num_x*cell_num_x + 
	  idx_y*cell_num_y + 
	  idx_z;
	counter++;

      }
    }
  } 
    
}

void get_cell_list_info(int *num_x, int *num_y, int *num_z,
			double *size_x, double *size_y, double *size_z){
  *num_x = cell_num_x;
  *num_y = cell_num_y;
  *num_z = cell_num_z;
  *size_x = cell_size_x;
  *size_y = cell_size_y;
  *size_z = cell_size_z;

}
