#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "init.h"
#include "read_input.h"
#include "cell_list.h"

static int neigh_num;
static int cell_num_x;
static int cell_num_y;
static int cell_num_z;
static double cell_size_x;
static double cell_size_y;
static double cell_size_z;

int *cell_list_alloc(int rows, int cols) {
  int (*arr)[cols] = malloc(sizeof(*arr) * rows);
  if (arr == NULL) {
    printf("ERROR: Failed cell list allocation\n");
    exit(EXIT_FAILURE);
  }
  return &arr[0][0];
}

void cell_list_init(int neigh_num, int (*neigh_mat)[neigh_num],
		    int max_part_cell, int (*part_cell)[max_part_cell]){

  // Total number of cells
  int cell_num_tot = cell_num_tot = cell_num_x * cell_num_y * cell_num_z;

  // Initialize cell list
  cell_list_new(cell_num_tot, in.neigh_max_part, part_cell);
  
  // Initialize neighbor matrix
  cell_neigh_init(cell_num_tot, neigh_num, neigh_mat);

}

void cell_list_new(int cell_num_tot, int max_part, int part_cell[cell_num_tot][max_part]){

  // Define number and size of the cells
  compute_cell_list_info();
 
  // Check for consistency
  if (cell_size_x < 1.0 || cell_size_y < 1.0 || cell_size_z < 1.0) {
    printf("ERROR: size of cells in cell list is smaller than the size of the particles\n");
    exit(EXIT_FAILURE);
  }

  // Initialize the cell lists
  for (int ii=0; ii<cell_num_tot; ii++){
    part_cell[ii][0] = 0;
    for (int jj=1; jj<max_part; jj++){
      part_cell[ii][jj] = -1;
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
    part_cell[idx_row][0] += 1;
    cell_list_check(idx_row, cell_num_tot, max_part, part_cell);
    // Update id of particles in cell
    part_cell[idx_row][idx_col[idx_row]] = ii;
    idx_col[idx_row] += 1;
  }
  free(idx_col);
  
}

void compute_cell_list_info(){

  // Number of cells
  cell_num_x = (int)floor(sim_box_info.lx/in.neigh_dr);
  cell_num_y = (int)floor(sim_box_info.ly/in.neigh_dr);
  cell_num_z = (int)floor(sim_box_info.lz/in.neigh_dr);
  int cell_num_tot = cell_num_x * cell_num_y * cell_num_z;
  
  // Do not allow less than 27 cells
  if (cell_num_tot < 27){
    printf("WARNING: The specified cell size has been reduced in order to accomodate 27 cells in the simulation box\n");
    cell_num_x = 3;
    cell_num_y = 3;
    cell_num_z = 3;
  }

  // Number of neighbor per cell
  neigh_num = 27;

  // Cell size
  cell_size_x = sim_box_info.lx/cell_num_x;
  cell_size_y = sim_box_info.ly/cell_num_y;
  cell_size_z = sim_box_info.lz/cell_num_z;

  
}

void cell_list_update(int cell_idx_del, int cell_idx_add, int part_idx,
		      int cell_num_tot, int max_part, int part_cell[cell_num_tot][max_part]){
  cell_list_del(cell_idx_del, part_idx, cell_num_tot, max_part, part_cell);
  cell_list_add(cell_idx_add, part_idx, cell_num_tot, max_part, part_cell);
}

void cell_list_add(int cell_idx, int part_idx,
		   int cell_num_tot, int max_part, int part_cell[cell_num_tot][max_part]){

  int n_part_cell = part_cell[cell_idx][0];
  part_cell[cell_idx][0] += 1;
  cell_list_check(cell_idx, cell_num_tot, max_part, part_cell);
  part_cell[cell_idx][n_part_cell+1] = part_idx;
  
}

void cell_list_del(int cell_idx, int part_idx,
		   int cell_num_tot, int max_part, int part_cell[cell_num_tot][max_part]){

  int n_part_cell = part_cell[cell_idx][0];
  int idx_remove = n_part_cell;
  bool shift_flag = false;


  part_cell[cell_idx][0] -= 1;
  cell_list_check(cell_idx, cell_num_tot, max_part, part_cell);
  for (int ii=1; ii<=n_part_cell; ii++){
    if (part_cell[cell_idx][ii] == part_idx){
      idx_remove = ii;
      part_cell[cell_idx][ii] = -1;
      shift_flag = true;
    }
    if (shift_flag && ii > idx_remove) {
      part_cell[cell_idx][ii-1] = part_cell[cell_idx][ii];
      part_cell[cell_idx][ii] = -1;
    }
  }
  if (!shift_flag){
    printf("ERROR: Attempt to delete particle %d from cell %d, but particle is not in cell\n",
	   part_idx, cell_idx);
    exit(EXIT_FAILURE);
  }
  
}


void cell_list_check(int cell_idx, 
		     int cell_num_tot, int max_part, int part_cell[cell_num_tot][max_part]){
  if (part_cell[cell_idx][0] < 0){
    printf("ERROR: Trying to remove one particle from an empty cell\n");
    exit(EXIT_FAILURE);
  }
  if (part_cell[cell_idx][0] > max_part){
    printf("ERROR: More than %d particles in one cell\n", max_part);
    exit(EXIT_FAILURE);
  }
}

int cell_part_idx(int id){

  return  (int)(part[id][1]/cell_size_x)*cell_num_x*cell_num_x
          + (int)(part[id][2]/cell_size_y)*cell_num_y
          + (int)(part[id][3]/cell_size_z);

}

void cell_neigh_init(int cell_num_tot, int neigh_num, int neigh_mat[cell_num_tot][neigh_num]){

  for (int ii=0; ii<cell_num_x; ii++){
    for (int jj=0; jj<cell_num_y; jj++){
      for (int kk=0; kk<cell_num_z; kk++){
  	neigh_id(ii,jj,kk,cell_num_tot,neigh_num,neigh_mat);
      }
    }
  }

}


void neigh_id(int ref_idx_x, int ref_idx_y, int ref_idx_z,
	      int cell_num_tot, int neigh_num, int neigh_mat[cell_num_tot][neigh_num]){

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
	neigh_mat[ref_idx][counter] =
	  idx_x*cell_num_x*cell_num_x +
	  idx_y*cell_num_y +
	  idx_z;
	counter++;

      }
    }
  }
    
}

void get_cell_list_info(int *num_neigh_cell, int *num_tot,
			int *num_x, int *num_y, int *num_z,
			double *size_x, double *size_y, double *size_z){
  
  if (num_neigh_cell != NULL) *num_neigh_cell = neigh_num;
  if (num_tot != NULL) *num_tot = cell_num_x * cell_num_y * cell_num_z;
  if (num_x != NULL) *num_x = cell_num_x;
  if (num_y != NULL) *num_y = cell_num_y;
  if (num_z != NULL) *num_z = cell_num_z;
  if (size_x != NULL) *size_x = cell_size_x;
  if (size_y != NULL) *size_y = cell_size_y;
  if (size_z != NULL) *size_z = cell_size_z;

}
