#include "cell_list.h"

int cl_cell_num;
int (*cl_neigh)[27];
int *cl_head, *cl_link;

void cell_list_init(){

  // Number of cells (use the number of cells employed to construct the simulation box)
  cl_cell_num = sim_box_info.cell_x * sim_box_info.cell_y * sim_box_info.cell_z;

  // Allocate linked list and neighbor matrix
  cell_list_alloc();

  // Initialize linked list
  cell_list_update();
  
  // Initialize neighbor matrix
  cell_neigh_init();

}


void cell_list_alloc(){

  // Arrays for linked lists
  cl_head = (int*)malloc(sizeof(int) * cl_cell_num);
  cl_link = (int*)malloc(sizeof(int) * part_info.NN);
  if (cl_head == NULL || cl_link == NULL){
    printf("ERROR: Failed linked list allocation\n");
    exit(EXIT_FAILURE);
  }

  // Allocate matrix to store cell neighbors
  cl_neigh = malloc(cl_cell_num * sizeof(*cl_neigh));
  if (cl_neigh == NULL){
    printf("ERROR: Failed cell list neighbor allocation\n");
    exit(EXIT_FAILURE);
  }

}


void cell_neigh_init(){

  for (int ii=0; ii<sim_box_info.cell_x; ii++){
    for (int jj=0; jj<sim_box_info.cell_y; jj++){
      for (int kk=0; kk<sim_box_info.cell_z; kk++){
  	neigh_id(ii,jj,kk);
      }
    }
  }

}

void cell_list_update(){

  int idx;
  for (int ii=0; ii<cl_cell_num; ii++){
    cl_head[ii] = 0;
  }
  for (int ii=1; ii<=part_info.NN; ii++){
    idx = cell_part_idx(ii-1);
    cl_link[ii] = cl_head[idx];
    cl_head[idx] = ii;
  }
  
}

int cell_part_idx(int id){

  return  (int)(part[id][1]/sim_box_info.cell_size)*sim_box_info.cell_x*sim_box_info.cell_x 
          + (int)(part[id][2]/sim_box_info.cell_size)*sim_box_info.cell_y 
          + (int)(part[id][3]/sim_box_info.cell_size);

}

void neigh_id(int ref_idx_x, int ref_idx_y, int ref_idx_z){

  int idx_x, idx_y, idx_z;
  int counter=0;
  int ref_idx = ref_idx_x*sim_box_info.cell_x*sim_box_info.cell_x + 
                ref_idx_y*sim_box_info.cell_y + 
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
	if (idx_x>sim_box_info.cell_x-1) idx_x -= sim_box_info.cell_x;
	else if (idx_x<0) idx_x += sim_box_info.cell_x;
	if (idx_y>sim_box_info.cell_y-1) idx_y -= sim_box_info.cell_y;
	else if (idx_y<0) idx_y += sim_box_info.cell_y;
	if (idx_z>sim_box_info.cell_z-1) idx_z -= sim_box_info.cell_z;
	else if (idx_z<0) idx_z += sim_box_info.cell_z;
	// Update neighbor matrix
	cl_neigh[ref_idx][counter] = 
	  idx_x*sim_box_info.cell_x*sim_box_info.cell_x + 
	  idx_y*sim_box_info.cell_y + 
	  idx_z;
	counter++;

      }
    }
  } 
    
}
