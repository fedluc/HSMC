#include "hsmc.h"

void cell_list_init(){

  // Variable declaration
  double cl_cell_size, cl_cell_num;
  int cell_idx;
  int *cl_head, *cl_list;

  // Number of cells (use the number of cells employed to construct the simulation box)
  cl_cell_num = sim_box_info.cell_x * sim_box_info.cell_y * sim_box_info.cell_z;

  /* // Arrays for linked lists */
  /* cl_head = (double*)malloc(sizeof(double) * cl_cell_num); */
  /* cl_list = (double*)malloc(sizeof(double) * in.NN); */

  /* // Initialize the values of the linked list */
  /* for (int ii=0; ii<cl_cell_num; ii++){ */
  /*   cell_head[ii] = 0; */
  /* } */
  
  /* for (int ii=0; ii<sim_box_info.cell_x; ii++){ */
  /*   for (int jj=0; jj<sim_box_info.cell_y; jj++){ */
  /*     for (int kk=0; kk<sim_box_info.cell_z; kk++){ */
  /* 	cell_idx = ii*sim_box_info.cell_x*sim_box_info.cell_x */
  /* 	           + jj*sim_box_info.cell_y */
  /* 	           + kk; */
  /* 	printf("%d\n",cell_idx); */
  /*     } */
  /*   } */
  /* } */
 

 
 /*  for (int ii=0; ii<sim_box_info.cell_x; ii++){ */
 /*    for (int jj=0; jj<sim_box_info.cell_y; jj++){ */
 /*      cell_idx = ii*sim_box_info.cell_x + jj; */
 /*      printf("%d\n",cell_idx); */
 /*    } */
 /* } */
  
  int cell_idx_1, cell_idx_2, cell_idx_3, cell_idx_4;
  int cell_idx_5, cell_idx_6, cell_idx_7, cell_idx_8;
  int iip1, iim1, jjp1, jjm1;
  for (int ii=0; ii<sim_box_info.cell_x; ii++){
    for (int jj=0; jj<sim_box_info.cell_y; jj++){
      cell_idx = ii*sim_box_info.cell_x + jj;
      iip1 = ii+1;
      iim1 = ii-1;
      jjp1 = jj+1;
      jjm1 = jj-1;
      if (iip1>sim_box_info.cell_x-1) iip1 -= sim_box_info.cell_x;
      if (iim1<0) iim1 += sim_box_info.cell_x;
      if (jjp1>sim_box_info.cell_y-1) jjp1 -= sim_box_info.cell_y;
      if (jjm1<0) jjm1 += sim_box_info.cell_y;
      cell_idx_1 = ii*sim_box_info.cell_x + jjp1;
      cell_idx_2 = ii*sim_box_info.cell_x + jjm1;
      cell_idx_3 = iim1*sim_box_info.cell_x + jj;
      cell_idx_4 = iim1*sim_box_info.cell_x + jjp1;
      cell_idx_5 = iim1*sim_box_info.cell_x + jjm1;
      cell_idx_6 = iip1*sim_box_info.cell_x + jj;
      cell_idx_7 = iip1*sim_box_info.cell_x + jjp1;
      cell_idx_8 = iip1*sim_box_info.cell_x + jjm1;
      printf("Cell %d interacts with cells %d %d %d %d %d %d %d %d\n",
	     cell_idx,cell_idx_1, cell_idx_2, cell_idx_3, cell_idx_4,
	     cell_idx_5, cell_idx_6, cell_idx_7, cell_idx_8);
    }
 }

  printf("----------------------------\n");
  int idx_1, idx_2;
  int ii_tmp[] = {0, 0, -1, -1, -1, 1, 1, 1};
  int jj_tmp[] = {1, -1, 0, 1, -1, 0, 1, -1};
  for (int ii=0; ii<sim_box_info.cell_x; ii++){
    for (int jj=0; jj<sim_box_info.cell_y; jj++){
      cell_idx = ii*sim_box_info.cell_x + jj;
      printf("Cell %d interacts with cells ", cell_idx);
      for (int nn=0; nn<8; nn++){
	idx_1 = ii + ii_tmp[nn];
	idx_2 = jj + jj_tmp[nn];
	if (idx_1>sim_box_info.cell_x-1) idx_1 -= sim_box_info.cell_x;
	else if (idx_1<0) idx_1 += sim_box_info.cell_x;
	if (idx_2>sim_box_info.cell_y-1) idx_2 -= sim_box_info.cell_y;
	else if (idx_2<0) idx_2 += sim_box_info.cell_y;
	cell_idx_1 = idx_1*sim_box_info.cell_x + idx_2;
	printf("%d ", cell_idx_1);
      }
      printf("\n");
    }
 }

}
