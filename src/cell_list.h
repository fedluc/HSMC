#ifndef CELL_LIST_H
#define CELL_LIST_H

#include "read_input.h"
#include "init.h"

extern int cl_cell_num;
extern int (*cl_neigh)[27];
extern int *cl_head, *cl_link;

void cell_list_init();
void cell_list_alloc();
void cell_list_update();
int cell_part_idx(int id);
void cell_neigh_init();
void neigh_id(int ref_idx_x, int ref_idx_y, int ref_idx_z);

#endif
