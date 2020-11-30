#ifndef CELL_LIST_H
#define CELL_LIST_H

extern int cl_cell_num, cl_neigh_num;
extern int *cl_head, *cl_link;
extern int (*cl_neigh)[27];

void cell_list_init();
void cell_list_alloc();
void cell_list_update();
int cell_part_idx(int id);
void cell_neigh_init();
void neigh_id(int ref_idx_x, int ref_idx_y, int ref_idx_z);

#endif
