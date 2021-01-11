#ifndef CELL_LIST_H
#define CELL_LIST_H

extern int cl_neigh_num;
extern int (*cl_neigh)[27];
extern int (*cl_part_cell)[11]; // The code stops if one cell has more than 10 particles

void cell_list_init();
void cell_list_alloc();
void cell_list_new();
void cell_list_update(int cell_idx_del, int cell_idx_add, int part_idx);
void cell_list_add(int cell_idx, int part_idx);
void cell_list_del(int cell_idx, int part_idx);
void cell_list_check();
int cell_part_idx(int id);
void cell_neigh_init();
void neigh_id(int ref_idx_x, int ref_idx_y, int ref_idx_z);

#endif
