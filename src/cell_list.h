#ifndef CELL_LIST_H
#define CELL_LIST_H

extern int cl_neigh_num;
extern int (*cl_neigh)[27];
extern int (*cl_part_cell)[11];

void cell_list_init();
int compute_cell_num(double sim_box_len, double dr);
void cell_list_alloc();
void cell_list_new();
void cell_list_update(int cell_idx_del, int cell_idx_add, int part_idx);
void cell_list_add(int cell_idx, int part_idx);
void cell_list_del(int cell_idx, int part_idx);
void cell_list_check();
int cell_part_idx(int id);
void cell_neigh_init();
void neigh_id(int ref_idx_x, int ref_idx_y, int ref_idx_z);
void get_cell_list_info(int *num_x, int *num_y, int *num_z,
                        double *size_x, double *size_y, double *size_z);

#endif
