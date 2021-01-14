#ifndef CELL_LIST_H
#define CELL_LIST_H

void cell_list_init();

void cell_list_free();

void cell_list_alloc_arr(int ***arr, int rows, int cols);

void cell_list_free_arr(int ***arr, int rows);

void compute_cell_list_info();

void get_cell_list_info(int ***pc, int ***nm,
			int *num_neigh_cell, int *num_tot,
			int *num_x, int *num_y, int *num_z,
			double *size_x, double *size_y, double *size_z);

void cell_list_new();

void cell_list_update(int cell_idx_del, int cell_idx_add, int part_idx);

void cell_list_add(int cell_idx, int part_idx);

void cell_list_del(int cell_idx, int part_idx);

void cell_list_check(int cell_idx);

int cell_part_idx(int id);

void cell_neigh_init();

void neigh_id(int ref_idx_x, int ref_idx_y, int ref_idx_z);

#endif
