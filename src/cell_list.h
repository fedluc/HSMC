#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <stdbool.h>

int *cell_list_alloc(int rows, int cols);

void cell_list_init(int neigh_num, int (*neigh_mat)[neigh_num],
                    int max_part_cell, int (*part_cell)[max_part_cell]);

void cell_list_new(int cell_num_tot, int max_part, int part_cell[cell_num_tot][max_part]);

void compute_cell_list_info();

int compute_cell_num(double sim_box_len, double dr);

void cell_list_update(int cell_idx_del, int cell_idx_add, int part_idx,
		      int cell_num_tot, int max_part, int part_cell[cell_num_tot][max_part]);

void cell_list_add(int cell_idx, int part_idx,
		   int cell_num_tot, int max_part, int part_cell[cell_num_tot][max_part]);

void cell_list_del(int cell_idx, int part_idx,
		   int cell_num_tot, int max_part, int part_cell[cell_num_tot][max_part]);

void cell_list_check(int cell_idx,
                     int cell_num_tot, int max_part, int part_cell[cell_num_tot][max_part]);

int cell_part_idx(int id);

void cell_neigh_init(int cell_num_tot, int neigh_num, int neigh_mat[cell_num_tot][neigh_num]);

void neigh_id(int ref_idx_x, int ref_idx_y, int ref_idx_z,
              int cell_num_tot, int neigh_num, int neigh_mat[cell_num_tot][neigh_num]);

void get_cell_list_info(int *num_neigh_cell, int *num_tot,
                        int *num_x, int *num_y, int *num_z,
                        double *size_x, double *size_y, double *size_z);

#endif
