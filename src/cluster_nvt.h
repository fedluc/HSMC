#ifndef CLUSTER_NVT_H
#define CLUSTER_NVT_H

void cluster_hs_nvt();

void cluster_run_nvt(bool prod_flag, int sweep_offset);

void cluster_sweep_nvt();

void cluster_move();

void cluster_move_part(int idx, double aa, double bb, double cc,
                       double uu, double vv, double ww);

void rotate(double xx, double yy, double zz,
            double aa, double bb, double cc,
            double uu, double vv, double ww,
            double theta, double *xx_new,
            double *yy_new, double *zz_new);

void apply_pbc(int idx);

void pocket_add(int idx_ref, int *pocket, int *part_pocket);

void pocket_del(int idx_ref, int *pocket, int *part_pocket);

#endif
