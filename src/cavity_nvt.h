#ifndef CAVITY_NVT_H
#define CAVITY_NVT_H

void cavity_hs_nvt();

void cavity_run_nvt(bool prod_flag, int sweep_offset);

void cavity_set_distance();

void cavity_sweep_nvt();

void cavity_psi_output();

void cavity_dist_output(bool init);

#endif
