#ifndef OPTIMIZER_H
#define OPTIMIZER_H

void opt_nvt();
void get_sample_nvt(double *dr, double *acc_ratio, 
		    int sample_iter);
void opt_npt();
void get_sample_npt(double *dr, double *acc_ratio_part,
                    double *dv, double *acc_ratio_vol,
                    int sample_iter);

void opt_cavity_nvt();
void get_sample_cavity_nvt(double *dr, double *acc_ratio, 
			   int sample_iter);


#endif
