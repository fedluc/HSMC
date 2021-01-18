#ifndef COMPUTE_PRESS_H
#define COMPUTE_PRESS_H

#include <stdbool.h>

void compute_pressv(bool init);

void pressv_hist_init();

void pressv_hist_alloc();

void pressv_hist_free();

void pressv_compute_hist();

void pressv_compute_rdf();

void pressv_output(bool init);

void compute_presst(bool init);

void presst_hist_init();

void presst_hist_alloc();

void presst_hist_free();

void presst_compute_hist();

void presst_output(bool init);

#endif
