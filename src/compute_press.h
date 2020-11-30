#ifndef PRESS_H
#define PRESS_H

#include <stdbool.h>

void compute_pressv(bool init);

void pressv_hist_init();

void pressv_compute_hist();

void pressv_compute_rdf();

void compute_presst(bool init);

void presst_hist_init();

void presst_compute_hist();


#endif
