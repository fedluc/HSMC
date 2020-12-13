#ifndef COMPUTE_ORDER_PARAMETER_H
#define COMPUTE_ORDER_PARAMETER_H

#include <stdbool.h>

void compute_op(bool init);

void global_ql_compute();

void ql_compute();

void qlm_compute(int ref_idx);

void global_ql_output(bool init);

#endif
