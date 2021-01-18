#ifndef COMPUTE_ORDER_PARAMETER_H
#define COMPUTE_ORDER_PARAMETER_H

#include <stdbool.h>

void compute_op(bool init);

void global_ql_compute();

void ql_alloc();

void ql_free();

void ql_compute();

void qlm2_compute(int ref_idx);

void global_ql_output(bool init);

#endif
