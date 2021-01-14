#ifndef RNG_H
#define RNG_H

#include <stdio.h>

void rng_init();

double rng_get_double();

int rng_get_int(int MM);

void rng_free();

void rng_write(FILE *fid);

void rng_read(FILE *fid);

#endif
