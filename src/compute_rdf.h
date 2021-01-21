#ifndef COMPUTE_RDF_H
#define COMPUTE_RDF_H

#include <stdbool.h>

void compute_rdf(bool init, int sweep);

void rdf_hist_init();

void rdf_hist_alloc();

void rdf_hist_free();

void rdf_hist_compute();

void rdf_hist_norm();

void rdf_output(bool init, int sweep);

#endif
