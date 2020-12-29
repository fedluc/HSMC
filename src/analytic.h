#ifndef ANALYTIC_H
#define ANALYTIC_H

double lny_gh(double xx, double eta, double sigma, bool init_coeff);

void coeff_gh(double eta);

double rdf_vw(double xx, double eta, double sigma, bool init_coeff);

void coeff_vw(double eta);

double rdf_py(double xx, double eta, double sigma, bool init_coeff);
  
void coeff_py(double eta);

#endif

