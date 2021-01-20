#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include "read_input.h"
#include "analytic.h"

// --------------------------------------------------------
// The module "analytic.c" provides a set of functions
// which can be employed to evaluate some known analytical
// results for the hard sphere system. The analytical 
// results computed with this module include:
//     --- The Grundke-Henderson parametrization for the 
//         cavity function
//     --- The Verlet-Weis parametrization for the  
//         radial distribution function (up to r = 2*sigma)
//     --- The Percus-Yevick solution for the
//         radial distribution function (up to r = 2*simga)
// --------------------------------------------------------

// ------ Global variables ------

// Coefficients for the Grundke-Henderson parametrization of the cavity function
static double ay0, ay1, ay2, ay3;
// Coefficients for the Verlet-Weiss parametrization of the rdf
static double AA, mu, etaW, dW;
// Coefficients for the Percus-Yevick solution of the rdf
static double complex t0, t1, t2;
static double complex a00, a01, a02;

// ------ Grundke-Henderson parametrization of the cavity function ------

double lny_gh(double xx, double eta, double sigma, bool init_coeff){

  // Initialize coefficients if necessary
  if (init_coeff) coeff_gh(eta);

  // Compute cavity function
  double lny;
  double xx_hs = xx/sigma;
  if (xx_hs < 1.0){
    lny = ay0 + ay1*xx_hs + ay2*xx_hs*xx_hs + ay3*xx_hs*xx_hs*xx_hs;
  }
  else {
    lny =  log(rdf_vw(xx, eta, sigma, false));
  }

  // Output
  return lny;

}


// ------ Coefficients for Grundke-Henderson parametrization of the cavity function ------

void coeff_gh(double eta){

  // Compute coefficients for the Verlet-Weiss paremetrization of the rdf  
  coeff_vw(eta);

  // Packing fraction and related quantities
  double eta2 = eta*eta, eta3 = eta2*eta;
  double ometa = 1. - eta, ometa3 = ometa*ometa*ometa;
  double ometaW = 1. - etaW, ometaW3 = ometaW*ometaW*ometaW;

  // Coefficients
  double g1 = rdf_py(1.0, etaW, dW, false) + AA;
  double dg1 = -(9./2.) * etaW * (1+etaW)/ometaW3 - AA*(1+mu);
  ay0 = (8.*eta -9.*eta2 + 3.*eta3)/ometa3;
  ay1 = -3. *eta * (2.-eta)/ ometa3;
  ay2 = -(dg1/g1) - 3.*ay0 - 2.*ay1 + 3.*log(g1);
  ay3 = 1./3.*( (dg1/g1) - ay1 - 2.*ay2);

}

// ------ Verlet-Weiss parametrization of the rdf ------

double rdf_vw(double xx, double eta, double sigma, bool init_coeff){
  
  // Initialize coefficients if necessary
  if (init_coeff) coeff_vw(eta);

  // Compute rdf
  double delta_g1, rdf;
  double xx_hs = xx/sigma;
  if (xx_hs < 1.0){
    rdf = 0.0;
  }
  else {
    delta_g1 = (AA/xx_hs) * exp(-mu*(xx_hs-1)) * cos(mu*(xx_hs-1));
    rdf = rdf_py(xx_hs, etaW, dW, false) + delta_g1;
  }

  // Output 
  return rdf;

}


// ------ Coefficients for Verlet-Weiss parametrization of the rdf ------

void coeff_vw(double eta){
  
  // Packing fraction and related quantities
  double eta2 = eta*eta;
  double eta3 = eta2*eta, eta4 = eta2*eta2;
  double ometa = 1. - eta, ometa3 = ometa*ometa*ometa;
  double ometa4 = ometa3*ometa;

  // Compressibility factor from Carnahan-Starling EOS
  double zz_cs = (1 + eta + eta2 - eta3)/ometa3;

  // Inverse compressibility from Carnahan-Starling EOS
  double mu_cs = ometa4/(1+ 4.*eta + 4.*eta2 - 4.*eta3 + eta4);
 
  // Effective hard-sphere radius and packing fraction for Percus-Yevick solution
  dW = pow(1-eta/16.,1./3.);
  double dW3 = dW*dW*dW, dW4 = dW3*dW;
  etaW = eta*dW3;
  double op2etaW2 = (1 + 2*etaW)*(1 + 2*etaW);
  double ometaW = 1. - etaW;
  double ometaW3 = ometaW*ometaW*ometaW;
  double ometaW4 = ometaW3*ometaW;

  // Inverse compressibility from the exact Percus-Yevick solution 
  double mu_py = ometaW4/op2etaW2;

  // Compute coefficients for Percus-Yevick solution for the rdf
  coeff_py(etaW);

  // Verlet-Weis parametrization for the radial distribution function
  AA = (zz_cs-1.0)/(4.0*eta) - rdf_py(1.0, etaW, dW, false);
  mu = 12 * AA * eta / (mu_cs - mu_py +
			       (etaW * (27.0*etaW*(1+etaW) - 8.0*dW*op2etaW2
					+ dW4*(8.0+ 5.0*etaW + 5.0*etaW*etaW)))
			       /(dW4 * (-ometaW3)));

}

// ------ Percus-Yevick solution of the rdf ------

double rdf_py(double xx, double eta, double sigma, bool init_coeff){
  
  // Initialize coefficients if necessary
  if (init_coeff) coeff_py(eta);

  // Compute rdf
  double complex g1;
  double complex xx_hs = xx/sigma;
  double rdf;
  if (creal(xx_hs) < 1.0){
    rdf = 0.0;
  }
  else {
    g1 = a00*cexp(t0*(xx_hs-1)) + a01*cexp(t1*(xx_hs-1))
          + a02*cexp(t2*(xx_hs-1));
    if (cimag(g1) > 1e-15){
      printf("WARNING: Imaginary part of PYHS solution is larger than threshold (%.8e)", cimag(g1));
    }
    rdf = creal(g1)/creal(xx_hs);
  }

  // Output 
  return rdf;

}

// ------ Coefficients for Percus-Yevick solution of the rdf ------

void coeff_py(double eta){

  // Quantities related to the packing fraction
  double complex eta2 = eta*eta;
  double complex ometa = 1.0 - eta, ometa2 = ometa*ometa;
  double complex opeta_2 = 1.0 + eta/2.0;
  double complex op2eta = 1.0 + 2*eta;

  // Note: Here we consider the solution only up to x = r/sigma = 2.0
  double complex ff = 3.0 + 3.*eta - eta2;
  double complex eta2_ff2 = (eta2/ff)*(eta2/ff);
  double complex tetaff = cpow(2*eta*ff,1./3.);
  double complex yp = cpow(sqrt(1.+ 2.*eta2_ff2) + 1.,1./3.);
  double complex ym = -cpow(sqrt(1.+ 2.*eta2_ff2) - 1.,1./3.); // WARNING: strange behavior if cpow has negative argument
  double complex jj = cexp((2./3.)* M_PI * I);
  double complex jj1 = jj;
  double complex jjm1 = 1./jj1;
  double complex jj2 = jj*jj;
  double complex jjm2 = 1./jj2;
  t0 = (-2*eta + tetaff*(yp+ym))/ometa;
  t1 = (-2*eta + tetaff*(yp*jj1+ym*jjm1))/ometa;
  t2 = (-2*eta + tetaff*(yp*jj2+ym*jjm2))/ometa;
  double complex L0 = opeta_2*t0 + op2eta;
  double complex L1 = opeta_2*t1 + op2eta;
  double complex L2 = opeta_2*t2 + op2eta;
  double complex S10 = 3*ometa2*t0*t0 + 12.*eta*ometa*t0 + 18.*eta2;
  double complex S11 = 3*ometa2*t1*t1 + 12.*eta*ometa*t1 + 18.*eta2;
  double complex S12 = 3*ometa2*t2*t2 + 12.*eta*ometa*t2 + 18.*eta2;
  a00 = t0 * L0 / S10;
  a01 = t1 * L1 / S11;
  a02 = t2 * L2 / S12;

}

