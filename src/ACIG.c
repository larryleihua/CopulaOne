#include "include/ACIG.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>

void psi1(double *s0, double *alp0, double *galp0, double *rst)
{ 
  double sroot, tem, alp, s, galp;
  alp = *alp0;
  galp = *galp0;
  s = *s0;  
  sroot = sqrt(s);
  tem = gsl_sf_bessel_Knu(fabs(alp-1), 2*sroot);
  //galp = gsl_sf_gamma(alp);
  rst[0] =  -2 * pow(s, ((alp-1)/2)) * tem / galp;
}

void psi2(double *s0, double *alp0, double *galp0, double *rst)
{ 
  double sroot, tem, alp, s, galp;
  alp = *alp0;
  galp = *galp0;
  s = *s0;  
  sroot = sqrt(s);
  tem = gsl_sf_bessel_Knu(fabs(alp-2), 2*sroot);
  //galp = gsl_sf_gamma(alp);
  rst[0] =  2 * pow(s, ((alp-2)/2)) * tem / galp;
}

double _psi(double s, double alp, double galp)
{ 
  double sroot, tem;
  sroot = sqrt(s);
  tem = gsl_sf_bessel_Knu(alp, 2*sroot);
  //galp = gsl_sf_gamma(alp);
  return(GSL_MIN(1, 2 * pow(s, (alp/2)) * tem / galp));
}

void lnK(double *nu, double *x, int *N, double *rst)
{
  int i;
  for(i=0; i < N[0]; i++)
  {
  if(x[i] < 1e-100){rst[i] = 1e99;}
  else
  {
   rst[i] = gsl_sf_bessel_lnKnu(fabs(nu[i]), x[i]);
  }
  }
}

void lnG(double *alp, int *N, double *rst)
{
  int i;
  for(i=0; i < N[0]; i++)
  {
    rst[i] = gsl_sf_lngamma(alp[i]);
  }  
}


// optimized
void psi(double *ss, double *alp, int *N, double *rst)
{ 
  double sroot, tem, galp;
  int i;
 for(i = 0; i < N[0]; i++)
 {
   sroot = sqrt(ss[i]);
   tem = gsl_sf_bessel_Knu(fabs(alp[i]), 2*sroot);
   galp = gsl_sf_gamma(alp[i]);
   rst[i] =  GSL_MIN(1, 2 * pow(ss[i], (alp[i]/2)) * tem / galp);
 }
}

void invpsi(double *tt, double *alp, int *N, double *rst)
{
  double sroot, diff, tol, s, g, gp, logK;
  int iter, maxiter, i;
  maxiter = 100;
  tol = pow(10, -8);
for(i = 0; i < N[0]; i++)
{
  iter = 0;
  diff = 1.0;
  s = 1/pow(tt[i], 0.3) - 1; // starting value
  while(fabs(diff) > tol && iter < maxiter)
  {
   sroot = sqrt(s); 
   logK = gsl_sf_bessel_lnKnu(fabs(alp[i]), 2*sroot);
   g = log(2) + alp[i]*log(s)/2 + logK - gsl_sf_lngamma(alp[i]) - log(tt[i]);
   gp = alp[i]/(2*s) - exp( gsl_sf_bessel_lnKnu(fabs(alp[i]-1), 2*sroot) - log(2*sroot) - logK) 
                  - exp( gsl_sf_bessel_lnKnu(fabs(alp[i]+1), 2*sroot) - log(2*sroot) - logK);
   diff = g / gp;
   s = s - diff;
   while(s <= 0) { diff = diff/2; s = s+diff; }
   iter = iter+1;
  }
  rst[i] = s;
}
}

// non-optimized
void invpsi0(double *tt0, double *alp0, double *galp0, double *rst)
{
  double alp, galp, tt, alp1, diff, tol, s, g, sroot, tem, gp;
  int iter, maxiter;
  maxiter = 100;
  tol = pow(10, -8);
  alp = *alp0;
  galp = *galp0;
  //galp = gsl_sf_gamma(alp);
  tt = *tt0;  

  if (tt >= 0.9999999 ) {rst[0]= 0.0000001; rst[1] = 9999;}
  else
  {
        if (tt <= 0.0000001 ) {rst[0] = pow(10,10); rst[1] = 9999;}
        else{
  alp1 = alp-1;
  iter = 0;
  diff = 1.0;
  s = 1/pow(tt, 0.3) - 1; // starting value
  while(fabs(diff) > tol && iter < maxiter)
  {
    g = _psi(s,alp,galp) - tt;
    sroot = sqrt(s);
    tem = gsl_sf_bessel_Knu(fabs(alp1), 2*sroot);
    gp = -2 * pow(s, (alp1/2)) * tem / galp;
    diff = g / gp;
    s = s - diff;
    while(s <= 0) { diff = diff/2; s = s+diff; }
    iter = iter+1;
  }
  rst[0] = s;
  rst[1] = alp;
  rst[2] = g;
  rst[3] = gp;
            }
  }          
}

////////////////////// Functions for the Gumbel copula  /////////////////

void invpsi_Gumbel(double *tt, double *alp, int *N, double *rst)
{
  int i;
  for(i=0; i < N[0]; i++)
  {
    rst[i] = pow( - log(tt[i]), alp[i] );
  }  
}

void logden_Gumbel(double *F1, double *F2, double *alp, double *lf1, double *lf2, int *N, double *rst)
{
  double s1, s2, s, alp1;
  int i;
  for(i=0; i < N[0]; i++)
  {
  alp1 = 1 / alp[i];
  s1 = pow( - log(F1[i]), alp[i] );
  s2 = pow( - log(F2[i]), alp[i] );
  s = s1 + s2;
  rst[i] = - pow(s, alp1) + log( pow(alp1,2)*pow(s, (2*alp1-2)) - (alp1*alp1 - alp1)*pow(s, (alp1 - 2)) ) + lf1[i] + lf2[i] - ( log(F1[i]) - log(alp[i]) + (1 - alp[i])*log(-log(F1[i])) ) - ( log(F2[i]) - log(alp[i]) + (1 - alp[i])*log(-log(F2[i])) );
  }
}




