#ifndef ACIG_H
#define ACIG_H

/*-----------------------------------------------------------*/
/* (A)rchimedean (C)opula based on LT of (I)nverse (G)amma   */
/* Lei Hua, Dec, 2016                                        */
/* Based on Example 4 of Hua and Joe, 2011, Tail order and intermediate tail dependence of multivariate copulas */
/*-----------------------------------------------------------*/

void psi1(double *s0, double *alp0, double *galp0, double *rst);
void psi2(double *s0, double *alp0, double *galp0, double *rst);
double _psi(double s, double alp, double galp);
void lnK(double *nu, double *x, int *N, double *rst);
void lnG(double *alp, int *N, double *rst);
void psi(double *ss, double *alp, int *N, double *rst);
void invpsi(double *tt, double *alp, int *N, double *rst);
void invpsi0(double *tt0, double *alp0, double *galp0, double *rst);
void invpsi_Gumbel(double *tt, double *alp, int *N, double *rst);
void logden_Gumbel(double *F1, double *F2, double *alp, double *lf1, double *lf2, int *N, double *rst);

#endif /* ACIG_H */
