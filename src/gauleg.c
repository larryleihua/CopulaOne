#include <stdio.h>
#include <math.h>
/* sample main program to compare with jacobi.f in Stroud and Secrest (1966) */
/* Gauss-Legendre quadrature */
// gcc -DMAIN -o gauleg gauleg.c -lm
// August 2012, changed to pointers for interface to R
#ifdef MAIN
#define NN 20
main(int argc, char *argv[])
{ int iq,nq;
  double wq[NN],xq[NN],x1,x2;
  void gauleg(double *x1,double *x2, int *n, double x[], double w[]);
  x1=0.; x2=1.;
  scanf("%d", &nq);
  while(nq>0)
  { gauleg(&x1,&x2,&nq,xq,wq);
    printf("\nnq=%d\n", nq);
    printf("iq    xq                wq\n");
    for(iq=1;iq<=nq;iq++)
      printf("%d %.15f %.15e\n", iq,xq[iq],wq[iq]);
    scanf("%d", &nq);
  }
}
#endif

#define EPS 3.0e-11
/* Inputs:
     x1 = lowerlimit, 
     x2 = upperlimit, 
     nq = #quadrature points, 
   Outputs:
     xq[] = quarature  nodes
     wq[] = quadrature weights 
*/
void gauleg(double *x1,double *x2, int *nq, double xq[], double wq[])
{ int m,j,i,n;
  double z1,z,xm,xl,pp,p3,p2,p1;
  n=*nq; m=(n+1)/2;
  xm=0.5*(*x2 + *x1);
  xl=0.5*(*x2 - *x1);
  for (i=1;i<=m;i++)
  { //z=cos(3.141592654*(i-0.25)/(n+0.5));
    z=cos(3.14159265358979323846*(i-0.25)/(n+0.5));
    do
    { p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++)
      { p3=p2;
        p2=p1;
        p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    }
    while (fabs(z-z1) > EPS);
    xq[i]=xm-xl*z;
    xq[n+1-i]=xm+xl*z;
    wq[i]=2.0*xl/((1.0-z*z)*pp*pp);
    wq[n+1-i]=wq[i];
  }
}
#undef EPS
