/*
   --------------------------------------------------------------
   File helm-vect-exam.C of ADOL-C vers. 1.8.0    as of Dec/01/98
   --------------------------------------------------------------
   based on helm-vect-exam.C version 1.7 

   Example: Helmholtz energy example 
            Computes gradient using AD driver reverse(..) 
            & using vector operations

   Last changes: 
     981201 olvo   new headers

   --------------------------------------------------------------
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "../../../adolc/adolc.h"

#include <math.h>


/****************************************************************************/
/*                                                    CONSTANTS & VARIABLES */
const double TE = 0.01; /* originally 0.0 */
const double R  = sqrt(2.0);
 
 
/****************************************************************************/
/*                                                         HELMHOLTZ ENERGY */
adouble energy(int n, const adoublev &x, const  adoublev &bv)
{ adouble he, xax, bx, tem;
  int i,j;
  xax = 0;
  bx  = 0;
  he  = 0;
  bx = bv*x;
  for (i=0; i<n; i++)
  { he += x[i]*log(x[i]);
    tem = (2.0/(1.0+i+i))*x[i];
    for (j=0; j<i; j++) 
      tem += (1.0/(1.0+i+j))*x[j];
    xax += x[i]*tem;
  }
  xax *= 0.5;
  he = 1.3625E-3*(he-TE*log(1.0-bx));
  he = he - log((1+bx*(1+R))/(1+bx*(1-R)))*xax/bx;
  return he;
}


/****************************************************************************/
/*                                                                     MAIN */
/* This program computes first order directional derivatives 
   for the helmholtz energy function. Uses vector operations */
int main() 
{ int nf, n, j, l;
  fprintf(stdout,"HELM-VECT-EXAM (ADOL-C Example)\n\n");
  fprintf(stdout," # of independents/10 =? \n ");
  scanf("%d",&nf);

/*--------------------------------------------------------------------------*/
  double result = 0.0;                                    /* Initilizations */
  n = 10 * nf; 
  double* grad = new double[n];

  adoublev bv(n);
  adoublev x(n);
  adouble he;

  double r = 1.0/n;
  for (j=0; j<n; j++) 
    bv[j] = 0.02*(1.0+fabs(sin(double(j))));
  
/*--------------------------------------------------------------------------*/
  int imd_rev = 1;                                     /* Tracing with keep */
  trace_on(1,imd_rev);
  for (j=0; j<n; j++) 
    x[j] <<= r*sqrt(1.0+j);
  he = energy(n,x,bv);
  he >>= result;
  trace_off();
  fprintf(stdout,"%14.6le -- energy\n",result);

/*--------------------------------------------------------------------------*/
  reverse(1,1,n,0,1.0,grad);             /* reverse computation of gradient */

/*--------------------------------------------------------------------------*/
  for (l=0; l<n; l++)                                            /* results */
    fprintf(stdout,"%3d: %14.6le,  \n",l,grad[l]);
  fprintf(stdout,"%14.6le -- energy\n",result);

  delete [] grad;

  return 1;
}


/****************************************************************************/
/*                                                               THAT'S ALL */

