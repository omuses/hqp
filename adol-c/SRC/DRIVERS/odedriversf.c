#define _ODEDRIVERSF_C_
#define _ADOLC_SRC_
/*
   ------------------------------------------------------------------------
   File odedriversf.c of ADOL-C version 1.8.0               as of Nov/30/98
   ------------------------------------------------------------------------
   Easy to use drivers for ordinary differential equations (ODE)
   (Implementation of the Fortran callable C interfaces).

   Last changes:
      981130  olvo:   newly created from driversc.c

  ------------------------------------------------------------------------
*/


/****************************************************************************/
/*                                                                 INCLUDES */
#include "dvlparms.h" /* Developers Parameters */
#include "usrparms.h" /* Users Parameters      */
#include "odedrivers.h"
#include "interfaces.h"
#include "adalloc.h"
#include "fortutils.h"

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
#include <malloc.h>
extern "C" {
#endif

/****************************************************************************/
/*                                                         DRIVERS FOR ODEs */

/*--------------------------------------------------------------------------*/
/*                                                                  forodec */
/* forodec(tag, n, tau, dold, dnew, X[n][d+1])                              */

fint forodec_(fint* ftag,    /* tape identifier */
             fint* fn,       /* space dimension */
             fdouble* ftau,  /* scaling defaults to 1.0 */
             fint* fdol,     /* previous degree defaults to zero */
             fint* fdeg,     /* New degree of consistency        */
	     fdouble* fy)    /* Taylor series                    */
{
  int rc= -1;
  int tag=*ftag, n=*fn, dol=*fdol, deg=*fdeg;
  int i;
  double tau=*ftau;
  double** Y = myalloc2(n,deg+1);
  for(i=0;i<n;i++)
    *Y[i] = fy[i];
  rc= forodec(tag,n,tau,dol,deg,Y);
  pack2(n,deg+1,Y,fy);
  free((char*)*Y); free((char*)Y);
  return rc;
}


/*--------------------------------------------------------------------------*/
/*                                                                  accodec */
/* accodec(n, tau, d, Z[n][n][d+1], B[n][n][d+1], nz[n][n])                 */

fint accodec_(fint* fn,             /* space dimension */
	      fdouble* ftau,        /* scaling defaults to 1.0 */
              fint* fdeg,           /* highest degree          */ 
	      fdouble* fa,          /* input tensor of "partial" Jacobians */
              fdouble* fb)          /* output tensor of "total" Jacobians  */
{ 
  int rc= 1;
  int n=*fn, deg=*fdeg;
  double tau=*ftau;
  double*** A = myalloc3(n,n,deg); 
  double*** B = myalloc3(n,n,deg); 
  spread3(n,n,deg,fa,A);
  accodec(n,tau,deg,A,B,0);
  pack3(n,n,deg,B,fb);
  free((char*)**A); free((char*)*A); free((char*)A);
  free((char*)**B); free((char*)*B); free((char*)B);
  return rc;
}


/****************************************************************************/
/*                                                               THAT'S ALL */

#ifdef __cplusplus
}
#endif

#undef _ADOLC_SRC_
#undef _ODEDRIVERSF_C_
