#define _FORTUTILSS_C_
#define _ADOLC_SRC_
/*
   ----------------------------------------------------------------
   File fortutils.c of ADOL-C version 1.8.0         as of Nov/30/98
   ----------------------------------------------------------------
   Internal tools to handle Fortran arrays

   Last changed : 
        981130 olvo   newly created from driversc.c

   ----------------------------------------------------------------
*/


/****************************************************************************/
/*                                                                 INCLUDES */
#include "dvlparms.h"
#include "usrparms.h"
#include "fortutils.h"

#ifdef __cplusplus
extern "C" {
#endif


/****************************************************************************/
/*                                              ROUTINES TO USE WITH ADOL-F */

/*--------------------------------------------------------------------------*/
void spread1(int m, fdouble* x, double* X)
{ int j;
  for (j=0; j<m; j++)
    X[j] = *x++;
}

/*--------------------------------------------------------------------------*/
void pack1(int m, double* X, fdouble* x)
{ int j;
  for (j=0; j<m; j++)
    *x++ = X[j];
}

/*--------------------------------------------------------------------------*/
void spread2(int m, int n, fdouble* x, double** X)
{ int i,j;
  for (j=0; j<n; j++)
    for (i=0; i<m; i++)
      X[i][j] = *x++;
}

/*--------------------------------------------------------------------------*/
void pack2(int m, int n, double** X, fdouble* x)
{ int i,j;
  for (j=0; j<n; j++)
    for (i=0; i<m; i++)
      *x++ = X[i][j];
}

/*--------------------------------------------------------------------------*/
void spread3(int m, int n, int p, fdouble* x, double*** X)
{ int i,j,k;
  for (k=0; k<p; k++)
    for (j=0; j<n; j++)
      for (i=0; i<m; i++)
	X[i][j][k] = *x++;
}

/*--------------------------------------------------------------------------*/
void pack3(int m, int n, int p, double*** X, fdouble* x)
{ int i,j,k;
  for (k=0; k<p; k++)
    for (j=0; j<n; j++)
      for (i=0; i<m; i++)
	*x++ = X[i][j][k];
}


/****************************************************************************/
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

#undef _ADOLC_SRC_
#undef _FORTUTILS_C_

