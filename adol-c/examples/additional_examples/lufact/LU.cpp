#define _LU_C
/*
   --------------------------------------------------------------
   File LU.C of ADOL-C version 1.8.5              as of Nov/22/99
   --------------------------------------------------------------

   Example:  'active' LU-decomposition and according solver

   Last changes:
     990922 olvo    first version 

   --------------------------------------------------------------
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "LU.h"              // check own interfaces


/****************************************************************************/
/*                                                          ADOUBLE ROUTINE */
/* Simple LU-factorization according to Crout's algorithm without pivoting */
void LUfact(int n, adouble **A)
{ int i, j, k;
  adouble dum;
  for (j=0; j<n; j++)
  { /* L-part */
    for (i=0; i<j; i++)
      for (k=0; k<i; k++) 
        A[i][j] -= A[i][k] * A[k][j];
    /* U-part */
    for (i=j; i<n; i++)
      for (k=0; k<j; k++)
        A[i][j] -= A[i][k] * A[k][j];
    if (A[j][j] != 0)   
    { dum = 1.0 / A[j][j];
      for (i=j+1; i<n; i++)
        A[i][j] *= dum; 
    }
    else
    { fprintf(stderr,"Error in LUfact(..): pivot is zero\n");
      exit(-99);
    } 
  }  
}


/****************************************************************************/
/*                                                          ADOUBLE ROUTINE */
/* Solution of A*x=b by forward and backward substitution */
void LUsolve(int n, adouble **A, adouble *bx)
{ int i, j;
  /* forward substitution */
  for (i=0; i<n; i++)
    for (j=0; j<i-1; j++)
      bx[i] -= A[i][j] * bx[j];
  /* backward substitution */
  for (i=n-1; i>=0; i--)
  { for (j=i+1; j<n; j++)
      bx[i] -= A[i][j] * bx[j];
    bx[i] /= A[i][i];
  }  
}


/****************************************************************************/
/*                                                              END OF FILE */
#undef _LU_C
