#define _ADALLOCC_C_
#define _ADOLC_SRC_
/*
   --------------------------------------------------------------
   File adallocc.c of ADOL-C version 1.8.6w     as of Feb/14/2000
   --------------------------------------------------------------
   C allocation of arrays of doubles in several dimensions 

   Last changes:
      20000217 olvo: Version Waechter
      20000214 olvo: The defininition of the macro USE_CALLOC
                     forces the ADOL-C allocation routines
                     to use 'calloc' instead of 'malloc'.
                     This may help in case of problems with
                     uninitialized memory as reported by Andreas
                     Waechter from CMU.
        990622 olvo: special identity allocations (2n-1-vectors)
                     routines for freeing memory
        981130 olvo: newly created.
   --------------------------------------------------------------
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "dvlparms.h" /* Developers Parameters */
#include "usrparms.h" /* Users Parameters      */
#include "adalloc.h"

#include <stdio.h>

#ifdef __cplusplus
#include <malloc.h>
extern "C" {
#endif

/****************************************************************************/
/*                                USE calloc INSTEAD OF malloc IF NECESSARY */
/* olvo 20000214: introduced after reported problems with uninitialized
                  memory calling ADOL-C drivers from FORTRAN on SGI */
#ifdef USE_CALLOC
#define ADOLC_MALLOC(n,m) calloc(n,m)
#else
#define ADOLC_MALLOC(n,m) malloc(n*m)
#endif


/****************************************************************************/
/*                                              MEMORY MANAGEMENT UTILITIES */

/*--------------------------------------------------------------------------*/
double* myalloc1(int m)
{ double* A = (double*)ADOLC_MALLOC(m,sizeof(double));
  if (A == NULL)
  { fprintf(DIAG_OUT,"ADOL-C error: myalloc1 cannot allocate %i bytes\n",
                     m*sizeof(double));
    exit (-1);
  }
  return A;
}

/*--------------------------------------------------------------------------*/
double** myalloc2(int m, int n)
{ double *Adum = (double*)ADOLC_MALLOC(m*n,sizeof(double));
  double   **A = (double**)ADOLC_MALLOC(m,sizeof(double*));
  int i;
  if (Adum == NULL)
  { fprintf(DIAG_OUT,"ADOL-C error: myalloc2 cannot allocate %i bytes\n",
                     m*n*sizeof(double));
    exit (-1);
  }
  if (A == NULL)
  { fprintf(DIAG_OUT,"ADOL-C error: myalloc2 cannot allocate %i bytes\n",
                     m*sizeof(double*));
    exit (-1);
  }
  for (i=0; i<m; i++)
  { A[i] = Adum;
    Adum += n;
  }
  return A;
}

/*--------------------------------------------------------------------------*/
double*** myalloc3(int m, int n, int p)
{ /* This function allocates 3-tensors contiguously */ 
  double *Adum = (double*) ADOLC_MALLOC(m*n*p,sizeof(double));
  double **Apt = (double**)malloc(m*n*sizeof(double*));
  double  ***A = (double***)malloc(m*sizeof(double**));
  int i,j;
  if (Adum == NULL)
  { fprintf(DIAG_OUT,"ADOL-C error: myalloc3 cannot allocate %i bytes\n",
                     m*n*p*sizeof(double));
    exit (-1);
  }
  if (Apt == NULL)
  { fprintf(DIAG_OUT,"ADOL-C error: myalloc3 cannot allocate %i bytes\n",
                     m*n*sizeof(double*));
    exit (-1);
  }
  if (A == NULL)
  { fprintf(DIAG_OUT,"ADOL-C error: myalloc3 cannot allocate %i bytes\n",
                     m*sizeof(double**));
    exit (-1);
  }
  for (i=0; i<m; i++)
  { A[i] = Apt;
    for (j=0; j<n; j++)
    { *Apt++ =  Adum;
      Adum += p;
    }
  }   
  return A;
}

/*--------------------------------------------------------------------------*/
void myfree1(double   *A)
{ free((char*) A);
}

/*--------------------------------------------------------------------------*/
void myfree2(double  **A)
{ free((char*)*A); free((char*) A);
}

/*--------------------------------------------------------------------------*/
void myfree3(double ***A)
{ free((char*)**A); free((char*)*A); free((char*) A);
}


/****************************************************************************/
/*                                          SPECIAL IDENTITY REPRESENTATION */

/*--------------------------------------------------------------------------*/
double   **myallocI2(int n)
{ double *Idum = (double*)ADOLC_MALLOC((2*n-1),sizeof(double));
  double   **I = (double**)malloc(n*sizeof(double*));
  int i;
  if (Idum == NULL)
  { fprintf(DIAG_OUT,"ADOL-C error: myallocI2 cannot allocate %i bytes\n",
                     (2*n-1)*sizeof(double));
    exit (-1);
  }
  if (I == NULL)
  { fprintf(DIAG_OUT,"ADOL-C error: myallocI2 cannot allocate %i bytes\n",
                     n*sizeof(double*));
    exit (-1);
  }
  I[0] = Idum+=(n-1); 
  *Idum = 1.0;
  for (i=1; i<n; i++)
  { I[i] = --Idum;
    *Idum = 0.0;
  }
  return I;
}

/*--------------------------------------------------------------------------*/
void myfreeI2(int n, double** I)
{ free((char*)I[n-1]); free((char*) I);
}


/****************************************************************************/
/*                                                                THAT'S ALL*/
#ifdef __cplusplus
}
#endif

#undef _ADOLC_SRC_
#undef _ADALLOCC_C_



