/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     convolute.c
 Revision: $Id: convolut.c,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: Convolution routines (used by ho_rev.mc)

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20040423 kowarz: adapted to configure - make - make install
          19981130 olvo:   last check 
          19980707 olvo:   changed index range
          19980616 olvo:   (1) void copyAndZeroset(..)
                               void inconv0(..)
                               void deconv0(..)
                           (2) code optimization

----------------------------------------------------------------------------*/

#include "convolut.h"

BEGIN_C_DECLS

/****************************************************************************/
/*                                                              CONVOLUTION */

/*--------------------------------------------------------------------------*/
/* Evaluates convolution of a and b to c */
void conv( int dim, double *a, revreal *b, double *c )
{ double tmpVal;
  int i,j;
  for (i=dim-1; i>=0; i--)
  { tmpVal = a[i]*b[0];
    for (j=1; j<=i; j++)
      tmpVal += a[i-j]*b[j];
    c[i] = tmpVal;
  }
}

/****************************************************************************/
/*                                                  INCREMENTAL CONVOLUTION */

/*--------------------------------------------------------------------------*/
/* Increments truncated convolution of a and b to c */
void inconv( int dim, double *a, revreal *b, double *c )
{ double tmpVal;
  int i,j;
  for (i=dim-1; i>=0; i--)
  { tmpVal = a[i]*b[0]; 
    for (j=1; j<=i; j++)
      tmpVal += a[i-j]*b[j];
    c[i] += tmpVal;
  }
}

/*--------------------------------------------------------------------------*/
/* olvo 980616 nf */
/* Increments truncated convolution of a and b to c and sets a to zero */
void inconv0( int dim, double *a, revreal *b, double *c )
{ double tmpVal;
  int i,j;
  for (i=dim-1; i>=0; i--)
  { tmpVal = a[i]*b[0];
    a[i] = 0;
    for (j=1; j<=i; j++)
      tmpVal += a[i-j]*b[j];
    c[i] += tmpVal;
  }
}

/****************************************************************************/
/*                                                  DECREMENTAL CONVOLUTION */

/*--------------------------------------------------------------------------*/
/* Decrements truncated convolution of a and b to c */
void deconv( int dim, double *a, double *b, double *c )
{ double tmpVal;
  int i,j;
  for (i=dim-1; i>=0; i--)
  { tmpVal = a[i]*b[0];
    for (j=1; j<=i; j++)
      tmpVal += a[i-j]*b[j];
    c[i] -= tmpVal;
  }
}

/*--------------------------------------------------------------------------*/
/* olvo 980616 nf */
/* Decrements truncated convolution of a and b to c and sets a to zero */
void deconv0( int dim, double *a, revreal *b, double *c )
{ double tmpVal;
  int i,j;
  for (i=dim-1; i>=0; i--)
  { tmpVal = a[i]*b[0];
    a[i] = 0;
    for (j=1; j<=i; j++)
      tmpVal += a[i-j]*b[j];
    c[i] -= tmpVal;
  }
}

/****************************************************************************/
/*                                                    OTHER USEFUL ROUTINES */

/*--------------------------------------------------------------------------*/
void divide( int dim, revreal *a, revreal *b, revreal *c )
{ int i,j;
  double rec = 1/b[0];
  for (i=0; i<dim; i++)
  { c[i] = a[i];
    for (j=0; j<i; j++)
      c[i] -= c[j]*b[i-j];
    c[i] *= rec;
  }
}

/*--------------------------------------------------------------------------*/
void recipr( int dim, double a, revreal *b, revreal *c )
{ int i,j;
  double rec = 1/b[0];
  c[0] = a*rec;
  for (i=1; i<dim; i++)
  { c[i] = 0;
    for (j=0; j<i; j++)
      c[i] -= c[j]*b[i-j];
    c[i] *= rec;
  }
}

/****************************************************************************/
/*                                                                  ZEROING */

/*--------------------------------------------------------------------------*/
/* Set a to zero */
void zeroset(int dim, double *a) 
{ 
  int i;
  for(i=0;i<dim;i++) 
    a[i] = 0; 
}

/*--------------------------------------------------------------------------*/
/* olvo 980616 nf */
/* Copies a to tmp and initializes a to zero */
void copyAndZeroset( int dim, double *a, double* tmp )
{ int i;
  for (i=0; i<dim; i++)
  { tmp[i] = a[i];
    a[i] = 0.0;
  } 
}

/****************************************************************************/
END_C_DECLS
