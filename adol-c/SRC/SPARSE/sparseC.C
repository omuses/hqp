#define _SPARSEC_CPP_
#define _ADOLC_SRC_
/*
   --------------------------------------------------------------
   File sparseC.C of ADOL-C version 1.8.2         as of Mar/09/99
   --------------------------------------------------------------
   All "Easy To Use" C++ interfaces of SPARSE package

   Last changes:
        990308 christo  bit patterns : unsigned int -> unsigned long int
        990308 christo  mode : short -> char
        981203 olvo: untransposing reverse
        981201 olvo: newly created from interfaces.C

   --------------------------------------------------------------
*/


/****************************************************************************/
/*                                                                 INCLUDES */
#include "dvlparms.h" /* Developers Parameters */
#include "usrparms.h" /* Users Parameters      */
#include "sparse.h"
#include "jacutils.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


/****************************************************************************/
/*                                    Bit pattern propagation; general call */
/*                                                                          */
int forward( short              tag,
             int                m,
             int                n,
             int                p,
             double             *x,
             unsigned long int  **X,
             double             *y,
             unsigned long int  **Y,
             char               mode)
/* forward(tag, m, n, p, x[n], X[n][p], y[m], Y[m][p], mode)                */
{
  int rc = -1;
  if (mode == 1) // tight version
    if (x != NULL) 
      rc = int_forward_tight(tag,m,n,p,x,X,y,Y);
    else
    { fprintf(DIAG_OUT,"ADOL-C error:  no basepoint for bit"
                         " pattern forward tight.\n");
        exit(-1);
    }
  else 
    if (mode == 0) // safe version
      rc = int_forward_safe(tag,m,n,p,X,Y);
    else
    { fprintf(DIAG_OUT,"ADOL-C error:  bad mode parameter to bit"
                       " pattern forward.\n");
      exit(-1);
    }
  return (rc);
}


/****************************************************************************/
/*                                    Bit pattern propagation; no basepoint */
/*                                                                          */
int forward( short              tag,
             int                m,
             int                n,
             int                p,
             unsigned long int  **X,
             unsigned long int  **Y,
             char               mode)
/* forward(tag, m, n, p, X[n][p], Y[m][p], mode)                            */
{ 
  if (mode != 0) // not safe
  { fprintf(DIAG_OUT,"ADOL-C error:  bad mode parameter to bit"
                     " pattern forward.\n");
    exit(-1);
  }
  return int_forward_safe(tag,m,n,p,X,Y);
}



/****************************************************************************/
/*                                                                          */
/*                                    Bit pattern propagation, general call */
/*                                                                          */
int reverse( short             tag, 
             int               m,
	     int               n,
             int               q,
	     unsigned long int **U,
	     unsigned long int **Z,
             char              mode)
/* reverse(tag, m, n, q, U[q][m], Z[q][n]) */
{ int rc=-1;

  /* ! use better the tight version, the safe version supports no subscripts*/

  if (mode == 0) // safe version
    rc = int_reverse_safe(tag,m,n,q,U,Z);
  else
    if (mode == 1)
      rc = int_reverse_tight(tag,m,n,q,U,Z);
    else
    { fprintf(DIAG_OUT,"ADOL-C error:  bad mode parameter"
                       " to bit pattern reverse.\n");
      exit(-1);
    }
  return rc;
}


/****************************************************************************/
/*                                                               THAT'S ALL */
#undef _ADOLC_SRC_
#undef _SPARSEC_CPP_
