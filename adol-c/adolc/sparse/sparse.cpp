/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     sparse/sparse.cpp
 Revision: $Id: sparse.cpp,v 1.2 2004/10/14 13:29:48 e_arnold Exp $
 Contents: All "Easy To Use" C++ interfaces of SPARSE package

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History: 20040414 kowarz:  adaption to configure - make - make install
          19990308 christo: bit patterns : unsigned int -> unsigned long int
          19990308 christo: mode : short -> char
          19981203 olvo: untransposing reverse
          19981201 olvo: newly created from interfaces.C

----------------------------------------------------------------------------*/


#include "../sparse/sparse.h"
#include "../sparse/jacutils.h"

#include <math.h>

#if defined(__cplusplus)

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

#endif
