/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     sparse/jacutils.h
 Revision: $Id: jacutils.h,v 1.2 2004/10/14 13:29:48 e_arnold Exp $
 Contents: This file containts declarations of jacobian utilities.

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History: 20040414 kowarz:  adaption to configure - make - make install
          19990308 christo: myalloc1_ushort -> myalloc1_uint
          19990308 christo: block_pattern : number of blocks : 
                            unsigned short -> unsigned int
          19990308 christo: bit patterns : unsigned int -> unsigned long int
          19981203 olvo:    untransposing reverse
          19981130 olvo:    includes changed
          18981126 olvo:    last check (p's & q's)
          19981125 christo: changed block_pattern() parameter list 
          19981118 christo: changed block_pattern() parameter list
          19981101 christo: changed block_pattern() parameter values

----------------------------------------------------------------------------*/
#if !defined(ADOLC_SPARSE_JACUTILS_H)
#define ADOLC_SPARSE_JACUTILS_H 1

#include "../common.h"

/* Max. number of unsigned ints to store the seed / jacobian matrix strips.  
   Reduce this value to x if your system happens to run out of memory. 
   x < 10 makes no sense. x = 50 or 100 is better
   x stays for ( x * sizeof(unsigned long int) * 8 ) 
   (block) variables at once                                            */

#define PQ_STRIPMINE_MAX 30

BEGIN_C_DECLS

/****************************************************************************/
/*                                                    BIT PATTERN UTILITIES */

/*--------------------------------------------------------------------------*/
/* int_forward_tight(tag, m, n, p, x[n], X[n][p], y[m], Y[m][p])            */

int int_forward_tight(short, int, int, int, double*, unsigned long int**,
                                            double*, unsigned long int**);


/*--------------------------------------------------------------------------*/
/* int_forward_safe(tag, m, n, p, X[n][p], Y[m][p])                         */

int int_forward_safe(short, int, int, int, unsigned long int**, 
                                           unsigned long int**);


/*--------------------------------------------------------------------------*/
/* int_reverse_tight(tag, m, n, q, U[q][m], Z[q][n])                        */

int int_reverse_tight(short, int, int, int, unsigned long int**, 
                                            unsigned long int**);


/*--------------------------------------------------------------------------*/
/* int_reverse_safe(tag, m, n, q, U[q][m], Z[q][n])                         */

int int_reverse_safe(short, int, int, int, unsigned long int**, 
                                           unsigned long int**);


/****************************************************************************/
/*                                                   JACOBIAN BLOCK PATTERN */

int block_pattern(short, int, int, double*, unsigned int*, unsigned int*,
                                            unsigned int**, int*);


/****************************************************************************/
/*                                              MEMORY MANAGEMENT UTILITIES */

unsigned int * myalloc1_uint(int);

unsigned long int *  myalloc1_ulong(int);
unsigned long int ** myalloc2_ulong(int, int);


/****************************************************************************/
/*                                                               THAT'S ALL */

END_C_DECLS

#endif
