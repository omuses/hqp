#ifndef _JACUTILS_H_
#define _JACUTILS_H_
/*
  ----------------------------------------------------------------- 
  File jacutils.h of ADOL-C version 1.8.2           as of Mar/09/99
  -----------------------------------------------------------------
  This file containts declarations of jacobian utilities. 

  Last changes:
    990308 christo: myalloc1_ushort -> myalloc1_uint
    990308 christo: block_pattern : number of blocks : 
                    unsigned short -> unsigned int
    990308 christo: bit patterns : unsigned int -> unsigned long int
    981203 olvo:    untransposing reverse
    981130 olvo:    includes changed
    981126 olvo:    last check (p's & q's)
    981125 christo: changed block_pattern() parameter list 
    981118 christo: changed block_pattern() parameter list
    981101 christo: changed block_pattern() parameter values
  
  -----------------------------------------------------------------
*/

/****************************************************************************/
/*                                                                   MAKROS */

/* Max. number of unsigned ints to store the seed / jacobian matrix strips.  
     Reduce this value to x if your system happens to run out of memory. 
     x < 10 makes no sense. x = 50 or 100 is better
     x stays for ( x * sizeof(unsigned long int) * 8 ) 
       (block) variables at once                                            */
#define PQ_STRIPMINE_MAX 30


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                      PUBLIC EXPORTS ONLY */

#ifdef __cplusplus
/****************************************************************************/
/****************************************************************************/
/*                                                            No C++ THINGS */


/****************************************************************************/
/****************************************************************************/
/*                                                         Now the C THINGS */
extern "C" {
#endif

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
#ifdef __cplusplus
}
#endif

#endif
