#ifndef _ADALLOC_H_
#define _ADALLOC_H_
/*
   ------------------------------------------------------------- 
   File adalloc.h of ADOL-C version 1.8.7      as of Mar/10/2000
   -------------------------------------------------------------
   Allocation of arrays of doubles in several dimensions 

   Last changes:
      20000310 olvo: removed superflous semicola
        990622 olvo: myfree routines & special identity 
                     allocations (2n-1-vectors) 
                     (MOSTLY INLINED)
        981130 olvo: newly created.
   --------------------------------------------------------------
*/

#include <malloc.h>
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                      PUBLIC EXPORTS ONLY */


#ifdef __cplusplus
extern "C" {
#endif
/****************************************************************************/
/****************************************************************************/
/*                                                         Now the C THINGS */

/****************************************************************************/
/*                                              MEMORY MANAGEMENT UTILITIES */
double    *myalloc1(int);
double   **myalloc2(int, int);
double  ***myalloc3(int, int, int);

void myfree1(double   *);
void myfree2(double  **);
void myfree3(double ***);


/****************************************************************************/
/*                                          SPECIAL IDENTITY REPRESENTATION */
double   **myallocI2(int);
void myfreeI2(int, double**);


/****************************************************************************/
/****************************************************************************/
/*                                                       Now the C++ THINGS */
#ifdef __cplusplus
}

/****************************************************************************/
/*                                              MEMORY MANAGEMENT UTILITIES */
inline double   * myalloc(int n) { return myalloc1(n); }
inline double  ** myalloc(int m, int n) { return myalloc2(m,n); }
inline double *** myalloc(int m, int n, int p) { return myalloc3(m,n,p); }

inline void myfree(double   *A) { myfree1(A); }
inline void myfree(double  **A) { myfree2(A); }
inline void myfree(double ***A) { myfree3(A); }

/****************************************************************************/
/*                                                               THAT'S ALL */
#endif

#endif

