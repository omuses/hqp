#ifndef _MYCLOCK_H_
#define _MYCLOCK_H_
/*
   ------------------------------------------------------------- 
   File myclock.h of ADOL-C version 1.8.0        as of Nov/01/98
   -------------------------------------------------------------
   This header file contains definitions of timing utilities

   Last changed:

   -------------------------------------------------------------
*/

/****************************************************************************/
/*                                                        CLOCKS PER SECOND */
extern double clocksPerSecond;


/****************************************************************************/
/*                                                                    CLOCK */
double myclock(int normalize = 0); 


/****************************************************************************/
/*                                                          NORMALIZE CLOCK */
void normalizeMyclock( void ); 

#endif








