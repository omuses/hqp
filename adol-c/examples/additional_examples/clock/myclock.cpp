/*
  ------------------------------------------------------------------------
  File myclock.C of ADOL-C version 1.8.0                   as of Nov/01/98
  ------------------------------------------------------------------------
  This file contains timing utilities

  Last changes:
     980806 olvo some little changes
                 now C++

  ------------------------------------------------------------------------
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include <sys/types.h>
#include <sys/times.h>
#include "myclock.h"


/****************************************************************************/
/*                                                        CLOCKS PER SECOND */
double clocksPerSecond = 100.0;


/****************************************************************************/
/*                                                          CLOCK UTILITIES */
double myclock( int normalize ) 
{ struct tms t;
  if (normalize)
    normalizeMyclock();
  times(&t);
  return ((double)t.tms_utime)/clocksPerSecond;
}


/****************************************************************************/
/*                                                          NORMALIZE CLOCK */
void normalizeMyclock( void ) 
{ long int utime;
  struct tms t;
  times(&t);
  utime = t.tms_utime;
  do
    times(&t);
  while (utime == t.tms_utime);
}


