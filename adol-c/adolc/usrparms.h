/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     usrparms.h
 Revision: $Id: usrparms.h,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
 Contents: User parameters:
           These parameters might affect the performance of the ADOL-C
           system; they are intended to be tweeked by users and local
           maintainence personal.

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20030306 andrea: introduce EPS
          20030305 andrea: introduce TBUFNUM
          20030304 andrea: change default allocation from malloc to calloc
          20000214 olvo:   The defininition of the macro USE_CALLOC
                           forces the ADOL-C allocation routines
                           to use 'calloc' instead of 'malloc'.
                           This may help in case of problems with
                           uninitialized memory as reported by Andreas
                           Waechter from CMU.
          19981130 olvo:   Fortran types
          19981030 olvo:   bufsize --> BUFSIZE & TBUFSIZE
          19980723 olvo:   new: DIAG_OUT as standard output  
                           FNAME3   as vs output
                           
----------------------------------------------------------------------------*/

#if !defined(ADOLC_USRPARMS_H)
#define ADOLC_USRPARMS_H 1

/*--------------------------------------------------------------------------*/
/* Buffer size for tapes */
#define BUFSIZE    65536 /* 16384 or  524288  */

/*--------------------------------------------------------------------------*/
/* Buffer size for temporary Taylor store */
#define TBUFSIZE   65536 /* 16384 or  524288  */

/*--------------------------------------------------------------------------*/
/* Number of temporary Taylor stores*/
#define TBUFNUM    32

/*--------------------------------------------------------------------------*/
/* ADOL-C data types */
#define locint     unsigned int   
#define revreal    double

/*--------------------------------------------------------------------------*/
/* Data types used by Fortran callable versions of functions */
#define fint       long 
#define fdouble    double

/*--------------------------------------------------------------------------*/
/* Definionion of inf and NaN */ 
#define inf_num    1.0     /* don't undefine these;  on non-IEEE machines */
#define inf_den    0.0     /* change the values to get small fractions    */
#define non_num    0.0     /* (inf_num/inf_den) and (non_num/non_den)     */
#define non_den    0.0     /* respectively, see the documentation         */
#define ADOLC_EPS  10E-20  /* for test on zero                            */

/*--------------------------------------------------------------------------*/
/* Enable/disable asinh, acosh,atanh, erf */ 
#undef ATRIG_ERF
/* #define ATRIG_ERF 1 */     

/****************************************************************************/
/* Standard output used for diagnostics by ADOL-C,                          */
/* e.g. stdout or stderr or whatever file identifier                        */
#define DIAG_OUT stdout

/*--------------------------------------------------------------------------*/
/* Use 'calloc' instead of 'malloc' in ADOL-C allocation routines. If you   */
/* have any trouble with uninitialized memory, then define ADOLC_USE_CALLOC.*/
#define ADOLC_USE_CALLOC 1

/*--------------------------------------------------------------------------*/
#endif
