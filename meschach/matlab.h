
/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/


/* matlab.h -- Header file for matlab.c and spmatlab.c for save/load formats */

#ifndef MATLAB_DEF

#define	MATLAB_DEF

/* structure required by MATLAB */
typedef struct {
	long    type;   /* matrix type */
	long    m;      /* # rows */
	long    n;      /* # cols */
	long    imag;   /* is complex? */
	long    namlen; /* length of variable name */
		} matlab;

/* macros for matrix storage type */
#define INTEL   0       /* for 80x87 format */
#define PC      INTEL
#define MOTOROLA        1       /* 6888x format */
#define SUN     MOTOROLA
#define APOLLO  MOTOROLA
#define MAC     MOTOROLA
#define VAX_D   2
#define VAX_G   3

#define COL_ORDER       0
#define ROW_ORDER       1

#define DOUBLE_PREC  0       /* double precision */
#define SINGLE_PREC  1       /* single precision */
#define INT_32  2       /* 32 bit integers (signed) */
#define INT_16  3       /* 16 bit integers (signed) */
#define INT_16u 4       /* 16 bit integers (unsigned) */
/* end of macros for matrix storage type */

#ifndef MACH_ID
#define MACH_ID         MOTOROLA
#endif

#define ORDER           ROW_ORDER

#if REAL == DOUBLE
#define PRECISION       DOUBLE_PREC
#elif REAL == FLOAT
#define PRECISION  	SINGLE_PREC
#endif


/* prototypes */

MESCH_API MAT *m_save(FILE *,MAT *,char *);
MESCH_API MAT *m_load(FILE *,char **);
MESCH_API VEC *v_save(FILE *,VEC *,char *);
MESCH_API double d_save(FILE *,double,char *);

/* complex variant */
#ifdef COMPLEX
#include "zmatrix.h"

MESCH_API ZMAT	*zm_save(FILE *fp,ZMAT *A,char *name);
MESCH_API ZVEC	*zv_save(FILE *fp,ZVEC *x,char *name);
MESCH_API complex	z_save(FILE *fp,complex z,char *name);
MESCH_API ZMAT	*zm_load(FILE *fp,char **name);

#endif

#endif
