/* machine.h.  Generated automatically by configure.  */
/* Any machine specific stuff goes here */
/* Add details necessary for your own installation here! */

/* RCS id: $Id: machine.h,v 1.4 2004/01/17 17:50:13 rfranke Exp $ */

/* This is for use with "configure" -- if you are not using configure
	then use machine.van for the "vanilla" version of machine.h */

/* Note special macros: ANSI_C (ANSI C syntax)
			SEGMENTED (segmented memory machine e.g. MS-DOS)
			MALLOCDECL (declared if malloc() etc have
					been declared) */

/* #undef const */

/* #undef MALLOCDECL */
#define NOT_SEGMENTED 1
#define HAVE_MEMORY_H 1
/* #undef HAVE_COMPLEX_H */
#define HAVE_MALLOC_H 1
/* #undef HAVE_SYS_CDEFS_H */
#define STDC_HEADERS 1
#define HAVE_BCOPY 1
#define HAVE_BZERO 1
#define CHAR0ISDBL0 1
/* #undef WORDS_BIGENDIAN */
/* #define U_INT_DEF 1 */
#define HAVE_SYS_TYPES_H 1
#define VARARGS 1
#define HAVE_PROTOTYPES 1
/* #undef HAVE_PROTOTYPES_IN_STRUCT */

/* for inclusion into C++ files */
#ifdef __cplusplus
#define ANSI_C 1
#ifndef HAVE_PROTOTYPES 
#define HAVE_PROTOTYPES 1
#endif
#ifndef HAVE_PROTOTYPES_IN_STRUCT
#define HAVE_PROTOTYPES_IN_STRUCT 1
#endif
#endif /* __cplusplus */

/* example usage: VEC *PROTO(v_get,(int dim)); */
#ifdef HAVE_PROTOTYPES
#define	PROTO(name,args)	name args
#else
#define PROTO(name,args)	name()
#endif /* HAVE_PROTOTYPES */
#ifdef HAVE_PROTOTYPES_IN_STRUCT
/* PROTO_() is to be used instead of PROTO() in struct's and typedef's */
#define	PROTO_(name,args)	name args
#else
#define PROTO_(name,args)	name()
#endif /* HAVE_PROTOTYPES_IN_STRUCT */

/* for basic or larger versions */
/* #define COMPLEX 1 */
#define SPARSE 1

/* for loop unrolling */
#define VUNROLL 1
#define MUNROLL 1

/* for segmented memory */
#ifndef NOT_SEGMENTED
#define	SEGMENTED
#endif

/* if the system has sys/types.h */
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

/* if the system has malloc.h */
#ifdef HAVE_MALLOC_H
#define	MALLOCDECL	1
#include	<malloc.h>
#endif

/* any compiler should have this header */
/* if not, change it */
#include        <stdio.h>


/* Check for ANSI C memmove and memset */
#ifdef STDC_HEADERS

/* standard copy & zero functions */
#define	MEM_COPY(from,to,size)	if ((size) > 0) memmove((to),(from),(size))
#define	MEM_ZERO(where,size)	if ((size) > 0) memset((where),'\0',(size))

#endif

/* standard headers */
#include	<stdlib.h>
#include	<stddef.h>
#include	<string.h>
#include	<float.h>
/*----------   E. Arnold   1999-04-29   ----------*/
#include	<math.h>
/*----------*/


/* if have bcopy & bzero and no alternatives yet known, use them */
#ifdef HAVE_BCOPY
#ifndef MEM_COPY
/* nonstandard copy function */
#define	MEM_COPY(from,to,size)	bcopy((char *)(from),(char *)(to),(int)(size))
#endif
#endif

#ifdef HAVE_BZERO
#ifndef MEM_ZERO
/* nonstandard zero function */
#define	MEM_ZERO(where,size)	bzero((char *)(where),(int)(size))
#endif
#endif

/* if the system has complex.h */
#ifdef HAVE_COMPLEX_H
#include	<complex.h>
#endif

/* floating point precision */

/* you can choose single, double or long double (if available) precision */

#define FLOAT 		1
#define DOUBLE 		2
#define LONG_DOUBLE 	3

/* #undef REAL_FLT */
/* #undef REAL_DBL */

/* if nothing is defined, choose double precision */
#ifndef REAL_DBL
#ifndef REAL_FLT
#define REAL_DBL 1
#endif
#endif

/* single precision */
#ifdef REAL_FLT
#define  Real float
#define  LongReal float
#define REAL FLOAT
#define LONGREAL FLOAT
#endif

/* double precision */
#ifdef REAL_DBL
#define Real double
#define LongReal double
#define REAL DOUBLE
#define LONGREAL DOUBLE
#endif


/* machine epsilon or unit roundoff error */
/* This is correct on most IEEE Real precision systems */
#ifdef DBL_EPSILON
#if REAL == DOUBLE
#define	MACHEPS	DBL_EPSILON
#elif REAL == FLOAT
#define	MACHEPS	FLT_EPSILON
#elif REAL == LONGDOUBLE
#define MACHEPS LDBL_EPSILON
#endif
#endif

#define F_MACHEPS 1.19209e-07
#define D_MACHEPS 2.22045e-16

#ifndef MACHEPS
#if REAL == DOUBLE
#define	MACHEPS	D_MACHEPS
#elif REAL == FLOAT  
#define MACHEPS F_MACHEPS
#elif REAL == LONGDOUBLE
#define MACHEPS D_MACHEPS
#endif
#endif

/* #undef M_MACHEPS */

/********************
#ifdef DBL_EPSILON
#define	MACHEPS	DBL_EPSILON
#endif
#ifdef M_MACHEPS
#ifndef MACHEPS
#define MACHEPS	M_MACHEPS
#endif
#endif
********************/

#define	M_MAX_INT 2147483647
#ifdef	M_MAX_INT
#ifndef MAX_RAND
#define	MAX_RAND ((double)(M_MAX_INT))
#endif
#endif

/* for non-ANSI systems */
#ifndef HUGE_VAL
#define HUGE_VAL HUGE
/*----------   E. Arnold   1999-04-29   ----------*/
/* #else */
/* #undef HUGE */
#endif
#ifndef HUGE
#define HUGE HUGE_VAL
/*----------*/
#endif


/* setup C environment */

#ifdef HAVE_SYS_CDEFS_H
#include <sys/cdefs.h>
#endif

#ifndef MESCH__P
#if defined(__STDC__) || defined(__cplusplus)
#define	MESCH__P(protos)	protos
#else
#define MESCH__P(protos)	()
#endif
#endif

#ifndef MESCH__BEGIN_DECLS
#ifdef __cplusplus
#define	MESCH__BEGIN_DECLS	extern "C" {
#define	MESCH__END_DECLS	};
#else
#define	MESCH__BEGIN_DECLS
#define	MESCH__END_DECLS
#endif
#endif

/* use _setjmp if defined */
#ifdef HAVE__SETJMP
#define setjmp(ENV)	_setjmp(ENV)
#endif

/** define SPARSE_COL_ACCESS for compilation of column accesses */
/* #define SPARSE_COL_ACCESS 1 */

/* modified files
 *  sparse.c:
 *    no function sp_col_access()
 */

/** define MESCH_API when compiling a Dynamic Link Library (DLL) */
#ifndef MESCH_API
#define MESCH_API 
#endif

