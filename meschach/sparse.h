
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


/*
	Header for sparse matrix stuff.
	Basic sparse routines to be held in sparse.c
*/

/* RCS id: $Id: sparse.h,v 1.1 2001/03/01 17:19:00 rfranke Exp $ */

#ifndef SPARSEH

#define SPARSEH 


#include        "matrix.h"


/* basic sparse types */

typedef struct row_elt	{
	Real	val;
	int	col;
#ifdef SPARSE_COL_ACCESS
	int nxt_row, nxt_idx;
#endif
		} row_elt;

typedef struct SPROW {
	int	len, maxlen, diag;
	row_elt	*elt;		/* elt[maxlen] */
		} SPROW;

typedef struct SPMAT {
	int	m, n, max_m, max_n;
	char	flag_col, flag_diag;
	SPROW	*row;		/* row[max_m] */
	int	*start_row;	/* start_row[max_n] */
	int	*start_idx;	/* start_idx[max_n] */
	      } SPMAT;

/* Note that the first allocated entry in column j is start_row[j];
	This starts the chain down the columns using the nxt_row and nxt_idx
	fields of each entry in each row. */

typedef struct SPPAIR { int pos; Real val; } SPPAIR;

typedef struct SPVEC {
	int	dim, max_dim;
	SPPAIR	*elt;		/* elt[max_dim] */
	       } SPVEC;

#define	SMNULL	((SPMAT*)NULL)
#define	SVNULL	((SPVEC*)NULL)

/* Macro for speedup */
#define	sprow_idx2(r,c,hint)	\
	( ( (hint) >= 0 && (hint) < (r)->len && \
	   (r)->elt[hint].col == (c)) ? (hint) : sprow_idx((r),(c)) )

MESCH__BEGIN_DECLS

/* memory functions */

int sp_get_vars MESCH__P((int m, int n, int deg, ...));
int sp_resize_vars MESCH__P((int m, int n, ...));
int sp_free_vars MESCH__P((SPMAT **, ...));

/* Sparse Matrix Operations and Utilities */

SPMAT	*sp_get MESCH__P((int, int, int));
SPMAT	*sp_copy MESCH__P((const SPMAT *));
SPMAT	*sp_copy2 MESCH__P((const SPMAT *, SPMAT *));
SPMAT	*sp_zero MESCH__P((SPMAT *));
SPMAT	*sp_resize MESCH__P((SPMAT *,int, int));
SPMAT	*sp_compact MESCH__P((SPMAT *, double));
double	sp_get_val MESCH__P((const SPMAT *, int, int));
double	sp_set_val MESCH__P((SPMAT *, int, int, double));
VEC	*sp_mv_mlt MESCH__P((const SPMAT *, const VEC *, VEC *));
VEC	*sp_vm_mlt MESCH__P((const SPMAT *, const VEC *,VEC *));
int	sp_free MESCH__P((SPMAT *));

/* Access path operations */
SPMAT	*sp_col_access MESCH__P((SPMAT *));
SPMAT	*sp_diag_access MESCH__P((SPMAT *));
int     chk_col_access MESCH__P((const SPMAT *));

/* Input/output operations */
SPMAT	*sp_finput MESCH__P((FILE *));
void	sp_foutput MESCH__P((FILE *, const SPMAT *));
void	sp_foutput2 MESCH__P((FILE * , const SPMAT *));

/* algebraic operations */
SPMAT 	*sp_smlt MESCH__P((const SPMAT *A, double alpha, SPMAT *B));
SPMAT	*sp_add MESCH__P((const SPMAT *A, const SPMAT *B, SPMAT *C));
SPMAT	*sp_sub MESCH__P((const SPMAT *A, const SPMAT *B, SPMAT *C));
SPMAT	*sp_mltadd MESCH__P((const SPMAT *A, const SPMAT *B, double alpha,
			SPMAT *C));

/* sparse row operations */
SPROW	*sprow_get MESCH__P((int)), *sprow_xpd MESCH__P((SPROW *r, int n, int type));
SPROW	*sprow_resize MESCH__P((SPROW *r, int n, int type));
SPROW	*sprow_merge MESCH__P((const SPROW *, const SPROW *, SPROW *, int type));
SPROW	*sprow_copy MESCH__P((const SPROW *, const SPROW *, SPROW *, int type));
SPROW	*sprow_mltadd MESCH__P((const SPROW *, const SPROW *, double, int,
			   SPROW *, int type));
SPROW 	*sprow_add MESCH__P((const SPROW *r1, const SPROW *r2, int j0,
			SPROW *r_out, int type));
SPROW 	*sprow_sub MESCH__P((const SPROW *r1, const SPROW *r2, int j0,
			SPROW *r_out, int type)); 
SPROW 	*sprow_smlt MESCH__P((const SPROW *r1, double alpha, int j0,
			 SPROW *r_out, int type));
double	sprow_set_val MESCH__P((SPROW *, int, double));
int     sprow_free MESCH__P((SPROW *));
int	sprow_idx MESCH__P((const SPROW *, int));
void	sprow_foutput MESCH__P((FILE *, const SPROW *));

/* dump */
void    sp_dump MESCH__P((FILE *fp, const SPMAT *A));
void    sprow_dump MESCH__P((FILE *fp, SPROW *r));
MAT	*sp_m2dense MESCH__P((const SPMAT *A, MAT *out));


/* MACROS */

#define	sp_input()	sp_finput(stdin)
#define	sp_output(A)	sp_foutput(stdout,(A))
#define	sp_output2(A)	sp_foutput2(stdout,(A))
#define	row_mltadd(r1,r2,alpha,out)	sprow_mltadd(r1,r2,alpha,0,out)
#define	out_row(r)	sprow_foutput(stdout,(r))

#define SP_FREE(A)    ( sp_free((A)),  (A)=(SPMAT *)NULL) 

/* utility for index computations -- ensures index returned >= 0 */
#define	fixindex(idx)	((idx) == -1 ? (error(E_BOUNDS,"fixindex"),0) : \
			 (idx) < 0 ? -((idx)+2) : (idx))


/*  NOT USED */

/* loop over the columns in a row */
/*
#define	loop_cols(r,e,code) \
    do { int _r_idx; row_elt *e; SPROW *_t_row;			\
	  _t_row = (r); e = &(_t_row->elt);				\
	  for ( _r_idx = 0; _r_idx < _t_row->len; _r_idx++, e++ )	\
	  {  code;  }  }  while ( 0 )
*/
/* loop over the rows in a column */
/*
#define	loop_cols(A,col,e,code) \
    do { int _r_num, _r_idx, _c; SPROW *_r; row_elt *e;		\
	  if ( ! (A)->flag_col )	sp_col_access((A));		\
	  col_num = (col);						\
	  if ( col_num < 0 || col_num >= A->n )				\
	      error(E_BOUNDS,"loop_cols");				\
          _r_num = (A)->start_row[_c]; _r_idx = (A)->start_idx[_c];	\
	  while ( _r_num >= 0 )  {					\
	      _r = &((A)->row[_r_num]);					\
              _r_idx = sprow_idx2(_r,_c,_r_idx);			\
              if ( _r_idx < 0 )  continue;				\
	      e = &(_r->elt[_r_idx]);	code;				\
	      _r_num = e->nxt_row;	_r_idx = e->nxt_idx;		\
	      } } while ( 0 )

*/

MESCH__END_DECLS

#endif

