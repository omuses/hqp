/*
	Header for sparse matrix stuff.
	Basic sparse routines to be held in sparse.c
*/

/* RCS id: $Header: /home/ruediger/hqp/hqp/meschach/sparsdef.h,v 1.1 2001/03/01 17:18:58 rfranke Exp $ */

#ifndef SPARSEDEFH

#define SPARSEDEFH 1

/*************************************************************************
* #ifndef MATRIXH
* #include	"matrix.h"
* #endif
*************************************************************************/

typedef struct row_elt	{
	int	col, nxt_row, nxt_idx;
	double	val;
		} row_elt;

typedef struct sp_row {
	int	len, maxlen, diag;
	row_elt	*elt;		/* elt[maxlen] */
		} sp_row;

typedef struct sp_mat {
	int	m, n, max_m, max_n;
	char	flag_col, flag_diag;
	sp_row	*row;		/* row[max_m] */
	int	*start_row;	/* start_row[max_n] */
	int	*start_idx;	/* start_idx[max_n] */
	      } sp_mat;

/* Note that the first allocated entry in column j is start_row[j];
	This starts the chain down the columns using the nxt_row and nxt_idx
	fields of each entry in each row. */

typedef struct SPPAIR { int pos; double val; } SPPAIR;

typedef struct sp_vec {
	int	dim, max_dim;
	SPPAIR	*elt;		/* elt[max_dim] */
	       } sp_vec;

#define	SMNULL	((sp_mat*)NULL)
#define	SVNULL	((sp_vec*)NULL)

/* Macro for speedup */
#define	sp_get_idx2(r,c,hint)	\
	( ( (hint) >= 0 && (hint) < (r)->len && \
	   (r)->elt[hint].col == (c)) ? (hint) : sp_get_idx((r),(c)) )

#endif
