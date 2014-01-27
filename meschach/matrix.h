
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
		Type definitions for general purpose maths package
*/

#ifndef	MATRIXH

/* RCS id: $Id: matrix.h,v 1.5 2004/01/17 17:50:13 rfranke Exp $ */

#define	MATRIXH	

#include	"machine.h"
#include        "err.h"
#include 	"meminfo.h"

/* unsigned integer type */
#ifndef U_INT_DEF
typedef	unsigned int	u_int;
#define U_INT_DEF
#endif

/* vector definition */
typedef	struct  {
		u_int	dim, max_dim;
		Real	*ve;
		} VEC;

/* matrix definition */
typedef	struct  {
		u_int	m, n;
		u_int	max_m, max_n, max_size;
		Real	**me,*base;	/* base is base of alloc'd mem */
		} MAT;

/* band matrix definition */
typedef struct {
               MAT   *mat;       /* matrix */
               int   lb,ub;    /* lower and upper bandwidth */
               } BAND;


/* permutation definition */
typedef	struct  {
		u_int	size, max_size, *pe;
		} PERM;

/* integer vector definition */
typedef struct  {
		u_int	dim, max_dim;
		int	*ive;
	        } IVEC;

MESCH__BEGIN_DECLS

/* Meschach version routine */
MESCH_API void	m_version MESCH__P((void));

/* allocate one object of given type */
#define	NEW(type)	((type *)calloc((size_t)1,(size_t)sizeof(type)))

/* allocate num objects of given type */
#define	NEW_A(num,type)	((type *)calloc((size_t)((num)>1?(num):1),(size_t)sizeof(type)))

/* re-allocate arry to have num objects of the given type */
#define	RENEW(var,num,type) \
    ((var)=(type *)((var) ? \
		    realloc((char *)(var),(size_t)(((num)>1?(num):1)*sizeof(type))) : \
		    calloc((size_t)((num)>1?(num):1),(size_t)sizeof(type))))

#define	MEMCOPY(from,to,n_items,type) \
 MEM_COPY((char *)(from),(char *)(to),(unsigned)(n_items)*sizeof(type))

/* type independent min and max operations */
#ifndef max
#define	max(a,b)	((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define	min(a,b)	((a) > (b) ? (b) : (a))
#endif


#undef TRUE
#define	TRUE	1
#undef FALSE
#define	FALSE	0


/* for input routines */
#define MAXLINE 81


/* Dynamic memory allocation */

/* Should use M_FREE/V_FREE/PX_FREE in programs instead of m/v/px_free()
   as this is considerably safer -- also provides a simple type check ! */

/* get/resize vector to given dimension */
MESCH_API VEC 	*v_get MESCH__P((int)), *v_resize MESCH__P((VEC *, int));

/* get/resize matrix to be m x n */
MESCH_API MAT 	*m_get MESCH__P((int, int)), *m_resize MESCH__P((MAT *, int, int));

/* get/resize permutation to have the given size */
MESCH_API PERM 	*px_get MESCH__P((int)), *px_resize MESCH__P((PERM *, int));

/* get/resize an integer vector to given dimension */
MESCH_API IVEC 	*iv_get MESCH__P((int)), *iv_resize MESCH__P((IVEC *, int));

/* get/resize a band matrix to given dimension */
MESCH_API BAND 	*bd_get MESCH__P((int, int, int));
MESCH_API BAND	*bd_resize MESCH__P((BAND *, int, int, int));

/* free (de-allocate) (band) matrices, vectors, permutations and 
   integer vectors */
MESCH_API int	iv_free MESCH__P((IVEC *));
MESCH_API int	m_free MESCH__P((MAT *)), v_free MESCH__P((VEC *)), 
    px_free MESCH__P((PERM *));
MESCH_API int	bd_free MESCH__P((BAND *));


/* MACROS */

/* macros that also check types and sets pointers to NULL */
#define	M_FREE(mat)	( m_free(mat),	(mat)=(MAT *)NULL )
#define V_FREE(vec)	( v_free(vec),	(vec)=(VEC *)NULL )
#define	PX_FREE(px)	( px_free(px),	(px)=(PERM *)NULL )
#define	IV_FREE(iv)	( iv_free(iv),	(iv)=(IVEC *)NULL )

#define MAXDIM  	M_MAX_INT


/* Entry level access to data structures */
#ifdef DEBUG

/* returns x[i] */
#define	v_entry(x,i)	(((i) < 0 || (i) >= (x)->dim) ? \
			 m_error(E_BOUNDS,"v_entry"), 0.0 : (x)->ve[i] )

/* x[i] <- val */
#define	v_set_val(x,i,val) ((x)->ve[i] = ((i) < 0 || (i) >= (x)->dim) ? \
			    m_error(E_BOUNDS,"v_set_val"), 0.0 : (val))

/* x[i] <- x[i] + val */
#define	v_add_val(x,i,val) ((x)->ve[i] += ((i) < 0 || (i) >= (x)->dim) ? \
			    m_error(E_BOUNDS,"v_add_val"), 0.0 : (val))

/* x[i] <- x[i] - val */
#define	v_sub_val(x,i,val) ((x)->ve[i] -= ((i) < 0 || (i) >= (x)->dim) ? \
			    m_error(E_BOUNDS,"v_sub_val"), 0.0 : (val))

/* returns A[i][j] */
#define	m_entry(A,i,j)	(((i) < 0 || (i) >= (A)->m || \
			  (j) < 0 || (j) >= (A)->n) ? \
			 m_error(E_BOUNDS,"m_entry"), 0.0 : (A)->me[i][j] )

/* A[i][j] <- val */
#define	m_set_val(A,i,j,val) ((A)->me[i][j] = ((i) < 0 || (i) >= (A)->m || \
					       (j) < 0 || (j) >= (A)->n) ? \
			      m_error(E_BOUNDS,"m_set_val"), 0.0 : (val) )

/* A[i][j] <- A[i][j] + val */
#define	m_add_val(A,i,j,val) ((A)->me[i][j] += ((i) < 0 || (i) >= (A)->m || \
						(j) < 0 || (j) >= (A)->n) ? \
			      m_error(E_BOUNDS,"m_add_val"), 0.0 : (val) )

/* A[i][j] <- A[i][j] - val */
#define	m_sub_val(A,i,j,val) ((A)->me[i][j] -= ((i) < 0 || (i) >= (A)->m || \
						(j) < 0 || (j) >= (A)->n) ? \
			      m_error(E_BOUNDS,"m_sub_val"), 0.0 : (val) )
#else

/* returns x[i] */
#define	v_entry(x,i)		((x)->ve[i])

/* x[i] <- val */
#define	v_set_val(x,i,val)	((x)->ve[i]  = (val))

/* x[i] <- x[i] + val */
#define	v_add_val(x,i,val)	((x)->ve[i] += (val))

 /* x[i] <- x[i] - val */
#define	v_sub_val(x,i,val)	((x)->ve[i] -= (val))

/* returns A[i][j] */
#define	m_entry(A,i,j)		((A)->me[i][j])

/* A[i][j] <- val */
#define	m_set_val(A,i,j,val)	((A)->me[i][j]  = (val) )

/* A[i][j] <- A[i][j] + val */
#define	m_add_val(A,i,j,val)	((A)->me[i][j] += (val) )

/* A[i][j] <- A[i][j] - val */
#define	m_sub_val(A,i,j,val)	((A)->me[i][j] -= (val) )

#endif


/* I/O routines */

/* print x on file fp */
MESCH_API void	v_foutput MESCH__P((FILE *fp, const VEC *x));

/* print A on file fp */
MESCH_API void	m_foutput MESCH__P((FILE *fp, const MAT *A));

/* print px on file fp */
MESCH_API void	px_foutput MESCH__P((FILE *fp, const PERM *px));

/* print ix on file fp */
MESCH_API void	iv_foutput MESCH__P((FILE *fp, const IVEC *ix));

/* Note: if out is NULL, then returned object is newly allocated;
        Also: if out is not NULL, then that size is assumed */

/* read in vector from fp */
MESCH_API VEC	*v_finput MESCH__P((FILE *fp, VEC *out));

/* read in matrix from fp */
MESCH_API MAT 	*m_finput MESCH__P((FILE *fp, MAT *out));

/* read in permutation from fp */
MESCH_API PERM 	*px_finput MESCH__P((FILE *fp, PERM *out));

/* read in int vector from fp */
MESCH_API IVEC 	*iv_finput MESCH__P((FILE *fp, IVEC *out));


/* fy_or_n -- yes-or-no to question in string s
        -- question written to stderr, input from fp 
        -- if fp is NOT a tty then return y_n_dflt */
MESCH_API int 	fy_or_n MESCH__P((FILE *fp, const char *s));

/* yn_dflt -- sets the value of y_n_dflt to val */
MESCH_API int 	yn_dflt MESCH__P((int val));

/* fin_int -- return integer read from file/stream fp
        -- prompt s on stderr if fp is a tty
        -- check that x lies between low and high: re-prompt if
                fp is a tty, error exit otherwise
        -- ignore check if low > high           */
MESCH_API int 	fin_int MESCH__P((FILE *fp, const char *s, int low, int high));

/* fin_double -- return double read from file/stream fp
        -- prompt s on stderr if fp is a tty
        -- check that x lies between low and high: re-prompt if
                fp is a tty, error exit otherwise
        -- ignore check if low > high           */
MESCH_API double 	fin_double MESCH__P((FILE *fp, const char *s, 
					     double low, double high));

/* it skips white spaces and strings of the form #....\n
   Here .... is a comment string */
MESCH_API int 	skipjunk MESCH__P((FILE *fp));


/* MACROS */

/* macros to use stdout and stdin instead of explicit fp */
#define	v_output(vec)	v_foutput(stdout,vec)
#define	v_input(vec)	v_finput(stdin,vec)
#define	m_output(mat)	m_foutput(stdout,mat)
#define	m_input(mat)	m_finput(stdin,mat)
#define	px_output(px)	px_foutput(stdout,px)
#define	px_input(px)	px_finput(stdin,px)
#define	iv_output(iv)	iv_foutput(stdout,iv)
#define	iv_input(iv)	iv_finput(stdin,iv)

/* general purpose input routine; skips comments # ... \n */
#define	finput(fp,prompt,fmt,var) \
	( ( isatty(fileno(fp)) ? fprintf(stderr,prompt) : skipjunk(fp) ), \
							fscanf(fp,fmt,var) )
#define	input(prompt,fmt,var)	finput(stdin,prompt,fmt,var)
#define	fprompter(fp,prompt) \
	( isatty(fileno(fp)) ? fprintf(stderr,prompt) : skipjunk(fp) )
#define	prompter(prompt)	fprompter(stdin,prompt)
#define	y_or_n(s)	fy_or_n(stdin,s)
#define	in_int(s,lo,hi)	fin_int(stdin,s,lo,hi)
#define	in_double(s,lo,hi)	fin_double(stdin,s,lo,hi)

/* Copying routines */

/* copy in to out starting at out[i0][j0] */
MESCH_API MAT	*_m_copy MESCH__P((const MAT *in, MAT *out, u_int i0, 
				   u_int j0));
MESCH_API MAT	* m_move MESCH__P((const MAT *in, int, int, int, int, 
				   MAT *out, int, int));
MESCH_API MAT	*vm_move MESCH__P((const VEC *in, int, MAT *out, int, int, 
				   int, int));

/* copy in to out starting at out[i0] */
MESCH_API VEC	*_v_copy MESCH__P((const VEC *in, VEC *out, u_int i0));
MESCH_API VEC	* v_move MESCH__P((const VEC *in, int, int, VEC *out, int));
MESCH_API VEC	*mv_move MESCH__P((const MAT *in, int, int, int, int, 
				   VEC *out, int));
MESCH_API PERM	*px_copy MESCH__P((const PERM *in, PERM *out));
MESCH_API IVEC	*iv_copy MESCH__P((const IVEC *in, IVEC *out));
MESCH_API IVEC	*iv_move MESCH__P((const IVEC *in, int, int, IVEC *out, int));
MESCH_API BAND	*bd_copy MESCH__P((const BAND *in, BAND *out));

/* MACROS */
#define	m_copy(in,out)	_m_copy(in,out,0,0)
#define	v_copy(in,out)	_v_copy(in,out,0)


/* Initialisation routines -- to be zero, ones, random or identity */
MESCH_API VEC	*v_zero MESCH__P((VEC *)); 
MESCH_API VEC	*v_rand MESCH__P((VEC *)); 
MESCH_API VEC	*v_ones MESCH__P((VEC *));
MESCH_API MAT 	*m_zero MESCH__P((MAT *));
MESCH_API MAT	*m_ident MESCH__P((MAT *));
MESCH_API MAT	*m_rand MESCH__P((MAT *));
MESCH_API MAT	*m_ones MESCH__P((MAT *));
MESCH_API PERM 	*px_ident MESCH__P((PERM *));
MESCH_API IVEC	*iv_zero MESCH__P((IVEC *));

/* Basic vector operations */

MESCH_API VEC	*sv_mlt MESCH__P((Real, const VEC *, VEC *));
MESCH_API VEC	*mv_mlt MESCH__P((const MAT *, const VEC *, VEC *));
MESCH_API VEC	*vm_mlt MESCH__P((const MAT *, const VEC *,VEC *));
MESCH_API VEC	*v_add MESCH__P((const VEC *, const VEC *, VEC *));
MESCH_API VEC	*v_sub MESCH__P((const VEC *, const VEC *, VEC *));
MESCH_API VEC	*px_vec MESCH__P((const PERM *, const VEC *, VEC *));
MESCH_API VEC	*pxinv_vec MESCH__P((const PERM *,const VEC *,VEC *));
MESCH_API VEC	*v_mltadd MESCH__P((const VEC *, const VEC *, Real, VEC *));

#ifdef HAVE_PROTOTYPES_IN_STRUCT
MESCH_API VEC	*v_map MESCH__P((Real (*f) MESCH__P((Real)), 
				 const VEC *, VEC *));
MESCH_API VEC	*_v_map MESCH__P((Real (*f) MESCH__P((void *, Real)), 
				  void *, const VEC *, VEC *));
#else
MESCH_API VEC	*v_map MESCH__P((Real (*f)(), const VEC *, VEC *));
MESCH_API VEC	*_v_map MESCH__P((Real (*f)(), void *, const VEC *, VEC *));
#endif

MESCH_API VEC	*v_lincomb MESCH__P((int, const VEC *[], const Real *, VEC *));
MESCH_API VEC	*v_linlist MESCH__P((VEC *out, const VEC *v1, Real a1, ...));

MESCH_API Real	v_min MESCH__P((const VEC *, int *));
MESCH_API Real	v_max MESCH__P((const VEC *, int *));
MESCH_API Real	v_sum MESCH__P((const VEC *));

/* Hadamard product: out[i] <- x[i].y[i] */
MESCH_API VEC	*v_star MESCH__P((const VEC *, const VEC *, VEC *));
MESCH_API VEC	*v_slash MESCH__P((const VEC *, const VEC *, VEC *));

/* sorts x, and sets order so that sorted x[i] = x[order[i]] */ 
MESCH_API VEC	*v_sort MESCH__P((VEC *, PERM *));

/* returns inner product starting at component i0 */
MESCH_API Real	_in_prod MESCH__P((const VEC *x, const VEC *y, 
					   u_int i0));

/* returns sum_{i=0}^{len-1} x[i].y[i] */
MESCH_API Real	__ip__ MESCH__P((const Real *, const Real *, int));

/* see v_mltadd(), v_add(), v_sub() and v_zero() */
MESCH_API void	__mltadd__ MESCH__P((const Real *, Real *, Real, int));
MESCH_API void	__add__ MESCH__P((const Real *, const Real *, Real *, int));
MESCH_API void	__sub__ MESCH__P((const Real *, const Real *, Real *, int));
MESCH_API void	__smlt__ MESCH__P((const Real *, Real, Real *, int));
MESCH_API void	__zero__ MESCH__P((Real *,int));


/* MACRO */
/* usual way of computing the inner product */
#define	in_prod(a,b)	_in_prod(a,b,0)

/* Norms */

/* scaled vector norms -- scale == NULL implies unscaled */
MESCH_API Real	_v_norm1 MESCH__P((const VEC *x, const VEC *scale));
MESCH_API Real	_v_norm2 MESCH__P((const VEC *x, const VEC *scale));
MESCH_API Real	_v_norm_inf MESCH__P((const VEC *x, const VEC *scale));

/* unscaled matrix norms */
MESCH_API Real	m_norm1 MESCH__P((const MAT *A));
MESCH_API Real	m_norm_inf MESCH__P((const MAT *A));
MESCH_API Real	m_norm_frob MESCH__P((const MAT *A));

/* MACROS */
/* unscaled vector norms */
#define	v_norm1(x)	_v_norm1(x,VNULL)
#define	v_norm2(x)	_v_norm2(x,VNULL)
#define	v_norm_inf(x)	_v_norm_inf(x,VNULL)

/* Basic matrix operations */
MESCH_API MAT	*sm_mlt MESCH__P((Real s, const MAT *A, MAT *out));
MESCH_API MAT	*m_mlt MESCH__P((const MAT *A, const MAT *B, MAT *out));
MESCH_API MAT	*mmtr_mlt MESCH__P((const MAT *A, const MAT *B, MAT *out));
MESCH_API MAT	*mtrm_mlt MESCH__P((const MAT *A, const MAT *B, MAT *out));
MESCH_API MAT	*m_add MESCH__P((const MAT *A, const MAT *B, MAT *out));
MESCH_API MAT	*m_sub MESCH__P((const MAT *A, const MAT *B, MAT *out));
MESCH_API MAT	*sub_mat MESCH__P((const MAT *A, u_int, u_int, u_int, u_int, 
				   MAT *out));
MESCH_API MAT	*m_transp MESCH__P((const MAT *A, MAT *out));
MESCH_API MAT	*ms_mltadd MESCH__P((const MAT *A, const MAT *B, Real s, 
				     MAT *out));

MESCH_API BAND	*bd_transp MESCH__P((const BAND *in, BAND *out));

MESCH_API MAT	*px_rows MESCH__P((const PERM *px, const MAT *A, MAT *out));
MESCH_API MAT	*px_cols MESCH__P((const PERM *px, const MAT *A, MAT *out));
MESCH_API MAT	*swap_rows MESCH__P((MAT *, int, int, int, int));
MESCH_API MAT	*swap_cols MESCH__P((MAT *, int, int, int, int));
MESCH_API MAT	*_set_col MESCH__P((MAT *A, u_int i, const VEC *out, 
				    u_int j0));
MESCH_API MAT	*_set_row MESCH__P((MAT *A, u_int j, const VEC *out, 
				    u_int i0));

MESCH_API VEC	*get_row MESCH__P((const MAT *, u_int, VEC *));
MESCH_API VEC	*get_col MESCH__P((const MAT *, u_int, VEC *));
MESCH_API VEC	*sub_vec MESCH__P((const VEC *, int, int, VEC *));
MESCH_API VEC	*mv_mltadd MESCH__P((const VEC *x, const VEC *y, 
				     const MAT *A, Real s, VEC *out));
MESCH_API VEC	*vm_mltadd MESCH__P((const VEC *x, const VEC *y, 
				     const MAT *A, Real s, VEC *out));


/* MACROS */
/* row i of A <- vec */
#define	set_row(mat,row,vec)	_set_row(mat,row,vec,0) 
/* col j of A <- vec */
#define	set_col(mat,col,vec)	_set_col(mat,col,vec,0)


/* Basic permutation operations */
MESCH_API PERM	*px_mlt MESCH__P((const PERM *px1, const PERM *px2, 
				  PERM *out));
MESCH_API PERM	*px_inv MESCH__P((const PERM *px, PERM *out));
MESCH_API PERM	*px_transp MESCH__P((PERM *px, u_int i, u_int j));

/* returns sign(px) = +1 if px product of even # transpositions
                      -1 if ps product of odd  # transpositions */
MESCH_API int	px_sign MESCH__P((const PERM *));


/* Basic integer vector operations */
MESCH_API IVEC	*iv_add MESCH__P((const IVEC *ix, const IVEC *iy, IVEC *out));
MESCH_API IVEC	*iv_sub MESCH__P((const IVEC *ix, const IVEC *iy, IVEC *out));
MESCH_API IVEC	*iv_sort MESCH__P((IVEC *ix, PERM *order));


/* miscellaneous functions */

MESCH_API Real	square MESCH__P((Real x));
/*MESCH_API Real	cube MESCH__P((Real x)); */
MESCH_API Real	mrand MESCH__P((void));
MESCH_API void	smrand MESCH__P((int seed));
MESCH_API void	mrandlist MESCH__P((Real *x, int len));

MESCH_API void	m_dump MESCH__P((FILE *fp, const MAT *a));
MESCH_API void	px_dump MESCH__P((FILE *fp, const PERM *px));
MESCH_API void	v_dump MESCH__P((FILE *fp, const VEC *x));
MESCH_API void	iv_dump MESCH__P((FILE *fp, const IVEC *ix));

MESCH_API MAT	*band2mat MESCH__P((const BAND *bA, MAT *A));
MESCH_API BAND	*mat2band MESCH__P((const MAT *A, int lb, int ub, BAND *bA));


/* miscellaneous constants */
#define	VNULL	((VEC *)NULL)
#define	MNULL	((MAT *)NULL)
#define	PNULL	((PERM *)NULL)
#define	IVNULL	((IVEC *)NULL)
#define BDNULL  ((BAND *)NULL)


/* varying number of arguments */

#include <stdarg.h>

MESCH_API int 	v_get_vars MESCH__P((int dim, ...));
MESCH_API int 	iv_get_vars MESCH__P((int dim, ...));
MESCH_API int 	m_get_vars MESCH__P((int m,int n, ...));
MESCH_API int 	px_get_vars MESCH__P((int dim, ...));

MESCH_API int 	v_resize_vars MESCH__P((int new_dim, ...));
MESCH_API int 	iv_resize_vars MESCH__P((int new_dim, ...));
MESCH_API int 	m_resize_vars MESCH__P((int m,int n, ...));
MESCH_API int 	px_resize_vars MESCH__P((int new_dim, ...));

MESCH_API int 	v_free_vars MESCH__P((VEC **, ...));
MESCH_API int 	iv_free_vars MESCH__P((IVEC **, ...));
MESCH_API int 	px_free_vars MESCH__P((PERM **, ...));
MESCH_API int 	m_free_vars MESCH__P((MAT **, ...));

MESCH__END_DECLS

#endif
