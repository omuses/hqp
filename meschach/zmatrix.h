
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


/* Main include file for zmeschach library -- complex vectors and matrices */

#ifndef ZMATRIXH
#define ZMATRIXH

#include "matrix.h"

MESCH__BEGIN_DECLS

          /*  Type definitions for complex vectors and matrices  */


/* complex definition */
typedef struct  {
                Real re,im;
        } complex;

/* complex vector definition */
typedef struct  {
                u_int   dim, max_dim;
                complex  *ve;
                } ZVEC;

/* complex matrix definition */
typedef struct  {
                u_int   m, n;
                u_int   max_m, max_n, max_size;
                complex *base;          /* base is base of alloc'd mem */
                complex **me;
                } ZMAT;

#define ZVNULL  ((ZVEC *)NULL)
#define ZMNULL  ((ZMAT *)NULL)

#define	Z_CONJ		1
#define	Z_NOCONJ	0


/* memory functions */

MESCH_API int zv_get_vars(int dim,...);
MESCH_API int zm_get_vars(int m,int n,...);
MESCH_API int zv_resize_vars(int new_dim,...);
MESCH_API int zm_resize_vars(int m,int n,...);
MESCH_API int zv_free_vars(ZVEC **,...);
MESCH_API int zm_free_vars(ZMAT **,...);

MESCH_API ZMAT	*_zm_copy(ZMAT *in,ZMAT *out,u_int i0,u_int j0);
MESCH_API ZMAT	* zm_move(ZMAT *, int, int, int, int, ZMAT *, int, int);
MESCH_API ZMAT	*zvm_move(ZVEC *, int, ZMAT *, int, int, int, int);
MESCH_API ZVEC	*_zv_copy(ZVEC *in,ZVEC *out,u_int i0);
MESCH_API ZVEC	* zv_move(ZVEC *, int, int, ZVEC *, int);
MESCH_API ZVEC	*zmv_move(ZMAT *, int, int, int, int, ZVEC *, int);
MESCH_API complex	z_finput(FILE *fp);
MESCH_API ZMAT	*zm_finput(FILE *fp,ZMAT *a);
MESCH_API ZVEC     *zv_finput(FILE *fp,ZVEC *x);
MESCH_API ZMAT	*zm_add(ZMAT *mat1,ZMAT *mat2,ZMAT *out);
MESCH_API ZMAT	*zm_sub(ZMAT *mat1,ZMAT *mat2,ZMAT *out);
MESCH_API ZMAT	*zm_mlt(ZMAT *A,ZMAT *B,ZMAT *OUT);
MESCH_API ZMAT	*zmma_mlt(ZMAT *A,ZMAT *B,ZMAT *OUT);
MESCH_API ZMAT	*zmam_mlt(ZMAT *A,ZMAT *B,ZMAT *OUT);
MESCH_API ZVEC	*zmv_mlt(ZMAT *A,ZVEC *b,ZVEC *out);
MESCH_API ZMAT	*zsm_mlt(complex scalar,ZMAT *matrix,ZMAT *out);
MESCH_API ZVEC	*zvm_mlt(ZMAT *A,ZVEC *b,ZVEC *out);
MESCH_API ZMAT	*zm_adjoint(ZMAT *in,ZMAT *out);
MESCH_API ZMAT	*zswap_rows(ZMAT *A,int i,int j,int lo,int hi);
MESCH_API ZMAT	*zswap_cols(ZMAT *A,int i,int j,int lo,int hi);
MESCH_API ZMAT	*mz_mltadd(ZMAT *A1,ZMAT *A2,complex s,ZMAT *out);
MESCH_API ZVEC	*zmv_mltadd(ZVEC *v1,ZVEC *v2,ZMAT *A,complex alpha,ZVEC *out);
MESCH_API ZVEC	*zvm_mltadd(ZVEC *v1,ZVEC *v2,ZMAT *A,complex alpha,ZVEC *out);
MESCH_API ZVEC	*zv_zero(ZVEC *x);
MESCH_API ZMAT	*zm_zero(ZMAT *A);
MESCH_API ZMAT	*zm_get(int m,int n);
MESCH_API ZVEC	*zv_get(int dim);
MESCH_API ZMAT	*zm_resize(ZMAT *A,int new_m,int new_n);
MESCH_API complex	_zin_prod(ZVEC *x,ZVEC *y,u_int i0,u_int flag);
MESCH_API ZVEC	*zv_resize(ZVEC *x,int new_dim);
MESCH_API ZVEC	*zv_mlt(complex scalar,ZVEC *vector,ZVEC *out);
MESCH_API ZVEC	*zv_add(ZVEC *vec1,ZVEC *vec2,ZVEC *out);
MESCH_API ZVEC	*zv_mltadd(ZVEC *v1,ZVEC *v2,complex scale,ZVEC *out);
MESCH_API ZVEC	*zv_sub(ZVEC *vec1,ZVEC *vec2,ZVEC *out);
#ifdef PROTOTYPES_IN_STRUCT
MESCH_API ZVEC	*zv_map(complex (*f)(),ZVEC *x,ZVEC *out);
MESCH_API ZVEC	*_zv_map(complex (*f)(),void *params,ZVEC *x,ZVEC *out);
#else
MESCH_API ZVEC	*zv_map(complex (*f)(complex),ZVEC *x,ZVEC *out);
MESCH_API ZVEC	*_zv_map(complex (*f)(void *,complex),void *params,ZVEC *x,ZVEC *out);
#endif
MESCH_API ZVEC	*zv_lincomb(int n,ZVEC *v[],complex a[],ZVEC *out);
MESCH_API ZVEC	*zv_linlist(ZVEC *out,ZVEC *v1,complex a1,...);
MESCH_API ZVEC	*zv_star(ZVEC *x1, ZVEC *x2, ZVEC *out);
MESCH_API ZVEC	*zv_slash(ZVEC *x1, ZVEC *x2, ZVEC *out);
MESCH_API int	zm_free(ZMAT *mat);
MESCH_API int	zv_free(ZVEC *vec);

MESCH_API ZVEC	*zv_rand(ZVEC *x);
MESCH_API ZMAT	*zm_rand(ZMAT *A);

MESCH_API ZVEC	*zget_row(ZMAT *A, int i, ZVEC *out);
MESCH_API ZVEC	*zget_col(ZMAT *A, int j, ZVEC *out);
MESCH_API ZMAT	*zset_row(ZMAT *A, int i, ZVEC *in);
MESCH_API ZMAT	*zset_col(ZMAT *A, int j, ZVEC *in);

MESCH_API ZVEC	*px_zvec(PERM *pi, ZVEC *in, ZVEC *out);
MESCH_API ZVEC	*pxinv_zvec(PERM *pi, ZVEC *in, ZVEC *out);

MESCH_API void	__zconj__(complex zp[], int len);
MESCH_API complex	__zip__(complex zp1[],complex zp2[],int len,int flag);
MESCH_API void	__zmltadd__(complex zp1[],complex zp2[],
			    complex s,int len,int flag);
MESCH_API void	__zmlt__(complex zp[],complex s,complex out[],int len);
MESCH_API void	__zadd__(complex zp1[],complex zp2[],complex out[],int len);
MESCH_API void	__zsub__(complex zp1[],complex zp2[],complex out[],int len);
MESCH_API void	__zzero__(complex zp[],int len);
MESCH_API void	z_foutput(FILE *fp,complex z);
MESCH_API void     zm_foutput(FILE *fp,ZMAT *a);
MESCH_API void     zv_foutput(FILE *fp,ZVEC *x);
MESCH_API void     zm_dump(FILE *fp,ZMAT *a);
MESCH_API void     zv_dump(FILE *fp,ZVEC *x);

MESCH_API double	_zv_norm1(ZVEC *x, VEC *scale);
MESCH_API double	_zv_norm2(ZVEC *x, VEC *scale);
MESCH_API double	_zv_norm_inf(ZVEC *x, VEC *scale);
MESCH_API double	zm_norm1(ZMAT *A);
MESCH_API double	zm_norm_inf(ZMAT *A);
MESCH_API double	zm_norm_frob(ZMAT *A);

MESCH_API complex	zmake(double real, double imag);
MESCH_API double	zabs(complex z);
MESCH_API complex zadd(complex z1,complex z2);
MESCH_API complex zsub(complex z1,complex z2);
MESCH_API complex	zmlt(complex z1,complex z2);
MESCH_API complex	zinv(complex z);
MESCH_API complex	zdiv(complex z1,complex z2);
MESCH_API complex	zsqrt(complex z);
MESCH_API complex	zexp(complex z);
MESCH_API complex	zlog(complex z);
MESCH_API complex	zconj(complex z);
MESCH_API complex	zneg(complex z);

#define	zv_copy(x,y)	_zv_copy(x,y,0)
#define	zm_copy(A,B)	_zm_copy(A,B,0,0)

#define	z_input()	z_finput(stdin)
#define	zv_input(x)	zv_finput(stdin,x)
#define	zm_input(A)	zm_finput(stdin,A)
#define	z_output(z)	z_foutput(stdout,z)
#define	zv_output(x)	zv_foutput(stdout,x)
#define	zm_output(A)	zm_foutput(stdout,A)

#define	ZV_FREE(x)	( zv_free(x), (x) = ZVNULL )
#define	ZM_FREE(A)	( zm_free(A), (A) = ZMNULL )

#define	zin_prod(x,y)	_zin_prod(x,y,0,Z_CONJ)

#define	zv_norm1(x)	_zv_norm1(x,VNULL)
#define	zv_norm2(x)	_zv_norm2(x,VNULL)
#define	zv_norm_inf(x)	_zv_norm_inf(x,VNULL)

MESCH__END_DECLS

#endif
