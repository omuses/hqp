/*
 * addon_hqp.c
 *   - some Meschach dense matrix add-ons and extensions
 *     developed for the Hqp project
 *
 * E. Arnold  03/07/97
 *
 */

/*
    Copyright (C) 1997--2002  Eckhard Arnold

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; 
    version 2 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library (file COPYING.LIB);
    if not, write to the Free Software Foundation, Inc.,
    59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <math.h>

#include "addon_hqp.h"
#include "matrix2.h"

/*
 * m_output_g -- Prints internal parameters of MAT * data structure
 * E. Arnold   10/08/96
 */
void m_output_g(const MAT *A)
{
    if ( ! A )
	printf("Matrix %p\n", A);
    else {
	printf("Matrix %p, m=%d, n=%d, max_m=%d, max_n=%d, ",
	       A, A->m, A->n, A->max_m, A->max_n);
	printf("max_size=%d, me=%p, base=%p\n",
	       A->max_size, A->me, A->base);
    }
}

/*
 * v_output_g -- Prints internal parameters of VEC * data structure
 * E. Arnold   10/08/96
 */
void v_output_g(const VEC *A)
{
    if ( ! A )
	printf("Vector %p\n", A);
    else
	printf("Vector %p, dim=%d, max_dim=%d, ve=%p\n",
	       A, A->dim, A->max_dim, A->ve);
}

/*
 * m_copy1 -- copies matrix into new area, resizes out to correct size
 *            same as m_copy(), but resizes out to correct size
 * E. Arnold   10/18/96
 */
MAT	*m_copy1(const MAT *in, MAT *out)
{
    u_int	i /* ,j */;

    if ( in==MNULL )
	m_error(E_NULL,"m_copy1");
    if ( in==out )
	return (out);
    if ( out==MNULL || out->m != in->m || out->n != in->n )
	out = m_resize(out,in->m,in->n);

    for ( i=0; i < in->m; i++ )
	MEM_COPY(&(in->me[i][0]),&(out->me[i][0]),
		 (in->n)*sizeof(Real));
    /* for ( j=0; j < in->n; j++ )
       out->me[i][j] = in->me[i][j]; */

    return (out);
}

/*
 * v_copy1 -- copies vector into new area, resizes out to correct size
 *            same as v_copy(), but resizes out to correct size.
 * E. Arnold   10/18/96
 */
VEC	*v_copy1(const VEC *in, VEC *out)
{
    /* u_int	i,j; */

    if ( in==VNULL )
	m_error(E_NULL,"_v_copy1");
    if ( in==out )
	return (out);
    if ( out==VNULL || out->dim != in->dim )
	out = v_resize(out,in->dim);

    MEM_COPY(&(in->ve[0]),&(out->ve[0]),(in->dim)*sizeof(Real));
    /* for ( i=0; i < in->dim; i++ )
		out->ve[i] = in->ve[i]; */

    return (out);
}

/*
 * v_diag -- get diagonal entries of a matrix
 * E. Arnold   02/22/97
 */
VEC *v_diag(MAT *A, VEC *C) {
    int i, n;
    if ( A == MNULL )
	m_error(E_NULL,"v_diag");
    n = A->m;
    if ( n >  (int) A->n ) 
	n = A->n;
    C = v_resize(C, n);
    for ( i = 0; i < n; i++ ) 
	v_set_val(C, i, m_entry(A, i, i));
    return C;
}

/*
 * dm_mlt -- diagonal matrix matrix multiplication
 * E. Arnold   02/22/97
 */
MAT *dm_mlt(MAT *A, VEC *B, MAT *C) {
    int i, j;
    double val;
    if ( ( A == MNULL ) || ( B == VNULL ) )
	m_error(E_NULL, "dm_mlt");
    if ( A->m !=  B->dim )
	m_error(E_SIZES, "dm_mlt");
    C = m_resize(C, A->m, A->n);
    for ( i = 0; i < (int) A->m; i++) {
	val = v_entry(B, i);
	for ( j = 0; j < (int)A->n; j++ )
	    m_set_val(C, i, j, val*m_entry(A, i, j));
    }
    return C;
}

/*
 * md_mlt -- matrix diagonal matrix multiplication
 * E. Arnold   02/22/97
 */
MAT *md_mlt(MAT *A, VEC *B, MAT *C) {
    int i, j;
    double val;
    if ( ( A == MNULL ) || ( B == VNULL ) )
	m_error(E_NULL, "md_mlt");
    if ( A->n !=  B->dim )
	m_error(E_SIZES, "md_mlt");
    C = m_resize(C, A->m, A->n);
    for ( i = 0; i < (int) A->n; i++ ) {
	val = v_entry(B, i);
	for ( j = 0; j < (int) A->m; j++ )
	    m_set_val(C, j, i, val*m_entry(A, j, i));
    }
    return C;
}

/*
 * m_symm -- make square matrix symmetric
 *           out = 0.5 * (in + in^T) 
 * E. Arnold   11/27/96
 *             07/03/97  bug: in->m != in->m  --> in->m != in->n
 */
MAT	*m_symm(const MAT *in, MAT *out)
{
    u_int	i , j;

    if ( in == MNULL )
	m_error(E_NULL,"m_symm");

    if ( in->m != in->n )
	m_error(E_SQUARE,"m_symm");

    if ( out == MNULL || out->m != in->m || out->n != in->n )
	out = m_resize(out, in->m, in->n);

    for ( i = 0; i < in->m; i++ )
	for ( j = i+1; j < in->m; j++ )
	    out->me[i][j] = out->me[j][i] = 0.5*(in->me[i][j]+in->me[j][i]);

    return (out);
}

/*
 * rel_symm -- check square matrix for symmetry
 *             rel_symm = max(max(abs(in-in')./(0.5*(abs(in)+abs(in')))
 * E. Arnold   11/27/96
 *             07/03/97  bug: in->m != in->m  --> in->m != in->n
 */
double	rel_symm(const MAT *in)
{
    u_int  i , j;
    double r;

    if ( in == MNULL )
	m_error(E_NULL,"rel_symm");

    if ( in->m != in->n )
	m_error(E_SQUARE,"rel_symm");

    for ( r = 0.0, i = 0; i < in->m; i++ )
	for ( j = i+1; j < in->m; j++ ) {
	    r = max(r, fabs(in->me[i][j]-in->me[j][i]));
	}
    return (r/m_norm_inf(in));
}

/*
 * CHsolveM -- CHsolve with matrix right hand side
 * E. Arnold   10/18/96
 *             07/03/97  bug: m_resize(C, A->m, B->m) --> 
 *                            m_resize(C, A->m, B->n)
 *                       bug: i < B->m --> i < B->n
 */
MAT *CHsolveM(MAT *A, MAT *B, MAT *C)
{
    static VEC *v1 = VNULL;
    u_int i;

    if ( ( A == MNULL ) || ( B == MNULL ) )
	m_error(E_NULL,"CHsolveM");
    if ( A->m != A->n )
	m_error(E_SQUARE,"CHsolveM");
    if ( A->m != B->m )
	m_error(E_SIZES,"CHsolveM");

    v1 = v_resize(v1, B->m);
    MEM_STAT_REG(v1, TYPE_VEC);
    C = m_resize(C, A->m, B->n);

    for (i = 0; i < B->n; i++) {
	v1 = get_col(B, i, v1);
	v1 = CHsolve(A, v1, v1);
	C = set_col(C, i, v1);
    }
    return C;
}

/*
 * CHsolveMT -- CHsolve with transposed matrix right hand side
 * E. Arnold   10/18/96
 */
MAT *CHsolveMT(MAT *A, MAT *B, MAT *C)
{
    static VEC *v1 = VNULL;
    u_int i;

    if ( ( A == MNULL ) || ( B == MNULL ) )
	m_error(E_NULL,"CHsolveMT");
    if ( A->m != A->n )
	m_error(E_SQUARE,"CHsolveMT");
    if ( A->m != B->n )
	m_error(E_SIZES,"CHsolveMT");
    if ( B == C)
	m_error(E_INSITU,"CHsolveMT");

    v1 = v_resize(v1, B->n);
    MEM_STAT_REG(v1, TYPE_VEC);
    C = m_resize(C, A->m, B->m);

    for (i = 0; i < B->m; i++) {
	v1 = get_row(B, i, v1);
	v1 = CHsolve(A, v1, v1);
	C = set_col(C, i, v1);
    }
    return C;
}

/*
 * BKPsolveM -- BKPsolve with matrix right hand side
 * E. Arnold   10/18/96
 */
MAT *BKPsolveM(MAT *A, PERM *pivot, PERM *blocks, MAT *B, MAT *C)
{
    static VEC *v1 = VNULL;
    u_int i;

    if ( ( A == MNULL ) || ( B == MNULL ) )
	m_error(E_NULL,"BKPsolveM");
    if ( A->m != A->n )
	m_error(E_SQUARE,"BKPsolveM");
    if ( A->m != B->m )
	m_error(E_SIZES,"BKPsolveM");

    v1 = v_resize(v1, B->m);
    MEM_STAT_REG(v1, TYPE_VEC);
    C = m_resize(C, A->m, B->n);

    for (i = 0; i < B->n; i++) {
	v1 = get_col(B, i, v1);
	v1 = BKPsolve(A, pivot, blocks, v1, v1);
	C = set_col(C, i, v1);
    }
    return C;
}

/*
 * BKPsolveMT -- BKPsolve with transposed matrix right hand side
 * E. Arnold   10/18/96
 */
MAT *BKPsolveMT(MAT *A, PERM *pivot, PERM *blocks, MAT *B, MAT *C)
{
    static VEC *v1 = VNULL;
    u_int i;

    if ( ( A == MNULL ) || ( B == MNULL ) )
	m_error(E_NULL,"BKPsolveMT");
    if ( A->m != A->n )
	m_error(E_SQUARE,"BKPsolveMT");
    if ( A->m != B->n )
	m_error(E_SIZES,"BKPsolveMT");
    if ( B == C)
	m_error(E_INSITU,"BKPsolveMT");

    v1 = v_resize(v1, B->n);
    MEM_STAT_REG(v1, TYPE_VEC);
    C = m_resize(C, A->m, B->m);

    for (i = 0; i < B->m; i++) {
	v1 = get_row(B, i, v1);
	v1 = BKPsolve(A, pivot, blocks, v1, v1);
	C = set_col(C, i, v1);
    }
    return C;
}

/*
 * LUsolveM -- LUsolve with matrix right hand side
 * E. Arnold   10/18/96
 */
MAT *LUsolveM(MAT *A, PERM *pivot, MAT *B, MAT *C)
{
    static VEC *v1 = VNULL;
    u_int i;

    if ( ( A == MNULL ) || ( B == MNULL ) )
	m_error(E_NULL,"LUsolveM");
    if ( A->m != A->n )
	m_error(E_SQUARE,"LUsolveM");
    if ( A->m != B->m )
	m_error(E_SIZES,"LUsolveM");

    v1 = v_resize(v1, B->m);
    MEM_STAT_REG(v1, TYPE_VEC);
    C = m_resize(C, A->m, B->n);

    for (i = 0; i < B->n; i++) {
	v1 = get_col(B, i, v1);
	v1 = LUsolve(A, pivot, v1, v1);
	C = set_col(C, i, v1);
    }
    return C;
}

/*
 * LUsolveMT -- LUsolve with transposed matrix right hand side
 * E. Arnold   10/18/96
 */
MAT *LUsolveMT(MAT *A, PERM *pivot, MAT *B, MAT *C)
{
    static VEC *v1 = VNULL;
    u_int i;

    if ( ( A == MNULL ) || ( B == MNULL ) )
	m_error(E_NULL,"LUsolveMT");
    if ( A->m != A->n )
	m_error(E_SQUARE,"LUsolveMT");
    if ( A->m != B->n )
	m_error(E_SIZES,"LUsolveMT");
    if ( B == C)
	m_error(E_INSITU,"LUsolveMT");

    v1 = v_resize(v1, B->n);
    MEM_STAT_REG(v1, TYPE_VEC);
    C = m_resize(C, A->m, B->m);

    for (i = 0; i < B->m; i++) {
	v1 = get_row(B, i, v1);
	v1 = LUsolve(A, pivot, v1, v1);
	C = set_col(C, i, v1);
    }
    return C;
}

/*
 * GE_QP -- Generalized elimination for quadratic programming
 * E. Arnold   10/08/96
 */
void GE_QP(MAT *A, MAT *S, MAT *Z, MAT *PQ, double eps)
{
    static VEC *diag = VNULL;
    static PERM *pivot = PNULL;
    static MAT *R1 = MNULL;
    static MAT *Q = MNULL;
    static MAT *R = MNULL;
    /*  double eps = 1.0e-8;  */
    u_int i, m, r;

    if ( ( A == MNULL ) || ( S == MNULL) || ( Z == MNULL ) || ( PQ == MNULL ) )
	m_error(E_NULL, "GE_QP");

    m = A->n;

    /*   Memory allocation   */

    diag = v_resize(diag, min(A->m, A->n));
    pivot = px_resize(pivot, A->m);
    R1 = m_resize(R1, m, m);
    Q = m_resize(Q, m, m);
    R = m_resize(R, m, m);

    /*   Registration of temporary variables   */

    MEM_STAT_REG(diag, TYPE_VEC);
    MEM_STAT_REG(pivot, TYPE_PERM);
    MEM_STAT_REG(R1, TYPE_MAT);
    MEM_STAT_REG(Q, TYPE_MAT);
    MEM_STAT_REG(R, TYPE_MAT);

    /*   1st QR factorization: A^T*P = Q*R   */

    PQ = m_transp(A, PQ);
    PQ = QRCPfactor(PQ, diag, pivot);
    Q = makeQ(PQ, diag, Q);
    R = m_resize(R, PQ->m, PQ->n);
    R = makeR(PQ, R);

    /*   find rank of A and form Z = Q(:,r+1:m)   */

    for (r = 0; r < min(R->m, R->n); r++)
	if (fabs(R->me[r][r]) < eps)
	    break;
    /*  printf("\n r = %d\n Q =", r); m_output(Q); */
    /*  printf("\n R ="); m_output(R); */
    Z = m_resize(Z, Q->m, Q->n-r);
    Z = m_move(Q, 0, r, Q->m, Q->n-r, Z, 0, 0);

    /*   2nd QR factorization: R^T=Q1*R1   */

    R1 = m_transp(R, R1);
    R1 = QRfactor(R1, diag);

    /*   form PQ = E*Qbar using S as temporary variable   */
  
    S = m_resize(S, R1->m, R1->m);
    S = makeQ(R1, diag, S);  
    pivot = px_inv(pivot, pivot);
    PQ = m_resize(PQ, S->m, S->n);
    PQ = px_rows(pivot, S, PQ);
  
    /*   form S = Q(:,1:r)*(R1)^-T   */
    R = m_resize(R, R1->m, R1->n);
    R = makeR(R1, R);
    R1 = m_resize(R1, r, r);
    R1 = m_move(R, 0, 0, r, r, R1, 0, 0);
    diag = v_resize(diag, r);
    S = m_resize(S, m, r);
    for (i = 0; i < m; i++) {
	diag = get_row(Q, i, diag);
	diag->dim = r;
	diag = UTsolve(R1, diag, diag, 0.0);
	S = set_row(S, i, diag);
    }
    /*  m_output(S);m_output(R1);printf("\nm = %d, r = %d\n",m,r); */
}

/*
 * m_concatc -- Matrix concatenation by rows C = [A; B]
 * E. Arnold   10/08/96
 */
MAT *m_concatc(MAT *A, MAT *B, MAT *C)
{
    int am;
    if ( ! A || ! B )
	m_error(E_NULL, "m_concatc");
    if ( A->n != B->n )
	m_error(E_SIZES, "m_concatc");
    if ( B == C ) 
	m_error(E_INSITU, "m_concatc"); 
    am = A->m;
    C = m_resize(C, am+B->m, A->n);
    if ( A != C )
	C = m_move(A, 0, 0, am, A->n, C, 0, 0);
    C = m_move(B, 0, 0, B->m, B->n, C, am, 0);
    return C;
}

/*
 * v_concat -- Vector concatenation C = [A; B]
 * E. Arnold   10/08/96
 */
VEC *v_concat(VEC *A, VEC *B, VEC *C)
{
    u_int adim;
    if ( ! A || ! B )
	m_error(E_NULL, "v_concat");
    if ( B == C ) 
	m_error(E_INSITU, "v_concat");
    adim = A->dim;
    C = v_resize(C, adim+B->dim);
    if ( A != C )
	C = v_move(A, 0, adim, C, 0);
    C = v_move(B, 0, B->dim, C, adim);
    return C;
}

/*
 * v_move_iv -- C(i) = A(iv(i))
 * E. Arnold   10/08/96
 */
VEC *v_move_iv(const VEC *A, const IVEC *iv, VEC *C)
{
    int i;
    if ( ! A || ! iv )
	m_error(E_NULL, "v_move_iv");
    if (  A == C )
	m_error(E_INSITU, "v_move_iv");
    C = v_resize(C, iv->dim);
    for ( i = 0; i < (int) iv->dim; i++)
	v_set_val(C, i, v_entry(A, iv_entry(iv, i)));
    return C;
}

/*
 * v_set_iv -- C(iv(i)) = A(i)
 * E. Arnold   10/08/96
 */
VEC *v_set_iv(VEC *A, IVEC *iv, VEC *C)
{
    int i;
    if ( ! A || ! iv || ! C)
	m_error(E_NULL, "v_set_iv");
    if ( A == C )
	m_error(E_INSITU, "v_set_iv");
    if ( ! ( A->dim == iv->dim ) )
	m_error(E_SIZES, "v_set_iv");
    for ( i = 0; i < (int) iv->dim; i++)
	v_set_val(C, iv_entry(iv, i), v_entry(A, i));
    return C;
}

/*
 * v_dist2 - 2-norm (Euclidean norm) of vector difference
 */
double v_dist2(const VEC *x, const VEC *y)
{
	int	i;
	double	sum;

	if ( x == (VEC *)NULL || y == (VEC *)NULL )
	    m_error(E_NULL, "v_dist_inf");
	if ( y->dim != x->dim )
	    m_error(E_SIZES, "v_dist_inf");

	for ( i = 0, sum = 0.0; i < x->dim; i++ ) 
	    sum += square(x->ve[i]-y->ve[i]);

	return sqrt(sum);
}

/* 
 * v_dist_inf -- infinity-norm (supremum norm) of vector difference 
 */
double v_dist_inf(const VEC *x, const VEC *y)
{
	int	i;
	double	maxval;

	if ( x == (VEC *)NULL || y == (VEC *)NULL )
	    m_error(E_NULL, "v_dist_inf");
	if ( y->dim != x->dim )
	    m_error(E_SIZES, "v_dist_inf");

	for ( i = 0, maxval = 0.0; i < x->dim; i++ ) 
	    maxval = max(maxval, fabs(x->ve[i]-y->ve[i]));

	return maxval;
}
