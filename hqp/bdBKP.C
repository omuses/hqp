/*
 * BKP factorization for band matrices
 *  - reference: J.R.Bunch, L.Kaufman, and B.N.Parlett:
 *      Decomposition of a Symmetric Matrix, Numer.Math. 27, 95--109 (1976)
 *  - factorize band matrix A in situ into P'AP = MDM'
 *  - M is unit lower triangular
 *  - D is block diagonal with 1x1 or 2x2 blocks
 *  - let strictly lower part of A untouched
 *  - P may only be used together with M by BKP-solve routine
 *
 * rf, 9/9/94
 *
 *  - restrict the pivot row to be within the band width
 * rf, 6/28/94
 */

/*
    Copyright (C) 1994--2002  Ruediger Franke

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

#include <stdio.h>
#include <math.h>

#include "Meschach.h"

#define alpha	0.6403882032022076 /* = (1+sqrt(17))/8 */

static PERM *bd_relief(const BAND *A, PERM *relief)
{
  int i, j, n, lb, ub;
  int offs;
  Real **A_me;

  if (!A)
    m_error(E_NULL,"bd_relief");
  n = A->mat->n;
  lb = A->lb;
  ub = A->ub;
  A_me = A->mat->me;
  if (!relief || (int)relief->size != n)
    relief = px_resize(relief, n);

  for (i = 0; i < n; i++) {
    j = min(i+ub, n-1);
    offs = lb - i + j;
    for (; j >= i; j--, offs--)
      if (A_me[offs][j] != 0.0)
	break;
    relief->pe[i] = j + 1;
  }

  return relief;
}


static void bdBKPchange(BAND *A, int i, int j, PERM *relief)
/*
 * - interchange row and col j of A and row and col i
 * - A is the reduced matrix of order n-i+1, i < j
 */
{
  Real	tmp, **A_me;
  int	lb, ub2, n;
  int	offs1, offs2;
  int	k, k_end;
  u_int	*re_ve;

  A_me = A->mat->me;
  re_ve = relief->pe;
  n = A->mat->n;
  lb = A->lb;
  ub2 = A->ub;
  offs1 = re_ve[i];
  offs2 = re_ve[j];
  re_ve[i] = offs2;
  re_ve[j] = offs1;
  k_end = max(offs1, offs2);
  // check storage limits of the band data structure
  if (k_end > i + 1 + ub2 || i >= j)
    m_error(E_INTERN, "bdBKPchange");
  k = j + 1;
  offs1 = lb + 1; /* - j + k */
  offs2 = lb - i + k;
  for (; k < k_end; k++, offs1++, offs2++) {
    tmp = A_me[offs1][k];
    A_me[offs1][k] = A_me[offs2][k];
    A_me[offs2][k] = tmp;
  }
  k = i + 1;
  offs1 = lb + 1; /* - i + k */
  offs2 = lb - k + j;
  for (; k < j; k++, offs1++, offs2--) {
    tmp = A_me[offs1][k];
    A_me[offs1][k] = A_me[offs2][j];
    A_me[offs2][j] = tmp;
    re_ve[k] = max((int)re_ve[k], j+1);
  }
  tmp = A_me[lb][i];
  A_me[lb][i] = A_me[lb][j];
  A_me[lb][j] = tmp;
}


BAND *bdBKPfactor(BAND *A, PERM *pivot, PERM *relief)
/*
 * - factorize A in situ into P'AP = MDM'
 * - P(i+1) == 0 for (i,i+1) is a 2x2 block
 * - A->ub is extended for pivoting
 * - use relief to store length of rows
 */
{
  int	i, ip1, ip2, j, k, n;
  int	k_end, ub, ub2, lb;
  u_int	*re_ve;
  Real	aii, aip1, aiip1, lambda, sigma, tmp;
  Real	det, s, t;
  Real	**A_me, *A_diag;
  int	offs1, offs2;

  if (!A || !pivot || !relief)
    m_error(E_NULL,"bdBKPfactor");
  n = A->mat->n;
  if (n != (int)pivot->size || n != (int)relief->size)
    m_error(E_SIZES,"bdBKPfactor");

  /*
   * initialize relief
   * extend band matrix for pivoting
   */

  bd_relief(A, relief);

  lb = A->lb;
  ub = A->ub;
  ub2 = ub * 2;
  ub2 = min(ub2, n-1);
  A = bd_resize(A, lb, ub2, n);
  A_me = A->mat->me;
  A_diag = A_me[lb];
  re_ve = relief->pe;

  px_ident(pivot);
  for (i = 0; i < n-1;) {
    ip1 = i+1;
    ip2 = i+2;
    aii = fabs(A_diag[i]);

    /*
     * find the maximum element in the first column of the reduced
     * matrix below the diagonal (go through first row as A is symmetric)
     */
    lambda = 0.0;
    j = i;
    // restrict the pivot row to be within the band width
    k_end = min((int)re_ve[i], i + ub);
    k = ip1;
    offs1 = lb + 1; /* - i + k */
    for (; k < k_end; k++, offs1++) {
      tmp = fabs(A_me[offs1][k]);
      if (tmp > lambda) {
	j = k;
	lambda = tmp;
      }
    }
    t = alpha * lambda;
    if (aii >= t)
      goto onebyone;
    
    /* 
     * determine the maximum element in the jth column of the 
     * reduced matrix off the diagonal
     */
    sigma = lambda;
    k = ip1;
    offs1 = lb - k + j;
    for (; k < j; k++, offs1--) {
      tmp = fabs(A_me[offs1][j]);
      if (tmp > sigma)
	sigma = tmp;
    }
    k_end = re_ve[j];
    k = j + 1;
    offs1 = lb + 1; /* - j + k */
    for (; k < k_end; k++, offs1++) {
      tmp = fabs(A_me[offs1][k]);
      if (tmp > sigma)
	sigma = tmp;
    }
    if (sigma * aii >= t * lambda)
      goto onebyone;
    if (fabs(A_diag[j]) >= alpha * sigma) {
      bdBKPchange(A, i, j, relief);
      pivot->pe[i] = j;
      goto onebyone;
    }

    /*
     * do 2x2 pivot
     */
    if (j > ip1) {
      bdBKPchange(A, ip1, j, relief);
      offs1 = lb - i + j;
      offs2 = lb + 1; /* - i + ip1 */
      tmp = A_me[offs1][j];
      A_me[offs1][j] = A_me[offs2][ip1];
      A_me[offs2][ip1] = tmp;
    }
    pivot->pe[i] = j;
    pivot->pe[ip1] = 0;
    tmp = A_me[lb+1][ip1]; /* lb - i + ip1 */
    det = A_diag[i] * A_diag[ip1] - tmp * tmp;
    aiip1 = tmp / det;
    aii = A_diag[i] / det;
    aip1 = A_diag[ip1] / det;
    k_end = max(re_ve[i], re_ve[ip1]);
    if (k_end > ip1 + ub2)
      m_error(E_INTERN, "bdBKPfactor");
    j = ip2;
    offs1 = lb + j;
    for (; j < k_end; j++, offs1++) {
      s = - aiip1 * A_me[offs1-ip1][j] + aip1 * A_me[offs1-i][j];
      t = - aiip1 * A_me[offs1-i][j] + aii * A_me[offs1-ip1][j];
      k = j;
      offs2 = lb + k;
      for (; k < k_end; k++, offs2++)
	A_me[offs2-j][k] -= s * A_me[offs2-i][k] + t * A_me[offs2-ip1][k];
      re_ve[j] = max((int)re_ve[j], k_end);
      A_me[offs1-i][j] = s;
      A_me[offs1-ip1][j] = t;
    }
    re_ve[i] = max((int)re_ve[i], k_end);
    re_ve[ip1] = max((int)re_ve[ip1], k_end);
    i = ip2;
    continue;

  onebyone:
    /*
     * do 1x1 pivot
     */
    aii = A_diag[i];
    if (aii != 0.0) {
      k_end = re_ve[i];
      j = ip1;
      offs1 = lb + 1; /* - i + j */
      for (; j < k_end; j++, offs1++) {
	s = A_me[offs1][j] / aii;
	k = j;
	offs2 = lb + k;
	for (; k < k_end; k++, offs2++)
	  A_me[offs2-j][k] -= s * A_me[offs2-i][k];
	re_ve[j] = max((int)re_ve[j], k_end);
	A_me[offs1][j] = s;
      }
      re_ve[i] = max((int)re_ve[i], k_end);
    }
    i = ip1;
  }

  return A;
}

VEC *bdBKPsolve(const BAND *A, const PERM *pivot, const PERM *relief,
		const VEC *b, VEC *x)
/*
 * - solve A*x = b after A has been factorized by bdBKPfactor
 * - raise an E_SING if A is singular
 */	
{
  int i, ii, j, k, ip1;
  int n, j_end;
  int lb, offs1;
  Real det, tmp, save;
  Real aiip1, aii, aip1;
  Real *x_ve, **A_me;
  u_int *p_pe, *re_ve;

  if (!A || !pivot || !b) 
    m_error(E_NULL, "bdBKPsolve");
  n = A->mat->n;
  if ((int)b->dim != n || (int)pivot->size != n || (int)relief->size != n)
    m_error(E_SIZES, "bdBKPsolve");
  if (!x || (int)x->dim != n)
    x = v_resize(x,n);

  p_pe = pivot->pe;
  re_ve = relief->pe;  
  A_me = A->mat->me;
  lb = A->lb;

  /*
   * solve MDy = b for y, where b is stored in x and store y in x
   */
  
  x = v_copy(b,x);
  x_ve = x->ve;
  for (i = 0; i < n-1;) {
    ip1 = i+1;
    save = x_ve[p_pe[i]];
    if (p_pe[ip1] > 0) {
      /*
       * 1x1 pivot
       */
      x_ve[p_pe[i]] = x_ve[i];
      aii = A_me[lb][i];
      if (aii == 0.0)
	m_error(E_SING, "bdBKPsolve");
      x_ve[i] = save / aii;
      j_end = re_ve[i];
      j = ip1;
      offs1 = lb + 1; /* - i + j */
      for (; j < j_end; j++, offs1++)
	x_ve[j] -= save * A_me[offs1][j];
      i = ip1;
    }
    else {
      /*
       * 2x2 pivot
       */
      tmp = x_ve[i];
      x_ve[p_pe[i]] = x_ve[ip1];
      aii = A_me[lb][i];
      aip1 = A_me[lb][ip1];
      aiip1 = A_me[lb+1][ip1]; /*lb - i + ip1 */
      det = aii * aip1 - aiip1 * aiip1;
      if (det == 0.0)
	m_error(E_SING, "bdBKPsolve");
      x_ve[i] = (tmp * aip1 - save * aiip1) / det;
      x_ve[ip1] = (save * aii - tmp * aiip1) / det;
      j_end = max(re_ve[i], re_ve[ip1]);
      j = i + 2;
      offs1 = lb + j;
      for (; j < j_end; j++, offs1++)
	x_ve[j] -= tmp * A_me[offs1-i][j] + save * A_me[offs1-ip1][j];
      i += 2;
    }
  }
  if (i == n-1) {
    aii = A_me[lb][i];
    if (aii == 0.0)
      m_error(E_SING, "bdBKPsolve");
    x_ve[i] /= aii;
    i = n - 2;
  }
  else
    i = n - 3;

  /*
   * solve M'x = y for x, where y is stored in x
   */
  while (i >= 0) {
    if (p_pe[i] > 0 || i == 0)
      ii = i;       /*  1x1 pivot */
    else
      ii = i-1;     /*  2x2 pivot */
    for (k = ii; k <= i; k++) {
      tmp = x_ve[k];
      j_end = re_ve[k];
      j = i + 1;
      offs1 = lb - k + j;
      for (; j < j_end; j++, offs1++)
	tmp -= A_me[offs1][j] * x_ve[j];
      x_ve[k] = tmp;
    }
    if (i != (int)p_pe[ii]) {
      tmp = x_ve[i];
      x_ve[i] = x_ve[p_pe[ii]];
      x_ve[p_pe[ii]] = tmp;
    }
    i = ii - 1;
  }

  return x;
}
