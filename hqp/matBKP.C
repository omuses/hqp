/*
 * BKP factorization -- sample implementation
 *  - reference: J.R.Bunch, L.Kaufman, and B.N.Parlett:
 *      Decomposition of a Symmetric Matrix, Numer.Math. 27, 95--109 (1976)
 *  - factorize band matrix A in situ into P'AP = MDM'
 *  - M is unit lower triangular
 *  - D is block diagonal with blocks 1x1 and 2x2
 *  - let strictly lower part of A untouched
 *  - P may only be used together with M by BKP-solve routine
 *
 * rf, 9/9/94
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

#define A(i,j) A->me[i][j]

static void interchange(MAT *A, int i, int j)
/*
 * - interchange row and col j of A and row and col i
 * - i < j
 * - A is the reduced matrix of order n-i+1
 */
{
  Real	tmp;
  int	n;
  int	k;

  n = A->n;
  for (k = j+1; k < n; k++) {
    tmp = A(j,k);
    A(j,k) = A(i,k);
    A(i,k) = tmp;
  }
  for (k = i+1; k < j; k++) {
    tmp = A(i,k);
    A(i,k) = A(k,j);
    A(k,j) = tmp;
  }
  tmp = A(i,i);
  A(i,i) = A(j,j);
  A(j,j) = tmp;
}


MAT *matBKPfactor(MAT *A, PERM *pivot)
/*
 * - factorize A in situ into P'AP = MDM'
 * - P(i+1) == 0 for (i,i+1) is a 2x2 block
 */
{
  int	i, ip1, ip2, j, k, n;
  Real	aii, aip1, aiip1, lambda, sigma, tmp;
  Real	det, s, t;

  if (!A || !pivot)
    m_error(E_NULL,"matBKPfactor");
  if (A->n != pivot->size)
    m_error(E_SIZES,"matBKPfactor");

  n = A->n;
  px_ident(pivot);

  for (i = 0; i < n-1;) {
    ip1 = i+1;
    ip2 = i+2;
    aii = fabs(A(i,i));

    /*
     * find the maximum element in the first column of the reduced
     * matrix below the diagonal (go through first row as A is symmetric)
     */
    lambda = 0.0;
    j = i;
    for (k = ip1; k < n; k++) {
      tmp = fabs(A(i,k));
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
    for (k = ip1; k < j; k++) {
      tmp = fabs(A(k,j));
      if (tmp > sigma)
	sigma = tmp;
    }
    for (k = j+1; k < n; k++) {
      tmp = fabs(A(j,k));
      if (tmp > sigma)
	sigma = tmp;
    }
    if (sigma * aii >= t * lambda)
      goto onebyone;
    if (fabs(A(j,j)) >= alpha * sigma) {
      interchange(A, i, j);
      pivot->pe[i] = j;
      goto onebyone;
    }

    /*
     * do 2x2 pivot
     */
    if (j > ip1) {
      interchange(A, ip1, j);
      tmp = A(i,j);
      A(i,j) = A(i,ip1);
      A(i,ip1) = tmp;
    }
    pivot->pe[i] = j;
    pivot->pe[ip1] = 0;
    tmp = A(i,ip1);
    det = A(i,i) * A(ip1,ip1) - tmp * tmp;
    aiip1 = tmp / det;
    aii = A(i,i) / det;
    aip1 = A(ip1,ip1) / det;
    for (j = ip2; j < n; j++) {
      s = - aiip1 * A(ip1,j) + aip1 * A(i,j);
      t = - aiip1 * A(i,j) + aii * A(ip1,j);
      for (k = j; k < n; k++)
	A(j,k) -= s * A(i,k) + t * A(ip1,k);
      A(i,j) = s;
      A(ip1,j) = t;
    }
    i = ip2;
    continue;

  onebyone:
    /*
     * do 1x1 pivot
     */
    aii = A(i,i);
    if (aii != 0.0) {
      for (j = ip1; j < n; j++) {
	s = A(i,j) / aii;
	for (k = j; k < n; k++)
	  A(j,k) -= s * A(i,k);
	A(i,j) = s;
      }
    }
    i = ip1;
  }

  return A;
}

VEC *matBKPsolve(const MAT *A, const PERM *pivot, const VEC *b, VEC *x)
/*
 * - solve A*x = b after A has been factorized by matBKPfactor
 * - raise an E_SING if A is singular
 */	
{
  int i, ii, j, k, ip1;
  int n;
  Real det, tmp, save;
  Real aiip1, aii, aip1;
  Real *x_ve;
  u_int *p_pe;

  if (!A || !pivot || !b) 
    m_error(E_NULL, "matBKPsolve");
  n = A->n;
  if ((int)b->dim != n || (int)pivot->size != n)
    m_error(E_SIZES, "matBKPsolve");
  if (!x || (int)x->dim != n)
    x = v_resize(x,n);

  p_pe = pivot->pe;

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
      aii = A(i,i);
      if (aii == 0.0)
	m_error(E_SING, "matBKPsolve");
      x_ve[i] = save / aii;
      for (j = ip1; j < n; j++)
	x_ve[j] -= save * A(i,j);
      i = ip1;
    }
    else {
      /*
       * 2x2 pivot
       */
      tmp = x_ve[i];
      x_ve[p_pe[i]] = x_ve[ip1];
      aii = A(i,i);
      aip1 = A(ip1,ip1);
      aiip1 = A(i,ip1);
      det = aii * aip1 - aiip1 * aiip1;
      if (det == 0.0)
	m_error(E_SING, "matBKPsolve");
      x_ve[i] = (tmp * aip1 - save * aiip1) / det;
      x_ve[ip1] = (save * aii - tmp * aiip1) / det;
      for (j = i+2; j < n; j++)
	x_ve[j] -= tmp * A(i,j) + save * A(ip1,j);
      i += 2;
    }
  }
  if (i == n-1) {
    aii = A(i,i);
    if (aii == 0.0)
      m_error(E_SING, "matBKPsolve");
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
      for (j = i+1; j < n; j++)
	tmp -= A(k,j) * x_ve[j];
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
