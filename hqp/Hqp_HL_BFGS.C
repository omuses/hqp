/*
 * Hqp_HL_BFGS.C -- 
 *   - class definition
 *
 * rf, 7/19/94
 *
 * rf, 2/13/97
 *  - consider an initial Hessian already given in init()
 *
 * rf, 4/13/99
 *  - change eigen_control default to false,
 *    (because Meschach's symmeig() may have endless iteration)
 * rf, 10/27/00
 *  - change eigen_control default to true again (+bug fix in Meschach)
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

#include <math.h>
extern "C" {
#include <meschach/matrix2.h>
}

#include <If_Int.h>
#include <If_Real.h>
#include <If_Bool.h>

#include "Hqp_HL_BFGS.h"
#include "Hqp_Program.h"
#include "Hqp_SqpProgram.h"

IF_CLASS_DEFINE("BFGS", Hqp_HL_BFGS, Hqp_HL);

//--------------------------------------------------------------------------
Hqp_HL_BFGS::Hqp_HL_BFGS()
{
  _v = VNULL;
  _sQ = VNULL;
  _Qs = VNULL;
  _b_Q = m_get(1, 1);

  _gamma = 0.1;
  _eigen_control = true;
  _bsize = -1;
  _max_bsize = -1;
  
  _ifList.append(new If_Real("sqp_hela_gamma", &_gamma));
  _ifList.append(new If_Bool("sqp_hela_eigen_control", &_eigen_control));
  _ifList.append(new If_Int("sqp_hela_bsize", &_bsize));
  _ifList.append(new If_Int("sqp_hela_max_bsize", &_max_bsize));
}

//--------------------------------------------------------------------------
Hqp_HL_BFGS::~Hqp_HL_BFGS()
{
  v_free(_v);
  v_free(_sQ);
  v_free(_Qs);
  m_free(_b_Q);
}

//--------------------------------------------------------------------------
void Hqp_HL_BFGS::setup(Hqp_SqpProgram *prg)
{
  SPMAT *Q = prg->qp()->Q;
  int	offs, size, old_size;
  MAT 	*b_Q_ones;

  // allocate entries for the block-diagonal structure

  _b_begin = 0;
  old_size = _b_Q->n;
  b_Q_ones = m_get(old_size, old_size);
  m_ones(b_Q_ones);
  _max_bsize = -1;

  while(next_block(Q, &offs, &size)) {
    _max_bsize = max(_max_bsize, size);
    if (size != old_size) {
      _b_Q = m_resize(_b_Q, size, size);
      b_Q_ones = m_resize(b_Q_ones, size, size);
      m_ones(b_Q_ones);
      old_size = size;
    }
    symsp_extract_mat(Q, offs, _b_Q);
    symsp_insert_symmat(Q, offs, b_Q_ones);	// allocate entries
    sp_update_mat(Q, offs, offs, _b_Q);
  }

  m_free(b_Q_ones);
}

#if 0
//--------------------------------------------------------------------------
static void restart_Q(MAT *Q, Real eps)
{
  Real m_norm;

  m_norm = m_norm_inf(Q);
  m_norm = max(m_norm, eps);
  m_norm = min(m_norm, 1.0/eps);
  m_ident(Q);
  sm_mlt(m_norm, Q, Q);
}
#else
//--------------------------------------------------------------------------
static void restart_Q(MAT *Q, Real eps)
{
  int i, j, dim;
  Real row_max, val;
  Real *Q_re;

  dim = Q->m;
  if (dim < 1)
    return;

  for (i = 0; i < dim; i++) {
    Q_re = Q->me[i];
    row_max = eps;
    for (j = 0; j < dim; j++) {
      val = fabs(Q_re[j]);
      row_max = max(row_max, val);
      Q_re[j] = 0.0;
    }
    row_max = min(row_max, 1.0/eps);
    Q_re[i] = row_max;
  }
}
#endif

//--------------------------------------------------------------------------
void Hqp_HL_BFGS::update_b_Q(const VEC *s, const VEC *u, Real alpha, MAT *Q)
{
  Real sv, sQs;
  Real theta, gamma;
  Real *v_ve, *Qs_ve, *sQ_ve, *Q_re;
  int i, i_end;
  int j, j_end;

  sv = in_prod(s, u);

  _sQ = vm_mlt(Q, s, _sQ);	// Q has lower and upper diagonal parts
  _Qs = mv_mlt(Q, s, _Qs);
  sQs = in_prod(_sQ, s);

  if (_gamma >= 0.0) {
    // standard damping
    gamma = _gamma;
  }
  else {
    // damping adapted to step length
    gamma = -_gamma;
    gamma = gamma + (1.0 - gamma) * (1.0 - alpha);
  }

  // Powells modification
  if (sv < gamma * sQs) {
    theta = (1.0 - gamma) * sQs / (sQs - sv);
    _v = sv_mlt(theta, u, _v);
    v_mltadd(_v, _Qs, 1.0 - theta, _v);
    sv = in_prod(s, _v);
  }
  else {
    _v = v_copy(u, _v);
  }

  if (!(sv != 0.0) || !(sQs != 0.0))
    return;

  v_ve = _v->ve;
  sQ_ve = _sQ->ve;
  Qs_ve = _Qs->ve;
  i_end = Q->m;
  j_end = Q->n;
  for (i = 0; i < i_end; i++) {
    Q_re = Q->me[i];
    for (j = i; j < j_end; j++) {	// update upper diagonal part
      Q_re[j] -= Qs_ve[i] * sQ_ve[j] / sQs;
      Q_re[j] += v_ve[i] * v_ve[j] / sv;
      if (_eigen_control)
	Q->me[j][i] = Q_re[j];		// update lower diagonal part
    }
  }

  if (_eigen_control) {
    theta = _eps * _eps;
    if (sQs < theta && sQs >= 0.0)
      theta = sQs;
    _sQ = symmeig(Q, MNULL, _sQ);
    sQs = v_min(_sQ, NULL);
    sQs -= theta;
    if (sQs < 0.0) {
      for (i = 0; i < i_end; i++)
	Q->me[i][i] -= sQs;
      if (_logging) {
	fprintf(stderr, "Eigenvalues: %g ... %g --> ",
		v_min(_sQ, NULL), v_max(_sQ, NULL));
	_sQ = symmeig(Q, MNULL, _sQ);
	fprintf(stderr, "%g ... %g\n",
		v_min(_sQ, NULL), v_max(_sQ, NULL));
      }
    }
  }
}

//--------------------------------------------------------------------------
void Hqp_HL_BFGS::update(const VEC *s, const VEC *u, Real alpha,
			 Hqp_SqpProgram *prg)
{
  SPMAT	*Q = prg->qp()->Q;
  VEC	b_s, b_u;
  int	offs, size;

  b_u.max_dim = 0;
  b_s.max_dim = 0;
  b_u.dim = _b_Q->n;
  b_s.dim = _b_Q->n;

  _b_begin = 0;

  while(next_block(Q, &offs, &size)) {
    if (size != (int)b_u.dim) {
      b_u.dim = size;
      b_s.dim = size;
      _b_Q = m_resize(_b_Q, size, size);
    }
    symsp_extract_mat(Q, offs, _b_Q); // build lower and upper diagonal parts
    b_u.ve = u->ve + offs;
    b_s.ve = s->ve + offs;
    update_b_Q(&b_s, &b_u, alpha, _b_Q);
    symsp_insert_symmat(Q, offs, _b_Q);	// insert upper diagonal part
  }
}

//--------------------------------------------------------------------------
// next_block:
//   find next diagonal block in Q, starting with index _b_begin
//--------------------------------------------------------------------------
int Hqp_HL_BFGS::next_block(const SPMAT *Q, int *offs, int *size)
{
  int	b_end, max_col;
  SPROW	*row;

  if (_b_begin >= Q->m)
    return 0;

  *offs = _b_begin;

  // _bsize < 0: automatic detection of the block size
  // _bsize = 0: full Hessian
  // _bsize > 0: block of dimension _bsize

  if (_bsize < 0) {
    b_end = _b_begin;
    while (_b_begin <= b_end) {
      row = &(Q->row[_b_begin]);
      if (row->len > 0)
	max_col = row->elt[row->len - 1].col;
      else
	max_col = _b_begin;
      b_end = max(b_end, max_col);
      _b_begin ++;
    }
    *size = b_end - *offs + 1;
  }
  else {
    *size = Q->n - _b_begin;
    if (_bsize > 0)
      *size = min(_bsize, *size);
    _b_begin += *size;
  }

  return 1;
}

//=========================================================================
