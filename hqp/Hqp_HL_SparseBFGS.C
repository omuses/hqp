/*
 * Hqp_HL_SparseBFGS.C -- 
 *   - class definition
 *
 * rf, 7/19/95
 */

/*
    Copyright (C) 1994--1998  Ruediger Franke

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
#include <If_Float.h>

#include "Hqp_HL_SparseBFGS.h"
#include "Hqp_Program.h"
#include "Hqp_SqpProgram.h"

#include "sprcm.h"

IF_CLASS_DEFINE("SparseBFGS", Hqp_HL_SparseBFGS, Hqp_HL);

//--------------------------------------------------------------------------
Hqp_HL_SparseBFGS::Hqp_HL_SparseBFGS()
{
  _gamma = 0.2;

  _Q2B = PNULL;
  _B2Q = PNULL;
  _v = VNULL;
  _sQ = VNULL;
  _Qs = VNULL;
  _b_Q = m_get(1, 1);
  _b_u = v_get(1);
  _b_s = v_get(1);
  
  _ifList.append(new If_Float("sqp_hela_gamma", &_gamma));
}

//--------------------------------------------------------------------------
Hqp_HL_SparseBFGS::~Hqp_HL_SparseBFGS()
{
  px_free(_Q2B);
  px_free(_B2Q);
  v_free(_v);
  v_free(_sQ);
  v_free(_Qs);
  m_free(_b_Q);
  v_free(_b_u);
  v_free(_b_s);
}

//--------------------------------------------------------------------------
void Hqp_HL_SparseBFGS::setup(Hqp_SqpProgram *prg)
{
  int i, idx, offs, size;
  SPMAT *Q = prg->qp()->Q;
  row_elt *elt;
  MAT 	*b_Q_ones;

  // find block diagonal structure

  _Q2B = sp_symrcm(Q, _Q2B);
  _B2Q = px_inv(_Q2B, _B2Q);

  // copy the upper diagonal part to the lower diagonal

  for (i = Q->m - 1; i >= 0; i--) {
    elt = Q->row[i].elt;
    size = Q->row[i].len;
    for (idx = 0; idx < size; idx++, elt++) {
      if (elt->col > i)
	sp_set_val(Q, elt->col, i, elt->val);
    }
  }

  // fill up the blocks

  pxinv_sprows(_Q2B, Q, Q);
  pxinv_spcols(_Q2B, Q, Q);
  _b_begin = 0;
  b_Q_ones = m_get(_b_Q->m, _b_Q->n);
  m_ones(b_Q_ones);
  while(next_block(Q, &offs, &size)) {
    if (!_b_Q || size != (int)_b_Q->m) {
      _b_Q = m_resize(_b_Q, size, size);
      b_Q_ones = m_resize(b_Q_ones, size, size);
      m_ones(b_Q_ones);
    }
    sp_extract_mat(Q, offs, offs, _b_Q);
    sp_insert_mat(Q, offs, offs, b_Q_ones);	// allocate entries
    sp_update_mat(Q, offs, offs, _b_Q);
  }
  pxinv_sprows(_B2Q, Q, Q);
  pxinv_spcols(_B2Q, Q, Q);

  m_free(b_Q_ones);
}

//--------------------------------------------------------------------------
void Hqp_HL_SparseBFGS::update_b_Q(const VEC *s, const VEC *u, Real, MAT *Q)
{
  Real sv, sQs;
  Real theta;
  Real *v_ve, *Qs_ve, *sQ_ve, *Q_re;
  int i, i_end;
  int j, j_end;

  sv = in_prod(s, u);

  /*
  // return if s or u are zero
  if (sv == 0.0) {
    m_ident(Q);
    sm_mlt(_eps, Q, Q);
    return;
  }
  */

  _sQ = vm_mlt(Q, s, _sQ);
  _Qs = mv_mlt(Q, s, _Qs);
  sQs = in_prod(_sQ, s);

  /*
  if (!(sQs > 0.0)) {
    m_ident(Q);
    sm_mlt(_eps, Q, Q);
    _sQ = vm_mlt(Q, s, _sQ);
    _Qs = mv_mlt(Q, s, _Qs);
    sQs = in_prod(_sQ, s);
  }
  */

  // Powells correcture
  if (sv < _gamma * sQs) {
    theta = (1.0 - _gamma) * sQs / (sQs - sv);
    if (!(theta <= 1.0)) {
      fprintf(stderr, "bingo %g\n", theta);
      theta = 1.0;
    }
    if (!(theta >= 0.0)) {
      fprintf(stderr, "bingo %g, %g, %g\n", theta, sQs, sv);
      theta = 0.0;
    }
    _v = sv_mlt(theta, u, _v);
    v_mltadd(_v, _Qs, 1.0 - theta, _v);
    sv = in_prod(s, _v);
  }
  else {
    _v = v_copy(u, _v);
  }

  if (!(sv != 0.0))
    return;

  v_ve = _v->ve;
  sQ_ve = _sQ->ve;
  Qs_ve = _Qs->ve;
  i_end = Q->m;
  j_end = Q->n;
  for (i = 0; i < i_end; i++) {
    Q_re = Q->me[i];
    for (j = 0; j < j_end; j++) {
      Q_re[j] -= Qs_ve[i] * sQ_ve[j] / sQs;
      Q_re[j] += v_ve[i] * v_ve[j] / sv;
    }
  }

  /*
  // positive definiteness against s
  _sQ = vm_mlt(Q, s, _sQ);
  sQs = in_prod(_sQ, s);
  if (!(sQs > 0.0)) {
    m_ident(Q);
    sm_mlt(_eps, Q, Q);
    return;
  }

  // positive definiteness against s
  _sQ = vm_mlt(Q, u, _sQ);
  sQs = in_prod(_sQ, u);
  if (!(sQs > 0.0)) {
    m_ident(Q);
    sm_mlt(_eps, Q, Q);
    return;
  }
  */

  /*
  // bounded
  sQs = m_norm_inf(Q);
  if (!(sQs >= _eps) || !(_eps*sQs < 1.0)) {
    m_ident(Q);
    sm_mlt(_eps, Q, Q);
    return;
  }
  */
}

//--------------------------------------------------------------------------
void Hqp_HL_SparseBFGS::update(const VEC *s, const VEC *u, Real alpha,
			       Hqp_SqpProgram *prg)
{
  SPMAT	*Q = prg->qp()->Q;
  int	i, offs, size;

  pxinv_sprows(_Q2B, Q, Q);
  pxinv_spcols(_Q2B, Q, Q);

  _b_begin = 0;

  while(next_block(Q, &offs, &size)) {
    _b_Q = m_resize(_b_Q, size, size);
    _b_u = v_resize(_b_u, size);
    _b_s = v_resize(_b_s, size);
    sp_extract_mat(Q, offs, offs, _b_Q);
    for (i = 0; i < size; i++) {
      _b_u->ve[i] = u->ve[_B2Q->pe[offs + i]];
      _b_s->ve[i] = s->ve[_B2Q->pe[offs + i]];
    }
    update_b_Q(_b_s, _b_u, alpha, _b_Q);
    sp_insert_mat(Q, offs, offs, _b_Q);
  }

  pxinv_sprows(_B2Q, Q, Q);
  pxinv_spcols(_B2Q, Q, Q);

/*
  FILE *fp;
  fp = fopen("hess.dat", "w");
//  sp_ones(Q);
  sp_foutput(fp, Q);
  fclose(fp);
  exit(0);
*/
}

//--------------------------------------------------------------------------
// next_block:
//   find next diagonal block in Q, starting with index _b_begin
//--------------------------------------------------------------------------
int Hqp_HL_SparseBFGS::next_block(const SPMAT *Q, int *offs, int *size)
{
  int	b_end, max_col;
  SPROW	*row;

  if (_b_begin >= Q->m)
    return 0;

  *offs = _b_begin;

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

  return 1;
}

//=========================================================================
