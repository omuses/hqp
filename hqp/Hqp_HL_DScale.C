/*
 * Hqp_HL_DScale.C -- 
 *   - class definition
 *
 * rf, 10/26/95
 *
 * rf, 2/13/97
 *  - consider an initial Hessian already given in init()
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

#include <If_Int.h>
#include <If_Float.h>

#include "Hqp_HL_DScale.h"
#include "Hqp_Program.h"
#include "Hqp_SqpProgram.h"

IF_CLASS_DEFINE("DScale", Hqp_HL_DScale, Hqp_HL);

//--------------------------------------------------------------------------
Hqp_HL_DScale::Hqp_HL_DScale()
{
  _b_q = VNULL;
  _sq = VNULL;
  _v = VNULL;
  _bsize = 0;
  _gamma = 0.2;

  _ifList.append(new If_Int("sqp_hela_bsize", &_bsize));
  _ifList.append(new If_Float("sqp_hela_gamma", &_gamma));
}

//--------------------------------------------------------------------------
Hqp_HL_DScale::~Hqp_HL_DScale()
{
  v_free(_b_q);
  v_free(_sq);
  v_free(_v);
}

//--------------------------------------------------------------------------
void Hqp_HL_DScale::setup(Hqp_SqpProgram *prg)
{
  SPMAT *Q = prg->qp()->Q;
  SPROW *row;
  int i, dim = Q->m;
  row_elt *elt;
  int j_idx, len;

  // clear out all off-diagonal entries
  for (i = 0; i < dim; i++) {
    row = Q->row + i;
    elt = row->elt;
    len = row->len;
    for (j_idx = 0; j_idx < len; j_idx++) {
      if (elt->col != i)
	elt->val = 0.0;
      else
	elt->val = max(elt->val, _eps);
      elt++;
    }
  }
  sp_compact(Q, 0.0);

  // ensure all diagonal entries
  for (i = 0; i < dim; i++) {
    if (sprow_idx(Q->row + i, i) < 0)
      sp_set_val(Q, i, i, _eps);
  }
}

//--------------------------------------------------------------------------
void Hqp_HL_DScale::update_b_Q(const VEC *s, const VEC *u, VEC *q)
{
  Real sv, sqs;
  Real theta;
  Real *v_ve, *sq_ve, *q_ve;
  int i, i_end;

  sv = in_prod(s, u);

  /*
  // special treatment if sv == 0.0
  if (sv == 0.0) {
    //v_set(q, _eps);
    return;
  }
  */

  _sq = v_star(s, q, _sq);
  sqs = in_prod(_sq, s);

  // Powells correcture
  if (sv < _gamma * sqs) {
    theta = (1.0 - _gamma) * sqs / (sqs - sv);
    _v = sv_mlt(theta, u, _v);
    v_mltadd(_v, _sq, 1.0 - theta, _v);
    sv = in_prod(s, _v);
  }
  else {
    _v = v_copy(u, _v);
  }

  if (!(sv != 0.0) || !(sqs != 0.0))
    return;

  v_ve = _v->ve;
  sq_ve = _sq->ve;
  q_ve = q->ve;
  i_end = q->dim;
  for (i=0; i<i_end; i++) {
    q_ve[i] -= sq_ve[i] * sq_ve[i] / sqs;
    q_ve[i] += v_ve[i] * v_ve[i] / sv;
  }
}

//--------------------------------------------------------------------------
void Hqp_HL_DScale::update(const VEC *s, const VEC *u, Real,
			   Hqp_SqpProgram *prg)
{
  SPMAT	*Q = prg->qp()->Q;
  SPROW *row;
  int 	i, m;
  VEC	b_s, b_u;
  int	offs, size, bsize;

  m = Q->m;
  offs = 0;
  bsize = _bsize > 0? _bsize: m;
  
  b_u.max_dim = 0;
  b_s.max_dim = 0;

  if (!Q->flag_diag)
    sp_diag_access(Q);
  
  while(offs < m) {
    size = min(bsize, m - offs);
    _b_q = v_resize(_b_q, size);
    b_s.dim = size;
    b_u.dim = size;

    row = Q->row + offs;
    for (i=0; i<size; i++, row++)
      _b_q->ve[i] = row->elt[row->diag].val;
    b_u.ve = u->ve + offs;
    b_s.ve = s->ve + offs;

    update_b_Q(&b_s, &b_u, _b_q);

    row = Q->row + offs;
    for (i=0; i<size; i++, row++)
      row->elt[row->diag].val = _b_q->ve[i];

    offs += size;
  }
}

//=========================================================================
