/*
 * Hqp_HL_Gangster.C -- 
 *   - class definition
 *
 * rf, 7/19/94
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

#include <If_Float.h>

#include "Hqp_HL_Gangster.h"
#include "Hqp_Program.h"
#include "Hqp_SqpProgram.h"

IF_CLASS_DEFINE("Gangster", Hqp_HL_Gangster, Hqp_HL);

//--------------------------------------------------------------------------
Hqp_HL_Gangster::Hqp_HL_Gangster()
{
  _v = VNULL;
  _sQ = VNULL;
  _Qs = VNULL;

  _gamma = 0.2;
  
  _ifList.append(new If_Float("sqp_hela_gamma", &_gamma));
}

//--------------------------------------------------------------------------
Hqp_HL_Gangster::~Hqp_HL_Gangster()
{
  v_free(_v);
  v_free(_sQ);
  v_free(_Qs);
}

//--------------------------------------------------------------------------
void Hqp_HL_Gangster::update_Q(const VEC *s, const VEC *u, Hqp_SqpProgram *prg)
{
  SPMAT *Q = prg->qp()->Q;
  Real	sQs, sv;
  Real	*v_ve, *Qs_ve, *sQ_ve;
  Real	theta;
  int	i, j, j_idx, len, m;
  SPROW	*row;
  row_elt *elt;

  sv = in_prod(s, u);

  _Qs = sp_mv_mlt(Q, s, _Qs);
  _sQ = sp_vm_mlt(Q, s, _sQ);
  sQs = in_prod(_sQ, s);

  // Powells correcture
  if (sv < _gamma * sQs) {
    theta = (1.0 - _gamma) * sQs / (sQs - sv);
    _v = sv_mlt(theta, u, _v);
    v_mltadd(_v, _Qs, 1.0 - theta, _v);
    sv = in_prod(s, _v);
  } else {
    _v = v_copy(u, _v);
  }

  // Hesse matrix update with Gangster operator
  v_ve = _v->ve;
  sQ_ve = _sQ->ve;
  Qs_ve = _Qs->ve;
  m = Q->m;
  row = Q->row;
  for (i=0; i<m; i++, row++) {
    elt = row->elt;
    len = row->len;
    for (j_idx=0; j_idx<len; j_idx++, elt++) {
      j = elt->col;
      elt->val -= Qs_ve[i] * sQ_ve[j] / sQs;
      elt->val += v_ve[i] * v_ve[j] / sv;
    }
  }
}

//=========================================================================
