/*
 * Hqp_HL_Gerschgorin.C -- 
 *   - class definition
 *
 * rf, 10/3/95
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

#include "Hqp_HL_Gerschgorin.h"
#include "Hqp_Program.h"
#include "Hqp_SqpProgram.h"

IF_CLASS_DEFINE("Gerschgorin", Hqp_HL_Gerschgorin, Hqp_HL);

//--------------------------------------------------------------------------
Hqp_HL_Gerschgorin::Hqp_HL_Gerschgorin()
{
  _rowsum = VNULL;
}

//--------------------------------------------------------------------------
Hqp_HL_Gerschgorin::~Hqp_HL_Gerschgorin()
{
  v_free(_rowsum);
}

//--------------------------------------------------------------------------
void Hqp_HL_Gerschgorin::setup(Hqp_SqpProgram *prg)
{
  SPMAT *Q = prg->qp()->Q;
  int i, dim = Q->m;

  _rowsum = v_resize(_rowsum, dim);

  // ensure all diagonal entries

  for (i = 0; i < dim; i++) {
    if (sprow_idx(Q->row + i, i) < 0)
      sp_set_val(Q, i, i, 0.0);
  }

  update(VNULL, VNULL, 0.0, prg);
}

//--------------------------------------------------------------------------
void Hqp_HL_Gerschgorin::update(const VEC *, const VEC *, Real,
				Hqp_SqpProgram *prg)
{
  SPMAT *Q;
  SPROW *row;
  row_elt *elt;
  Real *rs_ve, val;
  int i, i_end;
  int j, j_idx, j_end;

  // enshure positive definiteness
  Q = prg->qp()->Q;

  if (!Q->flag_diag)
    sp_diag_access(Q);

  v_zero(_rowsum);
  row = Q->row;
  rs_ve = _rowsum->ve;
  i_end = Q->m;
  for (i = 0; i < i_end; i++, row++) {
    j_end = row->len;
    j_idx = row->diag + 1;
    if (j_idx <= 0) {
      // there must always be a diagonal entry
      m_error(E_INTERN, "Hqp_HL_Gerschgorin::update_Q");
    }
    elt = row->elt + j_idx;
    for (; j_idx < j_end; j_idx++, elt++) {
      j = elt->col;
      val = fabs(elt->val);
      rs_ve[i] += val;
      rs_ve[j] += val;
    }
  }

  row = Q->row;
  for (i = 0; i < i_end; i++, row++, rs_ve++) {
    elt = row->elt + row->diag;
    elt->val = max(elt->val, *rs_ve + _eps);
  }

/*
  // enshure positive definitness
  Q = prg->qp()->Q;
  i_end = Q->m;
  for (i=0; i<i_end; i++) {
    rowsum = 0.0;
    j_end = Q->row[i].len;
    elt = Q->row[i].elt;
    diag_elt = NULL;
    for (j_idx=0; j_idx<j_end; j_idx++, elt++) {
      if (elt->col != i)
	rowsum += fabs(elt->val);
      else
	diag_elt = elt;
    }
    assert(diag_elt);
    diag_elt->val = max(diag_elt->val, rowsum + _eps);
  }    
*/
}

//=========================================================================
