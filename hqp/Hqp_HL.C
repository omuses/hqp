/*
 * Hqp_HL.C -- 
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

#include <math.h>

#include <If_Bool.h>
#include <If_Float.h>

#include "Hqp_HL.h"
#include "Hqp_Program.h"
#include "Hqp_SqpProgram.h"

IF_BASE_DEFINE(Hqp_HL);

//--------------------------------------------------------------------------
Hqp_HL::Hqp_HL()
{
  _scale = true;
  _eps = 1e-8;
  _logging = false;

  _rowsum = v_resize(v_get(1), 0);

  _ifList.append(new If_Bool("sqp_hela_scale", &_scale));
  _ifList.append(new If_Bool("sqp_hela_logging", &_logging));
  _ifList.append(new If_Float("sqp_hela_eps", &_eps));
}

//--------------------------------------------------------------------------
Hqp_HL::~Hqp_HL()
{
  v_free(_rowsum);
}

//--------------------------------------------------------------------------
static double dxi(double xi)
{
  return fabs(1e-4 * xi) + 1e-6;
}

//--------------------------------------------------------------------------
void Hqp_HL::init(const VEC *y, const VEC *z, Hqp_SqpProgram *prg)
{
  if (sp_norm_inf(prg->qp()->Q) > _eps) {
    // assume a user-specified Lagrangian Hessian
    posdef(prg);
    return;
  }

  Hqp_Program 	*qp = prg->qp();
  Real		f_bak;
  VEC	      	*c_bak=VNULL, *x_bak=VNULL;
  VEC	      	*b_bak=VNULL, *d_bak=VNULL;
  SPMAT		*A_bak, *C_bak;
  VEC		*dx, *x_mod, *dgL, *gL;
  int		i, n = qp->c->dim;
  Real		val;

  if (!_scale) {
    sp_ident(prg->qp()->Q);
  }
  else {
    f_bak = prg->f();
    x_bak = v_copy(prg->x(), x_bak);
    b_bak = v_copy(qp->b, b_bak);
    d_bak = v_copy(qp->d, d_bak);
    c_bak = v_copy(qp->c, c_bak);
    A_bak = qp->A;
    qp->A = sp_copy(A_bak);
    C_bak = qp->C;
    qp->C = sp_copy(C_bak);

    dx = v_get(n);
    x_mod = v_get(n);
    dgL = v_get(n);
    gL = v_get(n);

    gL = grd_L(y, z, qp, gL);

    dx = v_map(&dxi, prg->x(), dx);
    v_add(x_bak, dx, x_mod);
    prg->x(x_mod);
    prg->update((VEC *)y, (VEC *)z);
    dgL = grd_L(y, z, qp, dgL);
    dgL = v_sub(dgL, gL, dgL);

    sp_zero(qp->Q);
    for (i=0; i<n; i++) {
      val = dgL->ve[i] / dx->ve[i];
      sp_set_val(qp->Q, i, i, max(val, _eps));
    }

    prg->f(f_bak);
    prg->x(x_bak);
    v_copy(b_bak, qp->b);
    v_copy(d_bak, qp->d);
    v_copy(c_bak, qp->c);

    sp_free(qp->A);
    qp->A = A_bak;
    sp_free(qp->C);
    qp->C = C_bak;

    v_free(b_bak);
    v_free(d_bak);
    v_free(c_bak);
    v_free(x_bak);
    v_free(dgL);
    v_free(x_mod);
    v_free(dx);
    v_free(gL);
  }
}

//--------------------------------------------------------------------------
void Hqp_HL::posdef(Hqp_SqpProgram *prg)
{
  SPMAT *Q;
  SPROW *row;
  row_elt *elt;
  Real *rs_ve, val;
  int i, i_end;
  int j, j_idx, j_end;

  // ensure positive definiteness with diagonal offsets
  // (works for upper diagonal matrices)
  Q = prg->qp()->Q;

  if (!Q->flag_diag)
    sp_diag_access(Q);

  v_resize(_rowsum, Q->m);
  v_zero(_rowsum);
  row = Q->row;
  rs_ve = _rowsum->ve;
  i_end = Q->m;
  for (i = 0; i < i_end; i++, row++) {
    j_end = row->len;
    j_idx = row->diag + 1;
    if (j_idx <= 0) {
      // there must always be a diagonal entry
      error(E_INTERN, "Hqp_HL::posdef");
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
  // ensure positive definiteness
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

//--------------------------------------------------------------------------
VEC *Hqp_HL::grd_L(const VEC *y, const VEC *z, const Hqp_Program *qp,
		   VEC *out)
{
  int n = qp->Q->m;

  if (!out || (int)out->dim != n) {
    out = v_resize(out, n);
  }

  // calculate gradient of Lagrangian
  v_zero(out);
  sp_vm_mltadd(out, y, qp->A, 1.0, out);
  sp_vm_mltadd(out, z, qp->C, 1.0, out);
  v_sub(qp->c, out, out);

  return out;
}


//==========================================================================
