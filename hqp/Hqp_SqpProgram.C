/*
 * Hqp_SqpProgram.C -- class definition
 *
 * rf, 6/6/94
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

#include <assert.h>
#include <stdlib.h>

#include <If_Real.h>
#include <If_RealVec.h>
#include <If_Method.h>

#include "Hqp_SqpSolver.h"
#include "Hqp_SqpProgram.h"
#include "Hqp_Program.h"

// let currently created program be nodified by theSqpSolver
//----------------------------------------------------------
extern Hqp_SqpSolver *theSqpSolver;

typedef If_Method<Hqp_SqpProgram> If_Cmd;

IF_BASE_DEFINE(Hqp_SqpProgram);

#define GET_SET_CB(vartype, name) \
  "prg_"#name, \
  IF_GET_CB(vartype, Hqp_SqpProgram, name), \
  IF_SET_CB(vartype, Hqp_SqpProgram, set_##name)

//-------------------------------------------------------------------------
Hqp_SqpProgram::Hqp_SqpProgram()
{
  _qp = new Hqp_Program;
  _x = v_resize(v_get(1), 0);
  _f = 0.0;

  theSqpSolver->set_prg(this);

  _ifList.append(new If_Real(GET_SET_CB(Real, f)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, x)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, s)));
  _ifList.append(new If_Cmd("prg_setup", &Hqp_SqpProgram::setup, this));
  _ifList.append(new If_Cmd("prg_init_x", &Hqp_SqpProgram::init_x, this));
  _ifList.append(new If_Cmd("prg_qp_dump", &Hqp_SqpProgram::qp_dump, this));
  _ifList.append(new If_Cmd("prg_update_fbd",
			    &Hqp_SqpProgram::update_fbd, this));
  // note: need to cast as test method calls non-const member set_x 
  typedef double (Hqp_SqpProgram::*RealReadCb)() const;
  _ifList.append(new If_Real("prg_test",
			     new If_GetCb<double,Hqp_SqpProgram>
			     ((RealReadCb)&Hqp_SqpProgram::test, this)));
}

//-------------------------------------------------------------------------
Hqp_SqpProgram::~Hqp_SqpProgram()
{
  delete _qp;
  v_free(_x);
}

//-------------------------------------------------------------------------
void Hqp_SqpProgram::set_x(const VECP n_x)
{
  assert(n_x->dim == _x->dim);
  v_copy(n_x, _x);
}

//-------------------------------------------------------------------------
void Hqp_SqpProgram::reinit_bd()
{
  // do nothing per default
}

//-------------------------------------------------------------------------
// write difference between approximated and calculated gradients into qp
// -- _x, _f, and _qp must be initialized before
// -- results are stored in _qp->c, _qp->A, and _qp->C
//    (--> program should be reinitialzed afterwards)
// -- returns maximal found error
//
Real Hqp_SqpProgram::test()
{
  int i, j, idx;
  int xdim = _x->dim;
  int bdim = _qp->b->dim;
  int ddim = _qp->d->dim;
  Real xi_bak, dxi, f_bak;
  VECP dc, b_bak, d_bak;
  Real *ve;
  SPROW *row;
  Real val, maxval = 0.0;

  f_bak = _f;
  dc = v_get(xdim);
  b_bak = v_copy(_qp->b, VNULL);
  d_bak = v_copy(_qp->d, VNULL);
  for (i=0; i<xdim; i++) {
    xi_bak = _x[i];
    dxi = 1e-4 * fabs(xi_bak) + 1e-6;
    _x[i] += dxi;
    set_x(_x);	// write new _x with access method

    update_fbd();

    dc[i] = (_f - f_bak) / dxi;
    _x[i] = xi_bak;

    v_sub(_qp->b, b_bak, _qp->b);
    sv_mlt(1.0/dxi, _qp->b, _qp->b);
    ve = _qp->b->ve;
    row = _qp->A->row;
    for (j=0; j<bdim; j++, ve++, row++)
      if (*ve != 0.0) {
	idx = sprow_idx(row, i);
	if (idx >= 0)
	  val = row->elt[idx].val = *ve - row->elt[idx].val;
	else
	  val = sp_set_val(_qp->A, j, i, *ve);
	maxval = max(maxval, fabs(val));
      }
    
    v_sub(_qp->d, d_bak, _qp->d);
    sv_mlt(1.0/dxi, _qp->d, _qp->d);
    ve = _qp->d->ve;
    row = _qp->C->row;
    for (j=0; j<ddim; j++, ve++, row++)
      if (*ve != 0.0) {
	idx = sprow_idx(row, i);
	if (idx >= 0)
	  val = row->elt[idx].val = *ve - row->elt[idx].val;
	else
	  val = sp_set_val(_qp->C, j, i, *ve);
	maxval = max(maxval, fabs(val));
      }
  }
  _f = f_bak;

  v_sub(dc, _qp->c, _qp->c);
  maxval = max(maxval, v_max(_qp->c, &i));
  maxval = max(maxval, -v_min(_qp->c, &i));

  v_copy(b_bak, _qp->b);
  v_copy(d_bak, _qp->d);

  v_free(dc);
  v_free(b_bak);
  v_free(d_bak);

  return maxval;
}

//-------------------------------------------------------------------------
void Hqp_SqpProgram::qp_dump()
{
  const char *filename = "prg_qp_dump.out";
  FILE 	*fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    m_error(E_OVERWRITE, "Hqp_SqpProgram::qp_dump"
	    " that can't open file prg_qp_dump.out for writing");
  }
  _qp->foutput(fp);
  fclose(fp);
}


//=========================================================================
