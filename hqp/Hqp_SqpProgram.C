/*
 * Hqp_SqpProgram.C -- class definition
 *
 * rf, 6/6/94
 */

/*
    Copyright (C) 1994--2001  Ruediger Franke

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
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

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

//-------------------------------------------------------------------------
Hqp_SqpProgram::Hqp_SqpProgram()
{
  _qp = new Hqp_Program;
  _x = v_resize(v_get(1), 0);
  _f = 0.0;

  theSqpSolver->set_prg(this);

  _ifList.append(new If_Real("prg_f", &_f));
  _ifList.append(new If_RealVec("prg_x", &_x));
  _ifList.append(new If_Cmd("prg_setup", &Hqp_SqpProgram::setup_cmd, this));
  _ifList.append(new If_Cmd("prg_init_x", &Hqp_SqpProgram::init_x_cmd, this));
  _ifList.append(new If_Cmd("prg_test", &Hqp_SqpProgram::test_cmd, this));
  _ifList.append(new If_Cmd("prg_qp_dump",
			    &Hqp_SqpProgram::qp_dump_cmd, this));
  _ifList.append(new If_Cmd("prg_update_fbd",
			    &Hqp_SqpProgram::update_fbd_cmd, this));
}

//-------------------------------------------------------------------------
Hqp_SqpProgram::~Hqp_SqpProgram()
{
  delete _qp;
  v_free(_x);
}

//-------------------------------------------------------------------------
void Hqp_SqpProgram::x(const VECP n_x)
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
int Hqp_SqpProgram::setup_cmd(IF_CMD_ARGS)
{
  setup();

  return IF_OK;
}

//-------------------------------------------------------------------------
int Hqp_SqpProgram::init_x_cmd(IF_CMD_ARGS)
{
  init_x();

  return IF_OK;
}

//-------------------------------------------------------------------------
// write difference between approximated and calculated gradients into qp
// -- _x, _f, and _qp must be initialized before
// -- results are stored in _qp->c, _qp->A, and _qp->C
//    (--> program should be reinitialzed afterwards)
// -- returns maximal found error
//
int Hqp_SqpProgram::test_cmd(int, char *[], char **result)
{
  static char maxstr[15];
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
    x(_x);	// write new _x with access method

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

  if (result) {
    sprintf(maxstr, "%.5g", maxval);
    *result = maxstr;
  }

  v_copy(b_bak, _qp->b);
  v_copy(d_bak, _qp->d);

  v_free(dc);
  v_free(b_bak);
  v_free(d_bak);

  return IF_OK;
}

//-------------------------------------------------------------------------
int Hqp_SqpProgram::qp_dump_cmd(int argc, char *argv[], char **result)
{
  int	fd;
  char	*pos;
  FILE 	*fp;

  if (argc != 2) {
    *result = "wrong # args, should be: prg_qp_dump [<fd>|<filename>]";
    return IF_ERROR;
  }

  // try to interpret argv[1] as filedescriptor
  // if this fails, then create file with name argv[1]

  fd = (int)strtol(argv[1], &pos, 0);
  if (pos != argv[1]) {
    fp = fdopen(fd, "w");
    if (fp == NULL) {
      *result = strerror(errno);
      return IF_ERROR;
    }
    _qp->foutput(fp);
  }
  else {
    fp = fopen(argv[1], "w");
    if (fp == NULL) {
      *result = "can't open file for writing";
      return IF_ERROR;
    }
    _qp->foutput(fp);
    fclose(fp);
  }

  return IF_OK;
}

//-------------------------------------------------------------------------
int Hqp_SqpProgram::update_fbd_cmd(int argc, char *argv[], char **result)
{
  update_fbd();
  return IF_OK;
}


//=========================================================================
