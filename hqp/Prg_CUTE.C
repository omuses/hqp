/*
 * Prg_CUTE.C -- class definition
 *
 * rf, 10/27/95
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
#include <assert.h>
#include <string.h>

#include "Prg_CUTE.h"
#include "Hqp_Program.h"

#include <If_Int.h>
#include <If_Bool.h>
#include <If_Float.h>
#include <If_Method.h>
#include <If_Class.h>

extern "C" {
  csize_t csize_;
  cinit_t cinit_;
  cfn_t cfn_;
  csgr_t csgr_;
  csgrsh_t csgrsh_;
  csgreh_t csgreh_;
  cwrtsn_t cwrtsn_;
}

typedef If_Method<Prg_CUTE> If_Cmd;

IF_CLASS_DEFINE("CUTE", Prg_CUTE, Hqp_SqpProgram);

//-------------------------------------------------------------------------
Prg_CUTE::Prg_CUTE()
{
  _csize_p = NULL;
  _cinit_p = NULL;
  _cfn_p = NULL;
  _csgr_p = NULL;
  _csgrsh_p = NULL;
  _csgreh_p = NULL;
  _cwrtsn_p = NULL;

  _fscale = 1.0;
  _x0 = VNULL;
  _var_lb = VNULL;
  _var_ub = VNULL;
  _cns_lb = VNULL;
  _cns_ub = VNULL;
  _cns = VNULL;
  _J = VNULL;
  _H = VNULL;
  _v = VNULL;
  _var_ridx = IVNULL;
  _cns_ridx = IVNULL;

  _fbd_evals = 0;

  _INDVAR = NULL;
  _INDFUN = NULL;
  _IRNSH = NULL;
  _ICNSH = NULL;

  _hela = false;
  _hela_init = true;
  _ifList.append(new If_Float("prg_fscale", &_fscale));
  _ifList.append(new If_Bool("prg_hela", &_hela));
  _ifList.append(new If_Bool("prg_hela_init", &_hela_init));
  _ifList.append(new If_Int("prg_fbd_evals", &_fbd_evals));
  _ifList.append(new If_Cmd("prg_write_soln", &Prg_CUTE::write_soln, this));
}

//-------------------------------------------------------------------------
Prg_CUTE::~Prg_CUTE()
{
  v_free(_x0);
  v_free(_var_lb);
  v_free(_var_ub);
  v_free(_cns_lb);
  v_free(_cns_ub);
  v_free(_cns);
  v_free(_J);
  v_free(_H);
  v_free(_v);
  iv_free(_var_ridx);
  iv_free(_cns_ridx);

  delete [] _INDVAR;
  delete [] _INDFUN;
  delete [] _IRNSH;
  delete [] _ICNSH;
}

//-------------------------------------------------------------------------
void Prg_CUTE::parse_bounds(const Real *l_ve, const Real *u_ve, int n,
			    int *ridx_ive, int *me, int *m)
{
  int i;

  for (i=0; i<n; i++) {
    if (l_ve[i] == u_ve[i])
      ridx_ive[i] = (*me) ++;
    else {
      ridx_ive[i] = -1;
      if (l_ve[i] > -_Inf)
	ridx_ive[i] = (*m) ++;
      if (u_ve[i] < _Inf) {
	if (ridx_ive[i] == -1)
	  ridx_ive[i] = (*m) ++;
	else
	  (*m) ++;
      }
    }
  }
}

//-------------------------------------------------------------------------
void Prg_CUTE::setup()
{
  int me, m;	// number of equality and inequality constraints
  int i, offs;
  int *ridx_ive;
  Real *l_ve, *u_ve;
  fbool *EQUATN, *LINEAR;
  fbool EFIRST = 0;
  fbool LFIRST = 0;
  fbool NVFRST = 0;

  /*
   * bind the problem
   */

  _csize_p = (csize_t*)&csize_;
  _cinit_p = (cinit_t*)&cinit_;
  _cfn_p = (cfn_t*)&cfn_;
  _csgr_p = (csgr_t*)&csgr_;
  _csgrsh_p = (csgrsh_t*)&csgrsh_;
  _csgreh_p = (csgreh_t*)&csgreh_;
  _cwrtsn_p = (cwrtsn_t*)&cwrtsn_;

  /*
   * get the problem size; allocate and initialize variables
   */

  (*_csize_p)(&_N, &_M, &_NNZJMAX, &_NNZHMAX);

  _x0 = v_resize(_x0, _N);
  _var_lb = v_resize(_var_lb, _N);
  _var_ub = v_resize(_var_ub, _N);
  _v = v_resize(_v, _M);
  _cns_lb = v_resize(_cns_lb, _M);
  _cns_ub = v_resize(_cns_ub, _M);
  EQUATN = new fbool [_M];
  LINEAR = new fbool [_M];

  (*_cinit_p)(&_N, &_M, _x0->ve, _var_lb->ve, _var_ub->ve, &_Inf,
	      EQUATN, LINEAR, _v->ve, _cns_lb->ve, _cns_ub->ve,
	      &EFIRST, &LFIRST, &NVFRST);

  _x0->dim = _N;
  _var_lb->dim = _N;
  _var_ub->dim = _N;
  _var_ridx = iv_resize(_var_ridx, _N);
  _v->dim = _M;
  _cns_lb->dim = _M;
  _cns_ub->dim = _M;
  _cns_ridx = iv_resize(_cns_ridx, _M);
  _cns = v_resize(_cns, _M);

  me = 0;
  m = 0;
  parse_bounds(_var_lb->ve, _var_ub->ve, _N, _var_ridx->ive, &me, &m);
  _me_lin = me;
  _m_lin = m;
  parse_bounds(_cns_lb->ve, _cns_ub->ve, _M, _cns_ridx->ive, &me, &m);

  delete [] EQUATN;
  delete [] LINEAR;

  /*
   * allocate and initialize _x and _qp, including sparsity pattern
   */

  delete [] _INDVAR;
  delete [] _INDFUN;
  delete [] _IRNSH;
  delete [] _ICNSH;

  _J = v_resize(_J, _NNZJMAX);
  _INDVAR = new fint [_NNZJMAX];
  _INDFUN = new fint [_NNZJMAX];
  _H = v_resize(_H, _NNZHMAX);
  _IRNSH = new fint [_NNZHMAX];
  _ICNSH = new fint [_NNZHMAX];

  _qp->resize(_N, me, m);
  _x = v_resize(_x, _N);

  // variable bounds

  l_ve = _var_lb->ve;
  u_ve = _var_ub->ve;
  ridx_ive = _var_ridx->ive;
  for (i=0; i<_N; i++) {
    if (l_ve[i] == u_ve[i])
      sp_set_val(_qp->A, ridx_ive[i], i, 1.0);
    else {
      offs = 0;
      if (l_ve[i] > -_Inf) {
	sp_set_val(_qp->C, ridx_ive[i], i, 1.0);
	offs = 1;
      }
      if (u_ve[i] < _Inf)
	sp_set_val(_qp->C, ridx_ive[i] + offs, i, -1.0);
    }
  }

  /*
   * The structure resulting from general constraints can't be
   * initialized here, as there exists no appropriate procedure
   * call in the CUTE package.
   */

  _NNZJ = _NNZH = 0;

#if 0
  /*
   * evaluate csgreh_ and setup struct for Jacobians and Lagrangian Hessian
   */

  fbool GRLAGF = 0;
  fbool BYROWS = 1;
  fint NHI;
  fint *IRNHI = _IRNSH;
  fint LIRNHI = _NNZHMAX;
  fint LE = _NNZHMAX / 2;
  fint *IPRNHI = _ICNSH;
  fint *IPRHI = _ICNSH + LE;
  int j, k, idx, idx_end;
  int iidx, jidx, ijidx_end;
  
  sv_mlt(1.0/_fscale, _v, _v);
  if (_hela_init || _hela) {
    (*_csgreh_p)(&_N, &_M, _x0->ve, &GRLAGF, &_M, _v->ve,
		 &_NNZJ, &_NNZJMAX, _J->ve, _INDVAR, _INDFUN,
		 &NHI, IRNHI, &LIRNHI,
		 &LE, IPRNHI, _H->ve, &_NNZHMAX, IPRHI, &BYROWS);
  }
  else {
    (*_csgr_p)(&_N, &_M, &GRLAGF, &_M, _v->ve, _x0->ve,
	       &_NNZJ, &_NNZJMAX, _J->ve, _INDVAR, _INDFUN);
  }
  sv_mlt(_fscale, _v, _v);

  // general constraints

  v_zero(_qp->c);
  l_ve = _cns_lb->ve;
  u_ve = _cns_ub->ve;
  ridx_ive = _cns_ridx->ive;
  for (idx = 0; idx < _NNZJ; idx++) {
    i = _INDFUN[idx] - 1;
    j = _INDVAR[idx] - 1;
    if (i == -1)
      _qp->c[j] = _fscale * J_ve[idx];
    if (i >= 0) {
      if (l_ve[i] == u_ve[i])
	sp_set_val(_qp->A, ridx_ive[i], j, 0.0);
      else {
	offs = 0;
	if (l_ve[i] > -_Inf) {
	  sp_set_val(_qp->C, ridx_ive[i], j, 0.0);
	  offs = 1;
	}
	if (u_ve[i] < _Inf)
	  sp_set_val(_qp->C, ridx_ive[i] + offs, j, 0.0);
      }
    }
  }

  // Hessian

  if (_hela_init) {
    for (k = 0; k < NHI; k++) {
      ijidx_end = IPRNHI[k + 1] - 1;
      idx_end = IPRHI[k + 1] - 1;
      idx = IPRHI[k] - 1;
      for (iidx = IPRNHI[k] - 1; iidx < ijidx_end; iidx++) {
	i = IRNHI[iidx] - 1;
	for (jidx = IPRNHI[k] - 1; jidx < ijidx_end; jidx++) {
	  j = IRNHI[jidx] - 1;
	  if (i <= j) {
	    sp_add_val(_qp->Q, i, j, _H->ve[idx]);
	    idx++;
	  }
	}
      }
    }
  }

#else

  fbool GRLAGF = 0;
  Real *J_ve;
  int j, idx;

  sv_mlt(1.0/_fscale, _v, _v);
  if (_hela_init || _hela) {
    (*_csgrsh_p)(&_N, &_M, _x0->ve, &GRLAGF, &_M, _v->ve,
		 &_NNZJ, &_NNZJMAX, _J->ve, _INDVAR, _INDFUN,
		 &_NNZH, &_NNZHMAX, _H->ve, _IRNSH, _ICNSH);
  }
  else {
    (*_csgr_p)(&_N, &_M, &GRLAGF, &_M, _v->ve, _x0->ve,
	       &_NNZJ, &_NNZJMAX, _J->ve, _INDVAR, _INDFUN);
  }
  sv_mlt(_fscale, _v, _v);

  v_zero(_qp->c);
  l_ve = _cns_lb->ve;
  u_ve = _cns_ub->ve;
  ridx_ive = _cns_ridx->ive;
  J_ve = _J->ve;
  for (idx = 0; idx < _NNZJ; idx++) {
    i = _INDFUN[idx] - 1;
    j = _INDVAR[idx] - 1;
    if (i == -1)
      _qp->c[j] = _fscale * J_ve[idx];
    if (i >= 0) {
      if (l_ve[i] == u_ve[i])
	sp_set_val(_qp->A, ridx_ive[i], j, J_ve[idx]);
      else {
	offs = 0;
	if (l_ve[i] > -_Inf) {
	  sp_set_val(_qp->C, ridx_ive[i], j, J_ve[idx]);
	  offs = 1;
	}
	if (u_ve[i] < _Inf)
	  sp_set_val(_qp->C, ridx_ive[i] + offs, j, -J_ve[idx]);
      }
    }
  }

  if (_hela_init) {
    for (idx = 0; idx < _NNZH; idx++) {
      i = _IRNSH[idx] - 1;
      j = _ICNSH[idx] - 1;
      if (j >= i)
	sp_set_val(_qp->Q, i, j, _H->ve[idx]);
      else
	sp_set_val(_qp->Q, j, i, _H->ve[idx]);
    }
  }
#endif

  init_x();

  _fbd_evals = 1;
}

//-------------------------------------------------------------------------
void Prg_CUTE::init_x()
{
  v_copy(_x0, _x);
}

//-------------------------------------------------------------------------
void Prg_CUTE::update_bounds(const Real *val_ve,
			     const Real *l_ve, const Real *u_ve, int n,
			     const int *ridx_ive)
{
  int i, offs;
  Real *b_ve, *d_ve;

  b_ve = _qp->b->ve;
  d_ve = _qp->d->ve;

  for (i=0; i<n; i++) {
    if (l_ve[i] == u_ve[i])
      b_ve[ridx_ive[i]] = val_ve[i] - l_ve[i];
    else {
      offs = 0;
      if (l_ve[i] > -_Inf) {
	d_ve[ridx_ive[i]] = val_ve[i] - l_ve[i];
	offs = 1;
      }
      if (u_ve[i] < _Inf)
	d_ve[ridx_ive[i] + offs] = u_ve[i] - val_ve[i];
    }
  }
}

//-------------------------------------------------------------------------
void Prg_CUTE::update_fbd()
{
  if (!_cfn_p)
    error(E_NULL, "Prg_CUTE::update_fbd");
  
  (*_cfn_p)(&_N, &_M, _x->ve, &_f, &_M, _cns->ve);

  _f *= _fscale;
  update_bounds(_x->ve, _var_lb->ve, _var_ub->ve, _N, _var_ridx->ive);
  update_bounds(_cns->ve, _cns_lb->ve, _cns_ub->ve, _M, _cns_ridx->ive);

  _fbd_evals++;
}

//-------------------------------------------------------------------------
void Prg_CUTE::update(const VECP y, const VECP z)
{
  freal F;
  fbool GRLAGF = 0;
  int i, i_end, j, idx, offs;
  int *ridx_ive;
  Real *l_ve, *u_ve, *v_ve, *J_ve;
  SPROW *row;
  int newel_h, newel_j;

  update_fbd();
  _fbd_evals--;

  if (!(const VEC *)y || !(const VEC *)z)
    error(E_NULL, "Prg_CUTE::update");

  if (!_csgr_p || !_csgrsh_p)
    error(E_NULL, "Prg_CUTE::update");

  ridx_ive = _cns_ridx->ive;
  l_ve = _cns_lb->ve;
  u_ve = _cns_ub->ve;
  v_ve = _v->ve;
  for (i=0; i<_M; i++) {
    idx = ridx_ive[i];
    if (l_ve[i] == u_ve[i])
      v_ve[i] = y->ve[idx];
    else {
      offs = 0;
      v_ve[i] = 0.0;
      if (l_ve[i] > -_Inf) {
	v_ve[i] -= z->ve[idx];
	offs = 1;
      }
      if (u_ve[i] < _Inf)
	v_ve[i] += z->ve[idx + offs];
    }
  }

  sv_mlt(1.0/_fscale, _v, _v);
  if (_hela) {
    (*_csgrsh_p)(&_N, &_M, _x->ve, &GRLAGF, &_M, _v->ve,
		 &_NNZJ, &_NNZJMAX, _J->ve, _INDVAR, _INDFUN,
		 &_NNZH, &_NNZHMAX, _H->ve, _IRNSH, _ICNSH);
    sv_mlt(_fscale, _H, _H);
  }
  else {
    (*_csgr_p)(&_N, &_M, &GRLAGF, &_M, _v->ve, _x->ve,
	       &_NNZJ, &_NNZJMAX, _J->ve, _INDVAR, _INDFUN);
  }
  sv_mlt(_fscale, _v, _v);

  row = _qp->A->row + _me_lin;
  i_end = _qp->A->m;
  for (i = _me_lin; i < i_end; i++, row++)
    sprow_zero(row);

  row = _qp->C->row + _m_lin;
  i_end = _qp->C->m;
  for (i = _m_lin; i < i_end; i++, row++)
    sprow_zero(row);

  newel_j = 0;

  Real *c_ve;
  v_zero(_qp->c);
  c_ve = _qp->c->ve;
  l_ve = _cns_lb->ve;
  u_ve = _cns_ub->ve;
  ridx_ive = _cns_ridx->ive;
  J_ve = _J->ve;
  for (idx = 0; idx < _NNZJ; idx++) {
    i = _INDFUN[idx] - 1;
    j = _INDVAR[idx] - 1;
    if (i == -1)
      c_ve[j] = _fscale * J_ve[idx];
    if (i >= 0) {
      if (l_ve[i] == u_ve[i])
	newel_j |= sp_update_val(_qp->A, ridx_ive[i], j, J_ve[idx]);
      else {
	offs = 0;
	if (l_ve[i] > -_Inf) {
	  newel_j |= sp_update_val(_qp->C, ridx_ive[i], j, J_ve[idx]);
	  offs = 1;
	}
	if (u_ve[i] < _Inf)
	  newel_j |= sp_update_val(_qp->C, ridx_ive[i] + offs, j, -J_ve[idx]);
      }
    }
  }

  newel_h = 0;
  if (_hela) {
    sp_zero(_qp->Q);
    for (idx = 0; idx < _NNZH; idx++) {
      i = _IRNSH[idx] - 1;
      j = _ICNSH[idx] - 1;
      if (j >= i)
	newel_h |= sp_update_val(_qp->Q, i, j, _H->ve[idx]);
      else
	newel_h |= sp_update_val(_qp->Q, j, i, _H->ve[idx]);
    }
  }

  if (newel_j)
    warning(WARN_UNKNOWN,
	    "Prg_CUTE::update: Jacobian sparsity structure changed");
  if (newel_h)
    warning(WARN_UNKNOWN,
	    "Prg_CUTE::update: Hessian sparsity structure changed");

  extern Tcl_Interp *theInterp;
  if (newel_j || newel_h) {
    init_x();
    Tcl_Eval(theInterp, "sqp_init; sqp_qp_update");
  }
}

//-------------------------------------------------------------------------
int Prg_CUTE::write_soln(int argc, char *argv[], char **)
{
  char header[61];
  fint idum = 0;
  freal rdum = 0.0;
  extern const char *Hqp_Version;

  if (!_cwrtsn_p)
    error(E_NULL, "Prg_CUTE::write_soln");

  sprintf(header, "HQP %-4.4s: ", Hqp_Version);
  if (argc > 1) 
    sprintf(header + 10, "%-50.50s", argv[argc - 1]);
  else
    sprintf(header + 10, "%-50.50s", "");

  if (argc > 1 && strncmp(argv[1], "-nosol", 6) == 0)
    (*_cwrtsn_p)(&idum, &idum, header, &rdum, &rdum, &rdum);
  else
    (*_cwrtsn_p)(&_N, &_M, header, &_f, _x->ve, _v->ve);

  return IF_OK;
}


//=========================================================================
