/*
 * Prg_CUTE_ST.C -- class definition
 *
 * rf, 2/10/97
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

//#define RECALCULATE_CNS 1
//#define BLOW_OUTPUT 1

#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "Prg_CUTE_ST.h"
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
  cscifg_t cscifg_;
  csgreh_t csgreh_;
  cwrtsn_t cwrtsn_;
}

typedef If_Method<Prg_CUTE_ST> If_Cmd;

IF_CLASS_DEFINE("CUTE_ST", Prg_CUTE_ST, Hqp_SqpProgram);

//-------------------------------------------------------------------------
static int cmp_int(const int *i1, const int *i2)
{
  return *i1 - *i2;
}

//-------------------------------------------------------------------------
/*
 * sp_add_val:
 *  add val to an existing entry
 */
int sp_add_val(SPMAT *A, int i, int j, Real val)
{
  SPROW	*row;
  int	idx;
   
  if ( A == SMNULL )
    error(E_NULL,"sp_add_val");
  if ( i < 0 || i >= A->m || j < 0 || j >= A->n )
    error(E_SIZES,"sp_add_val");
   
  row = A->row+i;
  idx = sprow_idx(row, j);
  if ( idx >= 0 ) {
    row->elt[idx].val += val;
    return 0;
  }
  else {
    sp_set_val(A, i, j, val);
    return 1;
  }
}

//-------------------------------------------------------------------------
Prg_CUTE_ST::Prg_CUTE_ST()
{
  _csize_p = NULL;
  _cinit_p = NULL;
  _cfn_p = NULL;
  _csgr_p = NULL;
  _cscifg_p = NULL;
  _csgreh_p = NULL;
  _cwrtsn_p = NULL;

  _xscale = 1.0;
  _fscale = 1.0;
  _x0 = v_resize(v_get(1), 0);
  _xs = v_resize(v_get(1), 0);
  _var_lb = v_resize(v_get(1), 0);
  _var_ub = v_resize(v_get(1), 0);
  _cns_lb = v_resize(v_get(1), 0);
  _cns_ub = v_resize(v_get(1), 0);
  _cns = v_resize(v_get(1), 0);
  _J = v_resize(v_get(1), 0);
  _H = v_resize(v_get(1), 0);
  _v = v_resize(v_get(1), 0);
  _var_ridx = iv_resize(iv_get(1), 0);
  _cns_ridx = iv_resize(iv_get(1), 0);
  _h2s = iv_resize(iv_get(1), 0);
  _s2h0 = iv_resize(iv_get(1), 0);
  _s2h = iv_resize(iv_get(1), 0);
  _s2h_stage = iv_resize(iv_get(1), 0);
  _s2h_begin = iv_resize(iv_get(1), 0);
  _equ_idx0 = iv_resize(iv_get(1), 0);
  _equ_idx1 = iv_resize(iv_get(1), 0);
  _grp_stage = iv_resize(iv_get(1), 0);
  _cns_stage = iv_resize(iv_get(1), 0);

  _fbd_evals = 0;
  _nstages = 0;
  _NNZJ = 0;
  _NNZJMAX = 0;
  _INDVAR = NULL;
  _INDFUN = NULL;

  _NHI = 0;
  _NEL = 0;
  _NHIMAX = 0;
  _NELMAX = 0;
  _IPRHI = NULL;
  _IRNHI = NULL;
  _IPRNHI = NULL;
  _iprhi = iv_resize(iv_get(1), 0);
  _irnhi = iv_resize(iv_get(1), 0);
  _iprnhi = iv_resize(iv_get(1), 0);

  _hela = false;
  _hela_init = true;
  _stretch = true;
  _ifList.append(new If_Float("prg_xscale", &_xscale));
  _ifList.append(new If_Float("prg_fscale", &_fscale));
  _ifList.append(new If_Bool("prg_hela", &_hela));
  _ifList.append(new If_Bool("prg_hela_init", &_hela_init));
  _ifList.append(new If_Bool("prg_stretch", &_stretch));
  _ifList.append(new If_Int("prg_fbd_evals", &_fbd_evals));
  _ifList.append(new If_Int("prg_N", &_N));
  _ifList.append(new If_Int("prg_NEL", &_NEL));
  _ifList.append(new If_Int("prg_K", &_nstages));
  _ifList.append(new If_Cmd("prg_write_soln", &Prg_CUTE_ST::write_soln, this));
}

//-------------------------------------------------------------------------
Prg_CUTE_ST::~Prg_CUTE_ST()
{
  v_free(_x0);
  v_free(_xs);
  v_free(_var_lb);
  v_free(_var_ub);
  v_free(_cns_lb);
  v_free(_cns_ub);
  v_free(_cns);
  v_free(_J);
  v_free(_H);
  v_free(_v);
  iv_free(_grp_stage);
  iv_free(_cns_stage);
  iv_free(_var_ridx);
  iv_free(_cns_ridx);
  iv_free(_h2s);
  iv_free(_s2h0);
  iv_free(_s2h);
  iv_free(_s2h_stage);
  iv_free(_s2h_begin);
  iv_free(_equ_idx0);
  iv_free(_equ_idx1);
  iv_free(_iprnhi);
  iv_free(_irnhi);
  iv_free(_iprhi);

  delete [] _INDVAR;
  delete [] _INDFUN;
  delete [] _IPRHI;
  delete [] _IRNHI;
  delete [] _IPRNHI;
}

//-------------------------------------------------------------------------
void Prg_CUTE_ST::parse_bounds(const VECP lb, const VECP ub, int n,
			      IVECP ridx, int &me, int &m)
{
  int i;

  for (i=0; i<n; i++) {
    if (lb[i] == ub[i])
      ridx[i] = me++;
    else {
      ridx[i] = -1;
      if (lb[i] > -_Inf)
	ridx[i] = m++;
      if (ub[i] < _Inf) {
	if (ridx[i] == -1)
	  ridx[i] = m++;
	else
	  m++;
      }
    }
  }
}

//-------------------------------------------------------------------------
static int shared_vars(int *v1, int n1, int *v2, int n2)
{
  assert(n1 <= n2);

  int i1, i2;
  int shared = 0;

  for (i1 = 0, i2 = 0; i1 < n1 && i2 < n2; i1++) {
    while (v1[i1] > v2[i2] && i2 < n2) i2++;
    if (i2 < n2 && v1[i1] == v2[i2]) {
      shared++;
    }
  }
  return shared;
}

//-------------------------------------------------------------------------
static bool is_subset(int *v1, int n1, int *v2, int n2)
{
  return (shared_vars(v1, n1, v2, n2) == n1);
}

//-------------------------------------------------------------------------
static bool is_resembling(int *v1, int n1, int *v2, int n2)
{
  int sv = shared_vars(v1, n1, v2, n2);
  return (sv == n1 || sv >= n1 - 1 && n1 > 1);
}

//-------------------------------------------------------------------------
void Prg_CUTE_ST::setup()
{
  int n, me, m;	// number of variables, equality, and inequality constr.
  int i, offs;
  fbool *EQUATN, *LINEAR;
  fbool EFIRST = 0;
  fbool LFIRST = 0;
  fbool NVFRST = 0;

  /*
   * bind CUTE procedures
   */

  _csize_p = (csize_t*)&csize_;
  _cinit_p = (cinit_t*)&cinit_;
  _cfn_p = (cfn_t*)&cfn_;
  _csgr_p = (csgr_t*)&csgr_;
#ifdef RECALCULATE_CNS
  _cscifg_p = (cscifg_t*)&cscifg_;
#endif
  _csgreh_p = (csgreh_t*)&csgreh_;
  _cwrtsn_p = (cwrtsn_t*)&cwrtsn_;

  /*
   * get the problem size; allocate and initialize variables
   */

  (*_csize_p)(&_N, &_M, &_NNZJMAX, &_NHIMAX);

  _NELMAX = _NHIMAX / 2;	// number of Hessian elements (blocks)
  
  // allocate maximal dimensions
  delete [] _INDVAR;
  delete [] _INDFUN;
  delete [] _IPRHI;
  delete [] _IRNHI;
  delete [] _IPRNHI;

  v_resize(_J, _NNZJMAX);
  _INDVAR = new fint [_NNZJMAX];
  _INDFUN = new fint [_NNZJMAX];
  v_resize(_H, _NHIMAX);
  _IPRHI = new fint [_NELMAX];
  _IRNHI = new fint [_NHIMAX];
  _IPRNHI = new fint [_NELMAX];

  v_resize(_x0, _N);
  v_resize(_var_lb, _N);
  v_resize(_var_ub, _N);
  v_resize(_v, _M);
  v_resize(_cns_lb, _M);
  v_resize(_cns_ub, _M);
  EQUATN = new fbool [_M];
  LINEAR = new fbool [_M];

  (*_cinit_p)(&_N, &_M, _x0->ve, _var_lb->ve, _var_ub->ve, &_Inf,
	      EQUATN, LINEAR, _v->ve, _cns_lb->ve, _cns_ub->ve,
	      &EFIRST, &LFIRST, &NVFRST);

  // resize to actual dimensions
  _x0->dim = _N;
  _var_lb->dim = _N;
  _var_ub->dim = _N;
  _v->dim = _M;
  _cns_lb->dim = _M;
  _cns_ub->dim = _M;
  v_resize(_xs, _N);
  iv_resize(_var_ridx, _N);
  iv_resize(_cns_ridx, _M);
  v_resize(_cns, _M);

  sv_mlt(_xscale, _var_lb, _var_lb);
  sv_mlt(_xscale, _var_ub, _var_ub);
  if (_xscale > 1.0)
    _Inf *= _xscale;

  me = 0;
  m = 0;
  parse_bounds(_var_lb, _var_ub, _N, _var_ridx, me, m);
  _me_bnd = me;
  _m_lin = m;
  parse_bounds(_cns_lb, _cns_ub, _M, _cns_ridx, me, m);

  /*
   * evaluate csgreh_ to obtain Jacobian and Hessian structure
   */

  fbool GRLAGF = 0;
  fbool BYROWS = 1;
  int j, k, idx;
  int iidx, jidx, ijidx_end;
  int *ip, *jp;
  VECP qp_c;

  (*_csgreh_p)(&_N, &_M, _x0->ve, &GRLAGF, &_M, _v->ve,
	       &_NNZJ, &_NNZJMAX, _J->ve, _INDVAR, _INDFUN,
	       &_NEL, _IRNHI, &_NHIMAX,
	       &_NELMAX, _IPRNHI, _H->ve, &_NHIMAX, _IPRHI, &BYROWS);

  _J->dim = _NNZJ;
  _H->dim = _IPRHI[_NEL] - 1;
  n = _IPRNHI[_NEL] - 1;

  sv_mlt(1.0/_xscale, _J, _J);
  sv_mlt(1.0/_xscale/_xscale, _H, _H);

  /*
   * assign nonlinear groups to stages
   */

#ifdef DEBUG
  fprintf(stderr, "Sorting indices,\n");
  fflush(stderr);
#endif

  // initialize sorted pointers and indizes for finite-element Hessian
  iv_resize(_iprhi, _NEL);
  iv_resize(_iprnhi, _NEL + 1);
  iv_resize(_irnhi, n);
  iv_resize(_grp_stage, _NEL);

  // sort the groups with decreasing size
  // (use _grp_stage as storage and _iprhi for the size of each group)
  for (k = 0; k < _NEL; k++) {
    _iprhi[k] = _IPRNHI[k + 1] - _IPRNHI[k];
    _grp_stage[k] = k;
  }

  k = 1;
  while (k < _NEL) {
    if (_iprhi[k] > _iprhi[k - 1]) {
      idx = _iprhi[k];
      _iprhi[k] = _iprhi[k-1];
      _iprhi[k-1] = idx;
      idx = _grp_stage[k];
      _grp_stage[k] = _grp_stage[k-1];
      _grp_stage[k-1] = idx;
      k -= (k > 1)? 1: 0;
    }
    else
      k++;
  }

  idx = 0;
  for (k = 0; k < _NEL; k++) {

    // copy indizes
    _iprnhi[k] = idx;
    ijidx_end = _IPRNHI[_grp_stage[k] + 1] - 1;
    for (iidx = _IPRNHI[_grp_stage[k]] - 1; iidx < ijidx_end; iidx++)
      _irnhi[idx++] = _IRNHI[iidx] - 1;

    // store sorted _IPRHI in _iprhi
    _iprhi[k] = _IPRHI[_grp_stage[k]] - 1;
  }
  _iprnhi[_NEL] = idx;

  // sort the indizes of each group as required by shared_vars()

  for (k = 0; k < _NEL; k++) {
    qsort((void *)(_irnhi->ive + _iprnhi[k]), _iprnhi[k+1] - _iprnhi[k],
	  sizeof(int), (int (*)(const void *, const void *))cmp_int);
  }

#ifdef DEBUG
  fprintf(stderr, "Assigning groups to stages,\n");
  fflush(stderr);
#endif

  // assign the groups to stages

  int koffs, nk, l, loffs, nl;
  bool found;

  _nstages = 0;
  if (_stretch) {
    for (k = 0; k < _NEL; k++) {
      koffs = _iprnhi[k];
      nk = _iprnhi[k + 1] - _iprnhi[k];
      found = false;

      // check existing stages
      for (l = 0; l < k && !found; l++) {
	loffs = _iprnhi[l];
	nl = _iprnhi[l + 1] - _iprnhi[l];
	if (nl <= nk)
	  found = is_resembling(_irnhi->ive + loffs, nl,
				_irnhi->ive + koffs, nk);
	else
	  found = is_resembling(_irnhi->ive + koffs, nk,
				_irnhi->ive + loffs, nl);
	if (found)
	  _grp_stage[k] = _grp_stage[l];
      }
      if (!found) {
	// make up a new stage
	_grp_stage[k] = _nstages++;
      }
    }
  }
  else {
    _nstages = 1;
    iv_set(_grp_stage, 0);
  }

#ifdef DEBUG
  fprintf(stderr, "  Introduced %d stages for %d groups,\n", _nstages, _NEL);
  fflush(stderr);
#endif

#ifdef BLOW_OUTPUT
  // print out the stage assignments
  for (k = 0; k < _NEL; k++) {
    fprintf(stderr, "%3d(%3d):", k, _grp_stage[k]);
    ijidx_end = _iprnhi[k + 1];
    for (iidx = _iprnhi[k]; iidx < ijidx_end; iidx++)
      fprintf(stderr, "%4d", _irnhi[iidx]);
    fprintf(stderr, "\n");
  }
#endif

  IVECP var_stage = iv_get(_N);
  iv_resize(_s2h_begin, _N + 1);
  iv_resize(_cns_stage, _M);

#ifdef DEBUG
  fprintf(stderr, "Assigning constraints to stages,\n");
  fflush(stderr);
#endif

  // determine the stage of each constraint
  // linear constraints are 'stage' -1
  // mark groups of the objective

  IVECP indvar = iv_get(_NNZJ);
  IVECP obj_grp = iv_get(_NEL);
  iv_set(obj_grp, 1);
  for (idx = 0; idx < _NNZJ; idx++)
    indvar[idx] = _INDVAR[idx] - 1;
  iv_set(_cns_stage, -1);
  ijidx_end = 0;
  while (ijidx_end < _NNZJ) {
    idx = ijidx_end;
    j = _INDFUN[idx] - 1;
    while (ijidx_end < _NNZJ && _INDFUN[ijidx_end] - 1 == j)
      ijidx_end++;
    if (j < 0 || LINEAR[j]) {
      continue;
    }
    nl = ijidx_end - idx;
    qsort((void *)(indvar->ive + idx), nl, sizeof(int),
	  (int (*)(const void *, const void *))cmp_int);
    found = false;
    for (k = 0; k < _NEL && !found; k++) {
      nk = _iprnhi[k+1] - _iprnhi[k];
      if (nk <= nl)
	found = is_subset(_irnhi->ive + _iprnhi[k], nk,
			  indvar->ive + idx, nl);
      if (found) {
	_cns_stage[j] = _grp_stage[k];
	obj_grp[k] = 0;
      }
    }
  }

#ifdef BLOW_OUTPUT
  // print out constraint assignments
  ijidx_end = 0;
  while (ijidx_end < _NNZJ) {
    idx = ijidx_end;
    j = _INDFUN[idx] - 1;
    if (j == -1)
      fprintf(stderr, "%3d(-,-)%3d:", j, -1);
    else
      fprintf(stderr, "%3d(%d,%d)%3d:", j, EQUATN[j], LINEAR[j],
	      _cns_stage[j]);
    while (ijidx_end < _NNZJ && _INDFUN[ijidx_end] - 1 == j) {
      fprintf(stderr, "%4d", _INDVAR[ijidx_end] - 1);
      ijidx_end++;
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "obj_grp: ");
  iv_foutput(stderr, obj_grp);
#endif

  /*
   * allocate and initialize _x and _qp, including sparsity pattern
   */

#if 1
#ifdef DEBUG
  fprintf(stderr, "Ordering stages,\n");
  fflush(stderr);
#endif

  // order the stages
  PERM *order = px_get(_nstages);
#if 0
  SPMAT *SM = sp_get(_nstages, _nstages, 10);
  iv_set(var_stage, -1);
  for (k = 0; k < _nstages; k++) {
    sp_set_val(SM, k, k, 1.0);
    for (l = 0; l < _NEL; l++) {
      if (_grp_stage[l] != k)
	continue;
      ijidx_end = _iprnhi[l + 1];
      for (iidx = _iprnhi[l]; iidx < ijidx_end; iidx++) {
	i = _irnhi[iidx];
	if (var_stage[i] != k) {
	  if (var_stage[i] >= 0) {
	    sp_set_val(SM, min(k, var_stage[i]), max(k, var_stage[i]), 1.0);
	  }
	  var_stage[i] = k;
	}
      }
    }
  }
  sp_symrcm(SM, order);
  sp_free(SM);
#else
  for (i = 0; i < _nstages; i++)
    order->pe[i] = _nstages - i - 1;
#endif
  for (i = 0; i < _NEL; i++)
    _grp_stage[i] = order->pe[_grp_stage[i]];
  for (i = 0; i < _M; i++)
    if (_cns_stage[i] >= 0)
      _cns_stage[i] = order->pe[_cns_stage[i]];

#ifdef BLOW_OUTPUT
  fprintf(stderr, "order: ");
  px_foutput(stderr, order);
#endif

  px_free(order);
#endif

#ifdef DEBUG
  fprintf(stderr, "Allocating the QP problem,\n");
  fflush(stderr);
#endif

  // count the number of additional variables
  // store their numbers in _s2h_begin
  iv_set(_s2h_begin, 0);
  iv_set(var_stage, -1);
  n = _N;
  for (k = 0; k < _nstages; k++) {
    for (l = 0; l < _NEL; l++) {
      if (_grp_stage[l] != k)
	continue;
      ijidx_end = _iprnhi[l + 1];
      for (iidx = _iprnhi[l]; iidx < ijidx_end; iidx++) {
	i = _irnhi[iidx];
	if (var_stage[i] != k) {
	  if (var_stage[i] >= 0) {
	    n++;	// additional variable
	  }
	  _s2h_begin[i]++;
	  var_stage[i] = k;
	}
      }
    }
  }
  // sum up _s2h_begin
  // consider linear variables and additional bounds
  int m_offs;
  m_offs = 0;
  offs = 0;
  j = 0;
  for (i = 0; i < _N; i++) {
    nk = _s2h_begin[i];
    if (nk == 0) {
      nk = 1;
      _s2h_begin[i] = 1; // linear variable
      j++;
    }
    // additional variable bounds
    if (_var_lb[i] != _var_ub[i] && _var_ridx[i] >= 0) {
      if (_var_lb[i] > -_Inf)
	m_offs += nk - 1;
      if (_var_ub[i] < _Inf)
	m_offs += nk - 1;
    }
    offs += nk;
    _s2h_begin[i] = offs - nk;
  }
  _s2h_begin[_N] = offs;
  m += m_offs;
  _m_lin += m_offs;

#ifdef DEBUG
  fprintf(stderr, "  %d linear variables,\n", j);
  fflush(stderr);
#endif

  offs = n - _N; // number of additional variables and equ. constraints

  me += offs;
  _me_lin = _me_bnd + offs;

  // shift indices for general constraints
  for (i = 0; i < _M; i++) {
    if (_cns_lb[i] == _cns_ub[i])
      _cns_ridx[i] += offs;
    else
      _cns_ridx[i] += m_offs;
  }

  _qp->resize(n, me, m, 10, 10, 2);
  _x = v_resize(_x, n);

  iv_resize(_h2s, n);
  iv_resize(_equ_idx0, offs);
  iv_resize(_equ_idx1, offs);
  iv_resize(_s2h, n);
  iv_resize(_s2h_stage, n);
  iv_resize(_s2h0, _N);

#ifdef DEBUG
  fprintf(stderr, "Initializing additional variables,\n");
  fflush(stderr);
#endif

  iv_set(var_stage, -1);
  iv_set(_s2h0, -1);
  j = 0;
  offs = 0;
  for (k = 0; k < _nstages; k++) {
    for (l = 0; l < _NEL; l++) {
      if (_grp_stage[l] != k)
	continue;
      ijidx_end = _iprnhi[l + 1];
      for (iidx = _iprnhi[l]; iidx < ijidx_end; iidx++) {
	i = _irnhi[iidx];
	if (var_stage[i] != k) {
	  if (var_stage[i] >= 0) {
	    _equ_idx0[offs] = _s2h[_s2h_begin[i]-1];
	    _equ_idx1[offs] = j;
	    sp_set_val(_qp->A, _me_bnd + offs, _equ_idx0[offs], 1.0);
	    sp_set_val(_qp->A, _me_bnd + offs, _equ_idx1[offs], -1.0);
	    offs++;
	  }
	  _s2h[_s2h_begin[i]] = j;
	  _s2h_stage[_s2h_begin[i]] = k;
	  _h2s[j] = i;
	  if (obj_grp[l])
	    _s2h0[i] = j;
	  j++;
	  _s2h_begin[i]++;
	  var_stage[i] = k;
	}
	else {
	  if (obj_grp[l])
	    _s2h0[i] = _s2h[_s2h_begin[i]-1];
	}
      }
    }
  }

  // indices for variables that appear only linearly
  for (i = 0; i < _N; i++)
    if (var_stage[i] == -1) {
      _s2h[_s2h_begin[i]] = j;
      _s2h_stage[_s2h_begin[i]] = 0;
      _h2s[j] = i;
      j++;
      _s2h_begin[i]++;
    }

  // correct _s2h_begin by shifting it right
  for (i = _N; i > 0; i--)
    _s2h_begin[i] = _s2h_begin[i-1];
  _s2h_begin[0] = 0;

#ifdef DEBUG
  fprintf(stderr, "Initializing A and C,\n");
  fflush(stderr);
#endif

#ifdef BLOW_OUTPUT
  // print out variable assignments
  for (i = 0; i < _N; i++) {
    fprintf(stderr, "%3d:", i);
    for (j = _s2h_begin[i]; j < _s2h_begin[i+1]; j++)
      fprintf(stderr, "%3d(%3d) ", _s2h[j], _s2h_stage[j]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "_s2h0: ");
  iv_foutput(stderr, _s2h0);
#endif

  // fill up _s2h0 for variables that don't ingo into the objective
  // and check for overlapping objective groups
  for (i = 0; i < _N; i++)
    if (_s2h0[i] == -1)
      _s2h0[i] = _s2h[_s2h_begin[i]];

  // reinitialize _var_ridx, _var_lb, _var_ub

  IVECP var_ridx = iv_get(n);
  VECP var_lb = v_get(n);
  VECP var_ub = v_get(n);
  iv_set(var_ridx, -1);
  v_set(var_lb, -_Inf);
  v_set(var_ub, _Inf);
  m_offs = 0;
  for (i = 0; i < _N; i++) {
    jidx = _s2h_begin[i];
    var_lb[_s2h[jidx]] = _var_lb[i];
    var_ub[_s2h[jidx]] = _var_ub[i];
    var_ridx[_s2h[jidx]] = _var_ridx[i];
    if (_var_ridx[i] >= 0 && _var_lb[i] != _var_ub[i]) {
      ijidx_end = _s2h_begin[i+1];
      for (jidx = _s2h_begin[i]; jidx < ijidx_end; jidx++) {
	var_lb[_s2h[jidx]] = _var_lb[i];
	var_ub[_s2h[jidx]] = _var_ub[i];
	if (_var_lb[i] > -_Inf)
	  var_ridx[_s2h[jidx]] = m_offs++;
	if (_var_ub[i] < _Inf) {
	  if (var_ridx[_s2h[jidx]] == -1)
	    var_ridx[_s2h[jidx]] = m_offs++;
	  else
	    m_offs++;
	}
      }
    }
  }
  v_free(_var_lb);
  v_free(_var_ub);
  iv_free(_var_ridx);
  _var_lb = var_lb;
  _var_ub = var_ub;
  _var_ridx = var_ridx;

  // variable bounds

  for (i = 0; i < n; i++) {
    if (_var_lb[i] == _var_ub[i])
      sp_set_val(_qp->A, _var_ridx[i], i, 1.0);
    else {
      offs = 0;
      if (_var_lb[i] > -_Inf) {
	sp_update_val(_qp->C, _var_ridx[i], i, 1.0);
	offs = 1;
      }
      if (_var_ub[i] < _Inf)
	sp_update_val(_qp->C, _var_ridx[i] + offs, i, -1.0);
    }
  }

  /*
   * setup the sparsity structure using the starting values
   */

  // general constraints and objective gradient

  qp_c = _qp->c;
  v_zero(qp_c);

  for (idx = 0; idx < _NNZJ; idx++) {
    i = _INDFUN[idx] - 1;
    j = _INDVAR[idx] - 1;
    if (i >= 0 && _cns_stage[i] >= 0) {
#if 0
      ijidx_end = _s2h_begin[j+1];
      for (jidx = _s2h_begin[j]; ijidx_end; jidx++)
	if (_s2h_stage[jidx] == _cns_stage[i]) {
	  j = _s2h[jidx];
	  break;
	}
#else
      jp = (int *)
	bsearch(&_cns_stage[i], &_s2h_stage[_s2h_begin[j]],
		_s2h_begin[j+1] - _s2h_begin[j], sizeof(int),
		(int (*)(const void *, const void *))cmp_int);
      // there should always be an entry for the stage
      if (!jp)
	error(E_INTERN, "Prg_CUTE_ST::setup");
      j = _s2h[_s2h_begin[j] + (jp - &_s2h_stage[_s2h_begin[j]])];
#endif
    }
    else
      j = _s2h0[j];
    if (i == -1) {
      qp_c[j] = _fscale * _J[idx];
    }
    else {
      if (_cns_lb[i] == _cns_ub[i])
	sp_set_val(_qp->A, _cns_ridx[i], j, _J[idx]);
      else {
	offs = 0;
	if (_cns_lb[i] > -_Inf) {
	  sp_set_val(_qp->C, _cns_ridx[i], j, _J[idx]);
	  offs = 1;
	}
	if (_cns_ub[i] < _Inf)
	  sp_set_val(_qp->C, _cns_ridx[i] + offs, j, -_J[idx]);
      }
    }
  }

#ifdef DEBUG
  fprintf(stderr, "Initializing Q,\n");
  fflush(stderr);
#endif

  // Hessian

  if (_hela_init || _hela) {
    for (l = 0; l < _NEL; l++) {
      ijidx_end = _iprnhi[l + 1];
      idx = _iprhi[l];
      for (iidx = _iprnhi[l]; iidx < ijidx_end; iidx++) {
	for (jidx = iidx; jidx < ijidx_end; jidx++) {
	  i = _irnhi[iidx];
	  j = _irnhi[jidx];
	  ip = (int *)
	    bsearch(&_grp_stage[l], &_s2h_stage[_s2h_begin[i]],
		    _s2h_begin[i+1] - _s2h_begin[i], sizeof(int),
		    (int (*)(const void *, const void *))cmp_int);
	  jp = (int *)
	    bsearch(&_grp_stage[l], &_s2h_stage[_s2h_begin[j]],
		    _s2h_begin[j+1] - _s2h_begin[j], sizeof(int),
		    (int (*)(const void *, const void *))cmp_int);
	  // there should always be entries for the stage
	  if (!ip || !jp)
	    error(E_INTERN, "Prg_CUTE_ST::setup");
	  i = _s2h[_s2h_begin[i] + (ip - &_s2h_stage[_s2h_begin[i]])];
	  j = _s2h[_s2h_begin[j] + (jp - &_s2h_stage[_s2h_begin[j]])];
	  sp_add_val(_qp->Q, min(i, j), max(i, j), _H[idx]);
	  idx++;
	}
      }
    }
  }

  _fbd_evals = 1;

  delete [] EQUATN;
  delete [] LINEAR;
  iv_free(indvar);
  iv_free(var_stage);
  iv_free(obj_grp);

  init_x();

#ifdef DEBUG
  fprintf(stderr, "done.\n\n");
  fflush(stderr);
#endif
}

//-------------------------------------------------------------------------
void Prg_CUTE_ST::init_x()
{
  int i, n = _x->dim;
  for (i = 0; i < n; i++)
    _x[i] = _xscale * _x0[_h2s[i]];
}

//-------------------------------------------------------------------------
void Prg_CUTE_ST::update_bounds(const VECP val,
			       const VECP lb, const VECP ub, int n,
			       const IVECP ridx)
{
  int i, offs;
  VECP b, d;

  b = _qp->b;
  d = _qp->d;

  for (i = 0; i < n; i++) {
    if (lb[i] == ub[i])
      b[ridx[i]] = val[i] - lb[i];
    else {
      offs = 0;
      if (lb[i] > -_Inf) {
	d[ridx[i]] = val[i] - lb[i];
	offs = 1;
      }
      if (ub[i] < _Inf)
	d[ridx[i] + offs] = ub[i] - val[i];
    }
  }
}

//-------------------------------------------------------------------------
void Prg_CUTE_ST::update_fbd()
{
  int i, idx_end, idx, k;
  VECP qp_b = _qp->b;
  fint I, NNZSGC = 0;
  fbool GRAD = 0;

  if (!_cfn_p)
    error(E_NULL, "Prg_CUTE_ST::update_fbd");
  
#ifdef RECALCULATE_CNS
  if (!_cscifg_p)
    error(E_NULL, "Prg_CUTE_ST::update_fbd");
#endif

  for (i = 0; i < _N; i++)
    _xs[i] = 1.0/_xscale * _x[_s2h0[i]];

  (*_cfn_p)(&_N, &_M, _xs->ve, &_f, &_M, _cns->ve);

  _f *= _fscale;

  update_bounds(_x, _var_lb, _var_ub, _x->dim, _var_ridx);

#ifdef RECALCULATE_CNS
  // recalculate nonlinear constraints
  for (k = 0; k < _nstages; k++) {
    for (i = 0; i < _N; i++) {
      _xs[i] = 1.0/_xscale * _x[_s2h0[i]];
      idx_end = _s2h_begin[i+1];
      for (idx = _s2h_begin[i] + 1; idx < idx_end; idx++) {
	if (_s2h_stage[idx] == k) {
	  _xs[i] = 1.0/_xscale * _x[_s2h[idx]];
	  break;
	}
      }
    }
    for (i = 0; i < _M; i++) {
      if (_cns_stage[i] == k) {
	I = i + 1;
	(*_cscifg_p)(&_N, &I, _xs->ve, _cns->ve + i,
		     &NNZSGC, &NNZSGC, NULL, NULL, &GRAD);
      }
    }
  }
#endif

  update_bounds(_cns, _cns_lb, _cns_ub, _M, _cns_ridx);

  // update blow-up's
  idx_end = _me_lin - _me_bnd;
  for (idx = 0; idx < idx_end; idx++)
    qp_b[_me_bnd + idx] =
      1.0/_xscale * (_x[_equ_idx0[idx]] - _x[_equ_idx1[idx]]);
    
  _fbd_evals++;
}

//-------------------------------------------------------------------------
void Prg_CUTE_ST::update(const VECP y, const VECP z)
{
  fbool GRLAGF = 0;
  fbool BYROWS = 1;
  int i, i_end, offs;
  int j, k, l, idx, idx_end;
  int *ip, *jp;
  int iidx, jidx, ijidx_end;
  SPROW *row;
  int newel_h, newel_j;
  VECP qp_c = _qp->c;
  double val;
  fbool GRAD = 1;
  fint  I, NNZSGC, NNZSGCMAX;
  freal CI;

  update_fbd();
  _fbd_evals--;

  if (!(const VEC *)y || !(const VEC *)z)
    error(E_NULL, "Prg_CUTE_ST::update");

  if (!_csgreh_p || !_csgr_p)
    error(E_NULL, "Prg_CUTE_ST::update");

#ifdef RECALCULATE_CNS
  if (!_cscifg_p)
    error(E_NULL, "Prg_CUTE_ST::update");
#endif

  for (i = 0; i < _N; i++)
    _xs[i] = 1.0/_xscale * _x[_s2h0[i]];

  for (i = 0; i < _M; i++) {
    idx = _cns_ridx[i];
    if (_cns_lb[i] == _cns_ub[i])
      _v[i] = y->ve[idx];
    else {
      offs = 0;
      _v[i] = 0.0;
      if (_cns_lb[i] > -_Inf) {
	_v[i] -= z->ve[idx];
	offs = 1;
      }
      if (_cns_ub[i] < _Inf)
	_v[i] += z->ve[idx + offs];
    }
  }
  v_zero(_v);

  if (_hela) {
    sv_mlt(1.0/_fscale, _v, _v);
    (*_csgreh_p)(&_N, &_M, _xs->ve, &GRLAGF, &_M, _v->ve,
		 &_NNZJ, &_NNZJMAX, _J->ve, _INDVAR, _INDFUN,
		 &_NEL, _IRNHI, &_NHIMAX,
		 &_NELMAX, _IPRNHI, _H->ve, &_NHIMAX, _IPRHI, &BYROWS);
    _J->dim = _NNZJ;
    _H->dim = _IPRHI[_NEL] - 1;

    sv_mlt(_fscale, _v, _v);
    sv_mlt(_fscale, _H, _H);
    sv_mlt(1.0/_xscale/_xscale, _H, _H);
  }
  else {
    (*_csgr_p)(&_N, &_M, &GRLAGF, &_M, _v->ve, _xs->ve,
	       &_NNZJ, &_NNZJMAX, _J->ve, _INDVAR, _INDFUN);

    _J->dim = _NNZJ;
  }
  sv_mlt(1.0/_xscale, _J, _J);

#ifdef RECALCULATE_CNS
  // recalculate nonlinear constraints
  for (k = 0; k < _nstages; k++) {
    for (i = 0; i < _N; i++) {
      _xs[i] = 1.0/_xscale * _x[_s2h0[i]];
      idx_end = _s2h_begin[i+1];
      for (idx = _s2h_begin[i] + 1; idx < idx_end; idx++) {
	if (_s2h_stage[idx] == k) {
	  _xs[i] = 1.0/_xscale * _x[_s2h[idx]];
	  break;
	}
      }
    }
    idx_end = 0;
    for (i = 0; i < _M; i++) {
      idx = idx_end;
      while (idx_end < _NNZJ && _INDFUN[idx_end] - 1 == i)
	idx_end++;
      if (_cns_stage[i] == k) {
	I = i + 1;
	NNZSGCMAX = idx_end - idx;
	(*_cscifg_p)(&_N, &I, _xs->ve, &CI,
		     &NNZSGC, &NNZSGCMAX, _J->ve + idx, _INDVAR + idx, &GRAD);
      }
    }
  }
#endif
  // update Jacobians

  row = _qp->A->row + _me_lin;
  i_end = _qp->A->m;
  for (i = _me_lin; i < i_end; i++, row++)
    sprow_zero(row);

  row = _qp->C->row + _m_lin;
  i_end = _qp->C->m;
  for (i = _m_lin; i < i_end; i++, row++)
    sprow_zero(row);

  newel_j = 0;
  v_zero(qp_c);

  for (idx = 0; idx < _NNZJ; idx++) {
    val = _J[idx];
    if (val == 0.0)
      continue;
    i = _INDFUN[idx] - 1;
    j = _INDVAR[idx] - 1;
    if (i >= 0 && _cns_stage[i] >= 0) {
#if 0
      ijidx_end = _s2h_begin[j+1];
      for (jidx = _s2h_begin[j]; ijidx_end; jidx++)
	if (_s2h_stage[jidx] == _cns_stage[i]) {
	  j = _s2h[jidx];
	  break;
	}
#else
      jp = (int *)
	bsearch(&_cns_stage[i], &_s2h_stage[_s2h_begin[j]],
		_s2h_begin[j+1] - _s2h_begin[j], sizeof(int),
		(int (*)(const void *, const void *))cmp_int);
      // there should always be an entry for the stage
      if (!jp)
	error(E_INTERN, "Prg_CUTE_ST::setup");
      j = _s2h[_s2h_begin[j] + (jp - &_s2h_stage[_s2h_begin[j]])];
#endif
    }
    else
      j = _s2h0[j];
    if (i == -1) {
      qp_c[j] = _fscale * val;
    }
    else {
      if (_cns_lb[i] == _cns_ub[i])
	newel_j |= sp_update_val(_qp->A, _cns_ridx[i], j, val);
      else {
	offs = 0;
	if (_cns_lb[i] > -_Inf) {
	  newel_j |= sp_update_val(_qp->C, _cns_ridx[i], j, val);
	  offs = 1;
	}
	if (_cns_ub[i] < _Inf)
	  newel_j |= sp_update_val(_qp->C, _cns_ridx[i] + offs, j, -val);
      }
    }
  }

  // Hessian

  newel_h = 0;
  if (_hela) {
    sp_zero(_qp->Q);
    for (l = 0; l < _NEL; l++) {
      ijidx_end = _iprnhi[l + 1];
      idx = _iprhi[l];
      for (iidx = _iprnhi[l]; iidx < ijidx_end; iidx++) {
	for (jidx = iidx; jidx < ijidx_end; jidx++) {
	  i = _irnhi[iidx];
	  j = _irnhi[jidx];
	  ip = (int *)
	    bsearch(&_grp_stage[l], &_s2h_stage[_s2h_begin[i]],
		    _s2h_begin[i+1] - _s2h_begin[i], sizeof(int),
		    (int (*)(const void *, const void *))cmp_int);
	  jp = (int *)
	    bsearch(&_grp_stage[l], &_s2h_stage[_s2h_begin[j]],
		    _s2h_begin[j+1] - _s2h_begin[j], sizeof(int),
		    (int (*)(const void *, const void *))cmp_int);
	  // there should always be entries for the stage
	  if (!ip || !jp)
	    error(E_INTERN, "Prg_CUTE_ST::setup");
	  i = _s2h[_s2h_begin[i] + (ip - &_s2h_stage[_s2h_begin[i]])];
	  j = _s2h[_s2h_begin[j] + (jp - &_s2h_stage[_s2h_begin[j]])];
	  sp_add_val(_qp->Q, min(i, j), max(i, j), _H[idx]);
	  idx++;
	}
      }
    }
  }

  if (newel_j)
    warning(WARN_UNKNOWN,
	    "Prg_CUTE_ST::update: Jacobian sparsity structure changed");
  if (newel_h)
    warning(WARN_UNKNOWN,
	    "Prg_CUTE_ST::update: Hessian sparsity structure changed");

  extern Tcl_Interp *theInterp;
  if (newel_j || newel_h) {
    init_x();
    Tcl_Eval(theInterp, "sqp_init; sqp_qp_update");
  }
}

//-------------------------------------------------------------------------
int Prg_CUTE_ST::write_soln(int argc, char *argv[], char **)
{
  char header[61];
  fint idum = 0;
  freal rdum = 0.0;
  extern const char *Hqp_Version;

  if (!_cwrtsn_p)
    error(E_NULL, "Prg_CUTE_ST::write_soln");

  sprintf(header, "HQP %-4.4s: ", Hqp_Version);
  if (argc > 1) 
    sprintf(header + 10, "%-50.50s", argv[argc - 1]);
  else
    sprintf(header + 10, "%-50.50s", "");

  if (argc > 1 && strncmp(argv[1], "-nosol", 6) == 0)
    (*_cwrtsn_p)(&idum, &idum, header, &rdum, &rdum, &rdum);
  else
    (*_cwrtsn_p)(&_N, &_M, header, &_f, _xs->ve, _v->ve);

  return IF_OK;
}


//=========================================================================
