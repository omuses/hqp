/*
 * Omu_IntSDIRK.C --
 *   -- class for integrating an Ode over a stage
 *   -- using implicit Runge Kutta method
 *
 * hl, 31/07/00
 */

/*
    Copyright (C) 2000  Hartmut Linke

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
#include <math.h>
#include <iostream.h>

#ifdef OMU_WITH_ADOLC

#include <adutils.h>

#endif

#include <Hqp.h>
#include <If_Int.h>
#include <If_Bool.h>
#include <If_Real.h>
#include <If_Class.h>
#include <Meschach.h>

#include "Omu_IntSDIRK.h"

extern "C" {
#include <meschach/matrix2.h>
#include <meschach/sparse.h>
#include <meschach/sparse2.h>
}

SPMAT *sp_mmt2_mlt(const SPMAT *, const SPMAT *, SPMAT *);

IF_CLASS_DEFINE("SDIRK", Omu_IntSDIRK, Omu_Integrator);

//--------------------------------------------------------------------------

Omu_IntSDIRK::Omu_IntSDIRK()
{

  _sys = NULL;
  _cx_ptr = NULL;
  _cF_ptr = NULL;

  _recalc_jac = true;
  _stiffly_accurate = false;
  _sens_adolc = false;
  _sens_at_once = false;
  _lsqr_sol = false;
  _init_method = false; 

  _banded = true;
  _sparse_solver = true;

  _output = 0;
  _n_splitt_tape_eval = 30;
  _tag = 10;

  _jac_sbw = -1;

  _x_algebraic = iv_get(1);

  _ppivot= px_get(1);
  _Spivot= px_get(1);

  _cFh = v_get(1);

  _par = v_get(1);
  _y = v_get(1);
  _y0 = v_get(1);
  _yprime0 = v_get(1);
  _yprime1 = v_get(1);
  _yprime2 = v_get(1);
  _irk_delta = v_get(1);
  _irk_res = v_get(1);
  _irk_y = v_get(1);
  _irk_yprime = v_get(1);

  _yprime = m_get(1,1);

  _cFxh = m_get(1,1);
  _Sxd  = m_get(1,1);
  _Sxd0 = m_get(1,1);
  _Sxx  = m_get(1,1);
  _Sxx0 = m_get(1,1);
  _Sxu  = m_get(1,1);
  _Sxu0 = m_get(1,1);
  _SF   = m_get(1,1);
  _Sh   = m_get(1,1);

  _U = m_get(1,1);
  _Z = m_get(1,1);

  _irk_jac = m_get(1,1);
  _irk_jacf = m_get(1,1);
  _irk_jac_bd = bd_get(0,0,1);

  _irk_jac_sp = SMNULL;
  _smh = SMNULL;
  _smh1 = SMNULL;

  _na = 0;
  _nod = 0;

  _ml = 0;
  _mu = 0;

  _nsteps = 0;
  _modnewtonsteps = 5;
  _newtonsteps = 0;
  _maxiters = 100;

  _hinit = 0.0;
  _h = 0.0;
  _h_new = 0.0;

  _err = v_get(1);

  // irk parameters

  _irk_stages = 0;
  _a = m_get(1,1);
  _a_pred = m_get(1,1);
  _c = v_get(1);
  _b = v_resize(v_get(1),0);
  _b_err = v_resize(v_get(1),0);

  _vh = v_get(1);

  _gamma0 = 0.0;

#ifdef OMU_WITH_ADOLC

  _ifList.append(new If_Bool("prg_int_sensadolc", &_sens_adolc));
  _ifList.append(new If_Bool("prg_int_sensatonce", &_sens_at_once));

#endif

  _ifList.append(new If_Bool("prg_int_lsqrsol", &_lsqr_sol));
  _ifList.append(new If_Bool("prg_int_sparsol", &_sparse_solver));
  _ifList.append(new If_Bool("prg_int_banded", &_banded));
  _ifList.append(new If_Int("prg_int_out", &_output));
  _ifList.append(new If_Int("prg_int_method", &_method));
  _ifList.append(new If_Int("prg_int_nsteps", &_nsteps));
  _ifList.append(new If_Int("prg_int_nsplitt", &_n_splitt_tape_eval));
  _ifList.append(new If_Int("prg_int_modnewtonsteps", &_modnewtonsteps));
  _ifList.append(new If_Int("prg_int_maxiters", &_maxiters));
  _ifList.append(new If_Real("prg_int_hinit", &_hinit));

}

//--------------------------------------------------------------------------

Omu_IntSDIRK::~Omu_IntSDIRK()
{

  IV_FREE(_x_algebraic);

  PX_FREE(_ppivot);
  PX_FREE(_Spivot);

  V_FREE(_cFh);
  V_FREE(_par);
  V_FREE(_y);
  V_FREE(_y0);
  V_FREE(_yprime0);
  V_FREE(_yprime1);
  V_FREE(_yprime2);
  V_FREE(_irk_delta);
  V_FREE(_irk_res);
  V_FREE(_irk_y);
  V_FREE(_irk_yprime);

  V_FREE(_err);
  V_FREE(_vh);

  V_FREE(_b);
  V_FREE(_b_err);
  V_FREE(_c);
  M_FREE(_a);
  M_FREE(_a_pred);

  M_FREE(_yprime);
  M_FREE(_cFxh);

  M_FREE(_Sxd);
  M_FREE(_Sxd0);
  M_FREE(_Sxx);
  M_FREE(_Sxx0);
  M_FREE(_Sxu);
  M_FREE(_Sxu0);
  M_FREE(_SF);
  M_FREE(_Sh);

  M_FREE(_U);
  M_FREE(_Z);
  M_FREE(_irk_jac);
  M_FREE(_irk_jacf);

  bd_free(_irk_jac_bd);
  sp_free(_irk_jac_sp);
  sp_free(_smh);
  sp_free(_smh1);

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::init_method4()
{

  int   i;

  _irk_stages = 4;

  m_resize(_a,4,4);
  m_zero(_a);
  _a[0][0] = 1.0/4.0;
  _a[1][0] = 1.0/7.0;
  _a[1][1] = 1.0/4.0;
  _a[2][0] = 6.1e1/1.44e2;
  _a[2][1] = -4.9e1/1.44e2;
  _a[2][2] = 1.0/4.0;
  _a[3][0] = 0.0;
  _a[3][1] = 0.0;
  _a[3][2] = 3.0/4.0;
  _a[3][3] = 1.0/4.0;

  v_resize(_c,4);
  _c[0] = 1.0/4.0;
  _c[1] = 1.1e1/2.8e1;
  _c[2] = 1.0/3.0;
  _c[3] = 1.0;
 
  v_resize(_b,4);
  _b[0] = 0.0;
  _b[1] = 0.0;
  _b[2] = 3.0/4.0;
  _b[3] = 1.0/4.0;

  v_resize(_b_err,4);
  _b_err[0] = -6.1e1/7.04e2;
  _b_err[1] = 4.9e1/7.04e2;
  _b_err[2] = 6.9e1/8.8e1;
  _b_err[3] = 4.1e1/1.76e2;

  if(_c[(int)_c->dim-1] == 1.0) {
    _stiffly_accurate = true;
    for(i = 0; i < _irk_stages; i++)
      if(_b[i] != _a[_irk_stages-1][i]) {
	_stiffly_accurate = false;
	break;
      }
  }

  _gamma0 = 2.5e-1;
  // _gamma0 = 0.0;

  // simple linear interpolation

  m_resize(_a_pred,_irk_stages,_irk_stages+1);
  m_zero(_a_pred);

  _a_pred[0][0] = 1.0;
  _a_pred[1][0] = -1.0/7.0;
  _a_pred[1][1] = 1.25;
  _a_pred[2][1] = 5.0/1.2e1;
  _a_pred[2][2] = 7.0/1.2e1;
  _a_pred[3][1] = -1.7e1/4.0;
  _a_pred[3][2] = 2.1e1/4.0;

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::init_method5()
{

  int   i;

  _irk_stages = 5;

  m_resize(_a,5,5);
  m_zero(_a);
  _a[0][0] = 1.0/4.0;
  _a[1][0] = 1.0/2.0;
  _a[1][1] = 1.0/4.0;
  _a[2][0] = 1.7e1/5.0e1;
  _a[2][1] = -1.0/2.5e1;
  _a[2][2] = 1.0/4.0;
  _a[3][0] = 3.71e2/1.36e3;
  _a[3][1] = -1.37e2/2.72e3;
  _a[3][2] = 1.5e1/5.44e2;
  _a[3][3] = 1.0/4.0;
  _a[4][0] = 2.5e1/2.4e1;
  _a[4][1] = -4.9e1/4.8e1;
  _a[4][2] = 1.25e2/1.6e1;
  _a[4][3] = -8.5e1/1.2e1;
  _a[4][4] = 1.0/4.0;

  v_resize(_c,5);
  _c[0] = 1.0/4.0;
  _c[1] = 3.0/4.0;
  _c[2] = 11.0/20.0;
  _c[3] = 1.0/2.0;
  _c[4] = 1.0;
 
  v_resize(_b,5);
  _b[0] = 2.5e1/2.4e1;
  _b[1] = -4.9e1/4.8e1;
  _b[2] = 1.25e2/1.6e1;
  _b[3] = -8.5e1/1.2e1;
  _b[4] = 1.0/4.0;

  v_resize(_b_err,5);
  _b_err[0] = 5.9e1/4.8e1;
  _b_err[1] = -1.7e1/9.6e1;
  _b_err[2] = 2.25e2/3.2e1;
  _b_err[3] = -8.5e1/1.2e1;
  _b_err[4] = 0.0;

  if(_c[(int)_c->dim-1] == 1.0) {
    _stiffly_accurate = true;
    for(i = 0; i < _irk_stages; i++)
      if(_b[i] != _a[_irk_stages-1][i]) {
	_stiffly_accurate = false;
	break;
      }
  }

  // _gamma0 = 1.0/4.000264239115785;
  _gamma0 = 2.5e-1;

  // simple linear interpolation

  m_resize(_a_pred,_irk_stages,_irk_stages+1);
  m_zero(_a_pred);

  _a_pred[0][0] = 2.0;
  _a_pred[0][2] = -1.0;
  _a_pred[1][0] = -2.0;
  _a_pred[1][1] = 3.0;
  _a_pred[2][1] = 0.5;
  _a_pred[2][2] = 0.5;
  _a_pred[3][3] = 1.0;
  _a_pred[4][2] = -1.0;
  _a_pred[4][4] = 2.0;

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::init_method2()
{

  int     i;
  double  gamma;

  _irk_stages = 2;

  gamma = (2.0 - sqrt(2.0))/2.0;

  m_resize(_a,2,2);
  m_zero(_a);
  _a[0][0] = gamma;
  _a[1][0] = 1.0 - gamma;
  _a[1][1] = gamma;

  v_resize(_c,2);
  _c[0] = gamma;
  _c[1] = 1.0;
 
  v_resize(_b,2);
  _b[0] = 1.0 - gamma;
  _b[1] = gamma;

  v_resize(_b_err,2);
  _b_err[0] = -1;
  _b_err[1] = -1;

  if(_c[(int)_c->dim-1] == 1.0) {
    _stiffly_accurate = true;
    for(i = 0; i < _irk_stages; i++)
      if(_b[i] != _a[_irk_stages-1][i]) {
	_stiffly_accurate = false;
	break;
      }
  }

  // _gamma0 = 1.0/3.414213562373096;
  _gamma0 = 0.0;

  // simple linear interpolation

  m_resize(_a_pred,_irk_stages,_irk_stages+1);
  m_zero(_a_pred);

  _a_pred[0][0] = 1.0/(1.0 - gamma);
  _a_pred[0][1] = gamma/(gamma - 1.0);
  _a_pred[1][0] = (gamma - 1.0)/gamma;
  _a_pred[1][1] = 1.0/gamma;

}


//--------------------------------------------------------------------------
void Omu_IntSDIRK::resize()
{

  if ((int)_cxp->dim == _nd + _n && (int)_cu->dim == _nu)
    return;

  //
  // realloc variables for low level _sys->continuous callback
  //

  iv_resize(_x_algebraic,_n);

  v_resize(_cu, _nu);
 _cxp.resize(_nd + _n, _nx, _nu);

  v_resize(_cFh, _n);

  v_resize(_par, _nd+_nu);
  v_resize(_y, _n);
  v_resize(_y0, _n);
  v_resize(_yprime0, _n);
  v_resize(_yprime1, _n);
  v_resize(_yprime2, _n);

  m_resize(_yprime, _irk_stages, _n);

  m_resize(_cFxh, _n, _n);

  m_resize(_Sxd0,_n,_nd);
  m_resize(_Sxd,_n,_nd);
  m_resize(_Sxx,_n,_n);
  m_resize(_Sxx0,_n,_n);
  m_resize(_Sxu,_n,_nu);
  m_resize(_Sxu0,_n,_nu);
  m_resize(_SF,_irk_stages*_n,_irk_stages*_n);
  m_resize(_Sh,_irk_stages*_n,_nd+_n+_nu);
  px_resize(_Spivot,_irk_stages*_n);

  v_resize(_err,_n);

  px_resize(_ppivot,_n);
  v_resize(_irk_delta,_n);
  v_resize(_irk_res,_n);
  v_resize(_irk_y,_n);
  v_resize(_irk_yprime,_n);

  m_resize(_irk_jac,_n,_n);
  m_resize(_irk_jacf,_n,_n);

  if(_banded_solver)
    bd_resize(_irk_jac_bd,_ml,min(_n-1,_ml+_mu),_n);

  if(_sparse_solver) {
    sp_free(_irk_jac_sp);
    _irk_jac_sp = sp_get(_n,_n,_ml+_mu);
    sp_free(_smh);
    _smh = sp_get(_n,_n,_ml+_mu);
    sp_free(_smh1);
    _smh1 = sp_get(_n,_n,_ml+_mu);
  }
 
}


//--------------------------------------------------------------------------
void Omu_IntSDIRK::init_stage(int k,
			    const Omu_States &x, const Omu_Vector &u,
			    bool sa)
{

  int  i, j;

  if(!_init_method) {

    switch(_method) {

    case 2:

      if(_stepsize <= 0.0)
	m_error(E_INTERN,"Omu_IntSDIRK::init_stage - only fixed"
		  " stepsize allowed for this method!");
      init_method2();

    case 4:
      init_method4();

    default:
      init_method5();
    }
  }

  Omu_Integrator::init_stage(k, x, u, sa);
 
  if(x.na > 0) {
    _na = x.na;
    iv_resize(_x_algebraic,x->dim);
    iv_zero(_x_algebraic);
    for(i = 0; i < _n; i++) {
      if(x.flags[_nd+i] & Omu_States::Algebraic)
	_x_algebraic[i] = 1;	// indicate algebraic state
    }
  }

  _nv = x.nv;

  // check for sdirk method
  for(i = 0; i < _irk_stages-1; i++)
    for(j = i+1; j < _irk_stages; j++)
      if(_a[i][j] != 0.0)
	m_error(E_INTERN,"Omu_IntSDIRK::init_stage - " 
		"sdirk method not correctly defined!");

  for(i = 1; i < _irk_stages; i++)
      if(_a[i][i] != _a[0][0])
	m_error(E_INTERN,"Omu_IntSDIRK::init_stage - "
		"sdirk method not correctly defined!");

  if(_sparse_solver && !_banded_solver)
      _banded_solver = true;

  if (_jac_sbw >= 0) {
    // take user-defined value
    _ml = _mu = _jac_sbw;
    _banded_solver = _banded;
  }
  else {
    _ml = x.sbw_l;
    _mu = x.sbw_u;

    _banded_solver = _banded;
    // disable banded solver if full solver appears more efficient
    if (_ml+_ml+_mu >= _n)
      _banded_solver = false;
  }

  resize();

}

//--------------------------------------------------------------------------

void Omu_IntSDIRK::solve(int kk, Real tstart, Real tend,
			 const Omu_VariableVec &x, const Omu_VariableVec &u,
			 Omu_Program *sys, Omu_DependentVec &cF,
			 Omu_StateVec &cx)
{

  bool    ok, tend_ok;
  int     i, j, l, nsteps, steps;
  double  err;

  _kk = kk;	// propagate to sys->continuous
  _sys = sys;	// propagate to res() and jac()

  err = 0.0;

  _cx_ptr = &cx;
  _cF_ptr = &cF;

  // _stepsize overrides _nsteps
  nsteps = _nsteps;
  if (_stepsize > 0.0)
    nsteps = (int)ceil((tend - tstart) / _stepsize);

  if(_output > 1) {
    cout << endl;
    cout << "Omu_IntSDIRK::solve at stage " << kk;
    cout.precision(16);
    cout << "  tstart: " << tstart;
    cout << "  tend: " << tend << endl;
    cout.precision(8);
    v_output(cx);
  }

  m_zero(_Sxx0);
  m_zero(_Sxd0);
  m_zero(_Sxu0);

  if(_sa) {
    m_move(cx.Sx,_nd,0,_n-_nv,_nd,_Sxd0,0,0);
    m_move(cx.Sx,_nd,_nd,_n,_n-_nv,_Sxx0,0,0);
    m_move(cx.Su,_nd,0,_n,_nu,_Sxu0,0,0);
  }

  // compute initial values for yprime
  v_zero(_yprime0);
  for(i = 0; i < _n; i++)
    _yprime0[i] = 1.0e-5;
  init_yprime(kk,tstart,cx,u,_yprime0);

  if(nsteps > 0)
    _h_new = (tend - tstart)/nsteps;
  else
    if(_hinit > 0.0)
      _h_new = min(_hinit,tend-tstart);
    else
      if(_h_new == 0.0)
	_h_new = (tend - tstart)/1.0e3;

  v_zero(_y);

  for (i = 0; i < _nd; i++)
    _par[i] = cx[i];

  for (i = 0; i < _n; i++)
    _y[i] = cx[_nd + i];		// initial states

  for (i = 0; i < _nu; i++)
    _par[_nd + i] = u[i];

  v_zero(_irk_delta);

  _t = tstart;

  for(i = 0; i < _n; i++)
    if(_x_algebraic[i])
      _yprime0[i] =  _y[i];

  set_row(_yprime,1,_yprime0);

  tend_ok = false;
  steps = 0;

  _recalc_jac = true;

  if(_output > 1)
    cout << "h: " << _h_new << endl;
 
  while(_t < tend) {    

    steps++;
    _h = _h_new;

    if(_t + _h >= tend) {
      _h = tend - _t;
      tend_ok = true;
    }

    if(nsteps > 0 && steps == nsteps)
      tend_ok = true;

    v_copy(_y,_y0);
    v_copy(_yprime0,_yprime1);
    set_row(_yprime,0,_yprime0);
    ok = true;
    _newtonsteps = 0;

    do {

      // solve stages subsequently

      for(i = 0; i < _irk_stages; i++) {
      
	// prepare y and yprime for the current stage

	v_copy(_y0,_y);
	v_zero(_yprime0);

	for(j = 0; j < i; j++)
	  for(l = 0; l < _n; l++)
	    _yprime0[l] += _a_pred[i][0]*_yprime1[l];
	for(j = 0; j < i; j++)
	  for(l = 0; l < _n; l++)
	    _yprime0[l] += _a_pred[i][j+1]*_yprime[j][l];

	for(j = 0; j < i; j++)
	  for(l = 0; l < _n; l++)
	    if(!_x_algebraic[l])
	      _y[l] +=  _h*_a[i][j]*_yprime[j][l];

	if(_nod == 0 && !_lsqr_sol)
	  solve_stage(i+1,_y,_yprime0);
	else
	  solve_stage_lsqr(i+1,_y,_yprime0);

	for(l = 0; l < _n; l++)
	  if(_x_algebraic[l])
	    _yprime0[l] = _y[l];
	
	set_row(_yprime,i,_yprime0);
	
      }

      // now we have to check the error
      
      if(nsteps <= 0) {
	if((int)_b_err->dim == _irk_stages) {
	  v_copy(_y0,_err);
	  for(j = 0; j < _irk_stages; j++)
	    for(l = 0; l < _n; l++)
	      if(!_x_algebraic[l])
		_err[l] += _h*_b_err[j]*_yprime[j][l];
	      else
		_err[l] = _y[l];

	  /*
	  if(_na > 0)
	    solve_final(_err,_yprime2);                       // to do!!!
	  */
	}
      }
      
      if(!_stiffly_accurate) {
	v_copy(_y0,_y);
	for(j = 0; j < _irk_stages; j++)
	  for(l = 0; l < _n; l++)
	    if(!_x_algebraic[l])
	      _y[l] +=  _h*_b[j]*_yprime[j][l];
	    else
	      _y[l] = _yprime0[l];

	if(_na > 0)
	  solve_final(_y,_yprime0);
      }

      if(nsteps <= 0) {

	// algebraic variables are not considered in the error evaluation

	for(i = 0; i < _n; i++)
	  if(_x_algebraic[i])
	    _yprime1[i] = 0.0;

	if(_gamma0 > 0) {

	  v_sub(_err,_y,_irk_y);
	  jac(-3,_t,_y,_yprime0,_par,_irk_res);

	  if(_nod == 0 && !_lsqr_sol) {
	    
	    if(_banded_solver) {
	      if(_sparse_solver) {
		sp_resize(_irk_jac_sp,_n,_n);
		sp_insert_mat(_irk_jac_sp,0,0,_irk_jac);
		spLUfactor2(_irk_jac_sp, _ppivot);
		spLUsolve(_irk_jac_sp, _ppivot, _irk_y, _err);
	      }
	      else {
		bd_resize(_irk_jac_bd,_ml,_mu,_irk_jac->m);
		mat2band(_irk_jac,_ml,_mu,_irk_jac_bd);
		bdLUfactor(_irk_jac_bd, _ppivot);
		bdLUsolve(_irk_jac_bd, _ppivot, _irk_y, _err);
	      }
	    }
	    else {
	      LUfactor(_irk_jac, _ppivot);
	      LUsolve(_irk_jac, _ppivot, _irk_y, _err);
	    }
	  }
	  else {
	    m_resize(_irk_jacf,_irk_jac->n,_irk_jac->m);
	    m_resize(_Sh,_n,_n);
	    m_transp(_irk_jac,_irk_jacf);
	    m_mlt(_irk_jacf,_irk_jac,_Sh);
	    for(j = 0; j < _n; j++)
	      if(_x_algebraic[j] && _Sh[j][j] == 0.0)
		_Sh[j][j] = 1.0e-15;
	    
	    if(_banded_solver) {
	      if(_sparse_solver) {
		symsp_insert_symmat(_irk_jac_sp,0,_Sh);
		spBKPfactor(_irk_jac_sp, _ppivot,0.0);
		spBKPsolve(_irk_jac_sp, _ppivot, _irk_y, _err);	      
	      }
	      else {
		px_resize(_Spivot,_n);
		bd_resize(_irk_jac_bd,max(_ml,_mu),max(_ml,_mu),_Sh->m);
		mat2band(_Sh,max(_ml,_mu),max(_ml,_mu),_irk_jac_bd);
		bdBKPfactor(_irk_jac_bd, _ppivot, _Spivot);
		bdBKPsolve(_irk_jac_bd, _ppivot, _Spivot, _irk_y, _err);
	      }
	    }
	    else {
	      matBKPfactor(_Sh, _ppivot);
	      matBKPsolve(_Sh, _ppivot, _irk_y, _err);
	    }
	  }
	}
	else 
	  v_sub(_err,_y,_err);

	err = 0.0;

	for(i = 0; i < _n; i++)
	  if(!_x_algebraic[i]) {
	    _err[i] = _err[i]/(_atol+max(fabs(_y0[i]),fabs(_y[i]))*_rtol);
	    err = max(err,fabs(_err[i]));
	  }
	  else
	    _err[i] = 0.0;

	err = max(v_norm2(_err)/sqrt(double(_n-_na)),_atol);

	_h_new = 0.9*_h*pow(err,-0.25);
	
	if(err > 1.0) {
	  ok = false;
	  _h = _h_new;
	  tend_ok = false;
	}
	else
	  ok = true;

	if(_output > 1 && !ok) {
	  cout << "t: " << _t << "  err: " << err << "  h_new: ";
	  cout << _h_new << " accepted: " << ok;
	  cout << "  newtonsteps: " << _newtonsteps << endl;
	  v_output(_err);
	}
      }     
    }
    while(!ok);

    if(_output > 1) {
      cout.precision(16);
      cout << "t: " << _t;
      cout.precision(8);
      cout << "  err: " << err << "  h_new: ";
      cout << _h_new << " accepted: " << ok;
      cout << "  newtonsteps: " << _newtonsteps << endl;
    }

    if(_sa) {
      if(_sens_adolc) {

#ifdef OMU_WITH_ADOLC

	if(_nod == 0 && !_lsqr_sol)
	  sensitivity_adolc();
        else
	  sensitivity_lsqr_adolc();
#endif

      }
      else {
        if(_nod == 0 && !_lsqr_sol)
	  sensitivity();
        else
	  sensitivity_lsqr();
      }
    }

    _t += _h; 

    if(tend_ok)
      break;

  }

  _h_new = _h;

  // return results of integration process
  for(i = 0; i < _n; i++)
    cx[_nd+i] = _y[i];

  for(i = 0; i < _nd; i++)
    cx[i] = _par[i];
  
  if(_output > 1)
    v_output(_y);
    
  if(_sa) {
    m_move(_Sxd0,0,0,_n,_nd,cx.Sx,_nd,0);
    m_move(_Sxx0,0,0,_n,_n-_nv,cx.Sx,_nd,_nd);
    m_move(_Sxu0,0,0,_n,_nu,cx.Su,_nd,0);
  }

}

//--------------------------------------------------------------------------

void Omu_IntSDIRK::solve_stage(int stage, VECP y, VECP yprime)
{

  bool    ok;
  int     i, j;
  double  alpha, cur_err, old_err;

  v_ones(_irk_res);
  v_copy(y,_irk_y);
  v_copy(yprime,_irk_yprime);

  cur_err = 1.0e10;
  old_err = cur_err;
  ok = false;
  // _recalc_jac = false;
    
  if(_output > 2)
    cout << "   stage: " << stage << endl;
  
  // modified Newton's method    
  for(i = 0; i < _maxiters; i++ ) {
    
    // if(stage < 2 && (i % _modnewtonsteps) == 0)
    //   _recalc_jac = true;

    if(i > 0 && (i % _modnewtonsteps) == 0)
      _recalc_jac = true;

    alpha = 0.5;
    
    if(i > 0) {

      do {

	cur_err = 0.0;
	// check for convergence
	for(ok = true, j = 0 ; j < _n; j++) {
	  ok = ok  && 
	    (fabs(_irk_res[j]) < 1.0e-2*(_atol + _rtol*fabs(_irk_y[j])));
	  cur_err = max(cur_err,fabs(_irk_res[j]));
	}

	if(cur_err >= old_err) {

	  if(_output > 2) {
	    cout << "     alpha: " << alpha << "   cur_err: " << cur_err;
	    cout << "  old_err: " << old_err <<  endl;
	  }

	  // do damped step and recalculate jacobian
	  _recalc_jac = true;
	  for(j = 0; j < _n; j++)
	    if(!_x_algebraic[j])
	      _irk_yprime[j] -= alpha*_irk_delta[j];
	    else
	      _irk_y[j] -= alpha*_irk_delta[j];

	  alpha = 0.5*alpha;

	  v_copy(_irk_y,y);
	  for(j = 0; j < _n; j++)
	    if(!_x_algebraic[j]) {
	      y[j] += _h*_a[stage-1][stage-1]*_irk_yprime[j];
	      yprime[j] = _irk_yprime[j];
	    }
	  res(_t+_c[stage-1]*_h, y, yprime,_par,_irk_res);

	}
	else {
	  if(cur_err/old_err > 0.3)
	    _recalc_jac = true;
	  break;
	}
      } while(alpha > 1.0e-2);
    }
    
    if(i > 0 && _output > 2) {
      cout << "     iter: " << i << "  res: " << cur_err;
      cout << "  oldres: " << old_err;
      if(_recalc_jac)
	cout << "   new jacobian";
      cout << endl;
    }
    
    old_err = cur_err;
    
    v_copy(_irk_y,y);
    for(j = 0; j < _n; j++)
      if(!_x_algebraic[j]) {
	y[j] += _h*_a[stage-1][stage-1]*_irk_yprime[j];
	yprime[j] = _irk_yprime[j];
      }

    if(ok) 
      break;
    else if (i == _maxiters-1)
      m_error(E_CONV, 
	      "Omu_IntSDIRK::ode_solve Newton method failed to converge");
    
    // (re)calculate and factorize Jacobian
    if(_recalc_jac) {
      jac(stage,_t+_c[stage-1]*_h, y, yprime,_par,_irk_res);
      if(i == 0) {
	old_err = 0.0;
	// check for convergence
	for(ok = true, j = 0 ; j < _n; j++) {
	  ok = ok  && (fabs(_irk_res[j])<0.1*_atol);
	  old_err = max(old_err,fabs(_irk_res[j]));
	}
	old_err++;
      }

      if(_banded_solver) {
	if(_sparse_solver) {
	  sp_resize(_irk_jac_sp,_n,_n);
	  sp_insert_mat(_irk_jac_sp,0,0,_irk_jac);
	  spLUfactor2(_irk_jac_sp, _ppivot);
	}
	else {
	  bd_resize(_irk_jac_bd,_ml,_mu,_irk_jac->m);
	  mat2band(_irk_jac,_ml,_mu,_irk_jac_bd);
	  bdLUfactor(_irk_jac_bd, _ppivot);
	}
      }
      else
	LUfactor(_irk_jac, _ppivot);
      _recalc_jac = false;
    }
    else
      res(_t+_c[stage-1]*_h, y, yprime,_par,_irk_res);
    
    sv_mlt(-1.0,_irk_res,_irk_res);
    if(_banded_solver) {
      if(_sparse_solver)
	spLUsolve(_irk_jac_sp, _ppivot, _irk_res, _irk_delta);
      else
	bdLUsolve(_irk_jac_bd, _ppivot, _irk_res, _irk_delta);
    }
    else
      LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);

    // do step
    for(j = 0; j < _n; j++)
      if(!_x_algebraic[j])
	_irk_yprime[j] += _irk_delta[j];
      else
	_irk_y[j] += _irk_delta[j];
    
  }

  _newtonsteps += i;

}


//--------------------------------------------------------------------------

void Omu_IntSDIRK::solve_stage_lsqr(int stage, VECP y, VECP yprime)
{

  bool    ok;
  int     i, j;
  double  alpha, cur_err, old_err, min_err;

  v_resize(_vh,_irk_res->dim);
  v_ones(_irk_res);
  v_ones(_irk_delta);
  v_copy(y,_irk_y);
  v_copy(yprime,_irk_yprime);

  cur_err = 1.0e10;
  old_err = cur_err;
  min_err = cur_err;
  ok = false;
  // _recalc_jac = false;
    
  if(_output > 2)
    cout << "   stage: " << stage << endl;;

  // modified Newton's method    
  for(i = 0; i < _maxiters; i++ ) {
    
    // if(stage < 2 && (i % _modnewtonsteps) == 0)
    //   _recalc_jac = true;

    if(i > 0 && (i % _modnewtonsteps) == 0)
      _recalc_jac = true;

    alpha = 0.5;
    
    if(i > 0) {

      do {

	cur_err = 0.0;
	// check for convergence
	for(ok = true, j = 0 ; j < _n; j++) {
	  ok = ok  && 
	    (fabs(_irk_res[j]) < 1.0e-2*(_atol + _rtol*fabs(_irk_y[j])));
	  cur_err = max(cur_err,fabs(_irk_res[j]));
	}

	if(v_norm2(_irk_delta) < 1.0e-2*(_atol + _rtol*v_norm2(_irk_y)))
	  ok = true;

	if(cur_err >= old_err) {

	  if(_output > 2) {
	    cout << "     alpha: " << alpha << "   cur_err: " << cur_err;
	    cout << "  old_err: " << old_err <<  endl;
	  }

	  // do damped step and recalculate jacobian
	  _recalc_jac = true;
	  for(j = 0; j < _n; j++)
	    if(!_x_algebraic[j])
	      _irk_yprime[j] -= alpha*_irk_delta[j];
	    else
	      _irk_y[j] -= alpha*_irk_delta[j];

	  alpha = 0.5*alpha;

	  v_copy(_irk_y,y);
	  for(j = 0; j < _n; j++)
	    if(!_x_algebraic[j]) {
	      y[j] += _h*_a[stage-1][stage-1]*_irk_yprime[j];
	      yprime[j] = _irk_yprime[j];
	    }
	  res(_t+_c[stage-1]*_h, y, yprime,_par,_irk_res);

	}
	else {
	  if(cur_err/old_err > 0.3)
	    _recalc_jac = true;
	  break;
	}
      } while(alpha > 1.0e-2);
    }
    
    if(i > 0 && _output > 2) {
      cout << "lsqr     iter: " << i << "  res: " << cur_err;
      cout << "  oldres: " << old_err;
      cout << "  norm_dx: " <<  v_norm2(_irk_delta);
      cout << "  norm_x: " <<  v_norm2(_irk_res);
      if(_recalc_jac)
	cout << "   new jacobian";
      cout << endl;
    }
    
    old_err = cur_err;
    
    v_copy(_irk_y,y);
    for(j = 0; j < _n; j++)
      if(!_x_algebraic[j]) {
	y[j] += _h*_a[stage-1][stage-1]*_irk_yprime[j];
	yprime[j] = _irk_yprime[j];
      }

    if(ok) 
      break;
    else 
      if(i == _maxiters-1)
	m_error(E_CONV, 
		"Omu_IntSDIRK::ode_solve Newton method failed to converge");
    
    // (re)calculate and factorize Jacobian
    if(_recalc_jac) {
      jac(stage,_t+_c[stage-1]*_h, y, yprime,_par,_irk_res);
      if(i == 0) {
	old_err = 0.0;
	// check for convergence
	for(ok = true, j = 0 ; j < _n; j++) {
	  ok = ok  && (fabs(_irk_res[j])<0.1*_atol);
	  old_err = max(old_err,fabs(_irk_res[j]));
	}
	old_err++;
      }

      assert((int)_irk_jac->n == _n);

      if(_sparse_solver) {
	sp_resize(_irk_jac_sp,_n,_n);
	sp_resize(_smh,_irk_jac->m,_n);
	sp_resize(_smh1,_n,_irk_jac->m);
	sp_insert_mat(_smh,0,0,_irk_jac);
	sp_transp(_smh,_smh1);
	sp_mmt2_mlt(_smh,_smh1,_irk_jac_sp);
	spBKPfactor(_irk_jac_sp, _ppivot,0.0);
      }
      else {
	m_resize(_irk_jacf,_irk_jac->n,_irk_jac->m);
	m_resize(_Sh,_n,_n);
	m_transp(_irk_jac,_irk_jacf);
	m_mlt(_irk_jacf,_irk_jac,_Sh);
	for(j = 0; j < _n; j++)
	  if(_x_algebraic[j] && _Sh[j][j] == 0.0)
	    _Sh[j][j] = 1.0e-9;
	m_copy(_Sh,_irk_jac);

	if(_banded_solver) {
	  px_resize(_Spivot,_n);
	  bd_resize(_irk_jac_bd,max(_ml,_mu),max(_ml,_mu),_n);
	  mat2band(_irk_jac,max(_ml,_mu),max(_ml,_mu),_irk_jac_bd);
	  bdBKPfactor(_irk_jac_bd, _ppivot, _Spivot);
	}
	else
	  matBKPfactor(_irk_jac, _ppivot);
      }
      _recalc_jac = false;
    }
    else
      res(_t+_c[stage-1]*_h, y, yprime,_par,_irk_res);
    
    sv_mlt(-1.0,_irk_res,_irk_delta);

    if(_sparse_solver) {
      sp_mv_mlt(_smh1,_irk_delta,_vh);
      spBKPsolve(_irk_jac_sp, _ppivot, _vh, _irk_delta);
    }
    else {
      mv_mlt(_irk_jacf,_irk_delta,_vh);
      if(_banded_solver)
	bdBKPsolve(_irk_jac_bd, _ppivot, _Spivot, _vh, _irk_delta);
      else
	matBKPsolve(_irk_jac, _ppivot, _vh, _irk_delta);
    }

    // do step
    for(j = 0; j < _n; j++)
      if(!_x_algebraic[j])
	_irk_yprime[j] += _irk_delta[j];
      else
	_irk_y[j] += _irk_delta[j];

  }

  _newtonsteps += i;

}

//--------------------------------------------------------------------------

void Omu_IntSDIRK::solve_final(VECP y, VECP yprime)
{

  bool    ok;
  int     i, j;
  double  cur_err, old_err;

  v_ones(_irk_res);
  v_copy(y,_irk_y);
  v_copy(yprime,_irk_yprime);

  cur_err = 1.0e10;
  old_err = cur_err;
  ok = false;
  _recalc_jac = false;
    
  if(_output > 2)
    cout << "compute algebraic states and yprime at t+h" << endl;

  // modified Newton's method    
  for(i = 0; i < _maxiters; i++ ) {

    if((i % _modnewtonsteps) == 0)
      _recalc_jac = true;
    
    if(i > 0) {
      cur_err = 0.0;
      // check for convergence
      for(ok = true, j = 0 ; j < _n; j++) {
	ok = ok  && (fabs(_irk_res[j]) < 1.0e-2*_atol);
	cur_err = max(cur_err,fabs(_irk_res[j]));
      }
      
      if(cur_err > old_err) {
	// do damped step and recalculate jacobian
	_recalc_jac = true;
	for(j = 0; j < _n; j++)
	  if(_x_algebraic[j])
	    _irk_y[j] += 0.2*_irk_delta[j];
	  else
	    _irk_yprime[j] += 0.2*_irk_delta[j];
      }
      else {
	// do full step
	for(j = 0; j < _n; j++)
	  if(_x_algebraic[j])
	    _irk_y[j] += _irk_delta[j];
	  else
	    _irk_yprime[j] += _irk_delta[j];

	if(cur_err/old_err > 0.3)
	  _recalc_jac = true;   
      }
    }
    
    if(_output > 2) {
      cout << "   iter: " << i << "  res: " << cur_err;
      cout << "  oldres: " << old_err;
      if(_recalc_jac)
	cout << "   new jacobian";
      cout << endl;
    }
    
    old_err = cur_err;
    
    v_copy(_irk_y,y);
    v_copy(_irk_yprime,yprime);

    if(ok) 
      break;
    else if (i == _maxiters-1)
      m_error(E_CONV, 
	      "Omu_IntSDIRK::ode_solve Newton method failed to converge");
    
    // (re)calculate and factorize Jacobian
    if(_recalc_jac) {
      jac(0,_t+_h, y, yprime,_par,_irk_res);
      if(i == 0) {
	old_err = 0.0;
	// check for convergence
	for(ok = true, j = 0 ; j < _n; j++) {
	  ok = ok  && (fabs(_irk_res[j])<0.1*_atol);
	  old_err = max(old_err,fabs(_irk_res[j]));
	}
	old_err++;
      }

      if(_banded_solver) {
	if(_sparse_solver) {
	  sp_resize(_irk_jac_sp,_n,_n);
	  sp_insert_mat(_irk_jac_sp,0,0,_irk_jac);
	  spLUfactor2(_irk_jac_sp, _ppivot);
	}
	else {
	  bd_resize(_irk_jac_bd,_ml,_mu,_irk_jac->m);
	  mat2band(_irk_jac,_ml,_mu,_irk_jac_bd);
	  bdLUfactor(_irk_jac_bd, _ppivot);
	}
      }
      else
	LUfactor(_irk_jac, _ppivot);
      _recalc_jac = false;

    }
    else
      res(_t+_h, y, yprime,_par,_irk_res);
    
    sv_mlt(-1.0,_irk_res,_irk_res);

    if(_banded_solver) {
      if(_sparse_solver)
	spLUsolve(_irk_jac_sp, _ppivot, _irk_res, _irk_delta);
      else
	bdLUsolve(_irk_jac_bd, _ppivot, _irk_res, _irk_delta);
    }
    else
      LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);

  }

  _newtonsteps += i;

}


//--------------------------------------------------------------------------
// differentiate the whole integration scheme using adol-c

#ifdef OMU_WITH_ADOLC

void Omu_IntSDIRK::sensitivity_adolc()
{

  int i, j, l;

  double  *aindep, *f;

  adoublev au(_nu);
  adoublev ay(_nd+_n);
  adoublev ay0(_nd+_n);
  adoublev ayprime(_nd+_n);
  adoublev aF(_nd+_n);
  adoublev irk_ay(_irk_stages*(_nd+_n));
  adoublev irk_ayprime(_irk_stages*(_nd+_n));
  adoublev irk_aF(_irk_stages*(_nd+_n));

  f = new double[_irk_stages*_n];
  aindep = new double[_irk_stages*_n+_n+_nd+_nu];

  trace_on(_tag,1);

  // initialize variables

  for (i = 0; i < _nd+_n; i++)
    aF[i] = 0.0;

  for (i = 0; i < _irk_stages*(_nd+_n); i++) {
    irk_ay[i] = 0.0;
    irk_ayprime[i] = 0.0;
  }

  for(i = 0; i < _nd; i++) {
    aindep[i] = _par[i];
    ay0[i] <<= aindep[i];
  }
  
  for(i = 0; i < _n; i++) {
    aindep[_nd+i] = _y0[i];
    ay0[_nd+i] <<= aindep[_nd+i];
  }

  for(i = 0; i < _nu; i++) {
    aindep[_nd+_n+i] = _par[_nd+i];
    au[i] <<= aindep[_nd+_n+i];
  }
  
  for(i = 0; i < _irk_stages; i++)
    for(j = 0; j < _n; j++) 
      if(!_x_algebraic[j]) {
	aindep[_nd+(i+1)*_n+_nu+j] = _yprime[i][j];
	irk_ayprime[i*(_nd+_n)+_nd+j] <<= aindep[_nd+(i+1)*_n+_nu+j];
      }
      else {
	aindep[_nd+(i+1)*_n+_nu+j] = _yprime[i][j];
	irk_ay[i*(_nd+_n)+_nd+j] <<= aindep[_nd+(i+1)*_n+_nu+j];
      }	  
  
  for(i = 0; i < _irk_stages; i++) {

    for(j = 0; j < _nd; j++) {
      ay[j] = ay0[j];
      ayprime[j] = 0.0;
    }
    
    for(j = 0; j < _n; j++)
      if(!_x_algebraic[j]) {
     	ay[_nd+j] = ay0[_nd+j];
	for(l = 0; l <= i; l++)
	  ay[_nd+j] += _h*_a[i][l]*irk_ayprime[l*(_nd+_n)+_nd+j];
	ayprime[_nd+j] = irk_ayprime[i*(_nd+_n)+_nd+j];
      }
      else
	ay[_nd+j] = irk_ay[i*(_nd+_n)+_nd+j];

    _sys->continuous(_kk, _t+_c[i]*_h, ay, au, ayprime, aF);
    
    for(j = 0; j < _nd+_n; j++)
      irk_aF[i*(_nd+_n)+j] = aF[j];

  }

  for(i = 0; i < _irk_stages; i++)
    for(j = 0; j < _n; j++)
      irk_aF[i*(_nd+_n)+_nd+j] >>= f[i*_n+j];

  trace_off();

  if(_sens_at_once) {
    
    m_resize(_U,_irk_stages*_n,_irk_stages*_n);
    m_resize(_Z,_irk_stages*_n,(_irk_stages+1)*_n+_nd+_nu);
    m_ident(_U);
    m_zero(_Z);
  
    reverse(_tag, _irk_stages*_n, (_irk_stages+1)*_n+_nd+_nu, 0, 
	    _irk_stages*_n, _U->me, _Z->me);

    m_resize(_SF,_irk_stages*_n,_irk_stages*_n);
    m_resize(_irk_jacf,_irk_stages*_n,_irk_stages*_n);
    m_move(_Z,0,_n+_nd+_nu,_irk_stages*_n,_irk_stages*_n,_SF,0,0);
    m_resize(_Sh,_irk_stages*_n,_n+_nd+_nu);

    m_copy(_SF,_irk_jacf);
    LUfactor(_irk_jacf, _Spivot);
    
    v_resize(_irk_res,_irk_stages*_n);
    v_resize(_irk_delta,_irk_stages*_n);
    
    for(j = 0; j < _nd+_n+_nu; j++) {
      get_col(_Z,j,_irk_res);

      if(v_sum(_irk_res) != 0.0) {
	sv_mlt(-1.0,_irk_res,_irk_res);
	LUsolve(_irk_jacf, _Spivot, _irk_res, _irk_delta);
	set_col(_Sh,j,_irk_delta);
      }
      else
	set_col(_Sh,j,_irk_res);
    }
  }
  else {
    
    // compute sensitivities stagewise
    // treat the tape all at once
    if(_n < _n_splitt_tape_eval) {
      
      m_resize(_U,_irk_stages*_n,_irk_stages*_n);
      m_resize(_Z,_irk_stages*_n,(_irk_stages+1)*_n+_nd+_nu);
      m_ident(_U);
      m_zero(_Z);
      
      reverse(_tag, _irk_stages*_n, (_irk_stages+1)*_n+_nd+_nu, 0, 
	      _irk_stages*_n, _U->me, _Z->me);
      
      m_resize(_SF,_n,_n+_nd+_nu);
      m_resize(_Sxd,_n,_n+_nd+_nu);
      m_resize(_Sxu,_n,_n+_nd+_nu);
      m_resize(_Sh,_irk_stages*_n,_n+_nd+_nu);
      
      for(i = 0; i < _irk_stages; i++) {
	
	m_move(_Z,i*_n,(i+1)*_n+_nd+_nu,_n,_n,_irk_jac,0,0);
	m_move(_Z,i*_n,0,_n,_nd+_n+_nu,_Sxd,0,0);

	if(_banded_solver) {
	  bd_resize(_irk_jac_bd,_ml,_mu,_irk_jac->m);
	  mat2band(_irk_jac,_ml,_mu,_irk_jac_bd);
	  bdLUfactor(_irk_jac_bd, _ppivot);
	}
	else
	  LUfactor(_irk_jac, _ppivot);
	
	for(j = 0; j < _nd+_n+_nu; j++) {
	  get_col(_Sxd,j,_irk_res);
	  sv_mlt(-1.0,_irk_res,_irk_res);
	  if(_banded_solver)
	    bdLUsolve(_irk_jac_bd, _ppivot, _irk_res, _irk_delta);
	  else
	    LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);
	  set_col(_Sxd,j,_irk_delta);
	}
	
	for(j = 0; j < i; j++) {
	  
	  m_move(_Z,i*_n,(j+1)*_n+_nd+_nu,_n,_n,_Sxx,0,0);
	  m_move(_Sh,j*_n,0,_n,_nd+_n+_nu,_SF,0,0);
	  
	  for(l = 0; l < _n; l++) {
	    get_col(_Sxx,l,_irk_res);
	    sv_mlt(-1.0,_irk_res,_irk_res);
	    if(_banded_solver)
	      bdLUsolve(_irk_jac_bd, _ppivot, _irk_res, _irk_delta);
	    else
	      LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);
	    set_col(_Sxx,l,_irk_delta);
	  }
	  
	  m_mlt(_Sxx,_SF,_Sxu);
	  m_add(_Sxd,_Sxu,_Sxd);
	  
	}
	
	m_move(_Sxd,0,0,_n,_n+_nd+_nu,_Sh,i*_n,0);
	
      }
    }
    else {
      
      // compute sensitivities stagewise
      // treat the tape stagewise
      
      m_resize(_U,_n,_irk_stages*_n);
      m_resize(_Z,_n,(_irk_stages+1)*_n+_nd+_nu);
      m_resize(_SF,_n,_n+_nd+_nu);
      m_resize(_Sxd,_n,_n+_nd+_nu);
      m_resize(_Sxu,_n,_n+_nd+_nu);
      m_resize(_Sh,_irk_stages*_n,_n+_nd+_nu);
      
      m_zero(_U);
      
      // treat the tape stagewise
      for(i = 0; i < _irk_stages; i++) {
	
	for(j = 0; j < _n; j++)
	  _U[j][i*_n+j] = 1.0;
	
	reverse(_tag, _irk_stages*_n, (_irk_stages+1)*_n+_nd+_nu, 0, 
		_n, _U->me, _Z->me);

	m_move(_Z,0,(i+1)*_n+_nd+_nu,_n,_n,_irk_jac,0,0);
	m_move(_Z,0,0,_n,_nd+_n+_nu,_Sxd,0,0);

	if(_banded_solver) {
	  bd_resize(_irk_jac_bd,_ml,_mu,_irk_jac->m);
	  mat2band(_irk_jac,_ml,_mu,_irk_jac_bd);
	  bdLUfactor(_irk_jac_bd, _ppivot);
	}
	else
	  LUfactor(_irk_jac, _ppivot);

	for(j = 0; j < _nd+_n+_nu; j++) {
	  get_col(_Sxd,j,_irk_res);
	  sv_mlt(-1.0,_irk_res,_irk_res);
	  if(_banded_solver)
	    bdLUsolve(_irk_jac_bd, _ppivot, _irk_res, _irk_delta);
	  else
	    LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);
	  set_col(_Sxd,j,_irk_delta);
	}
	
	for(j = 0; j < i; j++) {
	  
	  m_move(_Z,0,(j+1)*_n+_nd+_nu,_n,_n,_Sxx,0,0);
	  m_move(_Sh,j*_n,0,_n,_nd+_n+_nu,_SF,0,0);
	  
	  for(l = 0; l < _n; l++) {
	    get_col(_Sxx,l,_irk_res);
	    sv_mlt(-1.0,_irk_res,_irk_res);
	    if(_banded_solver)
	      bdLUsolve(_irk_jac_bd, _ppivot, _irk_res, _irk_delta);
	    else
	      LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);
	    set_col(_Sxx,l,_irk_delta);
	  }
	  
	  m_mlt(_Sxx,_SF,_Sxu);
	  m_add(_Sxd,_Sxu,_Sxd);
	  
	}
	
	m_move(_Sxd,0,0,_n,_n+_nd+_nu,_Sh,i*_n,0);
	
	for(j = 0; j < _n; j++)
	  _U[j][i*_n+j] = 0.0;
	
      }
    } 
  }   

  m_resize(_Sxd,_n,_nd);
  m_resize(_Sxu,_n,_nu);
  v_resize(_irk_res,_n);
  v_resize(_irk_delta,_n);

  // Sxd
  m_zero(_Sxd);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _nd; j++)
	  _Sxd[i][j] += _h*_b[l]*_Sh[l*_n+i][j];
    else
      for(j = 0; j < _nd; j++)
	_Sxd[i][j] = _Sh[(_irk_stages-1)*_n+i][j];
  
  // Sxx
  m_zero(_Sxx);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _n; j++)
	  _Sxx[i][j] += _h*_b[l]*_Sh[l*_n+i][_nd+j];
    else
      for(j = 0; j < _n; j++)
	_Sxx[i][j] = _Sh[(_irk_stages-1)*_n+i][_nd+j];
  
  // Sxu
  m_zero(_Sxu);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _nu; j++)
	  _Sxu[i][j] += _h*_b[l]*_Sh[l*_n+i][_nd+_n+j];
    else
      for(j = 0; j < _nu; j++)
	_Sxu[i][j] = _Sh[(_irk_stages-1)*_n+i][_nd+_n+j];
  
  m_resize(_Sh,_n,_nd);
  m_mlt(_Sxx,_Sxd0,_Sh);

  for(j = 0; j < _n; j++)
    for(l = 0; l < _nd; l++)
      if(!_x_algebraic[j])
	_Sxd0[j][l] = _Sxd0[j][l] + _Sxd[j][l] +_Sh[j][l];
      else
	_Sxd0[j][l] = _Sxd[j][l] +_Sh[j][l];

  m_resize(_Sh,_n,_nu);
  m_mlt(_Sxx,_Sxu0,_Sh);

  for(j = 0; j < _n; j++)
    for(l = 0; l < _nu; l++)
      if(!_x_algebraic[j])
	_Sxu0[j][l] = _Sxu0[j][l] + _Sxu[j][l] +_Sh[j][l];
      else
	_Sxu0[j][l] = _Sxu[j][l] +_Sh[j][l];
  
  m_resize(_Sh,_n,_n);
  m_mlt(_Sxx,_Sxx0,_Sh);
  for(j = 0; j < _n; j++)
    for(l = 0; l < _n; l++)
      if(!_x_algebraic[j])
	_Sxx0[j][l] = _Sxx0[j][l] + _Sh[j][l];
      else
	_Sxx0[j][l] = _Sh[j][l];
  
  _sen_evals++;

  delete[] aindep;
  delete[] f;

}

#endif

//--------------------------------------------------------------------------

void Omu_IntSDIRK::sensitivity()
{

  int    i, j, k, l;

  Omu_StateVec &cx = *_cx_ptr;
  Omu_DependentVec &cF = *_cF_ptr;
  
  m_resize(_SF,_n,_n+_nd+_nu);
  m_resize(_Sxd,_n,_n+_nd+_nu);
  m_resize(_Sxu,_n,_n+_nd+_nu);
  m_resize(_Sh,_irk_stages*_n,_n+_nd+_nu);

  for (i = 0; i < _nd; i++) {
    cx[i] = _par[i];
    _cxp[i] = 0.0;
  }

  for (i = 0; i < _nu; i++) {
    _cu[i] = _par[_nd + i];
  }

  // evaluate residuals and Jacobians
  bool is_required_J_bak = cF.is_required_J();
  cF.set_required_J(true);
  
  for(i = 0; i < _irk_stages; i++) {

    for(j = 0; j < _n; j++)
      if(!_x_algebraic[j]) {
     	cx[_nd+j] = _y0[j];
	for(l = 0; l <= i; l++)
	  cx[_nd+j] += _h*_a[i][l]*_yprime[l][j];
	_cxp[_nd+j] = _yprime[i][j];
      }
      else
	cx[_nd+j] = _yprime[i][j];

    _sys->continuous(_kk, _t+_c[i]*_h, cx, _cu, _cxp, cF);

    for (j = 0; j < _n; j++)
      for (l = 0; l < _n; l++)
	if(!_x_algebraic[l]) {
	  _irk_jac[j][l] = cF.Jdx[_nd+ j][_nd+ l] + 
	    _h*_a[i][i]*cF.Jx[_nd+ j][_nd+ l];
	  _Sxd[j][_nd+l] = cF.Jx[_nd+ j][_nd+ l];
	}
	else
	  _irk_jac[j][l] = cF.Jx[_nd+ j][_nd+ l];
	
    m_move(cF.Jx,_nd,0,_n,_nd,_Sxd,0,0);
    m_move(cF.Ju,_nd,0,_n,_nu,_Sxd,0,_nd+_n);

    if(_banded_solver) {
      if(_sparse_solver) {
	sp_resize(_irk_jac_sp,_n,_n);
	sp_insert_mat(_irk_jac_sp,0,0,_irk_jac);
	spLUfactor2(_irk_jac_sp, _ppivot);
      }
      else {
	bd_resize(_irk_jac_bd,_ml,_mu,_irk_jac->m);
	mat2band(_irk_jac,_ml,_mu,_irk_jac_bd);
	bdLUfactor(_irk_jac_bd, _ppivot);
      }
    }
    else
      LUfactor(_irk_jac, _ppivot);
    
    for(j = 0; j < _nd+_n+_nu; j++) {
      get_col(_Sxd,j,_irk_res);
      sv_mlt(-1.0,_irk_res,_irk_res);
      if(_banded_solver) {
	if(_sparse_solver)
	  spLUsolve(_irk_jac_sp, _ppivot, _irk_res, _irk_delta);
	else
	  bdLUsolve(_irk_jac_bd, _ppivot, _irk_res, _irk_delta);
      }
      else
	LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);
      set_col(_Sxd,j,_irk_delta);
    }
    
    for(j = 0; j < i; j++) {

      for (k = 0; k < _n; k++)
	for (l = 0; l < _n; l++)
	  if(!_x_algebraic[l])
	    _Sxx[k][l] = _h*_a[i][j]*cF.Jx[_nd+ k][_nd+ l];

      m_move(_Sh,j*_n,0,_n,_nd+_n+_nu,_SF,0,0);
      
      for(l = 0; l < _n; l++) {
	get_col(_Sxx,l,_irk_res);
	sv_mlt(-1.0,_irk_res,_irk_res);
	if(_banded_solver) {
	  if(_sparse_solver)
	    spLUsolve(_irk_jac_sp, _ppivot, _irk_res, _irk_delta);
	  else
	    bdLUsolve(_irk_jac_bd, _ppivot, _irk_res, _irk_delta);
	}
	else
	  LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);
	set_col(_Sxx,l,_irk_delta);
      }
      
      m_mlt(_Sxx,_SF,_Sxu);
      m_add(_Sxd,_Sxu,_Sxd);
      
    }
    
    m_move(_Sxd,0,0,_n,_n+_nd+_nu,_Sh,i*_n,0);
    
  }

  cF.set_required_J(is_required_J_bak);

  m_resize(_Sxd,_n,_nd);
  m_resize(_Sxu,_n,_nu);
  v_resize(_irk_res,_n);
  v_resize(_irk_delta,_n);

  // Sxd
  m_zero(_Sxd);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _nd; j++)
	  _Sxd[i][j] += _h*_b[l]*_Sh[l*_n+i][j];
    else
      for(j = 0; j < _nd; j++)
	_Sxd[i][j] = _Sh[(_irk_stages-1)*_n+i][j];
  
  // Sxx
  m_zero(_Sxx);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _n; j++)
	  _Sxx[i][j] += _h*_b[l]*_Sh[l*_n+i][_nd+j];
    else
      for(j = 0; j < _n; j++)
	_Sxx[i][j] = _Sh[(_irk_stages-1)*_n+i][_nd+j];
  
  // Sxu
  m_zero(_Sxu);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _nu; j++)
	  _Sxu[i][j] += _h*_b[l]*_Sh[l*_n+i][_nd+_n+j];
    else
      for(j = 0; j < _nu; j++)
	_Sxu[i][j] = _Sh[(_irk_stages-1)*_n+i][_nd+_n+j];
  
  m_resize(_Sh,_n,_nd);
  m_mlt(_Sxx,_Sxd0,_Sh);

  for(j = 0; j < _n; j++)
    for(l = 0; l < _nd; l++)
      if(!_x_algebraic[j])
	_Sxd0[j][l] = _Sxd0[j][l] + _Sxd[j][l] +_Sh[j][l];
      else
	_Sxd0[j][l] = _Sxd[j][l] +_Sh[j][l];

  m_resize(_Sh,_n,_nu);
  m_mlt(_Sxx,_Sxu0,_Sh);

  for(j = 0; j < _n; j++)
    for(l = 0; l < _nu; l++)
      if(!_x_algebraic[j])
	_Sxu0[j][l] = _Sxu0[j][l] + _Sxu[j][l] +_Sh[j][l];
      else
	_Sxu0[j][l] = _Sxu[j][l] +_Sh[j][l];
  
  m_resize(_Sh,_n,_n);
  m_mlt(_Sxx,_Sxx0,_Sh);
  for(j = 0; j < _n; j++)
    for(l = 0; l < _n; l++)
      if(!_x_algebraic[j])
	_Sxx0[j][l] = _Sxx0[j][l] + _Sh[j][l];
      else
	_Sxx0[j][l] = _Sh[j][l];
  
  _sen_evals++;

}

//--------------------------------------------------------------------------
// differentiate the whole integration scheme using adol-c

#ifdef OMU_WITH_ADOLC

void Omu_IntSDIRK::sensitivity_lsqr_adolc()
{

  int i, j, l;

  double  *aindep, *f;

  adoublev au(_nu);
  adoublev ay(_nd+_n);
  adoublev ay0(_nd+_n);
  adoublev ayprime(_nd+_n);
  adoublev aF(_nd+_n);
  adoublev irk_ay(_irk_stages*(_nd+_n));
  adoublev irk_ayprime(_irk_stages*(_nd+_n));
  adoublev irk_aF(_irk_stages*(_nd+_n));

  f = new double[_irk_stages*_n];
  aindep = new double[_irk_stages*_n+_n+_nd+_nu];

  trace_on(_tag,1);

  // initialize variables

  for (i = 0; i < _nd+_n; i++)
    aF[i] = 0.0;

  for (i = 0; i < _irk_stages*(_nd+_n); i++) {
    irk_ay[i] = 0.0;
    irk_ayprime[i] = 0.0;
  }

  for(i = 0; i < _nd; i++) {
    aindep[i] = _par[i];
    ay0[i] <<= aindep[i];
  }
  
  for(i = 0; i < _n; i++) {
    aindep[_nd+i] = _y0[i];
    ay0[_nd+i] <<= aindep[_nd+i];
  }

  for(i = 0; i < _nu; i++) {
    aindep[_nd+_n+i] = _par[_nd+i];
    au[i] <<= aindep[_nd+_n+i];
  }
  
  for(i = 0; i < _irk_stages; i++)
    for(j = 0; j < _n; j++) 
      if(!_x_algebraic[j]) {
	aindep[_nd+(i+1)*_n+_nu+j] = _yprime[i][j];
	irk_ayprime[i*(_nd+_n)+_nd+j] <<= aindep[_nd+(i+1)*_n+_nu+j];
      }
      else {
	aindep[_nd+(i+1)*_n+_nu+j] = _yprime[i][j];
	irk_ay[i*(_nd+_n)+_nd+j] <<= aindep[_nd+(i+1)*_n+_nu+j];
      }	  
  
  for(i = 0; i < _irk_stages; i++) {

    for(j = 0; j < _nd; j++) {
      ay[j] = ay0[j];
      ayprime[j] = 0.0;
    }
    
    for(j = 0; j < _n; j++)
      if(!_x_algebraic[j]) {
     	ay[_nd+j] = ay0[_nd+j];
	for(l = 0; l <= i; l++)
	  ay[_nd+j] += _h*_a[i][l]*irk_ayprime[l*(_nd+_n)+_nd+j];
	ayprime[_nd+j] = irk_ayprime[i*(_nd+_n)+_nd+j];
      }
      else
	ay[_nd+j] = irk_ay[i*(_nd+_n)+_nd+j];

    _sys->continuous(_kk, _t+_c[i]*_h, ay, au, ayprime, aF);
    
    for(j = 0; j < _nd+_n; j++)
      irk_aF[i*(_nd+_n)+j] = aF[j];

  }

  for(i = 0; i < _irk_stages; i++)
    for(j = 0; j < _n; j++)
      irk_aF[i*(_nd+_n)+_nd+j] >>= f[i*_n+j];

  trace_off();

  // compute sensitivities stagewise
  // treat the tape stagewise
  
  m_resize(_U,_n,_irk_stages*_n);
  m_resize(_Z,_n,(_irk_stages+1)*_n+_nd+_nu);
  m_resize(_SF,_n,_n+_nd+_nu);
  m_resize(_Sxd,_n,_n+_nd+_nu);
  m_resize(_Sxu,_n,_n+_nd+_nu);
  m_resize(_Sh,_irk_stages*_n,_n+_nd+_nu);
  
  m_zero(_U);
  
  // treat the tape stagewise
  for(i = 0; i < _irk_stages; i++) {
    
    for(j = 0; j < _n; j++)
      _U[j][i*_n+j] = 1.0;
    
    reverse(_tag, _irk_stages*_n, (_irk_stages+1)*_n+_nd+_nu, 0, 
	    _n, _U->me, _Z->me);
    
    m_move(_Z,0,(i+1)*_n+_nd+_nu,_n,_n,_irk_jac,0,0);
    m_move(_Z,0,0,_n,_nd+_n+_nu,_Sxd,0,0);
    
    m_resize(_irk_jacf,_irk_jac->n,_irk_jac->m);
    m_resize(_cFxh,_n,_n);
    m_transp(_irk_jac,_irk_jacf);
    m_mlt(_irk_jacf,_irk_jac,_cFxh);
    for(j = 0; j < _n; j++)
      if(_x_algebraic[j] && _cFxh[j][j] == 0.0)
	_cFxh[j][j] = 1.0e-15;

    if(_banded_solver) {
      px_resize(_Spivot,_n);
      bd_resize(_irk_jac_bd,max(_ml,_mu),max(_ml,_mu),_cFxh->m);
      mat2band(_cFxh,max(_ml,_mu),max(_ml,_mu),_irk_jac_bd);
      bdBKPfactor(_irk_jac_bd, _ppivot, _Spivot);
    }
    else
      matBKPfactor(_cFxh, _ppivot);

    for(j = 0; j < _nd+_n+_nu; j++) {
      get_col(_Sxd,j,_irk_res);
      sv_mlt(-1.0,_irk_res,_irk_delta);
      mv_mlt(_irk_jacf,_irk_delta,_irk_res);
      if(_banded_solver)
	bdBKPsolve(_irk_jac_bd, _ppivot, _Spivot, _irk_res, _irk_delta);
      else
	matBKPsolve(_cFxh, _ppivot, _irk_res, _irk_delta);
      set_col(_Sxd,j,_irk_delta);
    }
    
    for(j = 0; j < i; j++) {
      
      m_move(_Z,0,(j+1)*_n+_nd+_nu,_n,_n,_Sxx,0,0);
      m_move(_Sh,j*_n,0,_n,_nd+_n+_nu,_SF,0,0);
      
      for(l = 0; l < _n; l++) {
	get_col(_Sxx,l,_irk_res);
	sv_mlt(-1.0,_irk_res,_irk_delta);
	mv_mlt(_irk_jacf,_irk_delta,_irk_res);
	if(_banded_solver)
	  bdBKPsolve(_irk_jac_bd, _ppivot, _Spivot, _irk_res, _irk_delta);
	else
	  matBKPsolve(_cFxh, _ppivot, _irk_res, _irk_delta);
	set_col(_Sxx,l,_irk_delta);
      }
      
      m_mlt(_Sxx,_SF,_Sxu);
      m_add(_Sxd,_Sxu,_Sxd);
      
    }
    
    m_move(_Sxd,0,0,_n,_n+_nd+_nu,_Sh,i*_n,0);
    
    for(j = 0; j < _n; j++)
      _U[j][i*_n+j] = 0.0;
    
  }

  m_resize(_Sxd,_n,_nd);
  m_resize(_Sxu,_n,_nu);
  v_resize(_irk_res,_n);
  v_resize(_irk_delta,_n);

  // Sxd
  m_zero(_Sxd);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _nd; j++)
	  _Sxd[i][j] += _h*_b[l]*_Sh[l*_n+i][j];
    else
      for(j = 0; j < _nd; j++)
	_Sxd[i][j] = _Sh[(_irk_stages-1)*_n+i][j];
  
  // Sxx
  m_zero(_Sxx);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _n; j++)
	  _Sxx[i][j] += _h*_b[l]*_Sh[l*_n+i][_nd+j];
    else
      for(j = 0; j < _n; j++)
	_Sxx[i][j] = _Sh[(_irk_stages-1)*_n+i][_nd+j];
  
  // Sxu
  m_zero(_Sxu);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _nu; j++)
	  _Sxu[i][j] += _h*_b[l]*_Sh[l*_n+i][_nd+_n+j];
    else
      for(j = 0; j < _nu; j++)
	_Sxu[i][j] = _Sh[(_irk_stages-1)*_n+i][_nd+_n+j];
  
  m_resize(_Sh,_n,_nd);
  m_mlt(_Sxx,_Sxd0,_Sh);

  for(j = 0; j < _n; j++)
    for(l = 0; l < _nd; l++)
      if(!_x_algebraic[j])
	_Sxd0[j][l] = _Sxd0[j][l] + _Sxd[j][l] +_Sh[j][l];
      else
	_Sxd0[j][l] = _Sxd[j][l] +_Sh[j][l];

  m_resize(_Sh,_n,_nu);
  m_mlt(_Sxx,_Sxu0,_Sh);

  for(j = 0; j < _n; j++)
    for(l = 0; l < _nu; l++)
      if(!_x_algebraic[j])
	_Sxu0[j][l] = _Sxu0[j][l] + _Sxu[j][l] +_Sh[j][l];
      else
	_Sxu0[j][l] = _Sxu[j][l] +_Sh[j][l];
  
  m_resize(_Sh,_n,_n);
  m_mlt(_Sxx,_Sxx0,_Sh);
  for(j = 0; j < _n; j++)
    for(l = 0; l < _n; l++)
      if(!_x_algebraic[j])
	_Sxx0[j][l] = _Sxx0[j][l] + _Sh[j][l];
      else
	_Sxx0[j][l] = _Sh[j][l];

  _sen_evals++;

  delete[] aindep;
  delete[] f;
  
}

#endif

//--------------------------------------------------------------------------

void Omu_IntSDIRK::sensitivity_lsqr()
{

  int i, j, k, l;

  Omu_StateVec &cx = *_cx_ptr;
  Omu_DependentVec &cF = *_cF_ptr;

  m_resize(_SF,_n,_n+_nd+_nu);
  m_resize(_Sxd,_n,_n+_nd+_nu);
  m_resize(_Sxu,_n,_n+_nd+_nu);
  m_resize(_Sh,_irk_stages*_n,_n+_nd+_nu);

  if(_nod > 0)
    cout << " nod: " << _nod << "  n: " << _n << endl;

  for (i = 0; i < _nd; i++) {
    cx[i] = _par[i];
    _cxp[i] = 0.0;
  }

  for (i = 0; i < _nu; i++) {
    _cu[i] = _par[_nd + i];
  }

  // evaluate residuals and Jacobians
  bool is_required_J_bak = cF.is_required_J();
  cF.set_required_J(true);
  
  for(i = 0; i < _irk_stages; i++) {

    for(j = 0; j < _n; j++)
      if(!_x_algebraic[j]) {
     	cx[_nd+j] = _y0[j];
	for(l = 0; l <= i; l++)
	  cx[_nd+j] += _h*_a[i][l]*_yprime[l][j];
	_cxp[_nd+j] = _yprime[i][j];
      }
      else
	cx[_nd+j] = _yprime[i][j];

    _sys->continuous(_kk, _t+_c[i]*_h, cx, _cu, _cxp, cF);

    for (j = 0; j < _n; j++)
      for (l = 0; l < _n; l++)
	if(!_x_algebraic[l]) {
	  _irk_jac[j][l] = cF.Jdx[_nd+ j][_nd+ l] +
	    _h*_a[i][i]*cF.Jx[_nd+ j][_nd+ l];
	  _Sxd[j][_nd+l] = cF.Jx[_nd+ j][_nd+ l];
	}
	else
	  _irk_jac[j][l] = cF.Jx[_nd+ j][_nd+ l];

    m_move(cF.Jx,_nd,0,_n,_nd,_Sxd,0,0);
    m_move(cF.Ju,_nd,0,_n,_nu,_Sxd,0,_nd+_n);

    if(_sparse_solver) {
      sp_resize(_irk_jac_sp,_n,_n);
      sp_resize(_smh,_irk_jac->m,_n);
      sp_resize(_smh1,_n,_irk_jac->m);
      sp_insert_mat(_smh,0,0,_irk_jac);
      sp_transp(_smh,_smh1);
      sp_mmt2_mlt(_smh,_smh1,_irk_jac_sp);
      symsp_insert_symmat(_irk_jac_sp,0,_cFxh);
      spBKPfactor(_irk_jac_sp, _ppivot,0.0);
    }
    else {
      m_resize(_irk_jacf,_irk_jac->n,_irk_jac->m);
      m_resize(_cFxh,_n,_n);
      m_transp(_irk_jac,_irk_jacf);
      m_mlt(_irk_jacf,_irk_jac,_cFxh);
      for(j = 0; j < _n; j++)
	if(_x_algebraic[j] && _cFxh[j][j] == 0.0)
	  _cFxh[j][j] = 1.0e-15;
      if(_banded_solver) {
	px_resize(_Spivot,_n);
	bd_resize(_irk_jac_bd,max(_ml,_mu),max(_ml,_mu),_cFxh->m);
	mat2band(_cFxh,max(_ml,_mu),max(_ml,_mu),_irk_jac_bd);
	bdBKPfactor(_irk_jac_bd, _ppivot, _Spivot);
      }
      else
	matBKPfactor(_cFxh, _ppivot);
    }

    for(j = 0; j < _nd+_n+_nu; j++) {

      get_col(_Sxd,j,_irk_res);
      sv_mlt(-1.0,_irk_res,_irk_delta);

      if(_sparse_solver) {
	sp_mv_mlt(_smh1,_irk_delta,_irk_res);
	spBKPsolve(_irk_jac_sp, _ppivot, _irk_res, _irk_delta);
      }
      else {
	mv_mlt(_irk_jacf,_irk_delta,_irk_res);
	if(_banded_solver)
	  bdBKPsolve(_irk_jac_bd, _ppivot, _Spivot, _irk_res, _irk_delta);
	else
	  matBKPsolve(_cFxh, _ppivot, _irk_res, _irk_delta);
      }
      set_col(_Sxd,j,_irk_delta);
    }
    
    for(j = 0; j < i; j++) {

      m_zero(_Sxx);
      
      for (k = 0; k < _n; k++)
	for (l = 0; l < _n; l++)
	  if(!_x_algebraic[l])
	    _Sxx[k][l] = _h*_a[i][j]*cF.Jx[_nd+ k][_nd+ l];

      m_move(_Sh,j*_n,0,_n,_nd+_n+_nu,_SF,0,0);
      
      for(l = 0; l < _n; l++) {
	get_col(_Sxx,l,_irk_res);
	sv_mlt(-1.0,_irk_res,_irk_delta);
	if(_sparse_solver) {
	  sp_mv_mlt(_smh1,_irk_delta,_irk_res);
	  spBKPsolve(_irk_jac_sp, _ppivot, _irk_res, _irk_delta);
	}
	else {
	  mv_mlt(_irk_jacf,_irk_delta,_irk_res);
	  if(_banded_solver)
	    bdBKPsolve(_irk_jac_bd, _ppivot, _Spivot, _irk_res, _irk_delta);
	  else
	    matBKPsolve(_cFxh, _ppivot, _irk_res, _irk_delta);
	}
	set_col(_Sxx,l,_irk_delta);
      }
      
      m_mlt(_Sxx,_SF,_Sxu);
      m_add(_Sxd,_Sxu,_Sxd);
      
    }
    
    m_move(_Sxd,0,0,_n,_n+_nd+_nu,_Sh,i*_n,0);
    
  }

  cF.set_required_J(is_required_J_bak);

  m_resize(_Sxd,_n,_nd);
  m_resize(_Sxu,_n,_nu);
  v_resize(_irk_res,_n);
  v_resize(_irk_delta,_n);

  // Sxd
  m_zero(_Sxd);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _nd; j++)
	  _Sxd[i][j] += _h*_b[l]*_Sh[l*_n+i][j];
    else
      for(j = 0; j < _nd; j++)
	_Sxd[i][j] = _Sh[(_irk_stages-1)*_n+i][j];
  
  // Sxx
  m_zero(_Sxx);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _n; j++)
	  _Sxx[i][j] += _h*_b[l]*_Sh[l*_n+i][_nd+j];
    else
      for(j = 0; j < _n; j++)
	_Sxx[i][j] = _Sh[(_irk_stages-1)*_n+i][_nd+j];
  
  // Sxu
  m_zero(_Sxu);
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      for(l = 0; l < _irk_stages; l++)
	for(j = 0; j < _nu; j++)
	  _Sxu[i][j] += _h*_b[l]*_Sh[l*_n+i][_nd+_n+j];
    else
      for(j = 0; j < _nu; j++)
	_Sxu[i][j] = _Sh[(_irk_stages-1)*_n+i][_nd+_n+j];
  
  m_resize(_Sh,_n,_nd);
  m_mlt(_Sxx,_Sxd0,_Sh);

  for(j = 0; j < _n; j++)
    for(l = 0; l < _nd; l++)
      if(!_x_algebraic[j])
	_Sxd0[j][l] = _Sxd0[j][l] + _Sxd[j][l] +_Sh[j][l];
      else
	_Sxd0[j][l] = _Sxd[j][l] +_Sh[j][l];

  m_resize(_Sh,_n,_nu);
  m_mlt(_Sxx,_Sxu0,_Sh);

  for(j = 0; j < _n; j++)
    for(l = 0; l < _nu; l++)
      if(!_x_algebraic[j])
	_Sxu0[j][l] = _Sxu0[j][l] + _Sxu[j][l] +_Sh[j][l];
      else
	_Sxu0[j][l] = _Sxu[j][l] +_Sh[j][l];
  
  m_resize(_Sh,_n,_n);
  m_mlt(_Sxx,_Sxx0,_Sh);
  for(j = 0; j < _n; j++)
    for(l = 0; l < _n; l++)
      if(!_x_algebraic[j])
	_Sxx0[j][l] = _Sxx0[j][l] + _Sh[j][l];
      else
	_Sxx0[j][l] = _Sh[j][l];

  _sen_evals++;
  
}


//--------------------------------------------------------------------------
void Omu_IntSDIRK::res(double t, const VECP  y, const VECP yprime, 
		       const VECP par, VECP res)
{

  int i;

  Omu_StateVec &cx = *_cx_ptr;
  Omu_DependentVec &cF = *_cF_ptr;

  // form vectors of independent variables

  for (i = 0; i < _nd; i++) {
    cx[i] = par[i];
    _cxp[i] = 0.0;
  }
  for (i = 0; i < _n; i++) {
    cx[_nd + i] = y[i];
    _cxp[_nd + i] = yprime[i];
  }
  for (i = 0; i < _nu; i++) {
    _cu[i] = par[_nd + i];
  }

  // evaluate residual

  _sys->continuous(_kk, t, cx, _cu, _cxp, cF);

  for (i = 0; i < _n; i++)
    res[i] = cF[_nd+i];

  _res_evals++;

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::jac(int stage, double t, const VECP y, const VECP yprime,
		       const VECP par, VECP res)
{

  int    i, j;

  Omu_StateVec &cx = *_cx_ptr;
  Omu_DependentVec &cF = *_cF_ptr;

  m_zero(_irk_jac);

  // form vectors of independent variables

  for (i = 0; i < _nd; i++) {
    cx[i] = par[i];
    _cxp[i] = 0.0;
  }
  for (i = 0; i < _n; i++) {
    cx[_nd + i] = y[i];
    _cxp[_nd + i] = yprime[i];
  }
  for (i = 0; i < _nu; i++) {
    _cu[i] = par[_nd + i];
  }

  // evaluate residuals and Jacobians
  bool is_required_J_bak = cF.is_required_J();
  cF.set_required_J(true);

  _sys->continuous(_kk, t, cx, _cu, _cxp, cF);

  cF.set_required_J(is_required_J_bak);

  for (i = 0; i < _n; i++)
    res[i] = cF[_nd+i];

  if(stage > 0) {

    v_set(_cFh,_h*_a[stage-1][stage-1]);
      
    if(_na > 0)    
      for(i = 0; i < _n; i++)
	if(_x_algebraic[i]) 
	  _cFh[i] = 1.0;
    
    for(i = 0; i < _n; i++)
      for(j = 0; j < _n; j++)
	_irk_jac[i][j] = cF.Jdx[_nd+i][_nd+j] + _cFh[j]*cF.Jx[_nd+i][_nd+j];
  }
  else {
    
    v_zero(_cFh);

    // yprime and algebraic states
    if(stage == 0) {

      if(_na > 0)    
	for(i = 0; i < _n; i++)
	  if(_x_algebraic[i]) 
	    _cFh[i] = 1.0;
      
      for(i = 0; i < _n; i++)
	for(j = 0; j < _n; j++)
	  _irk_jac[i][j] = cF.Jdx[_nd+i][_nd+j] + _cFh[j]*cF.Jx[_nd+i][_nd+j];
    }

    // only algebraic states
    if(stage == -1) {
      if(_na > 0)    
	for(i = 0; i < _n; i++)
	  if(_x_algebraic[i]) 
	    _cFh[i] = 1.0;
      for(i = 0; i < _n; i++)
	for(j = 0; j < _n; j++)
	  _irk_jac[i][j] = _cFh[j]*cF.Jx[_nd+i][_nd+j];
      for(i = 0; i < _n; i++)
	_irk_jac[i][i] += 1.0e-8;
    }

    // only yprime
    if(stage == -2) {
      if(_na > 0)    
	for(i = 0; i < _n; i++)
	  if(_x_algebraic[i]) 
	    _cFh[i] = 1.0;
      for(i = 0; i < _n; i++)
	for(j = 0; j < _n; j++)
	  _irk_jac[i][j] = cF.Jdx[_nd+i][_nd+j];
      for(i = 0; i < _n; i++)
	_irk_jac[i][i] += 1.0e-8;
    }

    // error evaluation
    if(stage == -3) {

      v_set(_cFh,_h*_gamma0);

      if(_na > 0)    
	for(i = 0; i < _n; i++)
	  if(_x_algebraic[i]) 
	    _cFh[i] = 0.0;
      
      for(i = 0; i < _n; i++)
	if(!_x_algebraic[i]) 
	  for(j = 0; j < _n; j++)
	    _irk_jac[i][j] = cF.Jdx[_nd+i][_nd+j] 
		+ _cFh[j]*cF.Jx[_nd+i][_nd+j];

      if(_na > 0)    
	for(i = 0; i < _n; i++)
	  if(_x_algebraic[i]) 
	    _irk_jac[i][i] = 1.0-8;
    }

  }
    
  _jac_evals++;

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::init_yprime(int k, double t, const Omu_StateVec &y, 
			       const Omu_Vector &u, VECP yprime)
{

  bool    ok;
  int     i, j;
  double  cur_err, old_err;
  
  Omu_DependentVec &cF = *_cF_ptr;

  cur_err = 0.0;
  old_err = 0.0;

  m_resize(_Sh,_n,_n);

  for(i = 0; i < _nd; i++)
    _cxp[i] = 0.0;
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      _cxp[_nd+i] = yprime[i];
    else
      _cxp[_nd+i] = 0.0;

  if(_output == -1) {
    cout << "initial values for yprime at time: " << _t << endl; 
    cout << "newton method" << endl;
  };

  // modified Newton's method    
  for(i = 0; i < _maxiters; i++ ) {

    _recalc_jac = false;
    ok = false;

    if((i % _modnewtonsteps) == 0)
      _recalc_jac = true;

    if(i > 0) {
      cur_err = 0.0;
      // check for convergence
      for(ok = true, j = 0 ; j < _n; j++)
	if(!_x_algebraic[j]) {
	  ok = ok  && (fabs(cF[_nd+j])<0.1*_atol);
	  cur_err = max(cur_err,fabs(cF[_nd+j]));
	  cur_err = max(cur_err,fabs(cF[_nd+j]));
	}

      if(cur_err > old_err) {
	// damped step and recalculate jacobian
	_recalc_jac = true;
	for(j = 0; j < _n; j++)
	  if(!_x_algebraic[j]) 
	    _cxp[_nd+j] += 0.2*_irk_delta[j];
      }
      else {
	// full step
	for(j = 0; j < _n; j++)
	  if(!_x_algebraic[j])
	    _cxp[_nd+j] += _irk_delta[j];
	if(cur_err/old_err > 0.3)
	  _recalc_jac = true;  
      }
      old_err = cur_err;
    }

    if(_output == -1)
      cout << "   iter: " << i << "  res: " << cur_err << endl;

    if(ok) 
      break;
    else if (i == _maxiters-1)
      m_error(E_CONV, 
	      "Omu_IntSDIRK::init_yprime Newton method failed to converge");
      

    // (re)calculate and factorize Jacobian
    if(_recalc_jac) {

      // evaluate residuals and Jacobians
      bool is_required_J_bak = cF.is_required_J();
      cF.set_required_J(true);

      _sys->continuous(k, t, y, u, _cxp, cF);

      cF.set_required_J(is_required_J_bak);

      if(i == 0) {
	old_err = 0.0;
	// check for convergence
	for(ok = true, j = 0 ; j < _n; j++)
	  if(!_x_algebraic[j]) {
	    ok = ok  && (fabs(cF[_nd+j])<0.1*_atol);
	    old_err = max(old_err,fabs(cF[_nd+j]));
	  }
	old_err++;
      }

      m_move(cF.Jdx,_nd,_nd,_n,_n,_Sh,0,0);
      for(j = 0; j < _n; j++)
	if(_x_algebraic[j]) 
	  _Sh[j][j] = 1.0;
      if(_banded_solver) {
	bd_resize(_irk_jac_bd,_ml,_mu,_Sh->m);
	mat2band(_Sh,_ml,_mu,_irk_jac_bd);
	bdLUfactor(_irk_jac_bd, _ppivot);
      }
      else
	LUfactor(_Sh, _ppivot);
    }
    else
      _sys->continuous(k, t, y, u, _cxp, cF);

    v_move(cF,_nd,_n,_irk_res,0);
    sv_mlt(-1.0,_irk_res,_irk_res);
    if(_banded_solver)
      bdLUsolve(_irk_jac_bd, _ppivot, _irk_res, _irk_delta);
    else
      LUsolve(_Sh, _ppivot, _irk_res, _irk_delta);
  }

  for(j = 0; j < _n; j++)
    yprime[j] = _cxp[_nd+j];

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::mat2bandf(const MATP mfull, int ml, int mu, MATP mband)
{

  int i, i1, i2, j, k, lda, m, n;

  lda = 2*ml+mu+1;
  m = ml+mu+1;
  n = mfull->n;

  m_resize(mband,n,lda);

  for(j = 0; j < n; j++) {

    i1 = max(0,j-mu);
    i2 = min(n-1,j+ml);
    
    for(i = i1; i < i2+1; i++) {
      k = i-j+m-1;
      mband[j][k] = mfull[i][j];
    }
  }
}
       












