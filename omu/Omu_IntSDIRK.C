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
#include <iostream>

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
  _xc_ptr = NULL;
  _Fc_ptr = NULL;

  _recalc_jac = true;
  _stiffly_accurate = false;
  _ode = false;
  _init_method = false; 

  _banded = false;
  _sparse_solver = false;

  _kappa = 0.1;

  _output = 0;

  _jac_sbw = -1;

  _x_algebraic = iv_get(1);

  _ppivot= px_get(1);
  _Spivot= px_get(1);

  _Fch = v_get(1);

  _senpar = v_get(1);
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

  _Fcxh = m_get(1,1);
  _Sh   = m_get(1,1);
  _Sh1   = m_get(1,1);
  _Sh2   = m_get(1,1);
  _Sh3   = m_get(1,1);

  _irk_jac = m_get(1,1);
  _irk_jacf = m_get(1,1);
  _irk_jac_bd = bd_get(0,0,1);

  _irk_jac_sp = SMNULL;
  _smh = SMNULL;
  _smh1 = SMNULL;

  _ml = 0;
  _mu = 0;

  _nsteps = 0;
  _modnewtonsteps = 7;
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

  _ifList.append(new If_Bool("prg_int_sparsol", &_sparse_solver));
  _ifList.append(new If_Bool("prg_int_banded", &_banded));
  _ifList.append(new If_Int("prg_int_out", &_output));
  _ifList.append(new If_Int("prg_int_method", &_method));
  _ifList.append(new If_Int("prg_int_nsteps", &_nsteps));
  _ifList.append(new If_Int("prg_int_modnewtonsteps", &_modnewtonsteps));
  _ifList.append(new If_Int("prg_int_maxiters", &_maxiters));
  _ifList.append(new If_Real("prg_int_hinit", &_hinit));
  _ifList.append(new If_Real("prg_int_kappa", &_kappa));

}

//--------------------------------------------------------------------------

Omu_IntSDIRK::~Omu_IntSDIRK()
{

  IV_FREE(_x_algebraic);

  PX_FREE(_ppivot);
  PX_FREE(_Spivot);

  V_FREE(_Fch);
  V_FREE(_senpar);
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
  M_FREE(_Fcxh);

  M_FREE(_Sh);
  M_FREE(_Sh1);
  M_FREE(_Sh2);
  M_FREE(_Sh3);

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

  _gamma0 = 1.0/4.0;

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

  _gamma0 = 1.0/4.0;

  // simple linear interpolation

  m_resize(_a_pred,_irk_stages,_irk_stages+1);
  m_zero(_a_pred);

  _a_pred[0][0] = 2.0;
  _a_pred[0][2] = -1.0;
  _a_pred[1][0] = -2.0;
  _a_pred[1][1] = 3.0;
  _a_pred[2][1] = 0.88;
  _a_pred[2][2] = 0.44;
  _a_pred[3][3] = 1.0;
  _a_pred[4][1] = _b_err[0];
  _a_pred[4][2] = _b_err[1];
  _a_pred[4][3] = _b_err[2];
  _a_pred[4][4] = _b_err[3];

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

  if((int)_xcp->dim == _n && (int)_q->dim == _nq && 
     (int)_senpar->dim == _nq)
    return;

  //
  // realloc variables for low level _sys->continuous callback
  //

  iv_resize(_x_algebraic,_n);

  _xcp.resize(_n, 0, 0, _nq);
  v_resize(_senpar, _nq);
  v_resize(_q, _nq);

  v_resize(_Fch, _n);
  v_resize(_y, _n);
  v_resize(_y0, _n);
  v_resize(_yprime0, _n);
  v_resize(_yprime1, _n);
  v_resize(_yprime2, _n);

  m_resize(_yprime, _irk_stages, _n);

  m_resize(_Fcxh, _n, _n);

  m_resize(_Sh,_irk_stages*_n,_nq);
  m_resize(_Sh1,_n,_nq);
  m_resize(_Sh2,_n,_nq);
  m_resize(_Sh3,_n,_nq);
  px_resize(_Spivot,_n);

  v_resize(_err,_n);

  px_resize(_ppivot,_n);
  v_resize(_irk_delta,_n);
  v_resize(_irk_res,_n);
  v_resize(_irk_y,_n);
  v_resize(_irk_yprime,_n);

  m_resize(_irk_jac,_n,_n);
  m_resize(_irk_jacf,_n,_n);

  if(_banded_solver && !_sparse_solver)
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
void Omu_IntSDIRK::init(int k,
			const Omu_StateVec &x, const Omu_Vec &u,
			const Omu_DependentVec &Fc, bool sa)
{

  int  i, j;

  if(!_init_method) {

    switch(_method) {

    case 2:

      if (_stepsize <= 0.0)
	m_error(E_INTERN,"Omu_IntSDIRK::init_stage - only fixed"
		  " stepsize allowed for this method!");
      init_method2();

    case 4:
      init_method4();

    default:
      init_method5();
    }
  }

  // check for ODE
  if (Fc.Jdx.is_scalar_constant())
      _ode = true;

  if (_na > 0) {
      iv_resize(_x_algebraic,x->dim);
      iv_zero(_x_algebraic);
      for (i = 0; i < _n; i++) {
	  if (Fc.Jdx.is_zero_column(i))
	      _x_algebraic[i] = 1;	// indicate algebraic state
      }
  }
  
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
  
  if(_sparse_solver)
    _banded_solver = true;
  else {
    if (_jac_sbw >= 0) {
      // take user-defined value
      _ml = _mu = _jac_sbw;
      _banded_solver = _banded;
    }
    else {
      _ml = Fc.Jx.sbw_lower();
      _mu = Fc.Jx.sbw_upper();
      
      _banded_solver = _banded;
      // disable banded solver if full solver appears more efficient
      if (_ml+_ml+_mu >= _n)
	_banded_solver = false;
    }
  }
  
  resize();

  assert(_maxiters > 5);

  _eta = Inf;
  
}

//--------------------------------------------------------------------------

void Omu_IntSDIRK::solve(int kk, double tstart, double tend,
			 Omu_StateVec &xc, Omu_StateVec &dxc, Omu_Vec &q,
			 Omu_DependentVec &Fc)
{

  bool    ok, tend_ok;
  int     i, j, l, nsteps, steps;
  double  err;

  _kk = kk;
  _xc_ptr = &xc;
  _Fc_ptr = &Fc;

  err = 0.0;

  // _stepsize overrides _nsteps
  nsteps = _nsteps;
  if (_stepsize > 0.0)
    nsteps = (int)ceil((tend - tstart) / _stepsize);

  if(_output > 1) {
    std::cout << std::endl;
    std::cout << "Omu_IntSDIRK::solve at stage " << kk;
    std::cout.precision(16);
    std::cout << "  tstart: " << tstart;
    std::cout << "  tend: " << tend << std::endl;
    std::cout.precision(8);
    v_output(xc);
  }

  if(nsteps > 0)
    _h_new = (tend - tstart)/nsteps;
  else
    if(_hinit > 0.0)
      _h_new = min(_hinit,tend-tstart);
    else
      if(_h_new == 0.0)
	_h_new = (tend - tstart)/1.0e3;

  _t = tstart;

  v_zero(_irk_delta);

  v_copy(xc,_y);		// initial states
  v_copy(q,_senpar);            // parameter vector

  // compute initial values for yprime

  if (!_ode) {
      v_rand(_yprime0);
      sv_mlt(1.0e-5,_yprime0,_yprime0);
      init_yprime(kk,tstart,xc,q,_yprime0);
  }
  else
    v_zero(_yprime0);

  for(i = 0; i < _n; i++)
    if(_x_algebraic[i])
      _yprime0[i] =  _y[i];

  set_row(_yprime,1,_yprime0);

  tend_ok = false;
  steps = 0;

  _recalc_jac = true;

  if(_output > 1)
    std::cout << "h: " << _h_new << std::endl;
 
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

	solve_stage(i+1,_y,_yprime0);

	if(_sa)
	  sensitivity_stage(i+1,_y,_yprime0,_senpar);

	for(l = 0; l < _n; l++)
	    if(_x_algebraic[l])
		_yprime0[l] = _y[l];
	
	set_row(_yprime,i,_yprime0);
	
      }

      // now we have to check the error
      if(nsteps <= 0) {

	  err = error_check();

	  _h_new = 0.9*_h*pow(err,-0.25);
    
	  if(err > 1.0) {
	      ok = false;
	      _h = _h_new;
	      tend_ok = false;
	  }
	  else
	      ok = true;
	  
	  if(_output > 1 && !ok) {
	      std::cout << "t: " << _t << "  err: " << err << "  h_new: ";
	      std::cout << _h_new << " accepted: " << ok;
	      std::cout << "  newtonsteps: " << _newtonsteps << std::endl;
	      v_output(_err);
	  }
      }
      else {
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
      }
    }
    while(!ok);

    if(_output > 1) {
      std::cout.precision(16);
      std::cout << "t: " << _t;
      std::cout.precision(8);
      std::cout << "  err: " << err << "  h_new: ";
      std::cout << _h_new << " accepted: " << ok;
      std::cout << "  newtonsteps: " << _newtonsteps << std::endl;
    }

    if(_sa)
	sensitivity_update();

    _t += _h; 

    if(tend_ok)
      break;

  }

  _h_new = _h;

  // return results of integration process
  v_copy(_y, xc);

  if(_output > 1)
    v_output(_y);
    
}

//--------------------------------------------------------------------------

void Omu_IntSDIRK::solve_stage(int stage, VECP y, VECP yprime)
{

    int i, j, result;
    double  del_z, del_z_m, norm_res, theta, tol;
    
    v_copy(y,_irk_y);
    v_copy(yprime,_irk_yprime);

    _eta = Inf;
    del_z_m = 1.0;

    for ( i = 0, tol = 0.0; i < _n; i++ )
	tol += square(y[i]);
    tol = _atol+_rtol*sqrt(tol);

    result = 1;
    _newtonsteps = 0;

    while (result > 0)  {

	// modified Newton's method    
	for(i = 0; i < _modnewtonsteps; i++ ) {
	
	    v_copy(_irk_y,y);
	    for(j = 0; j < _n; j++)
		if(!_x_algebraic[j]) {
		    y[j] += _h*_a[stage-1][stage-1]*_irk_yprime[j];
		    yprime[j] = _irk_yprime[j];
		}
	    
	    // (re)calculate and factorize Jacobian
	    if(_recalc_jac) {
		jac(stage,_t+_c[stage-1]*_h, y, yprime,_senpar,_irk_res);
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
		res(_t+_c[stage-1]*_h, y, yprime,_senpar,_irk_res);	    

	    norm_res = v_norm2(_irk_res);

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
	    
	    del_z = v_norm2(_irk_delta);
	    
	    // convergence measure
	    if ( i == 0 ) {
		theta = 0.0;
		_eta = pow(max(_eta, 1.0e-16), 0.8);
	    } else {
		theta = del_z/del_z_m;
		if ( theta < 1.0 ) 
		    _eta = theta/(1.0-theta);
	    }
	    
	    // check for convergence
	    if ( ( _eta*del_z <= _kappa*tol ) || 
		 ( norm_res < 1.0e-6*_kappa*tol ) ) {

		v_copy(_irk_y,y);
		for(j = 0; j < _n; j++)
		  if(!_x_algebraic[j]) {
		    y[j] += _h*_a[stage-1][stage-1]*_irk_yprime[j];
		    yprime[j] = _irk_yprime[j];
		  }

		result = 0;
		break;
	    }
	    
	    // check for divergence
	    if ( ( i >= 1 ) && ( ( theta >= 1.0 ) || 
		 ( pow(theta, _modnewtonsteps-1-i)/(1.0-theta)*del_z > 
		   _kappa*tol ) ) ) {
	      
	      if( (_nsteps > 0) || (_stepsize > 0.0) || (stage > 1)) {
		_recalc_jac = true;  // fixed stepsize
		if(_newtonsteps > _maxiters - 1)
		  m_error(E_CONV, "Omu_IntSDIRK::solve_stage");
	      }
	      else {
		if(_h < MACHEPS*_h_new)
		  m_error(E_CONV, "Omu_IntSDIRK::solve_stage");

		_h = 0.5*_h;
	      }
	    }

	    del_z_m = del_z;

	}

	_recalc_jac = true;
	_newtonsteps += i;	

    };

    _newtonsteps += i;

    return;

}

//--------------------------------------------------------------------------

void Omu_IntSDIRK::sensitivity_stage(int stage, const VECP y, 
				     const VECP yprime, const VECP senpar)
{

    int i, j, k;

    Omu_StateVec &xc = *_xc_ptr;
    Omu_DependentVec &Fc = *_Fc_ptr;

    // form vectors of independent variables

    for (i = 0; i < _n; i++) {
	xc[i] = y[i];
	_xcp[i] = yprime[i];
    }
    
    for (i = 0; i < _nq; i++) {
	_q[i] = senpar[i];
    }
    
    // evaluate residuals and Jacobians
    bool is_required_J_bak = Fc.is_required_J();
    Fc.set_required_J(true);
    
    residual(_kk, _t+_c[stage-1]*_h, xc, _xcp, _q, Fc);
    
    Fc.set_required_J(is_required_J_bak);

     if(_ode) {
	sm_mlt(_h*_a[stage-1][stage-1],Fc.Jx,_irk_jac);
	for(i = 0; i < _n; i++)	  
	    _irk_jac[i][i] += Fc.Jdx[i][i];
    }
    else {
	
	v_set(_Fch,_h*_a[stage-1][stage-1]);
	
	if(_na > 0)    
	    for(i = 0; i < _n; i++)
		if(_x_algebraic[i]) 
		    _Fch[i] = 1.0;
	
	for(i = 0; i < _n; i++)
	    for(j = 0; j < _n; j++)
		_irk_jac[i][j] = Fc.Jdx[i][j] + _Fch[j]*Fc.Jx[i][j];
    }

    m_copy(Fc.Jq,_Sh1);
    m_move(Fc.Jx,0,0,_n,_n,_Sh1,0,_nd);

    for(i = 0; i < stage-1; i++) {
	
	m_move(_Sh,_n*i,0,_n,_nq,_Sh2,0,0);
	
	if(_na > 0) {
	    v_set(_Fch,_h*_a[stage-1][i]);
	    for(j = 0; j < _n; j++)
		if(_x_algebraic[j]) 
		    _Fch[j] = 1.0;
	    for(j = 0; j < _n; j++)
		for(k = 0; k < _n; k++)
		    _irk_jacf[j][k] = _Fch[k]*Fc.Jx[j][k];	  
	}
	else
	    sm_mlt(_h*_a[stage-1][i],Fc.Jx,_irk_jacf);
	
	m_mlt(_irk_jacf,_Sh2,_Sh3);
	m_add(_Sh1,_Sh3,_Sh1);
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
    
    for(i = 0; i < _nq; i++) {
	get_col(_Sh1,i,_irk_res);
	sv_mlt(-1.0,_irk_res,_irk_res);
	
	if(_banded_solver) {
	    if(_sparse_solver)
		spLUsolve(_irk_jac_sp, _ppivot, _irk_res, _irk_delta);
	    else
		bdLUsolve(_irk_jac_bd, _ppivot, _irk_res, _irk_delta);
	}
	else
	    LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);
	
	set_col(_Sh1,i,_irk_delta);
    }
    
    m_move(_Sh1,0,0,_n,_nq,_Sh,(stage-1)*_n,0);
    
}

//--------------------------------------------------------------------------

void Omu_IntSDIRK::solve_final(VECP y, VECP yprime)
{

    int i, j, result;
    double  del_z, del_z_m, norm_res, theta, tol;
    
    v_copy(y,_irk_y);
    v_copy(yprime,_irk_yprime);

    _eta = Inf;
    del_z_m = 1.0;

    for ( i = 0, tol = 0.0; i < _n; i++ )
	tol += square(y[i]);
    tol = _atol+_rtol*sqrt(tol);

    result = 1;
    _newtonsteps = 0;

    while (result > 0)  {

	// modified Newton's method    
	for(i = 0; i < _modnewtonsteps; i++ ) {
	
	    v_copy(_irk_y,y);
	    v_copy(_irk_yprime,yprime);
	    
	    // (re)calculate and factorize Jacobian
	    if(_recalc_jac) {
		jac(0,_t+_h, y, yprime,_senpar,_irk_res);
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
		res(_t+_h, y, yprime,_senpar,_irk_res);	    

	    norm_res = v_norm2(_irk_res);

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
	    
	    del_z = v_norm2(_irk_delta);
	    
	    // convergence measure
	    if ( i == 0 ) {
		theta = 0.0;
		_eta = pow(max(_eta, 1.0e-16), 0.8);
	    } else {
		theta = del_z/del_z_m;
		if ( theta < 1.0 ) 
		    _eta = theta/(1.0-theta);
	    }
	    
	    // check for convergence
	    if ( ( _eta*del_z <= _kappa*tol ) || 
		 ( norm_res < 1.0e-6*_kappa*tol ) ) {

		v_copy(_irk_y,y);
		v_copy(_irk_yprime,yprime);

		result = 0;
		break;

	    }
	    
	    // check for divergence
	    if ( ( i >= 1 ) && ( ( theta >= 1.0 ) || 
		 ( pow(theta, _modnewtonsteps-1-i)/(1.0-theta)*del_z > 
		   _kappa*tol ) ) ) {

	      _recalc_jac = true;  // don't modify step size

	    }

	    del_z_m = del_z;

	}

	_recalc_jac = true;	
	_newtonsteps += i;

	if(_newtonsteps > _maxiters - 1)
	  m_error(E_CONV, "Omu_IntSDIRK::solve_final");

    };

    _newtonsteps += i;


    return;

}

//--------------------------------------------------------------------------

void Omu_IntSDIRK::sensitivity_update()
{

    int     i, j, k;
    double  val;
    Omu_StateVec &xc = *_xc_ptr;

    m_resize(_Sh1,_n,_nq);    
    m_resize(_Sh2,_n,_nq);    
    m_resize(_Sh3,_n,_n);    

    m_zero(_Sh1);
    m_zero(_Sh2);

    if(_na == 0) {

	for(i = 0; i < _irk_stages; i++) {
	    m_move(_Sh,i*_n,0,_n,_nq,_Sh2,0,0);
	    ms_mltadd(_Sh1,_Sh2,_h*_b[i],_Sh1);
	}

	m_move(_Sh1,0,_nd,_n,_n,_Sh3,0,0);
	m_mlt(_Sh3,xc.Sq,_Sh2);
	m_zero(_Sh3);
	m_move(_Sh3,0,0,_n,_n,_Sh1,0,_nd);
	m_add(xc.Sq,_Sh2,xc.Sq);
	m_add(xc.Sq,_Sh1,xc.Sq);

    }
    else {

	for(i = 0; i < _irk_stages; i++) {
	    val = _h*_b[i];
	    for (j = 0; j < _n; j++)
		if(!_x_algebraic[j])
		    for (k = 0; k < _nq; k++)
			_Sh1[j][k] += val*_Sh[i*_n+j][k];
	}
	
	for (j = 0; j < _n; j++)
	    if(!_x_algebraic[j])
		for (k = 0; k < _nq; k++)
		    _Sh1[j][k] += _Sh[(_irk_stages-1)*_n+j][k];

	m_move(_Sh1,0,_nd,_n,_n,_Sh3,0,0);
	m_mlt(_Sh3,xc.Sq,_Sh2);
	m_zero(_Sh3);
	m_move(_Sh3,0,0,_n,_n,_Sh1,0,_nd);

	for(i = 0; i < _n; i++) 
	    if(!_x_algebraic[i])
		for(j = 0; j < _nq; j++)
		    xc.Sq[i][j] += _Sh1[i][j] + _Sh2[i][j];
	    else
		for(j = 0; j < _nq; j++)
		    xc.Sq[i][j] = _Sh1[i][j] + _Sh2[i][j];

    }

    _sen_evals++;
    
}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::res(double t, const VECP  y, const VECP yprime, 
		       const VECP par, VECP res)
{

  int i;

  Omu_StateVec &xc = *_xc_ptr;
  Omu_DependentVec &Fc = *_Fc_ptr;

  // form vectors of independent variables

  for (i = 0; i < _n; i++) {
    xc[i] = y[i];
    _xcp[i] = yprime[i];
  }

  for (i = 0; i < _nq; i++) {
    _q[i] = par[i];
  }

  // evaluate residual
  residual(_kk, t, xc, _xcp, _q, Fc);

  v_copy(Fc,res);

  _res_evals++;

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::jac(int stage, double t, const VECP y, const VECP yprime,
		       const VECP par, VECP res)
{

  int i, j;

  Omu_StateVec &xc = *_xc_ptr;
  Omu_DependentVec &Fc = *_Fc_ptr;

  // form vectors of independent variables

  for (i = 0; i < _n; i++) {
    xc[i] = y[i];
    _xcp[i] = yprime[i];
  }

  for (i = 0; i < _nq; i++) {
    _q[i] = par[i];
  }

  m_zero(_irk_jac);

  // evaluate residuals and Jacobians
  bool is_required_J_bak = Fc.is_required_J();
  Fc.set_required_J(true);

  residual(_kk, t, xc, _xcp, _q, Fc);

  Fc.set_required_J(is_required_J_bak);

  v_copy(Fc, res);

  if(stage > 0) {

      if(_ode) {
	  sm_mlt(_h*_a[stage-1][stage-1],Fc.Jx,_irk_jac);
	  for(i = 0; i < _n; i++)	  
	      _irk_jac[i][i] += Fc.Jdx[i][i];
      }
      else {

	  v_set(_Fch,_h*_a[stage-1][stage-1]);
      
	  if(_na > 0)    
	      for(i = 0; i < _n; i++)
		  if(_x_algebraic[i]) 
		      _Fch[i] = 1.0;
	  
	  for(i = 0; i < _n; i++)
	      for(j = 0; j < _n; j++)
		  _irk_jac[i][j] = Fc.Jdx[i][j] + _Fch[j]*Fc.Jx[i][j];
      }
  }
  else {
    
      // yprime and algebraic states
      if(stage == 0) {
	  
	  v_zero(_Fch);
	  
	  if(_na > 0)    
	      for(i = 0; i < _n; i++)
		  if(_x_algebraic[i]) 
		      _Fch[i] = 1.0;
	  
	  for(i = 0; i < _n; i++)
	      for(j = 0; j < _n; j++)
		  _irk_jac[i][j] = Fc.Jdx[i][j] + _Fch[j]*Fc.Jx[i][j];
      }

      // only algebraic states
      if(stage == -1) {
	  if(_na > 0)    
	      for(i = 0; i < _n; i++)
		  if(_x_algebraic[i])
		      _Fch[i] = 1.0;
		  else
		      _irk_jac[i][i] += 1.0e-8;

	  for(i = 0; i < _n; i++)
	      for(j = 0; j < _n; j++)
		  _irk_jac[i][j] = _Fch[j]*Fc.Jx[i][j];
      }

      // only yprime
      if(stage == -2) {
	  for(i = 0; i < _n; i++)
	    if(!_x_algebraic[i]) 
	      for(j = 0; j < _n; j++)
		_irk_jac[i][j] = Fc.Jdx[i][j];
	    else
	      _irk_jac[i][i] = 1.0;

      }

    // error evaluation
    if(stage == -3) {

	v_set(_Fch,_h*_gamma0);

	if(_na > 0)    
	    for(i = 0; i < _n; i++)
		if(_x_algebraic[i]) {
		    _Fch[i] = 0.0;
		    _irk_jac[i][i] = 1.0-8;
		}
      
	for(i = 0; i < _n; i++)
	    if(!_x_algebraic[i]) 
		for(j = 0; j < _n; j++)
		    _irk_jac[i][j] = Fc.Jdx[i][j]+_Fch[j]*Fc.Jx[i][j];

    }
  }
    
  _jac_evals++;

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::init_yprime(int k, double t, const Omu_StateVec &y, 
			       const Omu_Vec &q, VECP yprime)
{

  int i, j, result;
  double  del_z, del_z_m, norm_res, theta, tol;

  _eta = Inf;
  del_z_m = 1.0;

  for ( i = 0, tol = 0.0; i < _n; i++ )
    tol += square(y[i]);
  tol = _atol+_rtol*sqrt(tol);

  result = 1;
  _newtonsteps = 0;

  _recalc_jac = true;

  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      _irk_yprime[i] = yprime[i];
    else
      _irk_yprime[i] = 0.0;

  if(_output == -1) {
    std::cout << "initial values for yprime at time: " << _t << std::endl; 
    std::cout << "newton method" << std::endl;
  };

  while (result > 0)  {

    // modified Newton's method    
    for(i = 0; i < _modnewtonsteps; i++ ) {

      v_copy(_irk_yprime,yprime);
	
      // (re)calculate and factorize Jacobian
      if(_recalc_jac) {
	jac(-2, t, y, yprime, q,_irk_res);
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
	res(t, y, yprime, q, _irk_res);	    

      norm_res = 0.0;
      for(j = 0; j < _n; j++)
	if(!_x_algebraic[j])
	  norm_res += _irk_res[j]*_irk_res[j];
      norm_res = sqrt(norm_res);

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
      
      del_z = 0.0;
      for(j = 0; j < _n; j++)
	if(!_x_algebraic[j])
	  del_z += _irk_res[j]*_irk_res[j];
      del_z = sqrt(norm_res);

      // convergence measure
      if ( i == 0 ) {
	theta = 0.0;
	_eta = pow(max(_eta, 1.0e-16), 0.8);
      } else {
	theta = del_z/del_z_m;
	if ( theta < 1.0 ) 
	  _eta = theta/(1.0-theta);
      }
      
      // check for convergence
      if ( ( _eta*del_z <= _kappa*tol ) || 
	   ( norm_res < 1.0e-6*_kappa*tol ) ) {
	
	v_copy(_irk_yprime,yprime);
	
	result = 0;
	break;
	
      }
      
      // check for divergence
      if ( ( i >= 1 ) && ( ( theta >= 1.0 ) || 
	   ( pow(theta, _modnewtonsteps-1-i)/(1.0-theta)*del_z > 
			     _kappa*tol ) ) ) {
	
	_recalc_jac = true;
	
      }
      
      del_z_m = del_z;
      
    }
    
    _recalc_jac = true;	
    _newtonsteps += i;

    if(_newtonsteps > _maxiters - 1)
      m_error(E_INTERN, "Omu_IntSDIRK::init_yprime");

  };
  
  _newtonsteps += i;

  return;
   
}
      
//--------------------------------------------------------------------------
double Omu_IntSDIRK::error_check()
{

    int     i, j, l;
    double  err;

    assert(_gamma0 > 0);

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
    
    // algebraic variables are not considered in the error evaluation
	
    for(i = 0; i < _n; i++)
	if(_x_algebraic[i])
	    _yprime1[i] = 0.0;
	
    v_sub(_err,_y,_irk_y);

    // factorized Jacobian for error control is availible
	
    if(_banded_solver) {
      if(_sparse_solver)
	spLUsolve(_irk_jac_sp, _ppivot, _irk_y, _err);
      else
	bdLUsolve(_irk_jac_bd, _ppivot, _irk_y, _err);
    }
    else
      LUsolve(_irk_jac, _ppivot, _irk_y, _err);

    err = 0.0;
    
    /*
    for(i = 0; i < _n; i++)
	if(!_x_algebraic[i]) {
	    _err[i] = _err[i]/(_atol+max(fabs(_y0[i]),fabs(_y[i]))*_rtol);
	    err = max(err,fabs(_err[i]));
	}
	else
	    _err[i] = 0.0;
    
    err = max(v_norm2(_err)/sqrt(double(_n-_na)),_atol);
    */

    err = v_norm2(_err)/(_atol+_rtol*v_norm2(_y));
    err = max(err,_atol);

    return err;

}     

    
    



