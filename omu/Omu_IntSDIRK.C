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
#include <fstream.h>

#include <adutils.h>

#include <Hqp.h>
#include <If_Int.h>
#include <If_Bool.h>
#include <If_Float.h>
#include <If_Class.h>
#include <Meschach.h>

#include "Omu_IntSDIRK.h"

extern "C" {
#include <matrix2.h>
}

IF_CLASS_DEFINE("SDIRK", Omu_IntSDIRK, Omu_Integrator);

//--------------------------------------------------------------------------

Omu_IntSDIRK::Omu_IntSDIRK()
{

  _recalc_jac = true;
  _stiffly_accurate = false;
  _sens_at_once = false;

  _output = 0;
  _n_splitt_tape_eval = 30;
  _tag = 10;

  _x_algebraic = iv_get(1);

  _ppivot= px_get(1);
  _Spivot= px_get(1);

  _cx = v_get(1);
  _cu = v_get(1);
  _cxp = v_get(1);
  _cF = v_get(1);
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

  _cFx  = m_get(1,1);
  _cFxh = m_get(1,1);
  _cFu  = m_get(1,1);
  _cFxp = m_get(1,1);
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

  _nz = NULL;

  _irk_jac = m_get(1,1);
  _irk_jacf = m_get(1,1);

  _na = 0;
  _nod = 0;

  _ml = 0;
  _mu = 0;

  _nsteps = 0;
  _modnewtonsteps = 5;
  _newtonsteps = 0;
  _maxiters = 20;

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

  _gamma0 = 1.0;

  _ifList.append(new If_Bool("prg_int_sens", &_sens_at_once));
  _ifList.append(new If_Int("prg_int_out", &_output));
  _ifList.append(new If_Int("prg_int_nsteps", &_nsteps));
  _ifList.append(new If_Int("prg_int_nsplitt", &_n_splitt_tape_eval));
  _ifList.append(new If_Int("prg_int_modnewtonsteps", &_modnewtonsteps));
  _ifList.append(new If_Int("prg_int_maxiters", &_maxiters));
  _ifList.append(new If_Float("prg_int_hinit", &_hinit));

  init_method();

}

//--------------------------------------------------------------------------

Omu_IntSDIRK::~Omu_IntSDIRK()
{

  int i;

  IV_FREE(_x_algebraic);

  PX_FREE(_ppivot);
  PX_FREE(_Spivot);

  V_FREE(_cx);
  V_FREE(_cu);
  V_FREE(_cxp);
  V_FREE(_cF);
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

  V_FREE(_b);
  V_FREE(_b_err);
  V_FREE(_c);
  M_FREE(_a);
  M_FREE(_a_pred);

  M_FREE(_yprime);
  M_FREE(_cFx);
  M_FREE(_cFxh);
  M_FREE(_cFu);
  M_FREE(_cFxp);

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

  if(_nz != NULL) {
    for(i = 0; i < _nd+_n+_m; i++)
      delete[] _nz[i];
    delete[] _nz;
  }

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::init_method()
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

  if(_c[_c->dim-1] == 1.0) {
    _stiffly_accurate = true;
    for(i = 0; i < _irk_stages; i++)
      if(_b[i] != _a[_irk_stages-1][i]) {
	_stiffly_accurate = false;
	break;
      }
  }

  init_yprime_pred(_a_pred);

  _gamma0 = 1.0/4.000271085;

}

//--------------------------------------------------------------------------

void Omu_IntSDIRK::init_yprime_pred(MATP pred)
{

  // simple linear interpolation

  m_resize(pred,_irk_stages,_irk_stages+1);

  pred[0][0] = 2.0;
  pred[0][2] = -1.0;
  pred[1][0] = -2.0;
  pred[1][1] = 3.0;
  pred[2][1] = 0.5;
  pred[2][2] = 0.5;
  pred[3][3] = 1.0;
  pred[4][2] = -1.0;
  pred[4][4] = 2.0;

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::realloc()
{

  int i;

  if (_cx->dim == _nd + _n && _cu->dim == _nu)
    return;

  //
  // realloc variables for low level _sys->continuous callback
  //

  iv_resize(_x_algebraic,_n);

  v_resize(_cx, _nd + _n);
  v_resize(_cu, _nu);
  v_resize(_cxp, _nd + _n);
  v_resize(_cF, _nd + _n);
  v_resize(_cFh, _n);

  v_resize(_par, _nd+_nu);
  v_resize(_y, _n);
  v_resize(_y0, _n);
  v_resize(_yprime0, _n);
  v_resize(_yprime1, _n);
  v_resize(_yprime2, _n);

  m_resize(_yprime, _irk_stages, _n);

  m_resize(_cFx, _nd + _n, _nd + _n);
  m_resize(_cFxh, _n, _n);
  m_resize(_cFu, _nd + _n, _nu);
  m_resize(_cFxp, _nd + _n, _nd + _n);

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

  if(_nz != NULL) {
    for(i = 0; i < _n; i++)
      delete[] _nz[i];
    delete[] _nz;
  }
  
  // set up sparsity array
  _nz = new short*[_n];
  for(i = 0; i < _n; i++)
    _nz[i] = new short[_n+_nd+_nu];

}


//--------------------------------------------------------------------------
void Omu_IntSDIRK::init_stage(int k,
			    const Omu_States &x, const Omu_Vector &u,
			    bool sa)
{

  int  i, j;

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

  // check for sdirk method
  for(i = 0; i < _irk_stages-1; i++)
    for(j = i+1; j < _irk_stages; j++)
      if(_a[i][j] != 0.0)
	error(E_INTERN,"Omu_IntSDIRK::init_stage - 
                        sdirk method not correctly defined!");

  for(i = 1; i < _irk_stages; i++)
      if(_a[i][i] != _a[0][0])
	error(E_INTERN,"Omu_IntSDIRK::init_stage - 
                        sdirk method not correctly defined!");

  realloc();

}

//--------------------------------------------------------------------------

void Omu_IntSDIRK::solve(int kk, Real tstart, Real tend,
			 const Omu_States &x, const Omu_Vector &u,
			 Omu_Program *sys, VECP xt,
			 MATP Sx, MATP Su)
{

  bool    ok, tend_ok;
  int     i, j, l, n, nv, m, steps;
  double  err;

  if(_output > 1) {
    cout << endl;
    cout << "Omu_IntSDIRK::solve at stage " << kk;
    cout.precision(16);
    cout << "  tstart: " << tstart;
    cout << "  tend: " << tend << endl;
    cout.precision(8);
    v_output(xt);
  }

  _sys = sys;	// propagate to res() and jac()

  m_zero(_cFx);
  m_zero(_cFu);
  m_zero(_cFxp);
  m_zero(_Sxx0);

  _nod = 0;

  if(_sa) {
    m_move(Sx,_nd,0,_n-_nod,_nd,_Sxd0,0,0);
    m_move(Sx,_nd,_nd,_n,_n-_nod,_Sxx0,0,0);
    m_move(Su,_nd,0,_n,_nu,_Sxu0,0,0);
  }

  // compute initial values for yprime
  v_zero(_yprime0);
  for(i = 0; i < _n; i++)
    _yprime0[i] = 1.0e-5;
  init_yprime(kk,tstart,xt,u,_yprime0);

  if(_nsteps > 0)
    _h_new = (tend - tstart)/_nsteps;
  else
    if(_hinit > 0.0)
      _h_new = min(_hinit,tend-tstart);
    else
      _h_new = (tend - tstart)/1.0e2;
  //      error(E_INTERN,"Omu_IntSDIRK::solve no step size provided");

  v_zero(_y);

  for (i = 0; i < _nd; i++)
    _par[i] = xt[i];

  for (i = 0; i < _n; i++)
    _y[i] = xt[_nd + i];		// initial states

  for (i = 0; i < _nu; i++)
    _par[_nd + i] = u[i];

  v_zero(_irk_delta);

  _h = _h_new;
  _t = tstart;

  for(i = 0; i < _n; i++)
    if(_x_algebraic[i])
      _yprime0[i] =  _y[i];

  set_row(_yprime,1,_yprime0);

  tend_ok = false;
  steps = 0;

  while(_t < tend) {    

    steps++;

    if(_t + _h >= tend) {
      _h = tend - _t;
      tend_ok = true;
    }

    if(_nsteps > 0 && steps == _nsteps)
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
	
	for(l = 0; l < _n; l++)
	  if(_x_algebraic[l])
	    _yprime0[l] = _y[l];
	
	set_row(_yprime,i,_yprime0);
	
      }

      // now we have to check the error
      
      if(_nsteps <= 0) {
	if(_b_err->dim == _irk_stages) {
	  v_copy(_y0,_err);
	  for(j = 0; j < _irk_stages; j++)
	    for(l = 0; l < _n; l++)
	      if(!_x_algebraic[l])
		_err[l] +=  _h*_b_err[j]*_yprime[j][l];
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

      if(_nsteps <= 0) {

	// algebraic variables are not considered in the error evaluation

	for(i = 0; i < _n; i++)
	  if(_x_algebraic[i])
	    _yprime1[i] = 0.0;

	v_sub(_err,_y,_irk_y);
	jac(-3,_t,_y0,_yprime1,_par,_irk_res);
	LUfactor(_irk_jac, _ppivot);
	LUsolve(_irk_jac, _ppivot, _irk_y, _err);

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

    if(_output > 2) {
      cout << "y_hat(t+h) - y(t+h):" << endl;
      v_output(_irk_y);
      v_output(_irk_res);
    }

    if(_sa)
      sensitivity();

    _t += _h; 
    _h = _h_new;

    if(tend_ok)
      break;

  }

  // return results of integration process
  for(i = 0; i < _n; i++)
    xt[_nd+i] = _y[i];

  if(_sa) {
    m_move(_Sxd0,0,0,_n,_nd,Sx,_nd,0);
    m_move(_Sxx0,0,0,_n,_n-_nod,Sx,_nd,_nd);
    m_move(_Sxu0,0,0,_n,_nu,Su,_nd,0);
  }

}

//--------------------------------------------------------------------------

void Omu_IntSDIRK::solve_stage(int stage, VECP y, VECP yprime)
{

  bool    ok;
  int     i, j, l, n, m, an, ap;
  double  alpha, cur_err, old_err;

  v_ones(_irk_res);
  v_copy(y,_irk_y);
  v_copy(yprime,_irk_yprime);

  cur_err = 1.0e10;
  ok = false;
  _recalc_jac = false;
    
  if(_output > 2)
    cout << "   stage: " << stage;
  
  // modified Newton's method    
  for(i = 0; i < _maxiters; i++ ) {
    
    if(stage < 2 && (i % _modnewtonsteps) == 0)
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

	if(cur_err > old_err) {

	  if(_output > 2) {
	    cout << "alpha: " << alpha << "   cur_err: " << cur_err;
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
      cout << "   iter: " << i << "  res: " << cur_err;
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
      error(E_INTERN, 
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

      LUfactor(_irk_jac, _ppivot);
      _recalc_jac = false;
    }
    else
      res(_t+_c[stage-1]*_h, y, yprime,_par,_irk_res);
    
    sv_mlt(-1.0,_irk_res,_irk_res);
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

void Omu_IntSDIRK::solve_final(VECP y, VECP yprime)
{

  bool    ok;
  int     i, j, l, n, m, an, ap;
  double  cur_err, old_err;

  v_ones(_irk_res);
  v_copy(y,_irk_y);
  v_copy(yprime,_irk_yprime);

  cur_err = 1.0e10;
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
      error(E_INTERN, 
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

      LUfactor(_irk_jac, _ppivot);
      _recalc_jac = false;

    }
    else
      res(_t+_h, y, yprime,_par,_irk_res);
    
    sv_mlt(-1.0,_irk_res,_irk_res);
    LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);

  }

  _newtonsteps += i;

}


//--------------------------------------------------------------------------
// differentiate the whole integration scheme using adol-c

void Omu_IntSDIRK::sensitivity()
{

  int i, j, l, m;

  double  aindep[_irk_stages*_n+_n+_nd+_nu];
  double  f[_irk_stages*_n];

  adoublev au(_nu);
  adoublev ay(_nd+_n);
  adoublev ay0(_nd+_n);
  adoublev ayprime(_nd+_n);
  adoublev aF(_nd+_n);
  adoublev irk_ay(_irk_stages*(_nd+_n));
  adoublev irk_ayprime(_irk_stages*(_nd+_n));
  adoublev irk_aF(_irk_stages*(_nd+_n));

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

	LUfactor(_irk_jac, _ppivot);
	
	for(j = 0; j < _nd+_n+_nu; j++) {
	  get_col(_Sxd,j,_irk_res);
	  sv_mlt(-1.0,_irk_res,_irk_res);
	  LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);
	  set_col(_Sxd,j,_irk_delta);
	}
	
	for(j = 0; j < i; j++) {
	  
	  m_move(_Z,i*_n,(j+1)*_n+_nd+_nu,_n,_n,_Sxx,0,0);
	  m_move(_Sh,j*_n,0,_n,_nd+_n+_nu,_SF,0,0);
	  
	  for(l = 0; l < _n; l++) {
	    get_col(_Sxx,l,_irk_res);
	    sv_mlt(-1.0,_irk_res,_irk_res);
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

	LUfactor(_irk_jac, _ppivot);

	for(j = 0; j < _nd+_n+_nu; j++) {
	  get_col(_Sxd,j,_irk_res);
	  sv_mlt(-1.0,_irk_res,_irk_res);
	  LUsolve(_irk_jac, _ppivot, _irk_res, _irk_delta);
	  set_col(_Sxd,j,_irk_delta);
	}
	
	for(j = 0; j < i; j++) {
	  
	  m_move(_Z,0,(j+1)*_n+_nd+_nu,_n,_n,_Sxx,0,0);
	  m_move(_Sh,j*_n,0,_n,_nd+_n+_nu,_SF,0,0);
	  
	  for(l = 0; l < _n; l++) {
	    get_col(_Sxx,l,_irk_res);
	    sv_mlt(-1.0,_irk_res,_irk_res);
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
  
}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::res(double t, const VECP  y, const VECP yprime, 
		       const VECP par, VECP res)
{

  int i, j;

  // form vectors of independent variables

  for (i = 0; i < _nd; i++) {
    _cx[i] = par[i];
    _cxp[i] = 0.0;
  }
  for (i = 0; i < _n; i++) {
    _cx[_nd + i] = y[i];
    _cxp[_nd + i] = yprime[i];
  }
  for (i = 0; i < _nu; i++) {
    _cu[i] = par[_nd + i];
  }

  // evaluate residual

  _sys->continuous(_kk, t, _cx, _cu, _cxp, _cF);

  for (i = 0; i < _n; i++)
    res[i] = _cF[_nd+i];

  _res_evals++;

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::jac(int stage, double t, const VECP y, const VECP yprime,
		       const VECP par, VECP res)
{

  int    i, j, n;
  double row_sum;

  // form vectors of independent variables

  for (i = 0; i < _nd; i++) {
    _cx[i] = par[i];
    _cxp[i] = 0.0;
  }
  for (i = 0; i < _n; i++) {
    _cx[_nd + i] = y[i];
    _cxp[_nd + i] = yprime[i];
  }
  for (i = 0; i < _nu; i++) {
    _cu[i] = par[_nd + i];
  }

  // evaluate residuals and Jacobians

  _sys->continuous(_kk, t, _cx, _cu, _cxp, _cF, _cFx, _cFu, _cFxp);

  for (i = 0; i < _n; i++)
    res[i] = _cF[_nd+i];

  if(stage > 0) {

    v_set(_cFh,_h*_a[stage-1][stage-1]);
      
    if(_na > 0)    
      for(i = 0; i < _n; i++)
	if(_x_algebraic[i]) 
	  _cFh[i] = 1.0;
    
    for(i = 0; i < _n; i++)
      for(j = 0; j < _n; j++)
	_irk_jac[i][j] = _cFxp[_nd+i][_nd+j] + _cFh[j]*_cFx[_nd+i][_nd+j];
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
	  _irk_jac[i][j] = _cFxp[_nd+i][_nd+j] + _cFh[j]*_cFx[_nd+i][_nd+j];
    }

    // only algebraic states
    if(stage == -1) {
      if(_na > 0)    
	for(i = 0; i < _n; i++)
	  if(_x_algebraic[i]) 
	    _cFh[i] = 1.0;
      for(i = 0; i < _n; i++)
	for(j = 0; j < _n; j++)
	  _irk_jac[i][j] = _cFh[j]*_cFx[_nd+i][_nd+j];
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
	  _irk_jac[i][j] = _cFxp[_nd+i][_nd+j];
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
	if(_x_algebraic[i]) 
	  for(j = 0; j < _n; j++)
	    _irk_jac[i][j] = _cFxp[_nd+i][_nd+j] + _cFh[j]*_cFx[_nd+i][_nd+j];

      if(_na > 0)    
	for(i = 0; i < _n; i++)
	  if(_x_algebraic[i]) 
	    _irk_jac[i][i] = 1.0-8;
    }
  }
    
  _res_evals++;
  _jac_evals++;

}

//--------------------------------------------------------------------------
void Omu_IntSDIRK::init_yprime(int k, double t, const VECP y, const VECP u, 
			     VECP yprime)
{

  bool    ok;
  int     i, j, l, m, n;
  double  cur_err, old_err;
  
  m_resize(_Sh,_n,_n);

  for(i = 0; i < _nd; i++)
    _cxp[i] = 0.0;
  for(i = 0; i < _n; i++)
    if(!_x_algebraic[i])
      _cxp[_nd+i] = yprime[i];
    else
      _cxp[_nd+i] = 0.0;

  if(_output > 2) {
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
	  ok = ok  && (fabs(_cF[_nd+j])<0.1*_atol);
	  cur_err = max(cur_err,fabs(_cF[_nd+j]));
	  cur_err = max(cur_err,fabs(_cF[_nd+j]));
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

    if(_output > 2)
      cout << "   iter: " << i << "  res: " << cur_err << endl;

    if(ok) 
      break;
    else if (i == _maxiters-1)
      error(E_INTERN, 
	    "Omu_IntSDIRK::init_yprime Newton method failed to converge");
      

    // (re)calculate and factorize Jacobian
    if(_recalc_jac) {
      _sys->continuous(k, t, y, u, _cxp, _cF, _cFx, _cFu, _cFxp);
      if(i == 0) {
	old_err = 0.0;
	// check for convergence
	for(ok = true, j = 0 ; j < _n; j++)
	  if(!_x_algebraic[j]) {
	    ok = ok  && (fabs(_cF[_nd+j])<0.1*_atol);
	    old_err = max(old_err,fabs(_cF[_nd+j]));
	  }
	old_err++;
      }

      m_move(_cFxp,_nd,_nd,_n,_n,_Sh,0,0);
      for(j = 0; j < _n; j++)
	if(_x_algebraic[j]) 
	  _Sh[j][j] = 1.0;
      LUfactor(_Sh, _ppivot);
    }
    else
      _sys->continuous(k, t, y, u, _cxp, _cF);

    v_move(_cF,_nd,_n,_irk_res,0);
    sv_mlt(-1.0,_irk_res,_irk_res);
    LUsolve(_Sh, _ppivot, _irk_res, _irk_delta);
  }

  for(j = 0; j < _n; j++)
    yprime[j] = _cxp[_nd+j];

}












