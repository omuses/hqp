/*
 * Omu_IntIMP.C --
 *   -- class implementation
 *
 * E. Arnold   1999-04-12
 *             2000-05-30 step size control
 *             2001-08-16 prg_int_nsteps --> _stepsize
 *
 */

/*
    Copyright (C) 1999--2001  Eckhard Arnold

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

#include <adutils.h>

#include <If_Int.h>
#include <If_Real.h>
#include <If_Class.h>

#include "Omu_IntIMP.h"

extern "C" {
#include <matrix2.h>
}

IF_CLASS_DEFINE("IMP", Omu_IntIMP, Omu_Integrator);

//--------------------------------------------------------------------------
Omu_IntIMP::Omu_IntIMP()
{

  _y     = v_get(1);
  _u     = v_get(1);
  _k1    = v_get(1);
  _y0    = v_get(1);
  _res   = v_get(1);
  _fh    = v_get(1);
  _yjac  = v_get(1);
  _yjacp = v_get(1);
  _yy    = m_get(1, 1);
  _yyn   = m_get(1, 1);
  _ppivot= px_get(1);
  _y1    = v_get(1);
  _y2    = v_get(1);

  _modnewtonsteps = 5;
  _maxiters = 20;
  _hinit = 0.0;
  _dt = 0.0;
  _IMP = 1;

  _ifList.append(new If_Int("prg_int_modnewtonsteps", &_modnewtonsteps));
  _ifList.append(new If_Int("prg_int_maxiters", &_maxiters));
  _ifList.append(new If_Real("prg_int_hinit", &_hinit));

  _res_evals = 0;
  _jac_evals = 0;

}

//--------------------------------------------------------------------------
Omu_IntIMP::~Omu_IntIMP()
{

  V_FREE(_y);
  V_FREE(_u);
  V_FREE(_k1);
  V_FREE(_y0);
  V_FREE(_res);
  V_FREE(_fh);
  V_FREE(_yjac);
  V_FREE(_yjacp);
  M_FREE(_yy);
  M_FREE(_yyn);
  PX_FREE(_ppivot);
  V_FREE(_y1);
  V_FREE(_y2);

}

//--------------------------------------------------------------------------
void Omu_IntIMP::realloc()
{

  v_resize(_u, _npar);
  v_resize(_y, _n);
  v_resize(_y0, _n);
  v_resize(_k1, _n);
  v_resize(_res, _n);
  v_resize(_fh, _n);

  v_resize(_yjac, _n*(1+_n+_npar));
  v_resize(_yjacp, _n*(1+_n+_npar));
  m_resize(_yy, _n, _n);
  px_resize(_ppivot, _n);

  if ( _sa ) {
    v_resize(_y1, _n*(1+_n+_npar));
    v_resize(_y2, _n*(1+_n+_npar));
  } else {
    v_resize(_y1, _n);
    v_resize(_y2, _n);
  }

}

//--------------------------------------------------------------------------
void eigenvalue(MATP A)
{

  MATP S, Q;
  VECP E_re, E_im;
  int i, imin, imax;
  
  S = m_get(A->m, A->m);
  m_copy(A, S);
  Q = m_get(A->m, A->m);
  schur(S, Q);
  E_re = v_get(A->m);
  E_im = v_get(A->m);
  schur_evals(S, E_re, E_im);
  imin = imax = 0;
  for ( i = 1; i < (int) A->m; i++ ) {
    if ( E_re[i] <  E_re[imin] )
      imin = i;
    if ( E_re[i] >  E_re[imax] )
      imax = i;
  }
  printf("l-min: %g+/-%g*j, l-max: %g+/-%g*j\n", 
	 E_re[imin], E_im[imin], E_re[imax], E_im[imax]);
  V_FREE(E_re);
  V_FREE(E_im);
  M_FREE(Q);
  M_FREE(S);

}

//--------------------------------------------------------------------------
void Omu_IntIMP::jac(double t, VECP y)
{

  bool sa_old;
  int  i, j;

  v_zero(_yjac);
  v_zero(_yjacp);

  for ( i = 0; i < _n; i++ ) {
    _yjac[i] = y[i];
    _yjac[(1+_nd+i)*_n+i] = 1.0;
  }

  sa_old = _sa;
  _sa = 1;
  syseq(t, _yjac, _u, _yjacp);
  _sa = sa_old;

  for ( i = 0; i < _n; i++ ) 
    for ( j = 0; j < _n; j++ )
      _yy[i][j] = _yjacp[(1+_nd+j)*_n+i];

  _jac_evals++;

  //  eigenvalue(_yy);

}

//--------------------------------------------------------------------------
// _IMP==0: BW Euler
// _IMP==1: IMP
void Omu_IntIMP::step(double tstep, double dt, VECP y)
{
  bool   ok, sa_old = _sa;
  int i, j, inewton;
  double t, cmeth;

  // IMP or BW Euler?
  cmeth = _IMP==1 ? 2.0 : 1.0;

  _sa = 0;
  for ( i = 0; i < _n; i++ )
    _y0[i] = y[i];
  v_zero(_k1);
  v_zero(_fh);
  t = tstep+dt/cmeth;

  // modified Newton's method    
  for ( inewton = 0; inewton < _maxiters; inewton++ ) {
    v_mltadd(_y0, _k1, dt/cmeth, _y);
    syseq(t, _y, _u, _fh);
    v_sub(_k1, _fh, _res);
    
    // check for convergence
    for ( ok = true, i = 0 ; i < _n; i++ )
      ok = ok && ( fabs(_res[i])<0.1*(_atol+fabs(_k1[i])) );
    if ( ok ) 
      break;
    else if ( inewton == _maxiters-1 )
      error(E_INTERN, 
	    "Omu_IntIMP::ode_solve Newton method failed to converge");
    
    // (re)calculate and factorize Jacobian
    if ( !(inewton % _modnewtonsteps) ) {
      jac(t, _y);
      sm_mlt(-dt/cmeth, _yy, _yyn);
      for ( i = 0; i < _n; i++ )
	_yyn[i][i] += 1.0;
      LUfactor(_yyn, _ppivot);
    }
    
    LUsolve(_yyn, _ppivot, _res, _fh);
    v_sub(_k1, _fh, _k1);
  }
  
  // calculate step 
  v_mltadd(_y0, _k1, dt, _y);
  for ( i = 0; i < _n; i++ )
    y[i] = _y[i];
  
  // sensitivities
  if ( sa_old ) {
    v_mltadd(_y0, _k1, dt/cmeth, _y);
    
    // (re)calculate and factorize Jacobian
    jac(t, _y);
    sm_mlt(-dt/cmeth, _yy, _yyn);
    for ( i = 0; i < _n; i++ )
      _yyn[i][i] += 1.0;
    LUfactor(_yyn, _ppivot);
    
    sm_mlt(dt/cmeth, _yy, _yy);
    for ( i = 0; i < _n; i++ )
      _yy[i][i] += 1.0;
    
    for ( i = 0; i < _n+_npar; i++ ) {
      for ( j = 0; j < _n; j++ )
	_res[j] = y[_n*(1+i)+j];
      mv_mlt(_yy, _res, _fh);
      if ( i >= _n )
	for ( j = 0; j < _n; j++ )
	  _fh[j] += dt*_yjacp[_n*(1+i)+j];
      
      LUsolve(_yyn, _ppivot, _fh, _fh);
      for ( j = 0; j < _n; j++ )
	y[_n*(1+i)+j] = _fh[j];
    }
  }

  _sa = sa_old; 
}

//--------------------------------------------------------------------------
void Omu_IntIMP::ode_solve(double tstart, VECP y, const VECP u, double tend)
{

  double t, dt, err, tol, ynorm, dtnew;
  int i;

  _npar = u->dim;
  realloc();
  v_copy(u, _u);

  if ( _stepsize > 0.0 ) {
      // with fixed step size 
      t = tstart;
      while ( t < tend ) {
	  dt = _stepsize;
	  if ( t+dt > tend ) 
	      dt = tend-t;
	  step(t, dt, y);
	  t += _stepsize;
      }
  } else {
    // with step size control
    if ( _hinit > 0.0 )
      dt = _hinit;
    else
      dt = _dt;
    if ( dt == 0.0 )
      dt = (tend-tstart)/10.0;

    t = tstart;
    while ( t < tend ) {
      _dt = dt;
      if ( dt < 10.0*MACHEPS*t ) 
	error(E_INTERN, 
	      "Omu_IntIMP::ode_solve step size too small");
      if ( t+dt > tend ) 
	dt = tend-t;

      // BW Euler step
      for ( i = 0; i < (int) _y1->dim; i++ )
	_y1[i] = y[i];
      _IMP = 0;
      step(t, dt, _y1);

      // IMP step
      for ( i = 0; i < (int) _y2->dim; i++ )
	_y2[i] = y[i];
      _IMP = 1;
      step(t, dt, _y2);
      
      // local error estimation
      for ( i = 0, err = 0.0, ynorm = 0.0; i < _n; i++ ) {
	err += square(_y1[i]-_y2[i]);
	ynorm += square(_y2[i]);
      }
      err = sqrt(err);
      ynorm = sqrt(ynorm);
      tol = _atol+_rtol*ynorm;
      dtnew = 0.9*dt*sqrt(tol/(MACHEPS+err));

      if ( err > tol ) {
	// reject step
	dt = max(0.25*dt, dtnew);
      } else {
	// accept step
	for ( i = 0; i < (int) _y2->dim; i++ )
	  y[i] = _y2[i];
	t += dt;
	dt = max(dt, min(dtnew, 4.0*dt));
      }
    }
  }

}

//========================================================================
