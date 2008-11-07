/*
 * Omu_IntIMP.C --
 *   -- class implementation
 *
 * E. Arnold   1999-04-12
 *             2000-05-12 _rtol, _atol -> Omu_Integrator 
 *             2000-05-30 step size control
 *             2001-08-16 prg_int_nsteps --> _stepsize
 *             2003-01-02 modified Newton's method from Hairer/Wanner
 *             2003-08-25 Omu_Integrator
 *             2003-09-06 step size control by Richardson extrapolation
 * R. Franke   2008-11-06 rework fixed step size
 *
 */

/*
    Copyright (C) 1999--2008  Eckhard Arnold, Ruediger Franke

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

#include <If_Int.h>
#include <If_Real.h>
#include <If_Bool.h>
#include <If_Class.h>

#include "Omu_IntIMP.h"

extern "C" {
#include <meschach/addon_hqp.h>
#include <meschach/matrix2.h>
}

IF_CLASS_DEFINE("IMP", Omu_IntIMP, Omu_Integrator);

//--------------------------------------------------------------------------
Omu_IntIMP::Omu_IntIMP()
{
    _x     = v_get(1);
    _y     = v_get(1);
    _k1    = v_get(1);
    _y0    = v_get(1);
    _res   = v_get(1);
    _fh    = v_get(1);
    _z     = v_get(1);
    _zp    = v_get(1);
    _z_old = v_get(1);
    _xs    = m_get(1, 1);
    _yy    = m_get(1, 1);
    _yyn   = m_get(1, 1);
    _yq    = m_get(1, 1);
    _yq1   = m_get(1, 1);
    _yq2   = m_get(1, 1);
    _yqq   = m_get(1, 1);
    _ppivot= px_get(1);
    _y1    = v_get(1);
    _y2    = v_get(1);

    _maxiters = 20;
    _hinit = 0.0;
    _dt = 0.0;
    _ixgz = -1;
    _kappa = 0.1;
    _max_modnewtonsteps = 10;
    _correct_der = 1;  // 0 much faster for many problems!
    _max_sing = 1;

    _ifList.append(new If_Int("prg_int_maxiters", &_maxiters));
    _ifList.append(new If_Real("prg_int_hinit", &_hinit));
    _ifList.append(new If_Int("prg_int_ixgz", &_ixgz));
    _ifList.append(new If_Real("prg_int_kappa", &_kappa));
    _ifList.append(new If_Int("prg_int_modnewtonsteps", &_max_modnewtonsteps));
    _ifList.append(new If_Bool("prg_int_correctder", &_correct_der)); 

    _res_evals = 0;
    _jac_evals = 0;
    _sen_evals = 0;
}

//--------------------------------------------------------------------------
Omu_IntIMP::~Omu_IntIMP()
{
    V_FREE(_x);
    V_FREE(_y);
    V_FREE(_k1);
    V_FREE(_y0);
    V_FREE(_res);
    V_FREE(_fh);
    V_FREE(_z);
    V_FREE(_zp);
    V_FREE(_z_old);
    M_FREE(_xs);
    M_FREE(_yy);
    M_FREE(_yyn);
    M_FREE(_yq);
    M_FREE(_yq1);
    M_FREE(_yq2);
    M_FREE(_yqq);
    PX_FREE(_ppivot);
    V_FREE(_y1);
    V_FREE(_y2);
}

//--------------------------------------------------------------------------
void Omu_IntIMP::resize()
{
    if ( (int) _dxc->dim != _n )
	_dxc.resize(_n, 0, 0, _nq);

    if ( (int)_y->dim != _n || (int)_yq->n != _nq ) {
	v_resize(_x, _n);
	v_resize(_y, _n);
	v_resize(_y0, _n);
	v_resize(_k1, _n);
	v_resize(_res, _n);
	v_resize(_fh, _n);
	v_resize(_z, _n);
	v_resize(_zp, _n);
	v_resize(_z_old, _n);
	m_resize(_xs, _n, _nq);
	m_resize(_yy, _n, _n);
	m_resize(_yyn, _n, _n);
	m_resize(_yq, _n, _nq);
	m_resize(_yq1, _n, _nq);
	m_resize(_yq2, _n, _nq);
	m_resize(_yqq, _n, _nq);
	px_resize(_ppivot, _n);
	v_resize(_y1, _n);
	v_resize(_y2, _n);

	v_zero(_z); 
	v_zero(_k1);
	_modnewtonsteps = 0;
	_nsing = 0;
    } 
//    _modnewtonsteps = 0;
//    v_zero(_k1);
}

//--------------------------------------------------------------------------
// calculate and print eigenvalues of matrix A
static void eigenvalue(MATP A)
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

//-----------------------------------------------------------------------------
// calculate ODE rhs
void Omu_IntIMP::sys(double t, VECP x, VECP xp)
{
    int i;
    Omu_StateVec &xc = *_xc_ptr;
    Omu_DependentVec &Fc = *_Fc_ptr;
    Omu_Vec &q = *_q_ptr;

    // non-negative state variables
    if ( _ixgz >= 0 )
	for ( i = _ixgz; i < _n; i++ )
	    x[i] = max(x[i], 0.0);

    // prepare call arguments
    for ( i = 0; i < _n; i++ ) {
	xc[i] = x[i];
	_dxc[i] = 0.0;
    }
    Fc.set_required_J(false);

    // evaluate residual
    residual(_kk, t, xc, _dxc, q, Fc);
    _res_evals++;

    // read and return result
    for ( i = 0; i < _n; i++ )
	xp[i] = Fc[i];
}

//-----------------------------------------------------------------------------
// calculate Jacobian of ODE rhs wrt x
void Omu_IntIMP::sys_jac(double t, VECP x, VECP xp, MATP fx)
{
    sys_jac(t, x, xp, fx, MNULL);
}

//-----------------------------------------------------------------------------
// calculate Jacobian of ODE rhs wrt x and q
void Omu_IntIMP::sys_jac(double t, VECP x, VECP xp, MATP fx, MATP fq)
{
    int i;
    Omu_StateVec &xc = *_xc_ptr;
    Omu_DependentVec &Fc = *_Fc_ptr;
    Omu_Vec &q = *_q_ptr;

    // non-negative state variables
    if ( _ixgz >= 0 )
	for ( i = _ixgz; i < _n; i++ )
	    x[i] = max(x[i], 0.0);

    // prepare call arguments
    for ( i = 0; i < _n; i++ ) {
	xc[i] = x[i];
	_dxc[i] = 0.0;
    }
    Fc.set_required_J(true);

    // evaluate residual
    residual(_kk, t, xc, _dxc, q, Fc);
    _jac_evals++;

    // read and return result
    if ( xp != VNULL )
  	for ( i = 0; i < _n; i++ )
  	    xp[i] = Fc[i];

    m_copy(Fc.Jx, fx);

    if ( fq != MNULL ) {
	m_copy(Fc.Jq, fq);
	_sen_evals++;
    }
}

//-----------------------------------------------------------------------------
// LU factorize I*gamma-delta*fx
// modifies _ppivot
int Omu_IntIMP::lufac_jac(double gamma, double delta, MATP fx)
{
    int i, result = 0;

    sm_mlt(-delta, fx, fx);
    for ( i = 0; i < _n; i++ ) 
	    fx[i][i] += gamma;

    m_catchall(LUfactor(fx, _ppivot);,
	       result = -1;
	       _nsing++; 
	       if ( _nsing >= _max_sing ) {
		   m_error(E_CONV, 
			   "Omu_IntIMP::lufac_jac singular Jacobian");
	       }
	);

    if ( result == 0 )
	_nsing = 0;

    return result;
}

//--------------------------------------------------------------------------
void Omu_IntIMP::init(int k,
		       const Omu_StateVec &xc, const Omu_Vec &q,
		       const Omu_DependentVec &Fc, bool sa)
{
    assert( _maxiters >= 3 );
    resize();
    _eta = Inf;
}

//--------------------------------------------------------------------------
void Omu_IntIMP::solve(int kk, double tstart, double tend,
			Omu_StateVec &xc, Omu_StateVec &dxc, Omu_Vec &q,
			Omu_DependentVec &Fc)
{
    double t, dt, err, tol, ynorm, dtnew, dtold, facmax;
    int i;
    double dtmin = _min_stepsize > 0.0? _min_stepsize: 10.0*MACHEPS*(tend - tstart);

    _kk = kk;
    _xc_ptr = &xc; 
    _Fc_ptr = &Fc;
    _q_ptr  = &q;

    v_copy(xc, _x);
    m_copy(xc.Sq, _xs);

    if ( _stepsize > 0.0 ) {    // with fixed step size
        if (_min_stepsize <= 0.0)
            dtmin = _stepsize; // use _stepsize if no _min_stepsize specified
	_max_sing = 5;
        _modnewtonsteps = -1; // force recalculation of Jacobian initially
	t = tstart;
        dt = _stepsize;
	while ( t < tend ) {
	    if ( t+dt > tend ) 
		dt = tend-t;
            m_catchall( fixedstep(t, dt, _x, _xs);,
                // reduce step size as fixedstep has failed
                dt *= 0.5;
	        if ( dt < dtmin )
                    m_error(E_CONV, 
		            "Omu_IntIMP::solve step size below min_stepsize");
                _modnewtonsteps = -1; // force recalculation of Jacobian
            );
            if (_modnewtonsteps >= 0) {
                // accept successful step
	        t += dt;
                // readjust dt to mesh defined by _stepsize if dt was reduced
                if (dt != _stepsize) {
                    for (dt = tstart; dt - t < dtmin; dt += _stepsize);
                    dt -= t;
                }
            }    
	}
    } else {    // with step size control
	// allow _max_sing singular Jacobians
	_max_sing = 5;

	// initial step size
	if ( _hinit > 0.0 )
	    dt = _hinit;
	else
	    dt = _dt;
	if ( dt == 0.0 )
	    dt = (tend-tstart)/10.0;

	// begin main loop
	t = tstart;
	_modnewtonsteps = -1;
	facmax = 4.0;

	while ( t < tend ) {
	    _dt = dt;  // keep last 'regular' step size
	    if ( dt < dtmin ) {
		m_error(E_CONV, 
			"Omu_IntIMP::solve step size too small");
	    }
	    if ( t+dt > tend ) 
		dt = tend-t;

	    // calculate ODE r.h.s. _k1 and Jacobian _yy
	    if ( _modnewtonsteps < 0 || 
		 _modnewtonsteps >= _max_modnewtonsteps ) {
		sys_jac(t, _x, _k1, _yy, _yq);
		_modnewtonsteps = 0;
	    } else
		_modnewtonsteps++;

	    // BW Euler step and IMP step
	    v_copy(_x, _y1);
	    v_copy(_x, _y2);
	    m_copy(_xs, _yq1);
	    m_copy(_xs, _yq2);

	    // tolerance threshold for Newton's method
	    tol = _atol+_rtol*v_norm2(_x);
	    // step size control using Richardson extrapolation
	    if ( ( step(t, dt, _y1, _yq1, tol) > 0 ) || 
		 ( step(t, dt/2.0, _y2, _yq2, tol) > 0 ) || 
		 ( step(t, dt/2.0, _y2, _yq2, tol) > 0 ) ) {
		if ( _modnewtonsteps > 0 ) // Newton's method did not converge:
		    _modnewtonsteps = -1;  // force recalculation of Jacobian
		else
		    dt *= 0.5;             // reduce step size
		continue;    
	    }

	    // local error estimation
	    for ( i = 0, err = 0.0; i < _n; i++ ) {
		err = max(err, fabs(_y2[i]-_y1[i])/
		    max(fabs(_y2[i]), max(fabs(_x[i]), 1.0e-6)));
	    }
	    err = 1.0/3.0*err; // 1/(2^p-1)
	    // new step size
	    ynorm = v_norm2(_y2);
	    tol = _atol+_rtol*ynorm;
	    dtnew = dt*min(facmax, max(0.25, 0.9*pow(tol/err, 1.0/3.0)));
	    if ( err > tol ) {
		// reject step
		dt = dtnew;
		facmax = 1.0;
	    } else {
		// accept step
		t += dt;
		dt = dtnew;
		facmax = 4.0;
		// step
		for ( i = 0; i < _n; i++ ) 
		    _x[i] = _y2[i]+(_y2[i]-_y1[i])/3.0; // 1/(2^p-1)

		// derivative of step function
		if ( _sa ) {
  		    m_sub(_yq2, _yq1, _xs);
		    ms_mltadd(_yq2, _xs, (1.0/3.0), _xs);
		}
	    }
	} // end main loop
    }

    v_copy(_x, xc);
    m_copy(_xs, xc.Sq);
}

//--------------------------------------------------------------------------
int Omu_IntIMP::step(double tstep, double dt, VECP y, MATP yq, double tol)
{
    int i, inewton, result;
    double cmeth, del_z, del_z_m = 0.0, norm_res, theta, dd;

    // IMP or BW Euler? cmeth = IMP==1 ? 2.0 : 1.0;
    cmeth = 2.0;

    v_copy(y, _y0);
    
    // factorize I-dt/cmeth*_yy
    m_copy(_yy, _yyn);
    lufac_jac(1.0, dt/cmeth, _yyn);

    // initialize _z
    sv_mlt(dt/cmeth, _k1, _z);

    for ( inewton = 0, result = 1; inewton < _maxiters; inewton++ ) {

	v_add(_y0, _z, _y);

	// ODE r.h.s.
	sys(tstep+dt/cmeth, _y, _fh);

	// residuum
	v_mltadd(_z, _fh, -dt/cmeth, _res);
	sv_mlt(-1.0, _res, _res);
	norm_res = v_norm2(_res);

	// solve linear system for Newton step
	LUsolve(_yyn, _ppivot, _res, _zp);

	// Newton step
	del_z = v_norm2(_zp);
	v_add(_z, _zp, _z);
    
	// convergence measure
	if ( inewton == 0 ) {
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
	    // calculate step 
	    v_mltadd(_y0, _z, cmeth, _y);
	    for ( i = 0; i < _n; i++ ) {
		dd = (y[i]+_y[i])/2.0;
		y[i] = _y[i];
		_y[i] = dd;
	    }
	    result = 0;
	    break;
	}

	// check for divergence
	if ( ( inewton >= 1 ) && ( ( theta >= 1.0 ) || 
	     ( pow(theta, _maxiters-1-inewton)/(1.0-theta)*del_z > 
	       _kappa*tol ) ) ) {
	    result = 1;
	    break;
	}
	del_z_m = del_z;
    }

    // sensitivities
    if ( _sa && ( result == 0 ) ) {
	if ( _correct_der ) {
	    // calculate ODE r.h.s. _k1 and Jacobian _yy
	    sys_jac(tstep+dt/2.0, _y, _k1, _yy, _yq);
	    _modnewtonsteps = 0;
	    // factorize I-dt/2*_yy
	    m_copy(_yy, _yyn);
	    lufac_jac(1.0, dt/2.0, _yyn);
	}
	m_mlt(_yy, yq, _yqq);
	m_add(_yqq, _yq, _yqq);
	LUsolveM(_yyn, _ppivot, _yqq, _yqq);
	ms_mltadd(yq, _yqq, dt, yq);
    }

    // return convergence status (0: yes, 1: no)
    return result;
}

//--------------------------------------------------------------------------
void Omu_IntIMP::fixedstep(double tstep, double dt, VECP y, MATP yq)
{
    int i, inewton;
    double t;

    for ( i = 0; i < _n; i++ )
	_y0[i] = y[i];

    // reuse _z
    v_zero(_fh);
    t = tstep+dt/2.0;

    // modified Newton's method    
    for ( inewton = 0; inewton < _maxiters; inewton++ ) {
	v_add(_y0, _z, _y);
	sys(t, _y, _fh);
	v_mltadd(_z, _fh, -dt/2.0, _res);

	// check for convergence
	if ( _modnewtonsteps >= 0 
             && v_norm2(_res) < 0.1*(_atol+_rtol*v_norm2(_z)) )
	    break;
	else if ( inewton == _maxiters-1 ) {
	    m_error(E_CONV, 
		    "Omu_IntIMP::fixedstep Newton method failed to converge");
	}

	// (re)calculate and factorize Jacobian
	if ( _modnewtonsteps < 0 || _modnewtonsteps >= _max_modnewtonsteps ) {
	    sys_jac(t, _y, _fh, _yy);
	    m_copy(_yy, _yyn);
	    lufac_jac(1.0, dt/2.0, _yyn);
	    _modnewtonsteps = 0;
	} else
	    _modnewtonsteps++;

	LUsolve(_yyn, _ppivot, _res, _fh);
	v_sub(_z, _fh, _z);
    }

    // sensitivities
    if ( _sa ) {
	v_add(_y0, _z, _y);
	// (re)calculate and factorize Jacobian
	sys_jac(t, _y, _fh, _yy, _yq);
	m_copy(_yy, _yyn);
	lufac_jac(1.0, dt/2.0, _yyn);
	_modnewtonsteps = 0;

	m_mlt(_yy, yq, _yqq);
	m_add(_yqq, _yq, _yqq);
	LUsolveM(_yyn, _ppivot, _yqq, _yqq);
	ms_mltadd(yq, _yqq, dt, yq);
    }

    // calculate step 
    v_mltadd(_y0, _z, 2.0, _y);
    v_copy(_y, y);
}

//========================================================================
