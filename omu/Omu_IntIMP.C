/*
 * Omu_IntIMP.C --
 *   -- class implementation
 *
 * E. Arnold   1999-04-12
 *             2000-05-12 _rtol, _atol -> Omu_Integrator 
 *             2000-05-30 step size control
 *             2001-08-16 prg_int_nsteps --> _stepsize
 *             2003-01-02 modified Newton's method from Hairer/Wanner
 *
 */

/*
    Copyright (C) 1999--2003  Eckhard Arnold

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
#include <meschach/matrix2.h>
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
    _z     = v_get(1);
    _zp    = v_get(1);
    _yjac  = v_get(1);
    _yjacp = v_get(1);
    _yy    = m_get(1, 1);
    _yyn   = m_get(1, 1);
    _ppivot= px_get(1);
    _y1    = v_get(1);
    _y2    = v_get(1);

    _maxiters = 10;
    _hinit = 0.0;
    _dt = 0.0;
    _ixgz = -1;
    _correct_der = 0;  // 0 much faster!
    _kappa = 0.1;
    _max_modnewtonsteps = 10;

    _ifList.append(new If_Int("prg_int_maxiters", &_maxiters));
    _ifList.append(new If_Real("prg_int_hinit", &_hinit));
    _ifList.append(new If_Int("prg_int_ixgz", &_ixgz));
    _ifList.append(new If_Bool("prg_int_correctder", &_correct_der));
    _ifList.append(new If_Real("prg_int_kappa", &_kappa));
    _ifList.append(new If_Int("prg_int_modnewtonsteps", &_max_modnewtonsteps));

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
    V_FREE(_z);
    V_FREE(_zp);
    V_FREE(_yjac);
    V_FREE(_yjacp);
    M_FREE(_yy);
    M_FREE(_yyn);
    PX_FREE(_ppivot);
    V_FREE(_y1);
    V_FREE(_y2);
}

//--------------------------------------------------------------------------
void Omu_IntIMP::resize()
{
    int nn = _n*(1+_n+_npar);

    if ( (int)_u->dim != _npar || (int)_y->dim != _n || 
	 (int)_yjac->dim != nn || (int)_y1->dim != _n ) {
	v_resize(_u, _npar);
	v_resize(_y, _n);
	v_resize(_y0, _n);
	v_resize(_k1, _n);
	v_zero(_k1);
	v_resize(_res, _n);
	v_resize(_fh, _n);
	v_resize(_z, _n);
	v_resize(_zp, _n);
	v_resize(_yjac, nn);
	v_resize(_yjacp, nn);
	m_resize(_yy, _n, _n);
	px_resize(_ppivot, _n);
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
// _k1 - ODE r.h.s 
// _yy - ODE r.h.s. Jacobian
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
    // non-negative state variables
    if ( _ixgz >= 0 )
	for ( i = _ixgz; i < _n; i++ )
	    _yjac[i] = max(_yjac[i], 0.0);

    sa_old = _sa;
    _sa = 1;
    syseq(t, _yjac, _u, _yjacp);
    _sa = sa_old;

    // _yjac, _yjacp: derivatives w.r.t.
    //                discrete states   _n ... _n*(1+_nd)-1 
    //                continuous states _n*(1+_nd)..._n*(1+_nd+_n)-1 
    //                controls          _n*(1+_nd+_n)...(end)
    for ( i = 0; i < _n; i++ ) {
	_k1[i] = _yjacp[i];
	for ( j = 0; j < _n; j++ )
	    _yy[i][j] = _yjacp[(1+_nd+j)*_n+i];
    }
    _jac_evals++;

    //  eigenvalue(_yy);
}

//--------------------------------------------------------------------------
// IMP==0: BW Euler (for step size control only!)
// IMP==1: IMP
int Omu_IntIMP::step(int IMP, double tstep, double dt, VECP y, double tol)
{
    int i, inewton, result;
    double cmeth, del_z, del_z_m, norm_res, theta;

    // IMP or BW Euler?
    cmeth = IMP==1 ? 2.0 : 1.0;

    for ( i = 0; i < _n; i++ )
	_y0[i] = y[i];
    
    // factorize I-dt/cmeth*_yy
    sm_mlt(-dt/cmeth, _yy, _yyn);
    for ( i = 0; i < _n; i++ )
	_yyn[i][i] += 1.0;
    
    m_catchall(LUfactor(_yyn, _ppivot);,
	       fprintf(stderr, 
		       "Omu_IntIMP::ode_solve singular Jacobian");
	       m_error(E_SING, 
		       "Omu_IntIMP::ode_solve singular Jacobian");
	);

    // initialize _z
    sv_mlt(dt/cmeth, _k1, _z);

    for ( inewton = 0, result = 1; inewton < _maxiters; inewton++ ) {

	v_add(_y0, _z, _y);
	if ( _ixgz >= 0 )
	    for ( i = _ixgz; i < _n; i++ )
		_y[i] = max(_y[i], 0.0);

	// ODE r.h.s.
	syseq(tstep+dt/cmeth, _y, _u, _fh);
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
	    for ( i = 0; i < _n; i++ )
		y[i] = _y[i];
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

    // return convergence status (0: yes, 1: no)
    return result;
}

//--------------------------------------------------------------------------
void Omu_IntIMP::init_stage(int k,
			    const Omu_States &x, const Omu_Vector &u,
			    bool sa)
{
    Omu_IntODE::init_stage(k, x, u, sa);

    if ( _nv ) 
	m_error(E_SIZES, "Omu_IntIMP::init_stage");
    
    assert( _maxiters >= 3 );

    _npar = _nd + _nu;
    resize();

    _eta = Inf;
}

//--------------------------------------------------------------------------
void Omu_IntIMP::ode_solve(double tstart, VECP y, const VECP u, double tend)
{
    double t, dt, err, tol, ynorm, dtnew;
    bool sa_old;
    int i, j;
    int modnewtonsteps;

    v_copy(u, _u);

    // non-negative state variables
    if ( _ixgz >= 0 )
	for ( i = _ixgz; i < _n; i++ )
	    y[i] = max(y[i], 0.0);

    // initial step size
    if ( _hinit > 0.0 )
	dt = _hinit;
    else
	dt = _dt;
    if ( dt == 0.0 )
	dt = (tend-tstart)/10.0;

    // calculate derivatives?
    sa_old = _sa;
    _sa = 0;

    // begin main loop
    t = tstart;
    modnewtonsteps = -1;
    while ( t < tend ) {
	_dt = dt;  // keep last 'regular' step size
	if ( dt < 10.0*MACHEPS*t ) {
	    fprintf(stderr, "Omu_IntIMP::ode_solve step size too small");
	    m_error(E_INTERN, 
		    "Omu_IntIMP::ode_solve step size too small");
	}
	if ( t+dt > tend ) 
	    dt = tend-t;

	// calculate ODE r.h.s. _k1 and Jacobian _yy
	if ( modnewtonsteps < 0 || modnewtonsteps >= _max_modnewtonsteps ) {
	    jac(t, y);
	    modnewtonsteps = 0;
	} else
	    modnewtonsteps++;

	// BW Euler step and IMP step
	for ( i = 0; i < _n; i++ ) {
	    _y1[i] = y[i];
	    _y2[i] = y[i];
	}
	// tolerance threshold for Newton's method
	for ( i = 0, tol = 0.0; i < _n; i++ )
	    tol += square(y[i]);
	tol = _atol+_rtol*sqrt(tol);

	if ( ( step(0, t, dt, _y1, tol) > 0 ) || 
	     ( step(1, t, dt, _y2, tol) > 0 ) ) {
	    if ( modnewtonsteps > 0 ) // Newton's method did not converge: 
		modnewtonsteps = -1;  // force recalculation of Jacobian
	    else
		dt *= 0.5;            // reduce step size
	    continue;    
	}
      
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
	    t += dt;
	    dt = min(dtnew, 4.0*dt);
	    for ( i = 0; i < _n; i++ ) {
		// interval mid-point
		_y[i] = (y[i]+_y2[i])/2.0;
		// step
		y[i] = _y2[i];
	    }
	    // non-negative state variables
	    if ( _ixgz >= 0 )
		for ( i = _ixgz; i < _n; i++ )
		    y[i] = max(y[i], 0.0);

	    // derivative of step function
	    if ( sa_old ) {
		if ( _correct_der ) { 
		    // calculate ODE r.h.s. _k1 and Jacobian _yy
		    jac(t-dt/2.0, _y);
		    // factorize I-dt/2*_yy
		    sm_mlt(-dt/2.0, _yy, _yyn);
		    for ( i = 0; i < _n; i++ )
			_yyn[i][i] += 1.0;
		    
		    m_catchall(LUfactor(_yyn, _ppivot);,
			       fprintf(stderr, 
				       "Omu_IntIMP::ode_solve singular Jacobian");
			       m_error(E_SING, 
				       "Omu_IntIMP::ode_solve singular Jacobian");
			);
		}    
		// derivatives
		for ( i = 0; i < _n+_npar; i++ ) {
		    for ( j = 0; j < _n; j++ )
			_res[j] = y[_n*(1+i)+j];
		    mv_mlt(_yy, _res, _fh);
		    if ( i < _nd || i >= _nd+_n )
			for ( j = 0; j < _n; j++ )
			    _fh[j] += _yjacp[_n*(1+i)+j];
		    
		    LUsolve(_yyn, _ppivot, _fh, _fh);
		    for ( j = 0; j < _n; j++ )
			y[_n*(1+i)+j] += dt*_fh[j];
		}
	    }
	}
    } // end main loop
    
    _sa = sa_old; 

}

//========================================================================
