/*
 * Omu_IntODE.C --
 *   -- class implementation
 *
 * rf, 10/2/96
 */

/*
    Copyright (C) 1997--2002  Ruediger Franke

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

#include <adutils.h>

#include <If_Class.h>

#include "Omu_IntODE.h"

//--------------------------------------------------------------------------
Omu_IntODE::Omu_IntODE()
{
  _sys = NULL;
  _xc_ptr = NULL;
  _Fc_ptr = NULL;

  _y = v_get(1);
  _u = v_get(1);

  _x = v_get(1);
  _v = v_get(1);
  _X = m_get(1, 1);
  _Y = m_get(1, 1);

  _Yx = m_get(1, 1);
  _Yu = m_get(1, 1);
}

//--------------------------------------------------------------------------
Omu_IntODE::~Omu_IntODE()
{
  m_free(_Y);
  m_free(_X);
  v_free(_v);
  v_free(_x);

  v_free(_u);
  v_free(_y);

  m_free(_Yx);
  m_free(_Yu);
}

//--------------------------------------------------------------------------
void Omu_IntODE::init_stage(int k,
			    const Omu_States &x, const Omu_Vector &u,
			    bool sa)
{
  if (!x.D_is_const) {
    m_error(E_FORMAT,
	    "Omu_IntODE::init_stage, which was called for non const dF/dxp");
  }
  if (x.na > 0) {
    m_error(E_FORMAT,
	    "Omu_IntODE::init_stage, which was called for algebraic states");
  }
  Omu_Integrator::init_stage(k, x, u, sa);
  realloc();
}

//--------------------------------------------------------------------------
void Omu_IntODE::realloc()
{
  if (_xcp->dim == _nx + _nd && _uc->dim == _nu && _u->dim == _nd + _nu)
    return;

  int neq = _n * (1 + _nx + _nu);
  v_resize(_y, neq);
  v_resize(_u, _nd + _nu);

  //
  // variables for ADOL-C
  //

  v_resize(_x, _nd + 2 * _n + _nu);
  v_resize(_v, _nd + 2 * _n + _nu);
  m_resize(_X, _nd + 2 * _n + _nu, _nx + _nu);
  m_resize(_Y, _n, _nx + _nu);

  //
  // variables for low level _sys->continuous callback
  //
  v_resize(_uc, _nu);
  _xcp.realloc(_nd + _n, _nx, _nu);
  m_resize(_Yx, _nd + _n, _nx);
  m_resize(_Yu, _nd + _n, _nu);
}

//--------------------------------------------------------------------------
void Omu_IntODE::solve(int kk, Real tstart, Real tend,
		       const Omu_States &x, const Omu_Vector &u,
		       Omu_Program *sys, Omu_DepVec &Fc, Omu_SVec &xc)
{
  int i, j;

  _sys = sys;	// propagate to syseq()
  _xc_ptr = &xc;
  _Fc_ptr = &Fc;

  v_zero(_y);
    
  for (i = 0; i < _nd; i++) {
    _u[i] = xc[i];
  }
  for (i = 0; i < _n; i++) {
    _y[i] = xc[_nd + i];		// initial states
  }
  for (i = 0; i < _nu; i++) {
    _u[_nd + i] = u[i];
  }

  v_zero(_xcp);	// time derivatives passed to continuous

  if (_sa) {
    for (i = 0; i < _n; i++) {
      for (j = 0; j < _nx; j++) {
	_y[(1 + j) * _n + i] = xc.Sx[_nd + i][j];
      }
      for (j = 0; j < _nu; j++) {
	_y[(1 + _nx + j) * _n + i] = xc.Su[_nd + i][j];
      }
    }
    m_zero(_xcp.Sx);
    m_zero(_xcp.Su);
  }

  _kk = kk;	// propagate to syseq()

  ode_solve(tstart, _y, _u, tend);

  for (i = 0; i < _n; i++) {
    xc[_nd + i] = _y[i];
  }

  if (_sa) {
    for (i = 0; i < _n; i++) {
      for (j = 0; j < _nx; j++) {
	xc.Sx[_nd + i][j] = _y[(1 + j) * _n + i];
      }
      for (j = 0; j < _nu; j++) {
	xc.Su[_nd + i][j] = _y[(1 + _nx + j) * _n + i];
      }
    }
  }
}

// default implementation calling low-level _sys->continuous()
//--------------------------------------------------------------------------
void Omu_IntODE::syseq(Real t, const VECP y, const VECP u,
		       VECP f)
{
  if (!_sys->has_low_level_continuous()) {
    // call faster version if no user defined low-level continuous
    syseq_forward(t, y, u, f);
    return;
  }

  int i, j;
  Omu_SVec &xc = *_xc_ptr;
  Omu_DepVec &Fc = *_Fc_ptr;

  //
  // prepare call arguments
  //

  for (i = 0; i < _nd; i++) {
    xc[i] = u[i];
  }
  for (i = 0; i < _n; i++) {
    xc[_nd + i] = y[i];
  }
  for (i = 0; i < _nu; i++) {
    _uc[i] = u[_nd + i];
  }
      
  if (_sa) {
    for (i = _nd; i < _nxt; i++) {
      for (j = 0; j < _nx; j++) {
	xc.Sx[i][j] = y[(1 + j) * _n + i - _nd];
      }
      for (j = 0; j < _nu; j++) {
	xc.Su[i][j] = y[(1 + _nx + j) * _n + i - _nd];
      }
    }
  }

  //
  // evaluate residual
  //

  Fc.set_required_J(_sa); // an integrator may request the Jacobian

  _sys->continuous(_kk, t, xc, _uc, _xcp, Fc);

  for (i = _nd; i < _nxt; i++) {
    // f = F * -(dF/dxp)^(-1)
    f[i - _nd] = Fc[i] / -Fc.Jxp[i][i];
  }
      
  if (_sa) {
    m_mlt(Fc.Jx, xc.Sx, _Yx);
    m_mlt(Fc.Jx, xc.Su, _Yu);
    m_add(_Yu, Fc.Ju, _Yu);

    for (i = _nd; i < _nxt; i++) {
      for (j = 0; j < _nx; j++) {
	f[(1 + j) * _n + i - _nd] = _Yx[i][j] / -Fc.Jxp[i][i];
      }
      for (j = 0; j < _nu; j++) {
	f[(1 + _nx + j) * _n + i - _nd] = _Yu[i][j] / -Fc.Jxp[i][i];
      }
    }
  }

  _res_evals++;
  if (_sa)
    _sen_evals++;
}

// alternative implementation calling high-level _sys->continuous
//--------------------------------------------------------------------------
void Omu_IntODE::syseq_forward(Real t, const VECP y, const VECP u,
			       VECP f)
{
  int i, j;
  Omu_DepVec &Fc = *_Fc_ptr;

  //
  // form a vector of independent variables
  //

  for (i = 0; i < _nd; i++) {
    _x[i] = u[i];
  }
  for (i = 0; i < _n; i++) {
    _x[_nd + i] = y[i];
    _x[_nd + _n + i] = 0.0;	// yprime[i]
  }
  for (i = 0; i < _nu; i++) {
    _x[_nd + 2 * _n + i] = u[_nd + i];
  }
      
  //
  // evaluate residual
  //

  adoublev ax(_nd + _n);
  adoublev axp(_nd + _n);
  adoublev au(_nu);
  adoublev aF(_nd + _n);

  for (i = 0; i < _nd; i++)
    axp[i] = 0.0;
  for (i = _nd; i < _nxt; i++)
    aF[i] = 0.0;

  if (_sa)
    trace_on(3);	// tape 3
  ax <<= _x->ve;
  for (i = 0; i < _n; i++)
    axp[_nd + i] <<= _x->ve[_nd + _n + i];
  au <<= _x->ve + _nd + 2 * _n;

  _sys->continuous(_kk, t, ax, au, axp, aF);

  for (i = _nd; i < _nxt; i++) {
    aF[i] >>= f[i - _nd];
    f[i - _nd] /= -Fc.Jxp[i][i];
  }
      
  if (_sa) {
    trace_off();

    int nindep = _nd + 2 * _n + _nu;
    int npar = _nx + _nu;

    m_zero(_X);
    for (i = 0; i < _nd; i++) {
      _X[i][i] = 1.0;
    }
    for (i = 0; i < _n; i++) {
      for (j = 0; j < npar; j++) {
	_X[_nd + i][j] = y[(1 + j) * _n + i];
	_X[_nd + _n + i][j] = 0.0; // yprime[(1 + j) * _n + i];
      }
    }
    for (i = 0; i < _nu; i++) {
      _X[_nd + 2 * _n + i][_nd + _n + i] = 1.0;
    }
      
    forward(3, _n, nindep, npar, _x->ve, _X->me, f->ve, _Y->me);

    for (i = _nd; i < _nxt; i++) {
      f[i - _nd] /= -Fc.Jxp[i][i];
      for (j = 0; j < npar; j++) {
	f[(1 + j) * _n + i - _nd] = _Y[i - _nd][j] / -Fc.Jxp[i][i];
      }
    }
  }

  _res_evals++;
  if (_sa)
    _sen_evals++;
}


//========================================================================
