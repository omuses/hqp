/*
 * Omu_Integrator.C --
 *   -- class definition
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

#include "Omu_Integrator.h"

#include <If_Real.h>
#include <If_Int.h>
#include <If_Bool.h>

IF_BASE_DEFINE(Omu_Integrator);

//--------------------------------------------------------------------------
Omu_Integrator::Omu_Integrator()
{
  _K = 0;
  _k = 0;
  _kk = 0;
  _nxt = 0;
  _nd = 0;
  _nv = 0;
  _na = 0;
  _nx = 0;
  _nu = 0;
  _sa = false;
  _serr = false;
  _n = 0;
  _m = 0;
  _stepsize = 0.0;
  _rtol = 1e-8;
  _atol = 1e-8;
  _res_evals = 0;
  _sen_evals = 0;
  _jac_evals = 0;

  // initialize references to data required by residual()
  _sys = NULL;
  _x_ptr = NULL;
  _u_ptr = NULL;
  _Ft_ptr = NULL;
  _xt_ptr = NULL;

  _Fcs = NULL;

  _ifList.append(new If_Bool("prg_int_serr", &_serr));
  _ifList.append(new If_Real("prg_int_stepsize", &_stepsize));
  _ifList.append(new If_Real("prg_int_rtol", &_rtol));
  _ifList.append(new If_Real("prg_int_atol", &_atol));
  _ifList.append(new If_Int("prg_int_res_evals", &_res_evals));
  _ifList.append(new If_Int("prg_int_sen_evals", &_sen_evals));
  _ifList.append(new If_Int("prg_int_jac_evals", &_jac_evals));
}

//--------------------------------------------------------------------------
Omu_Integrator::~Omu_Integrator()
{
  delete [] _Fcs;
}

//--------------------------------------------------------------------------
void Omu_Integrator::setup_stages(const Omu_Program *sys)
{
  delete [] _Fcs;
  _K = sys->K();
  _Fcs = new Omu_DepVec [_K];
}

//--------------------------------------------------------------------------
void Omu_Integrator::setup_struct(int k,
				  const Omu_VariableVec &x,
				  const Omu_VariableVec &u,
				  const Omu_DependentVec &Ft)
{
  // initialize Jacobians for high-level integrator interface

  if (k >= _K) {
    m_error(E_INTERN, "Omu_Integrator::setup_struct"
	    " that was called with wrong integrator setup");
  }

  Omu_SVarVec &sx = (Omu_SVarVec &)x;
  int i;
  int nd = sx.nd;
  int n = x->dim - nd;
  int nx = x->dim - sx.nv;
  int nu = u->dim;
  Omu_DepVec &Fc = _Fcs[k];

  Fc.size(n, n, 0, n, 0, nx+nu);
  m_move(Ft.Jx, nd, nd, n, n, Fc.Jx, 0, 0);
  m_move(Ft.Jdx, nd, nd, n, n, Fc.Jdx, 0, 0);
  m_zero(Fc.Jq); // zero Jq wrt. continuous states as Jx gets chained with Sx
  m_move(Ft.Jx, nd, 0, n, nd, Fc.Jq, 0, 0);
  m_move(Ft.Ju, nd, 0, n, nu, Fc.Jq, 0, nx);

  Fc.c_setup = true;
  for (i = 0; i < n; i++) {
    int wrt = 0;
    if (Ft.is_linear_element(nd + i, Omu_Dependent::WRT_x))
      wrt |= Omu_Dependent::WRT_x;
    if (Ft.is_linear_element(nd + i, Omu_Dependent::WRT_dx))
      wrt |= Omu_Dependent::WRT_dx;
    if ((nd == 0 || Ft.is_linear_element(nd + i, Omu_Dependent::WRT_x)) &&
	Ft.is_linear_element(nd + i, Omu_Dependent::WRT_u))
      wrt |= Omu_Dependent::WRT_q;
    Fc.set_linear_element(i, wrt);
  }
  Fc.analyze_struct();
  Fc.c_setup = false;
}

//--------------------------------------------------------------------------
void Omu_Integrator::init_stage(int k,
				const Omu_VariableVec &x,
				const Omu_VariableVec &u,
				const Omu_DependentVec &Ft,
				bool sa)
{
  if (k >= _K) {
    m_error(E_INTERN, "Omu_Integrator::init_stage"
	    " that was called with wrong integrator setup");
  }

  _k = k;

  // initialize dimensions
  Omu_SVarVec &sx = (Omu_SVarVec &)x;
  _nxt = x->dim;
  _nd = sx.nd;
  _nv = sx.nv;
  _na = sx.na;
  _nx = _nxt - _nv;
  _nu = u->dim;
  _nq = _nx + _nu;
  _sa = sa;

  _n = _nxt - _nd;	// number of states for integration
  _m = _nd + _nu;	// number of control parameters for integration

  if ((int)(_Fcs[k]->dim) != _n) {
    m_error(E_INTERN, "Omu_Integrator::solve"
	    " that was called with wrong integrator setup of stage");
  }

  // initialize call arguments for sys->continuous
  _ut.resize(_nu);
  _dxt.resize(_nxt, _nx, _nu);
  v_zero(_dxt); // zero time derivative of discrete states

  // initialize call arguments for high-level integrator interface
  _xc.resize(_n, 0, 0, _nq);
  _dxc.resize(_n, 0, 0, _nq);
  _q.resize(_nq);

  // call high-level init
  init(k, _xc, _q, _Fcs[k], sa);

  // call depreciated init_stage
  init_stage(k, sx, u, sa);
}

//--------------------------------------------------------------------------
void Omu_Integrator::solve(int kk, double tstart, double tend,
			   const Omu_VariableVec &x, const Omu_VariableVec &u,
			   Omu_Program *sys, Omu_DependentVec &Ft,
			   Omu_StateVec &xt)
{
  _kk = kk;

  // store references to data required by residual() for sys->continuous calls
  _sys = sys;
  _x_ptr = &x;
  _u_ptr = &u;
  _Ft_ptr = &Ft;
  _xt_ptr = &xt;

  // prepare call arguments for high-level solve
  v_move(xt, _nd, _n, _xc, 0);
  if (_sa) {
    m_move(xt.Sx, _nd, 0, _n, _nx, _xc.Sq, 0, 0);
    m_move(xt.Su, _nd, 0, _n, _nu, _xc.Sq, 0, _nx);
  }
  v_move(xt, 0, _nx, _q, 0);
  v_move(u, 0, _nu, _q, _nx);

  // call high-level solve
  solve(kk, tstart, tend, _xc, _dxc, _q, _Fcs[_k]);

  // return results
  v_move(_xc, 0, _n, xt, _nd);
  if (_sa) {
    m_move(_xc.Sq, 0, 0, _n, _nx, xt.Sx, _nd, 0);
    m_move(_xc.Sq, 0, _nx, _n, _nu, xt.Su, _nd, 0);
  }

  // release references to data required by residual()
  _sys = NULL;
  _x_ptr = NULL;
  _u_ptr = NULL;
  _Ft_ptr = NULL;
  _xt_ptr = NULL;
}

//--------------------------------------------------------------------------
void Omu_Integrator::solve(int kk, double tstart, double tend,
			   Omu_StateVec &xc, Omu_StateVec &dxc, Omu_Vec &q,
			   Omu_DependentVec &Fc)
{
  // call depreciated solve and store results in _xc
  Omu_StateVec &xt = *_xt_ptr;
  Omu_SVarVec &sx = (Omu_SVarVec &)*_x_ptr;
  if (_sa) {
    solve(kk, tstart, tend, sx, *_u_ptr, _sys, xt, xt.Sx, xt.Su);
    m_move(xt.Sx, _nd, 0, _n, _nx, _xc.Sq, 0, 0);
    m_move(xt.Su, _nd, 0, _n, _nu, _xc.Sq, 0, _nx);
  }
  else {
    solve(kk, tstart, tend, sx, *_u_ptr, _sys, xt, MNULL, MNULL);
  }
  v_move(xt, _nd, _n, _xc, 0);
}

//--------------------------------------------------------------------------
void Omu_Integrator::residual(int kk, double t,
			      const Omu_StateVec &xc, const Omu_StateVec &dxc,
			      const Omu_Vec &q, Omu_DependentVec &Fc)
{
  if (!_sys) {
    m_error(E_NULL, "Omu_Integrator::residual");
  }

  // prepare call arguments
  Omu_DependentVec &Ft = *_Ft_ptr;
  Omu_StateVec &xt = *_xt_ptr;

  v_move(q, 0, _nd, xt, 0);
  v_move(xc, 0, _n, xt, _nd);
  v_move(dxc, 0, _n, _dxt, _nd);
  v_move(q, _nx, _nu, _ut, 0);

  Ft.set_required_J(Fc.is_required_J());
  // note: treat seed derivatives xc.Sq and dxc.Sq to support forward mode

  // call continuous of Omu_Program
  _sys->continuous(kk, t, xt, _ut, _dxt, Ft);

  // return results
  v_move(Ft, _nd, _n, Fc, 0);

  if (Fc.is_required_J()) {
    if (!Fc.Jx.is_constant())
      m_move(Ft.Jx, _nd, _nd, _n, _n, Fc.Jx, 0, 0);
    if (!Fc.Jdx.is_constant())
      m_move(Ft.Jdx, _nd, _nd, _n, _n, Fc.Jdx, 0, 0);
    if (!Fc.Jq.is_constant()) {
      m_move(Ft.Jx, _nd, 0, _n, _nd, Fc.Jq, 0, 0);
      m_move(Ft.Ju, _nd, 0, _n, _nu, Fc.Jq, 0, _nx);
    }
  }
}

//==========================================================================
