/*
 * Hqp_Omuses.C --
 *   -- class implementation
 *
 * rf, 7/27/96
 *
 */

/*
    Copyright (C) 1997--2007  Ruediger Franke

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

#include "Hqp_Omuses.h"

#include <assert.h>

#include <If_Bool.h>
#include <If_Int.h>
#include <If_Real.h>
#include <If_Method.h>
#include <If_Module.h>

#include "Omu_Vars.h"
#include "Omu_Deps.h"
#include "Omu_Program.h"
#include "Omu_Integrator.h"

#include "Omu_IntDopri5.h"

typedef If_Method<Hqp_Omuses> If_Cmd;

//--------------------------------------------------------------------------
Hqp_Omuses::Hqp_Omuses()
{
  _prg = NULL;
  _integrator = new Omu_IntDopri5;
  _xs = NULL;
  _us = NULL;
  _css = NULL;

  _x0s = NULL;
  _xfs = NULL;

  _xts = NULL;
  _Fs = NULL;
  _fs = NULL;
  _f0s = NULL;
  _cs = NULL;

  _IS = m_get(1,1);

  _ad = true;
  _fscale = 1.0;
  _hela_setup = 1;
  _stages_ok = false;

  _ifList.append(new IF_MODULE("prg_name", &_prg, Omu_Program));
  _ifList.append(new IF_MODULE("prg_integrator", &_integrator,
			       Omu_Integrator));
  _ifList.append(new If_Cmd("prg_setup_stages",
			    &Hqp_Omuses::setup_stages, this));
  _ifList.append(new If_Bool("prg_ad", &_ad));
  _ifList.append(new If_Real("prg_fscale", &_fscale));
  _ifList.append(new If_Int("prg_hela_setup", &_hela_setup));
}

//--------------------------------------------------------------------------
Hqp_Omuses::~Hqp_Omuses()
{
  m_free(_IS);

  delete [] _cs;
  delete [] _f0s;
  delete [] _fs;
  delete [] _Fs;
  delete [] _xts;

  delete [] _xfs;
  delete [] _x0s;

  delete [] _css;
  delete [] _us;
  delete [] _xs;
  delete _integrator;
  delete _prg;
}

//--------------------------------------------------------------------------
void Hqp_Omuses::setup_stages()
{
  assert(_prg != NULL);

  int K, KK;

  _prg->setup_stages();

  K = _prg->ks()->dim - 1;
  KK = _prg->ts()->dim - 1;

  delete [] _css;
  delete [] _us;
  delete [] _xs;
  _xs = new Omu_SVarVec [K+1];
  _us = new Omu_VarVec [K+1];
  _css = new Omu_VarVec [K+1];

  delete [] _cs;
  delete [] _f0s;
  delete [] _fs;
  delete [] _Fs;
  delete [] _xts;
  delete [] _xfs;
  delete [] _x0s;
  _x0s = new Omu_SVec [K+1];
  _xfs = new Omu_SVec [K+1];
  _xts = new Omu_DepVec [K+1];
  _Fs = new Omu_DepVec [K+1];
  _fs = new Omu_DepVec [K+1];
  _f0s = new Omu_Dep [K+1];
  _cs = new Omu_DepVec [K+1];

  _stages_ok = true;

  // setup integrator
  _integrator->setup_stages(_prg);
}

//--------------------------------------------------------------------------
void Hqp_Omuses::setup_horizon(int &k0, int &kf)
{
  assert(_prg != NULL);

  // call setup_stages if this wasn't already done through the interface
  if (!_stages_ok)
    setup_stages();

  _stages_ok = false;	// for subsequent problem modifications

  //
  // setup the optimization horizon and related data
  //

  k0 = 0;
  kf = _prg->ks()->dim - 1;
}

//--------------------------------------------------------------------------
void Hqp_Omuses::setup_vars(int k,
			    VECP x, VECP xmin, VECP xmax,
			    VECP u, VECP umin, VECP umax,
			    VECP c, VECP cmin, VECP cmax)
{
  if (!_prg) {
    m_error(E_NULL, "Hqp_Omuses::setup_vars");
  }

  int i, nx;
  int K = _prg->ks()->dim - 1;
  Omu_SVarVec &xk = _xs[k];
  Omu_VarVec   &uk = _us[k];
  Omu_VarVec   &ck = _css[k];

  xk.c_setup = true;
  if (k < K) {
    xk.c_expand = true;
    uk.c_setup = true;
  }
  ck.c_setup = true;

  _prg->setup(k, xk, uk, ck);

  xk.c_setup = false;
  xk.c_expand = false;
  uk.c_setup = false;
  ck.c_setup = false;

  nx = xk->dim - xk.nv;
  alloc_vars(x, xmin, xmax, nx);
  for (i = 0; i < nx; i++) {
    x[i] = xk.initial[i];
    xmin[i] = xk.min[i];
    xmax[i] = xk.max[i];
  }

  alloc_vars(u, umin, umax, uk->dim);
  v_copy(uk.initial, u);
  v_copy(uk.min, umin);
  v_copy(uk.max, umax);

  alloc_vars(c, cmin, cmax, ck->dim);
  v_copy(ck.initial, c);
  v_copy(ck.min, cmin);
  v_copy(ck.max, cmax);
}

//--------------------------------------------------------------------------
void Hqp_Omuses::setup_struct(int k, const VECP, const VECP,
			      MATP fx, MATP fu, IVECP f_lin,
			      VECP f0x, VECP f0u, int &f0_lin,
			      MATP cx, MATP cu, IVECP c_lin,
			      MATP Lxx, MATP Luu, MATP Lxu)
{
  if (!_prg) {
    m_error(E_NULL, "Hqp_Omuses::setup_struct");
  }

  // obtain references to variables and dependents of stage k
  int K = _prg->ks()->dim - 1;
  Omu_SVarVec &xk = _xs[k];
  Omu_VarVec &uk = _us[k];
  Omu_VarVec &csk = _css[k];
  int nx = xk->dim - xk.nv;
  int nu = uk->dim;
  int nxt = xk->dim;
  int nf = k < K? _xs[k+1]->dim - _xs[k+1].nv: 0;
  int nc = csk->dim;
  Omu_SVec &x0k = _x0s[k];
  Omu_SVec &xfk = _xfs[k];
  Omu_DepVec &xtk = _xts[k];
  Omu_DepVec &Fk = _Fs[k];
  Omu_DepVec &fk = _fs[k];
  Omu_Dep &f0k = _f0s[k];
  Omu_DepVec &ck = _cs[k];
  int i, j;
  bool is_zero;

  // allocate dependents
  xtk.size(nxt, nxt, nu, 0, 0, 0);
  if (k < K) {
    Fk.size(nxt, nxt, nu, nxt, 0, 0);
    fk.size(max(nxt, nf), nxt, nu, 0, nxt, 0);
    f0k.size(nxt, nu, nxt);
    ck.size(nc, nxt, nu, 0, nxt, 0);
  }
  else {
    f0k.size(nxt, nu, 0);
    ck.size(nc, nxt, nu, 0, 0, 0);
  }
  xtk.c_setup = true;
  Fk.c_setup = true;
  fk.c_setup = true;
  f0k.c_setup = true;
  ck.c_setup = true;

  // initialize variables with defaults
  v_copy(xk.initial, xk);
  v_copy(uk.initial, uk);
  if (k < K)
    v_copy(xk.initial, fk);

  _prg->setup_struct(k, xk, uk, xtk, Fk, fk, f0k, ck);

  xtk.c_setup = false;
  Fk.c_setup = false;
  fk.c_setup = false;
  f0k.c_setup = false;
  ck.c_setup = false;

  // analyze Jacobian structures
  // (do this for all deps even at k==K to mark xtk, Fk, fk constant and zero)
  xtk.analyze_struct();
  Fk.analyze_struct();
  fk.analyze_struct();
  f0k.analyze_struct();
  ck.analyze_struct();

  if (k < K) {
    // check for discrete states
    xk.nd = 0;
    for (i = 0; i < nxt; i++) {
      is_zero = true;
      for (j = 0; j < nxt && is_zero; j++) {
	is_zero &= (Fk.Jx[i][j] == 0.0);
      }
      for (j = 0; j < nu && is_zero; j++) {
	is_zero &= (Fk.Ju[i][j] == 0.0);
      }
      for (j = 0; j < nxt && is_zero; j++) {
	is_zero &= (Fk.Jdx[i][j] == 0.0);
      }
      if (is_zero) {
	xk.flags[i] |= Omu_SVarVec::Discrete;
	xk.nd++;
      }
      else
	break;
    }

    // check for algebraic states
    xk.na = 0;
    for (j = xk.nd; j < nxt; j++) {
      is_zero = true;
      for (i = xk.nd; i < nxt; i++) {
	is_zero &= (Fk.Jdx[i][j] == 0.0);
      }
      if (is_zero && !(xk.flags[j] & Omu_SVarVec::Discrete)) {
	// state j does not appear differentiated in any equation
	xk.flags[j] |= Omu_SVarVec::Algebraic;
	xk.na++;
      }
    }

    // check if Fdx is diagonal and constant (can use ODE solver if also na=0)
    xk.D_is_const = Fk.Jdx.is_constant() && Fk.Jdx.sbw() < 1;
    for (i = 0; i < nxt; i++)
      xk.D[i] = Fk.Jdx[i][i];

    // propagate semi-bandwidths to integrator
    xk.sbw_l = max(Fk.Jx.sbw_lower(), Fk.Jdx.sbw_lower());
    xk.sbw_u = max(Fk.Jx.sbw_upper(), Fk.Jdx.sbw_upper());
  }
  else {
    // no continuous-time equations in last stage
    for (i = 0; i < nxt; i++)
      xk.flags[i] |= Omu_SVarVec::Discrete;
    xk.nd = nxt;
    v_zero(xk.D);
  }

  // allocate state variables
  x0k.size(nxt, nx, nu, 0);
  if (xk.nd < nxt)
    // final states do only exist if there are continuous equations
    xfk.size(nxt, nx, nu, 0);

  // communicate Jacobian struct of a pure discrete problem to Hqp_Docp
  if (k == K && xtk.Jx.is_ident() && xtk.Ju.is_zero()) {
    //  || xk.nv == 0 && xk.nd == nx && _prg->ks(k+1) - _prg->ks(k) == 1) {
    // We could also communicate the struct for a
    // discrete-time stage with one sample period,
    // though this may conflict with LQDOCP's Check_Structure!!!
    m_copy(fk.Jx, fx);
    m_copy(fk.Ju, fu);
    for (i = 0; i < nf; i++)
      f_lin[i] = (int)fk.is_linear_element(i);

    m_copy(ck.Jx, cx);
    m_copy(ck.Ju, cu);
    for (i = 0; i < nc; i++)
      c_lin[i] = (int)ck.is_linear_element(i);
  }

  // communicate struct of Hessian to Hqp_Docp
  // for now: make Hessian diagonal for all variables that are
  //   - discrete states or control parameters and
  //   - only appear linear in any constraint
  //   assuming that the objective does not introduce off-diagonal elements
  if (_hela_setup) {
    for (i = 0; i < nx; i++) {
      if ((xk.flags[i] & Omu_SVarVec::Discrete) != 0
          && (k == K || fk.is_linear_variable(Omu_Dependent::WRT_x, i))
          && ck.is_linear_variable(Omu_Dependent::WRT_x, i)) {
        // zero out off-diagonal elements (upper diagonal part is sufficient)
        for (j = 0; j < i; j++)
          Lxx[j][i] = 0.0;
        for (j = i + 1; j < nx; j++)
          Lxx[i][j] = 0.0;
        for (j = 0; j < nu; j++)
          Lxu[i][j] = 0.0;
      }
    }
    for (i = 0; i < nu; i++) {
      // check that control does not influence continuous-time equations
      is_zero = true;
      for (j = 0; j < nxt && is_zero; j++) {
        is_zero &= (Fk.Ju[j][i] == 0.0);
      }
      if (is_zero
          && fk.is_linear_variable(Omu_Dependent::WRT_u, i)
          && ck.is_linear_variable(Omu_Dependent::WRT_u, i)) {
        // zero out off-diagonal elements (upper diagonal part is sufficient)
        for (j = 0; j < i; j++)
          Luu[j][i] = 0.0;
        for (j = i + 1; j < nu; j++)
          Luu[i][j] = 0.0;
        for (j = 0; j < nx; j++)
          Lxu[j][i] = 0.0;
      }
    }
  }

  // setup integrator for stage
  if (k < K) {
    if (k == 0 && _integrator->K() != _prg->K()) {
      // _integrator needs to be set up completely, e.g. as it was exchanged
      _integrator->setup_stages(_prg);
    }
    _integrator->setup_struct(k, xk, uk, Fk);
  }
}

//--------------------------------------------------------------------------
void Hqp_Omuses::resetup_integrator()
{
  int k, K = _prg->K();

  _integrator->setup_stages(_prg);

  for (k = 0; k < K; k++)
    _integrator->setup_struct(k, _xs[k], _us[k], _Fs[k]);
}

//--------------------------------------------------------------------------
void Hqp_Omuses::init_simulation(int k, VECP x, VECP u)
{
  if (!_prg) {
    m_error(E_NULL, "Hqp_Omuses::init_simulation");
  }

  int i;
  int nx = x->dim;
  int nu = u->dim;
  Omu_SVarVec &xk = _xs[k];
  Omu_VarVec &uk = _us[k];
  int nxt = xk->dim;

  for (i = 0; i < nx; i++)
    xk[i] = x[i];
  for (i = nx; i < nxt; i++)
    xk[i] = xk.initial[i];
  v_copy(u, uk);

  _prg->init_simulation(k, xk, uk);
  if (xk->dim != nxt || uk->dim != nu) {
    m_error(E_FORMAT, "Hqp_Omuses::init_simulation; modified vector sizes");
  }

  for (i = 0; i < nx; i++)
    x[i] = xk[i];
  v_copy(uk, u);
}

//--------------------------------------------------------------------------
void Hqp_Omuses::update_vals(int k, const VECP x, const VECP u,
			     VECP f, Real &f0, VECP c)
{
  double val;
  int i, kk;
  int K = _prg->ks()->dim - 1;
  int kkend = (k < K)? _prg->ks(k+1): _prg->ks(k) + 1;
  Omu_SVarVec &xk = _xs[k];
  Omu_VarVec &uk = _us[k];
  Omu_SVec &x0k = _x0s[k];
  Omu_SVec &xfk = _xfs[k];
  Omu_DepVec &xtk = _xts[k];
  Omu_DepVec &Fk = _Fs[k];
  Omu_DepVec &fk = _fs[k];
  Omu_Dep &f0k = _f0s[k];
  Omu_DepVec &ck = _cs[k];
  int nxt = xk->dim;
  int nxf = xfk->dim;
  int nx = x->dim;
  int nu = u->dim;
  int nf = f->dim;
  int nc = c->dim;

  assert(xk.nd >= 0);	// setup_struct must have been called before
  if (!_prg) {
    m_error(E_NULL, "Hqp_Omuses::update_vals");
  }

  xtk.set_required_J(false);
  fk.set_required_J(false);
  f0k.set_required_g(false);
  ck.set_required_J(false);

  f0 = 0.0;
  v_zero(c);

  //
  // integrate system equations over each sample period of the stage
  // and update criterion and constraints after each sample period
  //

  if (nxf > 0) {
    if (k == 0 && _integrator->K() != _prg->K()) {
      // _integrator needs to be set up, e.g. as it was exchanged
      resetup_integrator();
    }
    _integrator->init_stage(k, xk, uk, Fk, false);
  }

  // initial states of stage
  for (i = 0; i < nx; i++)
    x0k[i] = x[i];
  for (i = nx; i < nxt; i++)
    x0k[i] = xk.initial[i];

  // control parameters for stage
  v_copy(u, uk);

  if (nxt != nf && kkend - _prg->ks(k) > 1) {
    // there is more than one sample period;
    // fk has size nxt until last sample period
    fk.adapt_size(nxt);
  }

  for (kk = _prg->ks(k); kk < kkend; kk++) {

    // initialize (expansion) states
    _prg->consistic(kk, _prg->ts(kk), x0k, uk, xtk);
    v_copy(xtk, x0k);

    // solve continuous equations over sample period
    if (nxf > 0) {
      v_copy(x0k, xfk);
      _integrator->solve(kk, _prg->ts(kk), _prg->ts(kk+1),
			 xk, uk, _prg, Fk, xfk);
    }

    f0k = 0.0;
    if (nxt != nf && kk == kkend - 1) {
      // final sample period of stage with more than one sample period;
      // fk gets its actual size nf
      fk.adapt_size(nf);
    }
    v_zero(ck);

    // call update
    _prg->update(kk, x0k, uk, xfk, fk, f0k, ck);

    f0 += _fscale * f0k;

    v_add(c, ck, c);

    if (kk < kkend - 1) {
      v_copy(fk, x0k);
    }
    else if (nf > 0) {
      v_copy(fk, f);
    }
  }
}

//--------------------------------------------------------------------------
void Hqp_Omuses::update_stage(int k, const VECP x, const VECP u,
			      VECP f, Real &f0, VECP c,
			      MATP fx, MATP fu, VECP f0x, VECP f0u,
			      MATP cx, MATP cu,
			      const VECP rf, const VECP rc,
			      MATP Lxx, MATP Luu, MATP Lxu)
{
  if (!_ad) {
    Hqp_DocpStub::update_stage(k, x, u, f, f0, c,
			       fx, fu, f0x, f0u, cx, cu,
			       rf, rc, Lxx, Luu, Lxu);
    return;
  }

  int i, j, kk;
  int K = _prg->ks()->dim - 1;
  int kkend = (k < K)? _prg->ks(k+1): _prg->ks(k) + 1;
  Omu_SVarVec &xk = _xs[k];
  Omu_VarVec &uk = _us[k];
  Omu_SVec &x0k = _x0s[k];
  Omu_SVec &xfk = _xfs[k];
  Omu_DepVec &xtk = _xts[k];
  Omu_DepVec &Fk = _Fs[k];
  Omu_DepVec &fk = _fs[k];
  Omu_Dep &f0k = _f0s[k];
  Omu_DepVec &ck = _cs[k];
  int nxt = xk->dim;
  int nxf = xfk->dim;
  int nx = x->dim;
  int nu = u->dim;
  int nf = f->dim;
  int nc = c->dim;
  bool is_ident_Sx_and_zero_Su;

  assert(xk.nd >= 0);	// init_vals() must have been called before
  if (!_prg) {
    m_error(E_NULL, "Hqp_Omuses::update_stage");
  }

  xtk.set_required_J(true);
  fk.set_required_J(true);
  f0k.set_required_g(true);
  ck.set_required_J(true);

  f0 = 0.0;
  v_zero(f0x);
  v_zero(f0u);
  v_zero(c);
  m_zero(cx);
  m_zero(cu);

  //
  // integrate system equations over each sample period of the stage
  // and update criterion and constraints after each sample period
  //

  if (nxf > 0) {
    if (k == 0 && _integrator->K() != _prg->K()) {
      // _integrator needs to be set up, e.g. as it was exchanged
      resetup_integrator();
    }
    _integrator->init_stage(k, xk, uk, Fk, true);

    // initial states of stage
    m_zero(x0k.Sx);
    m_zero(x0k.Su);
    for (i = 0; i < nx; i++) {
      x0k[i] = x[i];
      x0k.Sx[i][i] = 1.0;
    }
    for (i = nx; i < nxt; i++)
      x0k[i] = xk.initial[i];

    if (nxt == nx)
      is_ident_Sx_and_zero_Su = true;
    else
      is_ident_Sx_and_zero_Su = false;
  }
  else {
    // initialize stage without continuous equations
    v_copy(x, x0k);
    is_ident_Sx_and_zero_Su = true;
  }

  // control parameters for stage
  v_copy(u, uk);

  if (nxt != nf && kkend - _prg->ks(k) > 1) {
    // there is more than one sample period;
    // fk has size nxt until last sample period
    fk.adapt_size(nxt);
  }

  for (kk = _prg->ks(k); kk < kkend; kk++) {

    // initialize (expansion) states
    _prg->consistic(kk, _prg->ts(kk), x0k, uk, xtk);
    v_copy(xtk, x0k);
    if (!xtk.Jx.is_ident()) {
      m_copy(m_mlt(xtk.Jx, x0k.Sx, _IS), x0k.Sx);
      m_copy(m_mlt(xtk.Jx, x0k.Su, _IS), x0k.Su);
    }
    if (!xtk.Ju.is_zero())
      m_add(x0k.Su, xtk.Ju, x0k.Su);

    // solve continuous equations over sample period
    if (nxf > 0) {
      v_copy(x0k, xfk);
      m_copy(x0k.Sx, xfk.Sx);
      m_copy(x0k.Su, xfk.Su);
      _integrator->solve(kk, _prg->ts(kk), _prg->ts(kk+1),
                         xk, uk, _prg, Fk, xfk);
    }

    //
    // update criterion, constraints, and discrete states
    //

    // initialize dependent variables
    f0k = 0.0;
    if (nxt != nf && kk == kkend - 1) {
      // final sample period of stage with more than one sample period;
      // fk gets its actual size nf
      fk.adapt_size(nf);
    }
    v_zero(ck);

    // call update
    _prg->update(kk, x0k, uk, xfk, fk, f0k, ck);

    f0 += _fscale * f0k;

    if (is_ident_Sx_and_zero_Su) {
      for (i = 0; i < nx; i++)
	f0x[i] += _fscale * f0k.gx[i];
      for (i = 0; i < nu; i++)
	f0u[i] += _fscale * f0k.gu[i];
    }
    else {
      vm_mltadd(f0x, f0k.gx, x0k.Sx, _fscale, f0x);
      vm_mltadd(f0u, f0k.gx, x0k.Su, _fscale, f0u);
      for (i = 0; i < nu; i++)
	f0u[i] += _fscale * f0k.gu[i];
    }
    if (nxf > 0) {
      // add contribution of continuous equations
      vm_mltadd(f0x, f0k.gxf, xfk.Sx, _fscale, f0x);
      vm_mltadd(f0u, f0k.gxf, xfk.Su, _fscale, f0u);
    }

    v_add(c, ck, c);

    if (is_ident_Sx_and_zero_Su) {
      m_add(cx, ck.Jx, cx);
      m_add(cu, ck.Ju, cu);
    }
    else {
      m_mltadd(cx, ck.Jx, x0k.Sx, cx);
      m_mltadd(cu, ck.Jx, x0k.Su, cu);
      m_add(cu, ck.Ju, cu);
    }
    if (nxf > 0) {
      // add contribution of continuous equations
      m_mltadd(cx, ck.Jxf, xfk.Sx, cx);
      m_mltadd(cu, ck.Jxf, xfk.Su, cu);
    }

    if (kk < kkend - 1) {
      // intermediate sample period
      v_copy(fk, x0k);

      if (is_ident_Sx_and_zero_Su) {
	m_copy(fk.Jx, x0k.Sx);
	m_copy(fk.Ju, x0k.Su);
	is_ident_Sx_and_zero_Su = false;
      }
      else {
	m_copy(m_mlt(fk.Jx, x0k.Sx, _IS), x0k.Sx);
	m_copy(m_mlt(fk.Jx, x0k.Su, _IS), x0k.Su);
	m_add(x0k.Su, fk.Ju, x0k.Su);
      }
      if (nxf > 0) {
	// add contribution of continuous equations
	m_mltadd(x0k.Sx, fk.Jxf, xfk.Sx, x0k.Sx);
	m_mltadd(x0k.Su, fk.Jxf, xfk.Su, x0k.Su);
      }
    }
    else if (nf > 0) {
      // final sample period of stage
      v_copy(fk, f);

      if (is_ident_Sx_and_zero_Su) {
	m_copy(fk.Jx, fx);
	m_copy(fk.Ju, fu);
      }
      else {
	m_mlt(fk.Jx, x0k.Sx, fx);
	m_copy(fk.Ju, fu);
	m_mltadd(fu, fk.Jx, x0k.Su, fu);
      }
      if (nxf > 0) {
	// add contribution of continuous equations
	m_mltadd(fx, fk.Jxf, xfk.Sx, fx);
	m_mltadd(fu, fk.Jxf, xfk.Su, fu);
      }
    }
  }
}


//========================================================================
