/*
 * Hqp_Omuses.C --
 *   -- class implementation
 *
 * rf, 7/27/96
 *
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

#include "Hqp_Omuses.h"

#include <assert.h>

#ifdef OMU_WITH_ADOLC
#include <adutils.h>
#endif

#include <If_Bool.h>
#include <If_Real.h>
#include <If_Method.h>
#include <If_Module.h>

#include "Omu_Vars.h"
#include "Omu_Deps.h"
#include "Omu_Program.h"
#include "Omu_Integrator.h"

#include "Omu_IntDopri5.h"

typedef If_Method<Hqp_Omuses> If_Cmd;

#ifdef OMU_WITH_ADOLC
//--------------------------------------------------------------------------
static short** myalloc_short(int m, int n)
{
  short* Adum = (short*)malloc(m*n*sizeof(short));
  short**  A = (short**)malloc(m*sizeof(short*));
  int i;
  for(i=0;i<m;i++)
    {
      A[i] = Adum;
      Adum += n;
    }
  return A;

  /* To deallocate an array set up by   A = myalloc2(m,n)   */
  /*   use  free((char*)*A); free((char*)A); in that order  */
}
#endif

//--------------------------------------------------------------------------
Hqp_Omuses::Hqp_Omuses()
{
  _prg = NULL;
  _integrator = new Omu_IntDopri5;
  _integrator_setup = NULL;
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

  _xt = v_get(1);
  _xtxk = m_get(1,1);
  _xtuk = m_get(1,1);
  _Sx = m_get(1,1);
  _Sxk = m_get(1,1);
  _IS = m_get(1,1);
  _Su = m_get(1,1);
  _Suk = m_get(1,1);

  _fk = v_get(1);

  _fkxk = m_get(1, 1);
  _fkuk = m_get(1, 1);
  _fkfk = m_get(1, 1);
  _f0kxk = v_get(1);
  _f0kuk = v_get(1);
  _f0kfk = v_get(1);
  _ckxk = m_get(1, 1);
  _ckuk = m_get(1, 1);
  _ckfk = m_get(1, 1);

  _ad = true;
  _fscale = 1.0;
  _stages_ok = false;

  ad_alloc(1, 1, 1);

  _ifList.append(new IF_MODULE("prg_name", &_prg, Omu_Program));
  _ifList.append(new IF_MODULE("prg_integrator", &_integrator,
			       Omu_Integrator));
  _ifList.append(new If_Cmd("prg_setup_stages",
			    &Hqp_Omuses::setup_stages, this));
  _ifList.append(new If_Bool("prg_ad", &_ad));
  _ifList.append(new If_Real("prg_fscale", &_fscale));
}

//--------------------------------------------------------------------------
Hqp_Omuses::~Hqp_Omuses()
{
  ad_free();

  m_free(_ckfk);
  m_free(_ckuk);
  m_free(_ckxk);
  v_free(_f0kfk);
  v_free(_f0kuk);
  v_free(_f0kxk);
  m_free(_fkfk);
  m_free(_fkuk);
  m_free(_fkxk);

  v_free(_fk);
  m_free(_Suk);
  m_free(_Su);
  m_free(_IS);
  m_free(_Sxk);
  m_free(_Sx);
  m_free(_xtxk);
  m_free(_xtuk);
  v_free(_xt);

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
void Hqp_Omuses::ad_alloc(int m, int n, int p)
{
#ifdef OMU_WITH_ADOLC
  _Z3 = myalloc(m, n, 1);
  _nz = myalloc_short(m, n);

  _x = v_get(n);
  _y = v_get(m);
  _X = m_get(n, p);
  _Y = m_get(m, p);
  _U = m_get(m, m);
  m_ident(_U);

  _max_ndep = m;
  _max_nindep = n;
  _max_npar = p;
#else
  _Z3 = NULL;
  _nz = NULL;

  _x = NULL;
  _y = NULL;
  _X = NULL;
  _Y = NULL;
  _U = NULL;

  _max_ndep = 0;
  _max_nindep = 0;
  _max_npar = 0;
#endif
}

//--------------------------------------------------------------------------
void Hqp_Omuses::ad_free()
{
#ifdef OMU_WITH_ADOLC
  free(**_Z3);
  free(*_Z3);
  free(_Z3);
  free(*_nz);
  free(_nz);
  m_free(_U);
  m_free(_Y);
  m_free(_X);
  v_free(_y);
  v_free(_x);
#endif
}

//--------------------------------------------------------------------------
void Hqp_Omuses::ad_realloc(int ndep, int nindep, int npar)
{
#ifdef OMU_WITH_ADOLC
  if (ndep > _max_ndep
      || nindep > _max_nindep
      || npar > _max_npar) {
    ad_free();
    ad_alloc(ndep, nindep, npar);
  }
#endif
}

//--------------------------------------------------------------------------
int Hqp_Omuses::setup_stages(IF_CMD_ARGS)
{
  assert(_prg != NULL);

  int K, KK;

  _prg->setup_stages();

  K = _prg->ks()->dim - 1;
  KK = _prg->ts()->dim - 1;

  delete [] _css;
  delete [] _us;
  delete [] _xs;
  _xs = new Omu_DynVarVec [K+1];
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
  _integrator_setup = _integrator;

  return IF_OK;
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
  Omu_DynVarVec &xk = _xs[k];
  Omu_VarVec   &uk = _us[k];
  Omu_VarVec   &ck = _css[k];

  xk.c_alloc = true;
  if (k < K) {
    xk.c_expand = true;
    uk.c_alloc = true;
  }
  ck.c_alloc = true;

  _prg->setup(k, xk, uk, ck);

  xk.c_alloc = false;
  xk.c_expand = false;
  uk.c_alloc = false;
  ck.c_alloc = false;

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

  init_vars(k, x, u);
}

//--------------------------------------------------------------------------
void Hqp_Omuses::init_vars(int k, VECP x, VECP u)
{
  if (!_prg) {
    m_error(E_NULL, "Hqp_Omuses::init_vars");
  }

  int i;
  Omu_DynVarVec &xk = _xs[k];
  Omu_VarVec &uk = _us[k];
  int nx = x->dim;

  v_copy(_xs[k].initial, xk);
  nx = x->dim;
  for (i = 0; i < nx; i++)
    x[i] = xk[i];

  v_copy(_us[k].initial, uk);
  v_copy(uk, u);
}

#if 1 // new implementation calling low-level setup_struct
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
  Omu_DynVarVec &xk = _xs[k];
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
  if (k < K) {
    xtk.alloc(nxt, nxt, nu, 0, 0);
    Fk.alloc(nxt, nxt, nu, nxt, 0);
    fk.alloc(max(nxt, nf), nxt, nu, 0, nxt);
    f0k.alloc(nxt, nu, nxt);
    ck.alloc(nc, nxt, nu, 0, nxt);
  }
  else {
    x0k.alloc(nxt, nx, nu);
    f0k.alloc(nxt, nu, 0);
    ck.alloc(nc, nxt, nu, 0, 0);
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
	is_zero &= (Fk.Jxp[i][j] == 0.0);
      }
      if (is_zero) {
	xk.flags[i] |= Omu_DynVarVec::Discrete;
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
	is_zero &= (Fk.Jxp[i][j] == 0.0);
      }
      if (is_zero && !(xk.flags[j] & Omu_DynVarVec::Discrete)) {
	// state j does not appear differentiated in any equation
	xk.flags[j] |= Omu_DynVarVec::Algebraic;
	xk.na++;
      }
    }

    // check if Fxp is diagonal and constant (can use ODE solver if also na=0)
    xk.D_is_const = Fk.Jxp.is_constant() && Fk.Jxp.sbw() < 1;
    for (i = 0; i < nxt; i++)
      xk.D[i] = Fk.Jxp[i][i];

    // propagate semi-bandwidths to integrator
    xk.sbw_l = max(Fk.Jx.sbw_lower(), Fk.Jxp.sbw_lower());
    xk.sbw_u = max(Fk.Jx.sbw_upper(), Fk.Jxp.sbw_upper());
  }
  else {
    // no continuous-time equations in last stage
    xk.nd = nxt;
    v_zero(xk.D);
  }

  // allocate state variables
  x0k.alloc(nxt, nx, nu);
  if (xk.nd < nxt)
    // final states do only exist if there are continuous equations
    xfk.alloc(nxt, nx, nu);
  
  // communicate Jacobian struct of a pure discrete problem to Hqp_Docp
  if (k == K) {
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
}
#else // old implementation
//--------------------------------------------------------------------------
void Hqp_Omuses::setup_struct(int k, const VECP, const VECP,
			      MATP fx, MATP fu, IVECP f_lin,
			      VECP f0x, VECP f0u, int &f0_lin,
			      MATP cx, MATP cu, IVECP c_lin,
			      MATP Lxx, MATP Luu, MATP Lxu)
{
  // obtain discrete-time states and dF/dxt
  Omu_DynVarVec &xk = _xs[k];
  Omu_VarVec &uk = _us[k];
  obtain_structure(k, xk, uk);

  // for now don't analyze sparse structure as required by HQP
  // i.e. assume dense block matrices
}
#endif

//--------------------------------------------------------------------------
void Hqp_Omuses::obtain_structure(int k,
				  Omu_DynVarVec &xk, const Omu_VarVec &uk)
{
#ifdef OMU_WITH_ADOLC
  if (!_prg) {
    m_error(E_NULL, "Hqp_Omuses::obtain_structure");
  }

  //
  // determine the number of discrete-time states per sample period
  //

  int i, j, nd, kk, sbw;
  int nxt = xk->dim;
  int nx = nxt - xk.nv;
  int nu = uk->dim;
  int nxtu = nxt + nu;
  int K = _prg->ks()->dim - 1;

  if (k == K) {
    xk.nd = nxt;
    v_zero(xk.D);
    return;
  }

  adoublev ax(nxt);
  adoublev au(nu);
  adoublev axp(nxt);
  adoublev aF(nxt);
  int ndep, nindep;
  int depends, sumdepends;
  double dummy;

  ndep = nxt;
  nindep = 2*nxt + nu;
  ad_realloc(ndep, nindep, nx+nu);

  kk = _prg->ks(k);
  v_resize(_xt, nxt);
  _prg->consistic(kk, _prg->ts(kk), xk, uk, _xt);

  // record a residuum evaluation

  trace_on(2, 1);	// tape 2, keep = 1 for following reverse call
  ax <<= _xt->ve;
  au <<= uk->ve;
  for (i = 0; i < nxt; i++)
    axp[i] <<= 0.0;

  _prg->continuous(kk, _prg->ts(kk), ax, au, axp, aF);

  for (i = 0; i < nxt; i++)
    aF[i] >>= dummy;

  trace_off();

  // obtain nonzero pattern

  reverse(2, ndep, nindep, 0, _Z3, _nz);

  nd = 0;
  sumdepends = 0;
  xk.D_is_const = true;
  for (i = 0; i < nxt; i++) {
    depends = 0;
    for (j = 0; j < nxtu; j++) {
      sumdepends += _nz[i][j];
    }
    for (j = nxtu; j < nindep; j++) {
      sumdepends += _nz[i][j];
      depends += _nz[i][j];
    }
    if (sumdepends == 0) {
      xk.flags[i] |= Omu_DynVarVec::Discrete;
      xk.D[nd] = 0.0;
      nd++;
    }
    if (depends > 1 || _nz[i][nxtu + i] - depends != 0) {
      xk.D_is_const = false;
    }
    xk.D[i] = _Z3[i][nxtu + i][0];

  }
  xk.nd = nd;

  // obtain structure of continuous equations
  xk.na = 0;
  for (i = nd; i < nxt; i++) {
    // check for algebraic states in column i of dF/dxp
    depends = 0;
    for (j = nd; j < nxt; j++) {
      depends += _nz[j][nxtu+i];
    }
    if (depends == 0) {
      // state i does not appear differentiated in any equation
      xk.flags[i] |= Omu_DynVarVec::Algebraic;
      xk.na++;
    }

    // determine semi-bandwidths of dF/dx + dF/dxp
    sbw = i - nd;	// maximum lower sbw in this row
    for (j = nd; j < i; j++) {
      if (_nz[i][j] + _nz[i][nxtu+j] != 0)
	break;
      --sbw;
    }
    xk.sbw_l = max(sbw, xk.sbw_l);

    sbw = nxt - i - 1;	// maximum upper sbw in this row
    for (j = nxt - 1; j > i; j--) {
      if (_nz[i][j] + _nz[i][nxtu+j] != 0)
	break;
      --sbw;
    }
    xk.sbw_u = max(sbw, xk.sbw_u);
  }
#else
  m_error(E_NULL, "Hqp_Omuses::obtain_structure: was compiled without ADOL-C");
#endif
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
  Omu_DynVarVec &xk = _xs[k];
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

#if 1 // new implementation calling low-level update
//--------------------------------------------------------------------------
void Hqp_Omuses::update_vals(int k, const VECP x, const VECP u,
			     VECP f, Real &f0, VECP c)
{
  double val;
  int i, kk;
  int K = _prg->ks()->dim - 1;
  int kkend = (k < K)? _prg->ks(k+1): _prg->ks(k) + 1;
  Omu_DynVarVec &xk = _xs[k];
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
    if (_integrator != _integrator_setup) {
      // _integrator needs to be set up, e.g. as it was exchanged
      _integrator->setup_stages(_prg);
      _integrator_setup = _integrator;
    }
    _integrator->init_stage(k, xk, uk, Fk);
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

    // solve continuous equations over sample period
    if (nxf > 0) {

      _prg->consistic(kk, _prg->ts(kk), x0k, uk, xtk);
      v_copy(xtk, xfk);

      _integrator->init_sample(kk, _prg->ts(kk), _prg->ts(kk+1));
      m_catchall(// try
		 _integrator->solve(kk, _prg->ts(kk), _prg->ts(kk+1),
				    xk, uk, _prg, Fk, xfk),
		 // catch
		 f0 = Inf;
		 return);
    }

    f0k = 0.0;
    v_zero(ck);

    // call update
    _prg->update(kk, x0k, uk, xfk, fk, f0k, ck);

    f0 += _fscale * f0k;

    v_add(c, ck, c);

    if (kk < kkend - 1) {
      v_copy(fk, x0k);
    }
    else if (nf > 0) {
      fk.adapt_size(nf);
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
  Omu_DynVarVec &xk = _xs[k];
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
    if (_integrator != _integrator_setup) {
      // _integrator needs to be set up, e.g. as it was exchanged
      _integrator->setup_stages(_prg);
      _integrator_setup = _integrator;
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

    // solve continuous equations over sample period
    if (nxf > 0) {

      _prg->consistic(kk, _prg->ts(kk), x0k, uk, xtk);
      v_copy(xtk, xfk);
      m_mlt(xtk.Jx, x0k.Sx, xfk.Sx);
      m_mlt(xtk.Jx, x0k.Su, xfk.Su);
      m_add(xfk.Su, xtk.Ju, xfk.Su);

      _integrator->init_sample(kk, _prg->ts(kk), _prg->ts(kk+1));
      m_catchall(// try
		 _integrator->solve(kk, _prg->ts(kk), _prg->ts(kk+1),
				    xk, uk, _prg, Fk, xfk),
		 // catch
		 f0 = Inf;
		 return);
    }

    //
    // update criterion, constraints, and discrete states
    //

    // initialize dependent variables
    f0k = 0.0;
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
      fk.adapt_size(nf);
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

#else // old implementation calling high-level update

//--------------------------------------------------------------------------
void Hqp_Omuses::update_vals(int k, const VECP x, const VECP u,
			     VECP f, Real &f0, VECP c)
{
  double val;
  int i, kk, kkend;
  int K = _prg->ks()->dim - 1;
  Omu_DynVarVec &xk = _xs[k];
  Omu_VarVec &uk = _us[k];
  int nxt = xk->dim;
  int nx = x->dim;
  int nu = u->dim;
  int nf = f->dim;
  int nxtf = max(nxt, nf);
  int nc = c->dim;
  int nd = xk.nd;

  assert(nd >= 0);	// obtain_structure() must have been called before
  if (!_prg) {
    m_error(E_NULL, "Hqp_Omuses::update_vals");
  }

  adouble af0;
  adoublev ax(nxt);
  adoublev au(nu);
  adoublev af(nxtf);	// size may increase or decrease in the next stage
  adoublev ac(nc);

  v_resize(_xt, nxt);

  f0 = 0.0;
  v_zero(c);

  //
  // integrate the system equations over the stage
  // calculate the discrete part for each sample of the stage
  //

  if (k < K) {
    _integrator->init_stage(k, xk, uk);
  }

  // initial values for the integration
  for (i = 0; i < nx; i++)
    _xt[i] = x[i];
  for (i = nx; i < nxt; i++)
    _xt[i] = xk.initial[i];
  v_copy(u, uk);

  kkend = (k < K)? _prg->ks(k+1): _prg->ks(k) + 1;
  for (kk = _prg->ks(k); kk < kkend; kk++) {

    // backing store initial states
    v_copy(_xt, xk);

    if (nxt > nd) {
      _prg->consistic(kk, _prg->ts(kk), xk, uk, _xt);
      _integrator->init_sample(kk, _prg->ts(kk), _prg->ts(kk+1));
      m_catchall(// try
		 _integrator->solve(kk, _prg->ts(kk), _prg->ts(kk+1),
				    xk, uk, _prg, _xt),
		 // catch
		 f0 = Inf;
		 return);
    }

    af0 <<= 0.0;
    for (i = 0; i < nc; i++)
      ac[i] <<= 0.0;
    for (i = 0; i < nxt; i++) 
      af[i] = _xt[i];
    if (kk == kkend - 1) {
      // default values for possible additional fs
      for (i = nxt; i < nf; i++) 
	af[i] = 0.0;
    }
    ax <<= xk->ve;
    au <<= uk->ve;

    _prg->update(kk, ax, au, af, af0, ac);
    af0 *= _fscale;

    af0 >>= val;
    f0 += val;
    for (i = 0; i < nc; i++) {
      ac[i] >>= val;
      c[i] += val;
    }
    if (kk < kkend - 1) {
      for (i = 0; i < nxt; i++)
	af[i] >>= _xt[i];
    }
    else {
      for (i = 0; i < nf; i++)
	af[i] >>= f[i];
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
    Hqp_Docp_stub::update_stage(k, x, u, f, f0, c,
				fx, fu, f0x, f0u, cx, cu,
				rf, rc, Lxx, Luu, Lxu);
    return;
  }

  double val;
  int i, j;
  int K = _prg->ks()->dim - 1;
  Omu_DynVarVec &xk = _xs[k];
  Omu_VarVec &uk = _us[k];
  int nxt = xk->dim;
  int nx = x->dim;
  int nu = u->dim;
  int nf = f->dim;
  int nxtf = max(nxt, nf);
  int nc = c->dim;
  int nd = xk.nd;
  int kk, kkend;

  assert(nd >= 0);	// init_vals() must have been called before
  if (!_prg) {
    m_error(E_NULL, "Hqp_Omuses::update_stage");
  }

  adouble af0;
  adoublev ax(nxt);
  adoublev au(nu);
  adoublev af(nxtf);
  adoublev ac(nc);

  int nxtu = nxt + nu;
  int nindep = nxt + nu + nxtf;
  int ndep = 1 + nc + nxtf;

  ad_realloc(ndep, nindep, nx+nu);

  f0 = 0.0;
  v_zero(f0x);
  v_zero(f0u);
  v_zero(c);
  m_zero(cx);
  m_zero(cu);

  v_resize(_xt, nxt);
  m_resize(_xtxk, nxt, nxt);
  m_resize(_xtuk, nxt, nu);
  m_resize(_Sx, nxt, nx);
  m_resize(_Su, nxt, nu);

  //
  // integrate the system equations over the stage
  // calculate the criterion for each sample within the stage
  //

  if (k < K)
    _integrator->init_stage(k, xk, uk, true);

  // initial values for the integration
  m_zero(_Sx);
  m_zero(_Su);
  for (i = 0; i < nx; i++) {
    _xt[i] = x[i];
    _Sx[i][i] = 1.0;
  }
  for (i = nx; i < nxt; i++)
    _xt[i] = xk.initial[i];
  v_copy(u, uk);

  kkend = (k < K)? _prg->ks(k+1): _prg->ks(k) + 1;
  for (kk = _prg->ks(k); kk < kkend; kk++) {

    // backing store initial states of the sample period
    // initialize seed derivatives for the initial states

    v_copy(_xt, xk);
    m_zero(_X);

    if (nxt > nd) {
      // solve continuous-time equations
      m_zero(_xtxk);
      m_zero(_xtuk);
      _prg->consistic(kk, _prg->ts(kk), xk, uk, _xt, _xtxk, _xtuk);
      m_mlt(_xtxk, _Sx, _IS);
      m_copy(_IS, _Sx);
      m_mlt(_xtxk, _Su, _IS);
      m_add(_IS, _xtuk, _Su);
      for (i = 0; i < nxt; i++) {
	for (j = 0; j < nx; j++) {
	  _X[i][j] = _Sx[i][j];
	}
	for (j = 0; j < nu; j++) {
	  _X[i][nx + j] = _Su[i][j];
	}
      }
      _integrator->init_sample(kk, _prg->ts(kk), _prg->ts(kk+1));
      m_catchall(// try
		 _integrator->solve(kk, _prg->ts(kk), _prg->ts(kk+1),
				    xk, uk, _prg, _xt, _Sx, _Su),
		 // catch
		 f0 = Inf;
		 return);
    }
    else if (kk > _prg->ks(k)) {
      // second or further sample of stage without continuous-time equations
      for (i = 0; i < nxt; i++) {
	for (j = 0; j < nx; j++) {
	  _X[i][j] = _Sx[i][j];
	}
	for (j = 0; j < nu; j++) {
	  _X[i][nx + j] = _Su[i][j];
	}
      }
    }
    else {
      // first sample of stage without continuous-time equations
      for (i = 0; i < nx; i++) {
	_X[i][i] = 1.0;
      }
    }

    // init seed derivatives for control parameters
    for (i = 0; i < nu; i++) {
      _X[nxt + i][nx + i] = 1.0;
    }

    //
    // update criterion, constraints, and discrete states
    //

    // initialize dependent variables
    af0 = 0.0;
    ndep = 1;
    for (i = 0; i < nc; i++)
      ac[i] = 0.0;
    ndep += nc;
    if (kk < kkend - 1) {
      // intermediate sample
      for (i = 0; i < nxt; i++)
	af[i] = _xt[i];
      ndep += nxt;
    }
    else {
      // final sample
      for (i = 0; i < nxt; i++)
	af[i] = _xt[i];
      for (i = nxt; i < nf; i++)
	af[i] = 0.0;
      ndep += nf;
    }

    trace_on(2);

    // initialize variable vector
    nindep = 0;
    ax <<= xk->ve;
    for (i = 0; i < nxt; i++) {
      _x[i] = xk[i];
    }
    nindep += nxt;
    au <<= uk->ve;
    for (i = 0; i < nu; i++) {
      _x[nxt + i] = uk[i];
    }
    nindep += nu;
    for (i = nd; i < nxt; i++) {
      af[i] <<= _xt[i];
      _x[nxtu - nd + i] = _xt[i];
    }
    nindep += nxt - nd;
    // discrete-time states are constant, that is why
    // no extra independent variables necessary
    for (i = 0; i < nd; i++)
      af[i] = ax[i];

    _prg->update(kk, ax, au, af, af0, ac);
    af0 *= _fscale;

    af0 >>= val;
    f0 += val;
    for (i = 0; i < nc; i++) {
      ac[i] >>= val;
      c[i] += val;
    }
    if (kk < kkend - 1) {
      for (i = 0; i < nxt; i++)
	af[i] >>= _xt[i];
    }
    else {
      for (i = 0; i < nf; i++)
	af[i] >>= f[i];
    }

    trace_off();

    // initialize seed derivatives for final continuous-time states
    for (i = nd; i < nxt; i++) {
      for (j = 0; j < nx; j++) {
	_X[nxtu - nd + i][j] = _Sx[i][j];
      }
      for (j = 0; j < nu; j++) {
	_X[nxtu - nd + i][nx + j] = _Su[i][j];
      }
    }
      
    forward(2, ndep, nindep, nx+nu, _x->ve, _X->me, _y->ve, _Y->me);

    for (j = 0; j < nx; j++)
      f0x[j] += _Y[0][j];
    for (j = 0; j < nu; j++)
      f0u[j] += _Y[0][nx + j];

    for (i = 0; i < nc; i++) {
      for (j = 0; j < nx; j++)
	cx[i][j] += _Y[1 + i][j];
      for (j = 0; j < nu; j++)
	cu[i][j] += _Y[1 + i][nx + j];
    }

    if (kk < kkend - 1) {
      for (i = 0; i < nxt; i++) {
	for (j = 0; j < nx; j++)
	  _Sx[i][j] = _Y[1 + nc + i][j];
	for (j = 0; j < nu; j++)
	  _Su[i][j] = _Y[1 + nc + i][nx + j];
      }
    }
    else {
      for (i = 0; i < nf; i++) {
	for (j = 0; j < nx; j++)
	  fx[i][j] = _Y[1 + nc + i][j];
	for (j = 0; j < nu; j++)
	  fu[i][j] = _Y[1 + nc + i][nx + j];
      }
    }
  }
}
#endif


//========================================================================
