/*
 * Hqp_Omuses.C --
 *   -- class implementation
 *
 * rf, 7/27/96
 *
 * rf, 1/25/97
 *   -- discrete states are by default constant
 *   -- continuous states may be remapped in discrete()
 *
 * rf, 1/29/98
 *   -- bug fix in obtain_structure (as detected by H. Linke)
 */

/*
    Copyright (C) 1997--2001  Ruediger Franke

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
#include <adutils.h>
#include <If_Int.h>
#include <If_Bool.h>
#include <If_Float.h>
#include <If_FloatVec.h>
#include <If_IntVec.h>
#include <If_Method.h>
#include <If_Module.h>

#include "Omu_Vars.h"
#include "Omu_Program.h"
#include "Omu_Integrator.h"

#include "Omu_IntDopri5.h"

typedef If_Method<Hqp_Omuses> If_Cmd;

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

//--------------------------------------------------------------------------
Hqp_Omuses::Hqp_Omuses()
{
  _prg = NULL;
  _integrator = new Omu_IntDopri5;
  _xs = NULL;
  _us = NULL;
  _cs = NULL;

  _xt = v_get(1);
  _xtxk = m_get(1,1);
  _xtuk = m_get(1,1);
  _Sx = m_get(1,1);
  _IS = m_get(1,1);
  _Su = m_get(1,1);
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
  _ifList.append(new If_Float("prg_fscale", &_fscale));
}

//--------------------------------------------------------------------------
Hqp_Omuses::~Hqp_Omuses()
{
  ad_free();

  m_free(_Su);
  m_free(_IS);
  m_free(_Sx);
  m_free(_xtxk);
  m_free(_xtuk);
  v_free(_xt);

  delete [] _cs;
  delete [] _us;
  delete [] _xs;
  delete _integrator;
  delete _prg;
}

//--------------------------------------------------------------------------
void Hqp_Omuses::ad_alloc(int m, int n, int p)
{
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
}

//--------------------------------------------------------------------------
void Hqp_Omuses::ad_free()
{
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
}

//--------------------------------------------------------------------------
void Hqp_Omuses::ad_realloc(int ndep, int nindep, int npar)
{
  if (ndep > _max_ndep
      || nindep > _max_nindep
      || npar > _max_npar) {
    ad_free();
    ad_alloc(ndep, nindep, npar);
  }
}

//--------------------------------------------------------------------------
int Hqp_Omuses::setup_stages(IF_CMD_ARGS)
{
  assert(_prg != NULL);

  int K;

  _prg->setup_stages();

  K = _prg->ks()->dim - 1;

  delete [] _cs;
  delete [] _us;
  delete [] _xs;
  _xs = new Omu_States [K+1];
  _us = new Omu_Vars [K+1];
  _cs = new Omu_Vars [K+1];

  _stages_ok = true;

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
void Hqp_Omuses::obtain_structure(int k,
				  Omu_States &xk, const Omu_Vector &uk)
{
  if (!_prg) {
    error(E_NULL, "Hqp_Omuses::obtain_structure");
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
      xk.flags[i] |= Omu_States::Discrete;
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
      xk.flags[i] |= Omu_States::Algebraic;
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
}

//--------------------------------------------------------------------------
void Hqp_Omuses::setup_vars(int k,
			    VECP x, VECP xmin, VECP xmax,
			    VECP u, VECP umin, VECP umax,
			    VECP c, VECP cmin, VECP cmax)
{
  if (!_prg) {
    error(E_NULL, "Hqp_Omuses::setup_vars");
  }

  int i, nx;
  int K = _prg->ks()->dim - 1;
  Omu_States &xk = _xs[k];
  Omu_Vars   &uk = _us[k];
  Omu_Vars   &ck = _cs[k];

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
    error(E_NULL, "Hqp_Omuses::init_vars");
  }

  int i;
  Omu_States &xk = _xs[k];
  Omu_Vector &uk = _us[k];
  int nx = x->dim;

  v_copy(_xs[k].initial, xk);
  nx = x->dim;
  for (i = 0; i < nx; i++)
    x[i] = xk[i];

  v_copy(_us[k].initial, uk);
  v_copy(uk, u);
}

//--------------------------------------------------------------------------
void Hqp_Omuses::setup_struct(int k,
			      VECP f0x, VECP f0u, int &f0_lin,
			      MATP fx, MATP fu, IVECP f_lin,
			      MATP cx, MATP cu, IVECP c_lin,
			      MATP Lxx, MATP Luu, MATP Lxu)
{
  // obtain discrete-time states and dF/dxt
  Omu_States &xk = _xs[k];
  Omu_Vector &uk = _us[k];
  obtain_structure(k, xk, uk);

  // for now don't analyze sparse structure as required by HQP
  // i.e. assume dense block matrices
}

//--------------------------------------------------------------------------
void Hqp_Omuses::init_simulation(int k, VECP x, VECP u)
{
  if (!_prg) {
    error(E_NULL, "Hqp_Omuses::init_simulation");
  }

  int i;
  int nx = x->dim;
  int nu = u->dim;
  Omu_States &xk = _xs[k];
  Omu_Vector &uk = _us[k];
  int nxt = xk->dim;

  for (i = 0; i < nx; i++)
    xk[i] = x[i];
  for (i = nx; i < nxt; i++)
    xk[i] = xk.initial[i];
  v_copy(u, uk);

  _prg->init_simulation(k, xk, uk);
  if (xk->dim != nxt || uk->dim != nu) {
    error(E_FORMAT, "Hqp_Omuses::init_simulation; modified vector sizes");
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
  int i, kk, kkend;
  int K = _prg->ks()->dim - 1;
  Omu_States &xk = _xs[k];
  Omu_Vector &uk = _us[k];
  int nxt = xk->dim;
  int nx = x->dim;
  int nu = u->dim;
  int nf = f->dim;
  int nxtf = max(nxt, nf);
  int nc = c->dim;
  int nd = xk.nd;

  assert(nd >= 0);	// obtain_structure() must have been called before
  if (!_prg) {
    error(E_NULL, "Hqp_Omuses::update_vals");
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
      catchall(// try
	       _integrator->solve(kk, _prg->ts(kk), _prg->ts(kk+1), xk, uk, _prg, _xt),
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
  Omu_States &xk = _xs[k];
  Omu_Vector &uk = _us[k];
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
    error(E_NULL, "Hqp_Omuses::update_stage");
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
      catchall(// try
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
    // update the criterion, constraints, and discrete states
    //

    // initialize dependend variables
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
    // no extra independend variables necessary
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


//========================================================================
