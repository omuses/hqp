/*
 * Omu_IntDASPK.C --
 *   -- class implementation
 *
 * rf, 10/2/96
 */

/*
    Copyright (C) 1996--2001  Ruediger Franke and Hartmut Linke

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
#include <malloc.h>

#include <If_Int.h>
#include <If_Bool.h>
#include <If_Class.h>

#include "Omu_IntDASPK.h"

IF_CLASS_DEFINE("DASPK", Omu_IntDASPK, Omu_Integrator);

/*
 * functions for communication with DASPK
 */

// fuctions defined by Omu_IntDASPK
typedef void
(res_t)(freal *T, freal *Y, freal *XPRIME, freal *CJ,
	freal *DELTA, fint *IRES, freal *RPAR, fint *IPAR,
	freal *SENPAR);

typedef void
(jac_t)(freal *T, freal *Y, freal *YPRIME,
	freal *PD, freal *CJ, freal *RPAR, fint *IPAR,
	freal *SENPAR, fint *IJAC);

typedef  void
(kryjac_t)(res_t *RES, fint *IRES, fint *NEQ, freal *T, freal *Y,
	   freal *YPRIME, freal *REWT, freal *SAVR, freal *WK,
	   freal *H, freal *CJ, freal *WP, fint *IWP, fint *IERR,
	   freal *RPAR, fint *IPAR, freal *SENPAR);

typedef void
(psol_t)(fint *NEQ, freal *T, freal *Y, freal *YPRIME, freal *SAVR,
	 freal *WK, freal *CJ, freal *WGHT, freal *WP, fint *IWP,
	 freal *B, freal *EPLIN, fint *IER, freal *RPAR, fint *IPAR,
	 freal *SENPAR);

typedef void
(g_res_t)(/* currently not of interest */);

extern "C" {

  void dspsetup_(fint *NEQ, fint *LWP, fint *LIWP, freal *RPAR,
		 fint *IPAR, fint *IERR, fint *LWP_MIN,
		 fint *LIWP_MIN);

  void ddaspk_(res_t *RES, fint *NEQ, freal *T, freal *Y, freal *YPRIME,
	       freal *TOUT, fint *INFO, freal *RTOL, freal *ATOL,
	       fint *IDID, freal *RWORK, fint *LRW, fint *IWORK, fint *LIW,
	       freal *RPAR, fint *IPAR, jac_t *JAC, psol_t *PSOL,
	       freal *SENPAR, g_res_t *G_RES);

  kryjac_t dbanja_;	// iterative solver: generate banded preconditioner
  psol_t dbanps_;	// iterative solver: solve linear system
  kryjac_t djacilu_;	// iterative solver: generate ilu preconditioner
  psol_t dpsolilu_;	// iterative solver: solve linear system
}

Omu_IntDASPK *theOmu_IntDASPK;

//--------------------------------------------------------------------------
static void RES(freal *t, freal *y, freal *yprime, freal *,
		freal *delta, fint *ires, freal *rpar, fint *ipar,
		freal *senpar)
{
  theOmu_IntDASPK->res(t, y, yprime, delta, ires, rpar, ipar, senpar);
}

//--------------------------------------------------------------------------
static void PSOL()
{
  error(E_INTERN, "DASPK, which called PSOL");
}

//--------------------------------------------------------------------------
static void G_RES()
{
  error(E_INTERN, "DASPK, which called G_RES");
}

//--------------------------------------------------------------------------
static void JAC(freal *t, freal *y, freal *yprime,
                freal *pd, freal *cj, freal *rpar, fint *ipar,
		freal *senpar, fint *ijac)
{
  theOmu_IntDASPK->jac(t, y, yprime, pd, cj, rpar, ipar, senpar, ijac);
}

//--------------------------------------------------------------------------
Omu_IntDASPK::Omu_IntDASPK()
{
  theOmu_IntDASPK = this;

  _sys = NULL;
  _xc_ptr = NULL;
  _Fc_ptr = NULL;

  _Yx = m_get(1, 1);
  _Yu = m_get(1, 1);

  _ml = 0;
  _mu = 0;
  _y = v_get(1);
  _yprime = v_get(1);
  _info = iv_get(30);
  _rwork = v_get(1);
  _iwork = iv_get(1);
  _rpar = v_get(1);
  _ipar = iv_get(1);
  _senpar = v_get(1);

  _banded = true;
  _krylov = false;
  _krylov_prec = false;
  _with_jac = true;

  _nsteps = 0;

  _lwp = 0;
  _liwp = 0;
  _lwp_basic = 0;

  _jac_sbw = -1;	// i.e. automatic detection
  _ifList.append(new If_Int("prg_int_jac_sbw", &_jac_sbw));
  _ifList.append(new If_Bool("prg_int_banded", &_banded));
  _ifList.append(new If_Bool("prg_int_krylov", &_krylov));
  _ifList.append(new If_Bool("prg_int_krylov_prec", &_krylov_prec));
  _ifList.append(new If_Bool("prg_int_with_jac", &_with_jac));
  _ifList.append(new If_Int("prg_int_nsteps", &_nsteps));

}

//--------------------------------------------------------------------------
Omu_IntDASPK::~Omu_IntDASPK()
{
  m_free(_Yx);
  m_free(_Yu);

  v_free(_y);
  v_free(_yprime);
  iv_free(_info);
  v_free(_rwork);
  iv_free(_iwork);
  v_free(_rpar);
  iv_free(_ipar);
  v_free(_senpar);
}

//--------------------------------------------------------------------------
void Omu_IntDASPK::init_stage(int k,
			      const Omu_States &x, const Omu_Vector &u,
			      bool sa)
{

  Omu_Integrator::init_stage(k, x, u, sa);

  if (_jac_sbw >= 0) {
    // take user-defined value
    _ml = _mu = _jac_sbw;
    _banded_solver = _banded;
  }
  else {
    _ml = x.sbw_l;
    _mu = x.sbw_u;

    _banded_solver = _banded;
    // disable banded solver if full solver appears more efficient
    if (_ml+_ml+_mu >= _n)
      _banded_solver = false;
  }

  realloc();
  init_options(x);
}

//--------------------------------------------------------------------------
void Omu_IntDASPK::realloc()
{
  int kmp, maxl, maxord, neq, nrmax;

  if (_xcp->dim == _nd + _n && _uc->dim == _nu && _senpar->dim == _nd + _nu)
    return;

  //
  // realloc variables for low level _sys->continuous callback
  //
  v_resize(_uc, _nu);
  _xcp.realloc(_nd + _n, _nx, _nu);
  _xc_jac.realloc(_nd + _n, _nx, _nu);
  _xcp_jac.realloc(_nd + _n, _nx, _nu);
  m_resize(_Yx, _nd + _n, _nx);
  m_resize(_Yu, _nd + _n, _nu);

  //
  // re-allocate variables for DASPK
  //

  neq = _n * (1 + _nx + _nu);

  v_resize(_y, neq);
  v_resize(_yprime, neq);
  v_resize(_senpar, _nd + _nu);

  if (!_krylov) {
    //
    // direct matrix solver
    //

    int itmp1, itmp3;
    int lrw = 50 + (5 + 4)*neq;
    // if _info[16-1] = 1, i.e. no error test for algebraic variables:
    lrw += neq;
    if (_banded_solver) {
      lrw += (2*_ml + _mu + 1)*_n;
      // if _info[5-1] = 0, i.e. no analytical Jacobian:
      itmp1 = 2*(_n/(_ml+_mu+1)+1);
    }
    else {
      lrw += _n*_n;
      itmp1 = 0;
    }
    // if _info[20-1] < 2, i.e. no analytical sensitivity equations
    itmp3 = 4*_n;
    lrw += itmp1 > itmp3? itmp1: itmp3;

    v_resize(_rwork, lrw);
    iv_resize(_iwork, 40 + neq + _n);

    // store initial states and controls in RPAR/IPAR as proposed by DASPK
    v_resize(_rpar, _nd + _n + _nu);
    iv_resize(_ipar, _nd + _n + _nu);
  }
  else {
    //
    // iterative solver
    //

    v_resize(_rpar, _nd + _n + _nu);
    iv_resize(_ipar, 20 + _nd + _n + _nu);

    neq = _n * (1 + _nx + _nu);

    // daspk standard settings
    maxl = min(5,neq);
    kmp = maxl;
    nrmax = 5;
    maxord = 5;

    _lwp_basic = 50 + max(maxord+5,9)*neq + (maxl+3)*maxl+1;
    _lwp = 0;
    _liwp = 0;

    /*
    if(_sa)
      _lwp_basic += (maxl+3+max(0,min(1,maxl-kmp)))*neq/(_nd + _n + _nu + 1);
    else
    */

    _lwp_basic += (maxl+3+max(0,min(1,maxl-kmp)))*neq;

    if (!_krylov_prec) {
      v_resize(_rwork, _lwp_basic + _lwp);
      iv_resize(_iwork, 40 + 4*_n);
    }
    else {
      if (_banded_solver) {
	// banded preconditioner
	_lwp = (2*_ml + _mu + 1)*neq + 2*((neq/(_ml + _mu+1)) + 1);
	_liwp = neq;

	v_resize(_rwork, _lwp_basic + _lwp);
	iv_resize(_iwork, 40+2*neq+1+_liwp);
      }
      else {
	// parameters for incomplete LU factorization
	v_resize(_rwork, _lwp_basic + _lwp);
	iv_resize(_iwork, 40+2*neq+1+_liwp);
      }
    }
  }

  v_zero(_rwork);
  iv_zero(_iwork);
  v_zero(_rpar);
  iv_zero(_ipar);
}

//--------------------------------------------------------------------------
void Omu_IntDASPK::init_options(const Omu_States &x)
{

  int i, neq, n;

  iv_zero(_info);	// lots of defaults

  if (!_krylov) {
    _info[5-1] = _with_jac? 1: 0;    	// use analytical Jacobian

    if (_banded_solver) {
      _info[6-1] = 1;
      _iwork[1-1] = _ml;
      _iwork[2-1] = _mu;
    }
  }
  else {
    _info[12-1] = 1;	// use Krylov method

    neq = _n * (1 + _nx + _nu);

    if (!_krylov_prec) {
      _info[15-1] = 0;	// don't use preconditioner
    }
    else {
      _info[15-1] = 1;

      if (_banded_solver) {
	_ipar[1-1] = _ml;
	_ipar[2-1] = _mu;
      }
      else {
	// parameters for incomplete LU factorization

	fint lwp_act, liwp_act, IDID;

	n = (int) max(1.0,0.25*_n);

	_ipar[1-1] = min(_ml,n);
	_ipar[2-1] = min(_mu,n);
	_ipar[3-1] = 2;
	_ipar[4-1] = 2;
	_ipar[5-1] = 1;
	_ipar[6-1] = n;
	_ipar[7-1] = 0;
	_ipar[8-1] = 1;
	_ipar[9-1] = 0;
	_ipar[10-1] = 0;
	_ipar[11-1] = 1;

	_rpar[1-1] = 0.001;
	_rpar[2-1] = 0.01;

	lwp_act = _lwp;
	liwp_act = _liwp;

	// check
	dspsetup_(&neq, &lwp_act, &liwp_act, _rpar->ve, _ipar->ive, &IDID,
		  &_lwp, &_liwp);

	if(_lwp != lwp_act)
	  v_resize(_rwork, 50 + 10*neq + (neq+3)*neq + (neq+3)*neq + 1 + _lwp);
	v_zero(_rwork);

	if(_liwp != liwp_act)
	  iv_resize(_iwork, 40+2*neq+1+_liwp);
	iv_zero(_iwork);
      }
    }
    _iwork[27-1] = _lwp;
    _iwork[28-1] = _liwp;
  }

  if (x.na > 0) {

    // assume inconsistent initial algebraic vars and derivatives
    // if algebraic states are present
    // todo: check if this is valid and check other methods
    //       to treat DAE initialization problem
    _info[11-1] = 1;
    int LID = 40;
    for (i = 0; i < _n; i++) {
      if (x.flags[_nd+i] & Omu_States::Algebraic)
	_iwork[LID + i] = -1;	// indicate algebraic state
      else
	_iwork[LID + i] = 2;	// indicate fixed state variable
    }
  }

  // sensitivity options
  _info[20-1] = 2;	// res() provides sensitivity equations
  _info[22-1] = _nd + _nu; // number of parameters appearing in res()
  _info[23-1] = !_serr;	// exclude sensitivities from error test
  _info[25-1] = 2;	// staggert direct method
}

//--------------------------------------------------------------------------
void Omu_IntDASPK::solve(int kk, double tstart, double tend,
			 const Omu_States &x, const Omu_Vector &u,
			 Omu_Program *sys, Omu_DepVec &Fc, Omu_SVec &xc)
{
  int i, j;

  _sys = sys;	// propagate to res() and jac()
  _xc_ptr = &xc;
  _Fc_ptr = &Fc;

  // _stepsize overrides _nsteps
  int nsteps = _nsteps;
  if (_stepsize > 0.0)
    nsteps = (int)ceil((tend - tstart) / _stepsize);

  if (nsteps > 0) {
    // Fixed step size requires external Jacobian.
    _info[5-1] = 1;
  }

  //
  // init DASPK variables
  //

  fint NEQ;
  freal T = tstart;
  freal *Y = _y->ve;
  freal *YPRIME = _yprime->ve;
  freal TOUT = tend;
  fint  *INFO = _info->ive;
  freal RTOL = _rtol;
  freal ATOL = _atol;
  fint IDID;
  freal *RWORK = _rwork->ve;
  fint LRW = _rwork->dim;
  fint *IWORK = _iwork->ive;
  fint LIW = _iwork->dim;
  freal *RPAR = _rpar->ve;
  fint *IPAR = _ipar->ive;
  jac_t *JAC_P = &JAC;
  psol_t *PSOL_P = (psol_t*)&PSOL;
  freal *SENPAR = _senpar->ve;
  g_res_t *G_RES_P = &G_RES;

  if (_krylov) {
    if (!_krylov_prec) {
      PSOL_P = (psol_t *)&dpsolilu_;
    }
    else {
      if (_banded_solver) {
        JAC_P = (jac_t *)&dbanja_;
        PSOL_P = (psol_t *)&dbanps_;
      }
      else {
        JAC_P = (jac_t *)&djacilu_;
        PSOL_P = (psol_t *)&dpsolilu_;
      }
    }
  }

  v_zero(_y);
  v_zero(_yprime);

  for (i = 0; i < _nd; i++) {
    _senpar[i] = xc[i];
  }
  for (i = 0; i < _n; i++) {
    _y[i] = xc[_nd + i];		// initial states
  }
  for (i = 0; i < _nu; i++) {
    _senpar[_nd + i] = u[i];
  }

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
    m_zero(_xc_jac.Sx);
    m_zero(_xc_jac.Su);
    m_zero(_xcp_jac.Sx);
    m_zero(_xcp_jac.Su);
  }

  _info[1-1] = 0;	// restart new problem
  // Restarting a new problem might not be needed for multiple
  // sample periods in one stage.
  // But then the sensitivity error check should be enabled
  // (_serr, which goes to info[23-1])!
  // This check seems to be less important if problem is restarted.

  //_info[4-1] = 1;    	// don't go beyond TOUT
  //_rwork[1-1] = TOUT;

  // take user defined number of steps if nsteps > 0
  if (nsteps > 0) {
    //_info[3-1] = 1;			// return after each step
    RTOL = 1e37;			// disable error test
    ATOL = 1e37;
    _info[7-1] = 1;			// user defined max. step size
    _rwork[2-1] = (tend - tstart)/nsteps;	// maximal step size
    _info[8-1] = 1;			// user defined initial step size
    _rwork[3-1] = (tend - tstart)/nsteps;	// initial step size
    //_info[9-1] = 1;			// restrict maxord
    //_iwork[3-1] = 1;			// maxord=1
  }

  NEQ = _n;
  if (_sa) {
    _info[19-1] = _nx + _nu;
    NEQ += _n * _info[19-1];
  }
  else {
    _info[19-1] = 0;
  }

  //
  // integrate the system equations over the stage
  //
  while (true) {
    ddaspk_(&RES, &NEQ, &T, Y, YPRIME, &TOUT, INFO, &RTOL, &ATOL,
	    &IDID, RWORK, &LRW, IWORK, &LIW, RPAR, IPAR,
	    JAC_P, PSOL_P, SENPAR, G_RES_P);
    if (IDID != -1)
      break;
    _info[1-1] = 1;
  }

  if (IDID != 3 || T != TOUT) {
    /*
    fprintf(stderr, "equations: %d\n",NEQ);
    v_foutput(stderr, x);
    v_foutput(stderr, u);
    v_foutput(stderr, _yprime);
    v_foutput(stderr, _y);
    */
    error(E_UNKNOWN,
	  "Omu_IntDASPK::solve that failed in calling DASPK");
  }

  //
  // read and return results
  //

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

//--------------------------------------------------------------------------
void Omu_IntDASPK::res(freal *t, freal *y, freal *yprime,
		       freal *delta, fint *ires, freal *rpar, fint *,
		       freal *senpar)
{
  int i, j;
  Omu_SVec &xc = *_xc_ptr;
  Omu_DepVec &Fc = *_Fc_ptr;

  // prepare call arguments

  for (i = 0; i < _nd; i++) {
    xc[i] = senpar[i];
    _xcp[i] = 0.0;
  }
  for (i = 0; i < _n; i++) {
    xc[_nd + i] = y[i];
    _xcp[_nd + i] = yprime[i];
  }
  for (i = 0; i < _nu; i++) {
    _uc[i] = senpar[_nd + i];
  }

  if (*ires == 1) {
    // init seed derivatives
    for (i = 0; i < _n; i++) {
      for (j = 0; j < _nx; j++) {
	xc.Sx[_nd + i][j] = y[(1 + j) * _n + i];
	_xcp.Sx[_nd + i][j] = yprime[(1 + j) * _n + i];
      }
      for (j = 0; j < _nu; j++) {
	xc.Su[_nd + i][j] = y[(1 + _nx + j) * _n + i];
	_xcp.Su[_nd + i][j] = yprime[(1 + _nx + j) * _n + i];
      }
    }
  }

  // evaluate residual
  _sys->continuous(_kk, *t, xc, _uc, _xcp, Fc);

  // read and return result
  for (i = 0; i < _n; i++)
    delta[i] = Fc[_nd+i];

  if (*ires == 1) {
    // _Yx = dF/dxk = dF/dx * dx/dxk + dF/dxp * dxp/dxk
    m_mlt(Fc.Jx, xc.Sx, _Yx);
    m_mltadd(_Yx, Fc.Jxp, _xcp.Sx, _Yx);

    // _Yu = dF/duk = dF/dx * dx/duk + dF/dxp * dxp/duk + dF/duk
    m_mlt(Fc.Jx, xc.Su, _Yu);
    m_mltadd(_Yu, Fc.Jxp, _xcp.Su, _Yu);
    m_add(_Yu, Fc.Ju, _Yu);

    // write result
    for (i = 0; i < _n; i++) {
      for (j = 0; j < _nx; j++) {
	delta[(1 + j) * _n + i] = _Yx[_nd + i][j];
      }
      for (j = 0; j < _nu; j++) {
	delta[(1 + _nx + j) * _n + i] = _Yu[_nd + i][j];
      }
    }
  }

  _res_evals++;
  if (*ires == 1)
    _sen_evals++;
}

//--------------------------------------------------------------------------
void Omu_IntDASPK::jac(freal *t, freal *y, freal *yprime,
		       freal *pd, freal *cj, freal *rpar, fint *,
		       freal *senpar, fint *ijac)
{
  if (*ijac != 0)
    // don't know how to treat
    error(E_INTERN, "Omu_IntDASPK::jac");

  int i, j;
  Omu_DepVec &Fc = *_Fc_ptr;

  // prepare call arguments

  for (i = 0; i < _nd; i++) {
    _xc_jac[i] = senpar[i];
    _xcp_jac[i] = 0.0;
  }
  for (i = 0; i < _n; i++) {
    _xc_jac[_nd + i] = y[i];
    _xcp_jac[_nd + i] = yprime[i];
  }
  for (i = 0; i < _nu; i++) {
    _uc[i] = senpar[_nd + i];
  }

  for (i = _nd; i < _nxt; i++) {
    _xc_jac.Sx[i][i] = 1.0;
    _xcp_jac.Sx[i][i] = *cj;
  }

  // evaluate residuals and Jacobians
  bool is_required_J_bak = Fc.is_required_J();
  Fc.set_required_J(true);

  _sys->continuous(_kk, *t, _xc_jac, _uc, _xcp_jac, Fc);

  Fc.set_required_J(is_required_J_bak);

  // read and return results

  if (!_banded_solver)
    // full Jacobian
    for (i = 0; i < _n; i++) {
      for (j = 0; j < _n; j++) {
	pd[i + _n * j] = Fc.Jx[_nd+i][_nd+j] + *cj * Fc.Jxp[_nd+i][_nd+j];
      }
    }
  else {
    // banded Jacobian
    int j_start, j_end;
    int nb = 2*_ml + _mu + 1;
    int irow;
    for (i = 0; i < _n; i++) {
      j_start = i - _ml;
      j_start = max(j_start, 0);
      j_end = i + _mu + 1;
      j_end = min(j_end, _n);
      for (j = j_start; j < j_end; j++) {
	irow = i - j + _ml + _mu;
	pd[irow + nb * j] = Fc.Jx[_nd+i][_nd+j] + *cj * Fc.Jxp[_nd+i][_nd+j];
      }
    }
  }

  _jac_evals++;
}


//========================================================================
