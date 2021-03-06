/*
 * Omu_IntDASPK.C --
 *   -- class implementation
 *
 * rf, 10/2/96
 */

/*
    Copyright (C) 1996--2014  Ruediger Franke and Hartmut Linke

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
  m_error(E_INTERN, "DASPK, which called PSOL");
}

//--------------------------------------------------------------------------
static void G_RES()
{
  m_error(E_INTERN, "DASPK, which called G_RES");
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

  _xc_ptr = NULL;
  _Fc_ptr = NULL;

  _Yq = m_get(1, 1);

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
  m_free(_Yq);

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
void Omu_IntDASPK::init(int k,
			const Omu_StateVec &x, const Omu_Vec &u,
			const Omu_DependentVec &Fc, bool sa)
{
  if (_jac_sbw >= 0) {
    // take user-defined value
    _ml = _mu = _jac_sbw;
    _banded_solver = _banded;
  }
  else {
    _ml = max(Fc.Jx.sbw_lower(), Fc.Jdx.sbw_lower());
    _mu = max(Fc.Jx.sbw_upper(), Fc.Jdx.sbw_upper());

    _banded_solver = _banded;
    // disable banded solver if full solver appears more efficient
    if (_ml+_ml+_mu >= _n)
      _banded_solver = false;
  }

  resize();
  init_options(Fc);
}

//--------------------------------------------------------------------------
void Omu_IntDASPK::resize()
{
  int kmp, maxl, maxord, neq, nrmax;

  if (_dxc->dim == _n && _q->dim == _nq && _senpar->dim == _nq)
    return;

  //
  // realloc variables for low level residual callback
  //
  v_resize(_q, _nq);
  _dxc.resize(_n, 0, 0, _nq);
  _xc_jac.resize(_n, _n, 0);
  _dxc_jac.resize(_n, _n, 0);
  m_resize(_Yq, _n, _nq);

  //
  // re-allocate variables for DASPK
  //

  neq = _n * (1 + _nq);

  v_resize(_y, neq);
  v_resize(_yprime, neq);
  v_resize(_senpar, _nq);

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
    v_resize(_rpar, _n + _nq);
    iv_resize(_ipar, _n + _nq);
  }
  else {
    //
    // iterative solver
    //

    v_resize(_rpar, _n + _nq);
    iv_resize(_ipar, 20 + _n + _nq);

    neq = _n * (1 + _nq);

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
      _lwp_basic += (maxl+3+max(0,min(1,maxl-kmp)))*neq/(_n + _nq + 1);
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
void Omu_IntDASPK::init_options(const Omu_DependentVec &Fc)
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

    neq = _n * (1 + _nq);

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

  if (_na > 0) {

    // assume inconsistent initial algebraic vars and derivatives
    // if algebraic states are present
    // todo: check if this is valid and check other methods
    //       to treat DAE initialization problem
    _info[11-1] = 1;
    int LID = 40;
    for (i = 0; i < _n; i++) {
      if (Fc.Jdx.is_zero_column(i))
	_iwork[LID + i] = -1;	// indicate algebraic state
      else
	_iwork[LID + i] = 2;	// indicate fixed state variable
    }
  }

  // sensitivity options
  _info[20-1] = 2;	// res() provides sensitivity equations
  _info[22-1] = _nq; 	// number of parameters appearing in res()
  _info[23-1] = !_serr;	// exclude sensitivities from error test
  _info[25-1] = 2;	// staggert direct method
}

//--------------------------------------------------------------------------
void Omu_IntDASPK::solve(int kk, double tstart, double tend,
			 Omu_StateVec &xc, Omu_StateVec &dxc, Omu_Vec &q,
			 Omu_DependentVec &Fc)
{
  int i, j;

  _xc_ptr = &xc; // propagate to res() and jac()
  _Fc_ptr = &Fc;

  // _stepsize overrides _nsteps
  int nsteps = _nsteps;
  if (_stepsize > 0.0)
    nsteps = (int)ceil((tend - tstart) / _stepsize);

  if (nsteps > 0) {
    // Fixed step size requires external Jacobian.
    _info[5-1] = 1;
  }
  else if (!_sys->has_low_level_continuous()) {
    // Don't use analytical Jacobian if not explicitly given
    // and if variable step size is enabled. DASPK seems to adapt
    // its behavior assuming a cheap analytical Jacobian evaluation,
    // which is not true when using ADOL-C.
    _info[5-1] = 0;
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

  for (i = 0; i < _n; i++) {
    _y[i] = xc[i];		// initial states
  }
  for (i = 0; i < _nq; i++) {
    _senpar[i] = q[i];
  }

  if (_sa) {
    for (i = 0; i < _n; i++) {
      for (j = 0; j < _nq; j++) {
	_y[(1 + j) * _n + i] = xc.Sq[i][j];
      }
    }
    m_zero(_dxc.Sq);
    m_zero(_xc_jac.Sx);
    m_zero(_dxc_jac.Sx);
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
    _info[19-1] = _nq;
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
    m_error(E_CONV,
	    "Omu_IntDASPK::solve that failed in calling DASPK");
  }

  //
  // read and return results
  //

  for (i = 0; i < _n; i++) {
    xc[i] = _y[i];
  }

  if (_sa) {
    for (i = 0; i < _n; i++) {
      for (j = 0; j < _nq; j++) {
	xc.Sq[i][j] = _y[(1 + j) * _n + i];
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
  Omu_StateVec &xc = *_xc_ptr;
  Omu_DependentVec &Fc = *_Fc_ptr;

  // prepare call arguments

  for (i = 0; i < _n; i++) {
    xc[i] = y[i];
    _dxc[i] = yprime[i];
  }
  for (i = 0; i < _nq; i++) {
    _q[i] = senpar[i];
  }

  if (*ires == 1) {
    // init seed derivatives
    for (i = 0; i < _n; i++) {
      for (j = 0; j < _nq; j++) {
	xc.Sq[i][j] = y[(1 + j) * _n + i];
	_dxc.Sq[i][j] = yprime[(1 + j) * _n + i];
      }
    }
    Fc.set_required_J(true);
  }
  else
    Fc.set_required_J(false);

  // evaluate residual
  residual(_kk, *t, xc, _dxc, _q, Fc);

  // read and return result
  for (i = 0; i < _n; i++)
    delta[i] = Fc[i];

  if (*ires == 1) {
    // _Yq = dF/dq = pF/px*dx/dq + pF/pdx*ddx/dq + pF/pq
    m_mlt(Fc.Jx, xc.Sq, _Yq);
    m_mltadd(_Yq, Fc.Jdx, _dxc.Sq, _Yq);
    m_add(_Yq, Fc.Jq, _Yq);

    // write result
    for (i = 0; i < _n; i++) {
      for (j = 0; j < _nq; j++) {
	delta[(1 + j) * _n + i] = _Yq[i][j];
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
    m_error(E_INTERN, "Omu_IntDASPK::jac");

  int i, j;
  Omu_DependentVec &Fc = *_Fc_ptr;

  // prepare call arguments

  for (i = 0; i < _n; i++) {
    _xc_jac[i] = y[i];
    _dxc_jac[i] = yprime[i];
  }
  for (i = 0; i < _nq; i++) {
    _q[i] = senpar[i];
  }

  for (i = 0; i < _n; i++) {
    _xc_jac.Sx[i][i] = 1.0;
    _dxc_jac.Sx[i][i] = *cj;
  }

  // evaluate residuals and Jacobians
  Fc.set_required_J(true);

  residual(_kk, *t, _xc_jac, _dxc_jac, _q, Fc);

  // read and return results

  if (!_banded_solver)
    // full Jacobian
    for (i = 0; i < _n; i++) {
      for (j = 0; j < _n; j++) {
	pd[i + _n * j] = Fc.Jx[i][j] + *cj * Fc.Jdx[i][j];
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
	pd[irow + nb * j] = Fc.Jx[i][j] + *cj * Fc.Jdx[i][j];
      }
    }
  }

  _jac_evals++;
}


//========================================================================
