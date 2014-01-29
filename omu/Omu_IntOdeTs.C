/*
 * Omu_IntOdeTs.C --
 *   -- class for integrating an Ode over a stage
 *   -- using Taylor series expansion of the Ode
 *   -- derived from adolc driver routine forodec
 *
 * hl, 28/01/98
 */

/*
    Copyright (C) 1997--1999  Hartmut Linke

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

#ifdef OMU_WITH_ADOLC
#include <adolc/interfaces.h>
#include <adolc/adalloc.h>
#include <adolc/drivers/odedrivers.h>
#include <adolc/taping.h>
#include "adoublev.h"
#endif
#include "Omu_Program.h"

#include <If_Int.h>
#include <If_Bool.h>
#include <If_Real.h>
#include <If_Class.h>

#include "Omu_IntOdeTs.h"

IF_CLASS_DEFINE("OdeTs", Omu_IntOdeTs, Omu_Integrator);

#define mindec(a,b) if ((a) > (b)) (a) = (b)

//--------------------------------------------------------------------------

Omu_IntOdeTs::Omu_IntOdeTs()
{

  _adtaylor = false;
  _multiple_record = true;
  _check_autonomous = false;
  _output = 0;

  _max_deg = 31;
  _max_deg0 = 30;

  _kk = 0;

  _tau = 1.0;
  _rho = 0.0;

  _a = v_resize(v_get(1),0);

  _Sxd  = m_resize(m_get(1,1),0,0);
  _Sxd0 = m_resize(m_get(1,1),0,0);
  _Sxx  = m_resize(m_get(1,1),0,0);
  _Sxx0 = m_resize(m_get(1,1),0,0);
  _Sxu  = m_resize(m_get(1,1),0,0);
  _Sxu0 = m_resize(m_get(1,1),0,0);
  _Sh   = m_resize(m_get(1,1),0,0);

  _nnz = 0;
  _nz = NULL;

  _ifList.append(new If_Real("prg_int_rho", &_rho));
  _ifList.append(new If_Real("prg_int_tau", &_tau));
  _ifList.append(new If_Int("prg_int_maxdeg", &_max_deg0));
  _ifList.append(new If_Int("prg_int_out", &_output));
  _ifList.append(new If_Bool("prg_int_multrec", &_multiple_record));

}

//--------------------------------------------------------------------------

Omu_IntOdeTs::Omu_IntOdeTs(int max_deg)
{

  _adtaylor = false;
  _multiple_record = true;
  _output = 0;

  _max_deg = max_deg+1;
  _max_deg0 = max_deg;

  _kk = 0;

  _tau = 1.0;
  _rho = 0.0;

  _a = v_resize(v_get(1),0);

  _Sxd  = m_resize(m_get(1,1),0,0);
  _Sxd0 = m_resize(m_get(1,1),0,0);
  _Sxx  = m_resize(m_get(1,1),0,0);
  _Sxx0 = m_resize(m_get(1,1),0,0);
  _Sxu  = m_resize(m_get(1,1),0,0);
  _Sxu0 = m_resize(m_get(1,1),0,0);
  _Sh   = m_resize(m_get(1,1),0,0);

  _nnz = 0;
  _nz = NULL;

  _ifList.append(new If_Real("prg_int_rho", &_rho));
  _ifList.append(new If_Real("prg_int_tau", &_tau));
  _ifList.append(new If_Int("prg_int_out", &_output));
  _ifList.append(new If_Bool("prg_int_multrec", &_multiple_record));

}

//--------------------------------------------------------------------------

Omu_IntOdeTs::~Omu_IntOdeTs()
{

  int i;

  M_FREE(_Sxd);
  M_FREE(_Sxd0);
  M_FREE(_Sxx);
  M_FREE(_Sxx0);
  M_FREE(_Sxu);
  M_FREE(_Sxu0);
  M_FREE(_Sh);

  V_FREE(_a);


  if(_adtaylor)
    free_adtaylor();


  if(_nz != NULL) {
    for(i = 0; i < _nnz; i++)
      delete[] _nz[i];
    delete[] _nz;
  }
}

//--------------------------------------------------------------------------
void Omu_IntOdeTs::resize(int deg)
{
  if ((int)_Sxx->m == _nd && (int)_Sxd->n == _nd && (int)_Sxu->n == _nu 
      && (deg+1) <= _max_deg)
    return;

  if ((int)_Sxx->m == _n && (int)_Sxd->n == _nd &&  (int)_Sxu->n == _nu) {

    if(_adtaylor)
      free_adtaylor();

    _max_deg = deg+1;

    _aa1 = myalloc1(_nd+_n+_nu);
    _ab1 = myalloc1(_nd+_n+_nu);
    _aa = myalloc2(_nd+_n+_nu,_max_deg);
    _ab = myalloc2(_nd+_n+_nu,_max_deg);
    _ai = myalloc2(_nd+_n+_nu,_nd+_n+_nu);
    _aA = myalloc3(_nd+_n+_nu,_nd+_n+_nu,_max_deg);
    _aB = myalloc3(_nd+_n+_nu,_nd+_n+_nu,_max_deg);
    _adtaylor = true;

  }
  else {

    int i;

    if(_nz != NULL) {
      for(i = 0; i < _nnz; i++)
	delete[] _nz[i];
      delete[] _nz;
    }

    v_resize(_a,_nd+_n+_nu);

    m_resize(_Sxd0,_n,_nd);
    m_resize(_Sxd,_n,_nd);
    m_resize(_Sxx,_n,_n);
    m_resize(_Sxx0,_n,_n);
    m_resize(_Sxu,_n,_nu);
    m_resize(_Sxu0,_n,_nu);
    m_resize(_Sh,_n,_n);

    _cxp.resize(_nd + _n, _nx, _nu);

    _nnz = _nd+_n+_nu;

    _nz = new short*[_nnz];
    for(i = 0; i < _nnz; i++)
      _nz[i] = new short[_nnz];             // set up sparsity array

    if(_adtaylor)
      free_adtaylor();

    _max_deg = deg+1;

    _aa1 = myalloc1(_nd+_n+_nu);
    _ab1 = myalloc1(_nd+_n+_nu);
    _aa = myalloc2(_nd+_n+_nu,_max_deg);
    _ab = myalloc2(_nd+_n+_nu,_max_deg);
    _ai = myalloc2(_nd+_n+_nu,_nd+_n+_nu);
    _aA = myalloc3(_nd+_n+_nu,_nd+_n+_nu,_max_deg);
    _aB = myalloc3(_nd+_n+_nu,_nd+_n+_nu,_max_deg);
    _adtaylor = true;
  }

}
	    
//--------------------------------------------------------------------------
void Omu_IntOdeTs::init_stage(int k, const Omu_VariableVec &x, 
			      const Omu_VariableVec &u,
			      const Omu_DependentVec &Fc, bool sa)
{

  int  i;

  if (!(Fc.Jdx.is_diagonal())) {
    m_error(E_FORMAT,
	    "Omu_IntOdeTs::init_stage, only Ode allowed");
  }
  else {
    for(i = 0; i < _nx; i++)
      if (!((fabs(Fc.Jdx[i][i]) == 1) || (Fc.Jdx[i][i] == 0.0))) {
	std::cout << i << std::endl;
	m_output(Fc.Jdx);
	m_error(E_FORMAT,
		"Omu_IntOdeTs::init_stage, only Ode's allowed");
      }
  }

  if (_na > 0) {
    m_error(E_FORMAT,
	    "Omu_IntOdeTs::init_stage, which was called for algebraic states");
  }

  Omu_Integrator::init_stage(k, x, u, Fc, sa);
 
  if(_tau <= 0) {
    printf("time scaling must be greather then zero!\n");
    printf("use original time\n");
    _tau = 1.0;
  }

  resize(_max_deg0);

}

//--------------------------------------------------------------------------

void Omu_IntOdeTs::solve(int kk, double tstart, double tend,
	     const Omu_VariableVec &x, const Omu_VariableVec &u,
	     Omu_Program *sys, Omu_DependentVec &cF, Omu_StateVec &cx)
{

  bool   one_step, multiple_record;
  short  tag;
  int    deg, i, j, k, rc, tayl_finite;
  double err, h, hp, t, taut, tayl_max;
  adouble  au0;

  VECP   hpv;

  multiple_record = true;
  tag = 1;

  //  adoublev ax(_nd+_n);
  static adoublev ax; ax.alloc(_nd+_n);
  //  adoublev axp(_nd+_n);
  static adoublev axp; axp.alloc(_nd+_n);
  //  adoublev au(_m);
  static adoublev au; au.alloc(_m);
  //  adoublev aF(_nd+_n);
  static adoublev aF; aF.alloc(_nd+_n);

  double *f;

  f = new double[_nd+_n+_m];

  h = (tend - tstart);

  if(h <= 0.0)
    return;

  if(_output > 0) {
    printf("t0: %g tf: %g max_deg: %d\n",tstart,tend,_max_deg-1);
    v_output(cx);
  }

  hpv = v_get(_max_deg);

  rc = 3;
  one_step = false;
  t = tstart;
  au0 = 0.0;

  // check for autonomous system
  if(!_check_autonomous) {

    v_zero(_cxp);
    v_zero(cF);
    sys->continuous(kk, t, x, u, _cxp, cF);
    err = v_norm2(cF);
    v_zero(_cxp);
    v_zero(cF);
    sys->continuous(kk, t+0.1*h/sqrt(2.0), x, u, _cxp, cF);

    if(fabs(err - v_norm2(cF))/err > MACHEPS) {
      std::cout << kk << "   ";
      std::cout << err << "    " << v_norm2(cF) << "     ";
      std::cout << fabs(err - v_norm2(cF)) << std::endl;
      std::cout <<  "Only autonomous systems are allowed for OdeTs!";
      std::cout << std::endl;
      m_error(E_FORMAT,"Omu_IntOdeTs::solve");
    }

    if(kk < _kk)
      _check_autonomous = true;
    
    _kk++;

  }

  for(i = 0; i < _nd+_n; i++)
    _a[i] = cx[i];
  for(i = 0; i < _nu; i++)
    _a[i+_nd+_n] = u[i];

  if(_sa) {
    m_move(cx.Sx,_nd,0,_n-_nv,_nd,_Sxd0,0,0);
    m_move(cx.Sx,_nd,_nd,_n,_n-_nv,_Sxx0,0,0);
    m_move(cx.Su,_nd,0,_n,_nu,_Sxu0,0,0);
  }

  while(t < tend) {

    _res_evals++;
    if(_sa)
      _sen_evals++;

    err = _atol;
    hp = 1.0;
    deg = 0;
    tayl_finite = 1;

    if(multiple_record) {
      
      trace_on(tag);	
      
      for (i = 0; i < _nd+_n; i++)
	axp[i] = 0.0;
      for (i = 0; i < _nd+_n; i++)
	aF[i] = 0.0;

      for(i = 0; i < _nd+_n; i++)
	ax[i] <<= _a[i];
      for(i = 0; i < _nu; i++)
	au[i] <<= _a[i+_nd+_n];

      sys->continuous(kk, t, ax, au, axp, aF);

      for (i = 0; i < _nd+_n; i++)
	aF[i] >>= f[i];
      for(i = 0; i < _nu; i++)
	au0 >>= f[i+_nd+_n];

      trace_off();

      if(!_multiple_record)
	multiple_record = false;

    }

    for(i = 0; i < _nnz; i++) {
      _aa1[i] = _a[i];
      for(j = 0; j < _max_deg; j++)
	_aa[i][j] = 0.0;
    }	

    mindec(rc,zos_forward(tag,_nnz,_nnz,0,_aa1,_ab1)); 

    if( rc < 0)
      m_error(E_INTERN,"Omu_IntOdeTs::solve");

    taut = _tau/(1+deg);      
    tayl_max = 0.0;

    for (i = 0; i < _nnz; i++)
      _aa[i][0] = taut*_ab1[i];
    deg++;
 
    while(deg < (_max_deg-2)) {

      mindec(rc,hos_forward(tag,_nnz,_nnz,deg,0,_aa1,_aa,_ab1,_ab));	
      if( rc < 0)
        m_error(E_INTERN,"Omu_IntOdeTs::solve");

      taut = _tau/(1+deg);      
      tayl_max = 0.0;

      for(i = 0; i < _nnz; i++)
	if(!is_finite(_ab[i][deg-1])) {
	  tayl_finite = 0;
	  printf("Taylor-series contains infinite or NaN elements!\n");
	  m_error(E_INTERN,"Omu_IntOdeTs::solve");
	}

      for(i = 0; i < _nnz; i++) {
	_aa[i][deg] = taut*_ab[i][deg-1];
	tayl_max = max(tayl_max,fabs(_aa[i][deg]));
      }

      hp *= h/_tau;
      deg++;

      if((tayl_max*hp/(1.0-_rho)) < err) {
	one_step = true;
	break;
      }
    }

    // compute the next taylor coefficients and keep the
    // values in preparation

    mindec(rc,hos_forward(tag,_nnz,_nnz,deg,deg+1,_aa1,_aa,_ab1,_ab));

    if( rc < 0)
      m_error(E_INTERN,"Omu_IntOdeTs::solve");
    
    taut = _tau/(1+deg);
    tayl_max = 0.0;

    for(i = 0; i < _nnz; i++) {
      _aa[i][deg] = taut*_ab[i][deg-1];
      tayl_max = max(tayl_max,fabs(_aa[i][deg]));
    }

    for(i = 0; i < _nnz; i++) {
      for(k = deg+1; k > 0; k--)
	_aa[i][k] = _aa[i][k-1];
      _aa[i][0] = _aa1[i];
    }
    
    // determine h
    if(!one_step) {
      h = _tau*pow(err*(1.0-_rho)/tayl_max,1.0/(1.0*deg));
      if((tend - (t + h)) < h/3.0)
	h = 2.0*h/3.0;
      h = min(h,tend-t);
    }

    if(h == 0.0) {
      m_error(E_SING,"Omu_IntOdeTs::solve");
      printf("Taylor-coeficient: %g \n",tayl_max);
    }

    if(_output > 0) 
      printf("deg: %d time step: %g\n",deg+1,h);

    if(_output > 1) {
      printf("taylor-series:\n");
      for(j = 0; j < _nnz; j++) {
	for(i = 0; i < (deg+1); i++)
	  printf("%g  ",_aa[j][i]);
	printf("\n");
      }
    }

    for(i = _nd; i < _nd+_n; i++) {
      _a[i] = 0.0;
    }

    // compute continuous states after step
    hp = 1.0;

    for(j = 0; j < deg+1; j++) {
      for(i = _nd; i < _nd+_n; i++) {
	_a[i] += _aa[i][j]*hp;
      }
      hp *= h/_tau;
      hpv[j] = hp;
      
    }

    hp *= h/_tau;
    hpv[deg+1] = hp;

    if(_sa) {

      for(i = 0; i < _nnz; i++) {
	for(j = 0; j < _nnz; j++)
	  _ai[i][j] = 0.0;
	_ai[i][i] = 1.0;
      }

      hov_reverse(tag,_nnz,_nnz,deg,_nnz,_ai,_aA,_nz);
      
      accodec(_nnz,_tau,deg,_aA,_aB,_nz);

      if(_sa && (_output > 3)) {
	for(i=0;i<deg+1;i++){
	  printf("\n\t<-- _aB(%d)\n",i);
	  for(j = 0; j < _nnz; j++){
	    for(k = 0; k < _nnz; k++)
	      printf("%g  ",_aB[j][k][i]);
	    printf("\n");
	  }
	}
      }

      m_zero(_Sxd);
      m_ident(_Sxx);
      m_zero(_Sxu);

      for(i = _nd; i < _nd+_n ; i++)
	for(j = 0; j < _nd; j++)
	  for(k = 0; k < deg+1; k++)
	    _Sxd[i-_nd][j] += _aB[i][j][k]*hpv[k];

      for(i = _nd; i < _nd+_n ; i++)
	for(j = 0; j < _nu; j++)
	  for(k = 0; k < deg+1; k++)
	    _Sxu[i-_nd][j] += _aB[i][j+_nd+_n][k]*hpv[k];

      for(i = _nd; i < _nd+_n ; i++)
	for(j = 0; j < _n; j++)
	  for(k = 0; k < deg+1; k++)
	    _Sxx[i-_nd][j] += _aB[i][j+_nd][k]*hpv[k];

      m_resize(_Sh,_n,_nd);
      m_mlt(_Sxx0,_Sxd,_Sh);
      m_add(_Sxd0,_Sh,_Sxd0);

      m_resize(_Sh,_n,_nu);
      m_mlt(_Sxx0,_Sxu,_Sh);
      m_add(_Sxu0,_Sh,_Sxu0);

      m_resize(_Sh,_n,_n);
      m_copy(_Sxx0,_Sh);
      m_mlt(_Sxx,_Sh,_Sxx0);
            
    }

    t += h;
    h = (tend - t);

  }

  if(_sa) {

    m_move(_Sxd0,0,0,_n,_nd,cx.Sx,_nd,0);
    m_move(_Sxx0,0,0,_n,_n-_nv,cx.Sx,_nd,_nd);
    m_move(_Sxu0,0,0,_n,_nu,cx.Su,_nd,0);

      if(_output > 0) {
	m_output(cx.Sx);
	m_output(cx.Su);
      }
  }

  for(i = _nd; i < _nd+_n; i++)
    cx[i] = _a[i];

  if(_output > 0)
    v_output(cx);

  delete[] f;

  V_FREE(hpv);

}

//--------------------------------------------------------------------------


void Omu_IntOdeTs::free_adtaylor()
{

  if(_adtaylor) {
    free((char*) _aa1);
    free((char*) _ab1);
    free((char*) *_aa); free((char*) _aa);
    free((char*) *_ab); free((char*) _ab);
    free((char*) *_ai); free((char*) _ai);
    free((char*)**_aA); free((char*)*_aA); free((char*)_aA);
    free((char*)**_aB); free((char*)*_aB); free((char*)_aB);

    _adtaylor = false;

  }

}











