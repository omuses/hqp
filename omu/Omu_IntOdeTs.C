/*
 * Omu_IntOdeTs.C --
 *   -- class for integrating an Ode over a stage
 *   -- using Taylor series
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
#include <math.h>

#include <adutils.h>

#include <Hqp.h>
#include <If_Int.h>
#include <If_Bool.h>
#include <If_Float.h>
#include <If_Class.h>

#include "Omu_IntOdeTs.h"

IF_CLASS_DEFINE("OdeTs", Omu_IntOdeTs, Omu_Integrator);

#define mindec(a,b) if ((a) > (b)) (a) = (b)

//--------------------------------------------------------------------------

Omu_IntOdeTs::Omu_IntOdeTs()
{

  int i;

  _adtaylor = false;
  _multiple_record = true;
  _output = 0;

  _max_deg = 31;
  _max_deg0 = 30;

  _n = 0;
  _m = 0;
  _kk = 0;

  _tau = 1.0;
  _rho = 0.0;

//    _atol = 1.0e-6;
//    _rtol = 1.0e-3;

  _a = v_resize(v_get(1),0);
  _u = v_resize(v_get(1),0);
  _x = v_resize(v_get(1),0);

  _Sxd  = m_resize(m_get(1,1),0,0);
  _Sxd0 = m_resize(m_get(1,1),0,0);
  _Sxx  = m_resize(m_get(1,1),0,0);
  _Sxx0 = m_resize(m_get(1,1),0,0);
  _Sxu  = m_resize(m_get(1,1),0,0);
  _Sxu0 = m_resize(m_get(1,1),0,0);
  _Sh   = m_resize(m_get(1,1),0,0);

  _nz = NULL;

//    _ifList.append(new If_Float("prg_int_atol", &_atol));
//    _ifList.append(new If_Float("prg_int_rtol", &_rtol));
  _ifList.append(new If_Float("prg_int_rho", &_rho));
  _ifList.append(new If_Float("prg_int_tau", &_tau));
  _ifList.append(new If_Int("prg_int_maxdeg", &_max_deg0));
  _ifList.append(new If_Int("prg_int_out", &_output));
  _ifList.append(new If_Bool("prg_int_multrec", &_multiple_record));

}

//--------------------------------------------------------------------------

Omu_IntOdeTs::Omu_IntOdeTs(int max_deg)
{

  int i;

  _adtaylor = false;
  _multiple_record = true;
  _output = 0;

  _max_deg = max_deg+1;
  _max_deg0 = max_deg;

  _n = 0;
  _m = 0;
  _kk = 0;

  _tau = 1.0;
  _rho = 0.0;

//    _atol = 1.0e-3;
//    _rtol = 1.0e-6;

  _a = v_resize(v_get(1),0);
  _u = v_resize(v_get(1),0);
  _x = v_resize(v_get(1),0);

  _Sxd  = m_resize(m_get(1,1),0,0);
  _Sxd0 = m_resize(m_get(1,1),0,0);
  _Sxx  = m_resize(m_get(1,1),0,0);
  _Sxx0 = m_resize(m_get(1,1),0,0);
  _Sxu  = m_resize(m_get(1,1),0,0);
  _Sxu0 = m_resize(m_get(1,1),0,0);
  _Sh   = m_resize(m_get(1,1),0,0);

  _nz = NULL;

//    _ifList.append(new If_Float("prg_int_atol", &_atol));
//    _ifList.append(new If_Float("prg_int_rtol", &_rtol));
  _ifList.append(new If_Float("prg_int_rho", &_rho));
  _ifList.append(new If_Float("prg_int_tau", &_tau));
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

  V_FREE(_x);
  V_FREE(_u);
  V_FREE(_a);


  if(_adtaylor)
    free_adtaylor();


  if(_nz != NULL) {
    for(i = 0; i < _nd+_n+_m; i++)
      delete[] _nz[i];
    delete[] _nz;
  }
}

//--------------------------------------------------------------------------
void Omu_IntOdeTs::realloc(int nd, int n, int nv, int m, int deg)
{
  if (nd == _nd && n == _n && m == _m && nv == _nv && (deg+1) <= _max_deg)
    return;

  int i;

  if (nd == _nd && n == _n && m == _m && nv == _nv) {

    if(_adtaylor)
      free_adtaylor();

    _max_deg = deg+1;

    _aa1 = myalloc1(_nd+_n+_m);
    _ab1 = myalloc1(_nd+_n+_m);
    _aa = myalloc2(_nd+_n+_m,_max_deg);
    _ab = myalloc2(_nd+_n+_m,_max_deg);
    _ai = myalloc2(_nd+_n+_m,_nd+_n+_m);
    _aA = myalloc3(_nd+_n+_m,_nd+_n+_m,_max_deg);
    _aB = myalloc3(_nd+_n+_m,_nd+_n+_m,_max_deg);
    _adtaylor = true;

  }
  else {

    if(_nz != NULL) {
      for(i = 0; i < _nd+_n+_m; i++)
	delete[] _nz[i];
      delete[] _nz;
    }

    _n = n;
    _m = m;
    _nv = nv;

    assert(_nv == 0);	        	       // not supported by Ode


    v_resize(_a,_nd+_n+_m);
    v_resize(_x,_nd+_n+_m);

    m_resize(_Sxd0,_n,_nd);
    m_resize(_Sxd,_n,_nd);
    m_resize(_Sxx,_n,_n);
    m_resize(_Sxx0,_n,_n);
    m_resize(_Sxu,_n,_m);
    m_resize(_Sxu0,_n,_m);
    m_resize(_Sh,_n,_n);

    _nz = new short*[_nd+_n+_m];
    for(i = 0; i < _nd+_n+_m; i++)
      _nz[i] = new short[_nd+_n+_m];             // set up sparsity array

    if(_adtaylor)
      free_adtaylor();

    _max_deg = deg+1;

    _aa1 = myalloc1(_nd+_n+_m);
    _ab1 = myalloc1(_nd+_n+_m);
    _aa = myalloc2(_nd+_n+_m,_max_deg);
    _ab = myalloc2(_nd+_n+_m,_max_deg);
    _ai = myalloc2(_nd+_n+_m,_nd+_n+_m);
    _aA = myalloc3(_nd+_n+_m,_nd+_n+_m,_max_deg);
    _aB = myalloc3(_nd+_n+_m,_nd+_n+_m,_max_deg);
    _adtaylor = true;
  }

}
	    

//--------------------------------------------------------------------------

void Omu_IntOdeTs::solve(int kk, Real tstart, Real tend,
			 const Omu_States &x, const Omu_Vector &u,
			 Omu_Program *sys, VECP xt,
			 MATP Sx, MATP Su)
{

  bool   one_step, multiple_record;
  short  tag;
  int    deg, i, j, k, n, rc, tayl_finite;
  double err, h, hp, t, taut, tayl_max;
  adouble  au0;

  VECP   hpv;

  multiple_record = true;
  tag = 1;

  realloc(_nd, _nxt - _nd - _nv, _nv, _nu,_max_deg0);

  adoublev ax(_nd+_n);
  adoublev axp(_nd+_n);
  adoublev au(_m);
  adoublev aF(_nd+_n);

  double f[_nd+_n+_m];

  if(_tau <= 0) {
    printf("time scaling must be greather then zero!\n");
    printf("use original time\n");
    _tau = 1.0;
  }

  h = (tend - tstart);

  if(h <= 0.0)
    return;

  if(_output > 0) {
    printf("t0: %g tf: %g max_deg: %d\n",tstart,tend,_max_deg-1);
    v_output(xt);
  }

  hpv = v_get(_max_deg);

  rc = 3;
  one_step = false;
  t = tstart;
  n = _nd + _n + _m;
  au0 = 0.0;

  for(i = 0; i < _nd+_n; i++)
    _a[i] = xt[i];
  for(i = 0; i < _m; i++)
    _a[i+_nd+_n] = u[i];

  if(_sa) {
    m_move(Sx,_nd,0,_n,_nd,_Sxd0,0,0);
    m_move(Sx,_nd,_nd,_n,_n,_Sxx0,0,0);
    m_move(Su,_nd,0,_n,_m,_Sxu0,0,0);
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
      for(i = 0; i < _m; i++)
	au[i] <<= _a[i+_nd+_n];

      sys->continuous(_kk, t, ax, au, axp, aF);

      for (i = 0; i < _nd+_n; i++)
	aF[i] >>= f[i];
      for(i = 0; i < _m; i++)
	au0 >>= f[i+_nd+_n];

      trace_off();

      if(!_multiple_record)
	multiple_record = false;

    }

    for(i = 0; i < n; i++) {
      _aa1[i] = _a[i];
      for(j = 0; j < _max_deg; j++)
	_aa[i][j] = 0.0;
    }	

    mindec(rc,zos_forward(tag,n,n,0,_aa1,_ab1)); 

    if( rc < 0)
      error(E_INTERN,"Omu_IntOdeTs::solve");

    taut = _tau/(1+deg);      
    tayl_max = 0.0;

    for (i=0;i<n;i++)
      _aa[i][0] = taut*_ab1[i];
    deg++;
 
    while(deg < (_max_deg-2)) {

      mindec(rc,hos_forward(tag,n,n,deg,0,_aa1,_aa,_ab1,_ab));	
      if( rc < 0)
        error(E_INTERN,"Omu_IntOdeTs::solve");

      taut = _tau/(1+deg);      
      tayl_max = 0.0;

      for(i = 0; i < n; i++)
	if(!is_finite(_ab[i][deg-1])) {
	  tayl_finite = 0;
	  printf("Taylor-series contains infinite or NaN elements!\n");
	  error(E_INTERN,"Omu_IntOdeTs::solve");
	}

      /*
      for(i = 0; i < n; i++)
	if(!is_finite(_ab[i][deg-1])) {
	  tayl_finite = 0;
	  break;
	}

      if(!tayl_finite) {
	for(j = 0; j < n; j++)
	  _ab[j][deg-1] = 0.0;
	one_step = false;
	deg--;
	break;
      }
      */
      
      for(i = 0; i < n; i++) {
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

    mindec(rc,hos_forward(tag,n,n,deg,deg+1,_aa1,_aa,_ab1,_ab));

    if( rc < 0)
      error(E_INTERN,"Omu_IntOdeTs::solve");
    
    taut = _tau/(1+deg);
    tayl_max = 0.0;

    for(i = 0; i < n; i++) {
      _aa[i][deg] = taut*_ab[i][deg-1];
      tayl_max = max(tayl_max,fabs(_aa[i][deg]));
    }

    for(i = 0; i < n; i++) {
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
      error(E_SING,"Omu_IntOdeTs::solve");
      printf("Taylor-coeficient: %g \n",tayl_max);
    }

    if(_output > 0) 
      printf("deg: %d time step: %g\n",deg+1,h);

    if(_output > 1) {
      printf("taylor-series:\n");
      for(j = 0; j < n; j++) {
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

      for(i = 0; i < n; i++) {
	for(j = 0; j < n; j++)
	  _ai[i][j] = 0.0;
	_ai[i][i] = 1.0;
      }

      hov_reverse(tag,n,n,deg,n,_ai,_aA,_nz);
      
      accodec(n,_tau,deg,_aA,_aB,_nz);

      if(_sa && (_output > 3)) {
	for(i=0;i<deg+1;i++){
	  printf("\n\t<-- _aB(%d)\n",i);
	  for(j=0;j<n;j++){
	    for(k=0;k<n;k++)
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
	for(j = 0; j < _m; j++)
	  for(k = 0; k < deg+1; k++)
	    _Sxu[i-_nd][j] += _aB[i][j+_nd+_n][k]*hpv[k];

      for(i = _nd; i < _nd+_n ; i++)
	for(j = 0; j < _n; j++)
	  for(k = 0; k < deg+1; k++)
	    _Sxx[i-_nd][j] += _aB[i][j+_nd][k]*hpv[k];

      m_resize(_Sh,_n,_nd);
      m_mlt(_Sxx0,_Sxd,_Sh);
      m_add(_Sxd0,_Sh,_Sxd0);

      m_resize(_Sh,_n,_m);
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

    m_move(_Sxd0,0,0,_n,_nd,Sx,_nd,0);
    m_move(_Sxx0,0,0,_n,_n,Sx,_nd,_nd);
    m_move(_Sxu0,0,0,_n,_m,Su,_nd,0);

      if(_output > 0) {
	m_output(Sx);
	m_output(Su);
      }
  }



  for(i = _nd; i < _nd+_n; i++)
    xt[i] = _a[i];

  if(_output > 0)
    v_output(xt);

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











