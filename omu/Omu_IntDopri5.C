/*
 * Omu_IntDopri5.C --
 *   -- class implementation
 *
 * E. Arnold   1999-01-23 C++ version
 *                        continuous output not used!
 *             1999-04-23 egcs, g77, libg2c 
 *             2000-03-29 Real -> double
 *                        error estimation using first _n components only
 *                        _rtol, _atol -> Omu_Integrator
 *
 */

/*
    Copyright (C) 1997--1999  Eckhard Arnold

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

#include <adutils.h>

#include <If_Int.h>
#include <If_Real.h>
#include <If_Class.h>

#include "Omu_IntDopri5.h"

IF_CLASS_DEFINE("Dopri5", Omu_IntDopri5, Omu_Integrator);

//--------------------------------------------------------------------------

Omu_IntDopri5::Omu_IntDopri5()
{

  _y = v_get(1);
  k1 = v_get(1);
  k2 = v_get(1);
  k3 = v_get(1);
  k4 = v_get(1);
  k5 = v_get(1);
  k6 = v_get(1);
  y1 = v_get(1);
  ysti = v_get(1);
  cont = v_get(1);
  cont = v_resize(cont, 0);
  cont1 = v_get(1);
  cont2 = v_get(1);
  cont3 = v_get(1);
  cont4 = v_get(1);
  cont5 = v_get(1);
  _icont = px_get(1);
  _icont->pe[0] = 0;
  _icont = px_resize(_icont, 0);
  _u = v_get(1);

  //  _rtol = _atol = 1.0e-6;
  //  _rtol = 1e-6;                         // see Omu_IntRKsuite.C
  //  _atol = 1e-12;                        // 

  _nmax = 100000;
  _nstiff = 1000;
  _hinit = 0.0;
  _hmaxinit = 0.0;
  _uround = MACHEPS;
  _safe = 0.9;
  _beta = 0.04;
  _fac1 = 0.2;
  _fac2 = 10.0;

  //  _ifList.append(new If_Real("prg_int_rtol", &_rtol));
  //  _ifList.append(new If_Real("prg_int_atol", &_atol));
  _ifList.append(new If_Real("prg_int_hinit", &_hinit));
  _ifList.append(new If_Real("prg_int_hmax", &_hmaxinit));
  // additional interfaces for _nmax, _nstiff, _uround, _safe, _beta,
  //                           _fac1, _fac2, _nfcn, _nstep, _naccpt, _nrejct

}

//--------------------------------------------------------------------------

Omu_IntDopri5::~Omu_IntDopri5()
{

  V_FREE(_y);
  V_FREE(k1);
  V_FREE(k2);
  V_FREE(k3);
  V_FREE(k4);
  V_FREE(k5);
  V_FREE(k6);
  V_FREE(y1);
  V_FREE(ysti);
  V_FREE(cont);
  V_FREE(cont1);
  V_FREE(cont2);
  V_FREE(cont3);
  V_FREE(cont4);
  V_FREE(cont5);
  PX_FREE(_icont);
  V_FREE(_u);

}

//--------------------------------------------------------------------------

void Omu_IntDopri5::ode_solve(double tstart, VECP y, const VECP u, double tend)
{

  int i, res;

  res = 0;
  _x = tstart;
  _xend = tend;

  // resize _y: without sensitivity analysis _y->dim=_n
  if ( _sa ) 
    v_resize(_y, y->dim);
  else
    v_resize(_y, _n);
  for ( i=0; i<(int)_y->dim; i++ )
    _y[i] = y[i];
  v_resize(_u, u->dim);
  _u = v_copy(u, _u);

  // check parameters
  if ( _rtol <= 0.0 || _atol <= 0.0 ) {
    fprintf(stderr, "Curious input for tolerances rtol = %g, atol = %g\n", 
	   _rtol, _atol);
    res = -1;
  }
  if ( _nmax < 0 ) {
    fprintf(stderr, "Wrong input nmax = %ld\n", _nmax);
    res = -1;
  }
  if ( _nstiff < 0 )
    _nstiff = _nmax+10;
  if ( _uround <= 1e-35 || _uround >= 1.0 ) {
    fprintf(stderr, "Which machine do you have? uround = %g\n", _uround);
    res = -1;
  }
  if ( _safe <= 1e-4 || _safe >= 1.0 ) {
    fprintf(stderr, "Curious input for safety factor safe = %g\n", _safe);
    res = -1;
  }
  if ( _beta < 0.0 ) 
    _beta = 0.0;
  else if ( _beta > 0.2 ) {
    fprintf(stderr, "Curious input for stabilization parameter beta = %g\n", 
	    _beta);
    res = -1;
  }
  if ( _fac1 < 0.0 || _fac1 > 1.0 || _fac2 < 1.0 ) {
    fprintf(stderr, 
	    "Curious input for step size parameters fac1 = %g, fac2 = %g\n",
	    _fac1, _fac2);
    res = -1;
  } else {
    _facc1 = 1.0/_fac1; 
    _facc2 = 1.0/_fac2; 
  }

  // start simulation
  _h = fabs(_hinit);
  _hmax = fabs(_hmaxinit);
  if ( res >= 0 )
    res = simulation();
  if ( res < 0 )
    error(E_UNKNOWN, "Omu_IntDopri5::ode_solve");

  for ( i=0; i<(int)_y->dim; i++ )
    y[i] = _y[i];

}

//-----------------------------------------------------------------------------

void Omu_IntDopri5::set_cont(PERM *icont)
{

  _icont = px_resize(_icont, (int) icont->size);
  _icont = px_copy(icont, _icont);

}

//-----------------------------------------------------------------------------
// Computation of an initial stepsize guess 
double Omu_IntDopri5::hinit()
{
  // Compute a first guess for explicit EULER as
  // h = 0.01*norm(y0)/norm(f0)
  // the increment for explicit EULER is small compared to the solution 
  // uses _x, _y, k1, k2, k3, _hmax, _posneg, _atol, _rtol 
  double sk, dnf, dny, der2, der12, h, h1;
  int i;
  //  for ((dnf=dny=0.0,i=0); i<(int)_y->dim; i++) {
  for ( (dnf=dny=0.0,i=0); i<_n; i++ ) {
    sk = _atol+_rtol*fabs(_y[i]);
    dnf += square(k1[i]/sk); 
    dny += square(_y[i]/sk);
  }
  if ( dnf <= 1e-10 || dny <= 1e-10 )
    h = 1e-6;
  else
    h = sqrt(dny/dnf)*0.01;
  h = min(h,_hmax)*_posneg;
  // perform an explicit EULER step 
  for (i=0; i<(int)_y->dim; i++) 
    k3[i] = _y[i]+h*k1[i]; 
  f(_x+h, k3, k2);
  _nfcn++; 
  // estimate the 2nd derivative of the solution 
  //  for ((der2=0.0,i=0); i<(int)_y->dim; i++) 
  for ( (der2=0.0,i=0); i<_n; i++ ) 
    der2 += square((k2[i]-k1[i])/(_atol+_rtol*fabs(_y[i])));
  der2 = sqrt(der2)/h;
  // step-size is computed such that
  // h^iord*max(norm(f0),norm(der2)) = 0.01 
  der12 = max(fabs(der2),sqrt(dnf));  
  if ( der12 <= 1e-15 )
    h1 = max(1e-6,fabs(h)*1e-3);
  else
    h1 = pow(0.01/der12,1.0/iord);
  return min(100*fabs(h),min(h1,_hmax))*_posneg;       

} 

//-----------------------------------------------------------------------------
// Continuous output in connection with output function 
VECP Omu_IntDopri5::contd5(const double x, VECP cont)
{

  int i;
  double theta = (x-_xold)/_h, theta1 = 1.0-theta;
  for ( i=0; i<(int)_icont->size; i++ ) 
    cont[i] = cont1[i]+
      theta*(cont2[i]+theta1*(cont3[i]+theta*(cont4[i]+theta1*cont5[i]))); 
  return cont;

}

//-----------------------------------------------------------------------------

int Omu_IntDopri5::simulation()
{

  int i, j, last = 0, reject = 0, nonsti = 0, iasti = 0;
  double xph, err, facold = 1.0e-4, hlamb = 0.0, expo1,
    fac, fac11, h, hnew, stnum, stden;
  char format979[100] = "  Exit of DOPRI at x = %g, "; 

  // ensure that k1, ..., k6, ysti are of the correct size 
  v_resize(k1, _y->dim);
  v_resize(k2, _y->dim);
  v_resize(k3, _y->dim);
  v_resize(k4, _y->dim);
  v_resize(k5, _y->dim);
  v_resize(k6, _y->dim);
  v_resize(y1, _y->dim);
  v_resize(ysti, _y->dim);
  
  // ensure that cont1, ..., cont are of the correct size 
  v_resize(cont1, _icont->size);
  v_resize(cont2, _icont->size);
  v_resize(cont3, _icont->size);
  v_resize(cont4, _icont->size);
  v_resize(cont5, _icont->size);
  
  _nfcn = _nstep = _naccpt = _nrejct = 0;
  expo1 = 0.2-0.04*0.75; 
  
  h = _h;
  if ( _hmax == 0.0 ) 
    _hmax = fabs(_xend-_x);
  if ( (_xend > _x && h <= 0) || (_xend < _x && h > 0.0) ) 
    _h = h = -h;
  
  _posneg = (_xend-_x)/fabs(_xend-_x);
  
  f(_x, _y, k1);
  _nfcn++;
  if (h == 0.0) 
    _h = h = hinit();
  _xold = _x; 
  
  for ( j=0; j<(int)_icont->size; j++ ) 
    cont[j] = _y[_icont->pe[j]];
  if ( fout(_naccpt+1, _x, cont) ) {
    fprintf(stderr, format979, _x);
    return 2;
  }
  
  // basic integration step 
  while ( ! last ) {
    if ( _nstep > _nmax ) {
      fprintf(stderr, 
	      strcat(format979,"more than nmax = %d steps are needed.\n"), 
	      _x, _nmax); 
      return -2;
    } 
    if ( 0.1*fabs(h) <= _uround*fabs(_x) ) {
      fprintf(stderr, strcat(format979,"stepsize to small, h = %g\n"), _x, h);
      return -3;
    }     
    last = ((_x+1.01*h-_xend)*_posneg > 0.0); 
    if ( last ) 
      _h = h = _xend-_x;
    _nstep++;
    
    // the first 6 stages 
    for ( i=0; i<(int)_y->dim; i++ ) 
      y1[i] = _y[i]+h*a21*k1[i]; 
    f(_x+c2*h, y1, k2);       
    for ( i=0; i<(int)_y->dim; i++ ) 
      y1[i] = _y[i]+h*(a31*k1[i]+a32*k2[i]); 
    f(_x+c3*h, y1, k3);    
    for ( i=0; i<(int)_y->dim; i++ ) 
      y1[i] = _y[i]+h*(a41*k1[i]+a42*k2[i]+a43*k3[i]); 
    f(_x+c4*h, y1, k4);    
    for ( i=0; i<(int)_y->dim; i++ ) 
      y1[i] = _y[i]+h*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]); 
    f(_x+c5*h, y1, k5);    
    for ( i=0; i<(int)_y->dim; i++ ) 
      ysti[i] = _y[i]+h*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]); 
    xph = _x+h;
    f(xph, ysti, k6);
    for ( i=0; i<(int)_y->dim; i++ ) 
      y1[i] = _y[i]+h*(a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
    f(xph, y1, k2);
    _nfcn += 6;
    
    for ( j=0; j<(int)_icont->size; j++ ) {
      i = _icont->pe[j];
      cont5[i] = h*(d1*k1[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+d7*k2[i]);
    }     
    for ( i=0; i<(int)_y->dim; i++ ) 
      k4[i] = h*(e1*k1[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*k2[i]); 
    
    // error estimation 
    for ( (err=0.0,i=0); i<(!_serr?_n:(int)_y->dim); i++ ) 
      err += square(k4[i]/(_atol+_rtol*max(fabs(_y[i]),fabs(y1[i])))); 
    err = sqrt(err/(!_serr?_n:(int)_y->dim));
      
    // computation of hnew 
    fac11 = pow(err,expo1);
    // LUND stabilization 
    fac = fac11/pow(facold,_beta);
    // we require fac1 <= hnew/h <= fac2 
    fac = max(_facc2,min(_facc1,fac/_safe));
    hnew = h/fac;

    if ( err<=1.0 ) { // step is accepted 
      facold = max(err,1e-4);
      _naccpt++;
      // stiffness detection 
      if ( (_naccpt % _nstiff) == 0 || iasti > 0 ) {
	for ( (stnum=stden=0.0, i=0); i<(!_serr?_n:(int)_y->dim); i++ ) {
	  stnum += square(k2[i]-k6[i]);
	  stden += square(y1[i]-ysti[i]);
	}
	if ( stden > 0.0 ) 
	  hlamb = h*sqrt(stnum/stden);
	if ( hlamb > 3.25 ) {
	  nonsti = 0;
	  iasti++;
	  if (iasti == 15) {
	    fprintf(stderr,
		    "The problem seems to become stiff at x = %g\n", _x);
	    return -4;
	  }
	} else {
	  nonsti++;
	  if ( nonsti == 6 ) 
	    iasti = 0;
	}
      } 
      for ( j=0; j<(int)_icont->size; j++ ) {
	i = _icont->pe[j];
	cont1[i] = _y[i];
	cont2[i] = y1[i]-_y[i];
	cont3[i] = h*k1[i]-cont2[i];
	cont4[i] = -h*k2[i]+cont2[i]-cont3[i];
      }     
      v_copy(k2, k1); 
      v_copy(y1, _y);
      _xold = _x;
      _x = xph;
      if ( fout(_naccpt+1, _x, cont) ) {
	fprintf(stderr, format979, _x);
	return 2;
      }
      if ( last ) {
	_h = h = hnew;
	return 1;
      }
      if ( fabs(hnew) > _hmax ) 
	hnew = _posneg*_hmax;
      if ( reject ) 
	hnew = _posneg*min(fabs(hnew),fabs(h));
      reject = 0;
    } else {   // step is rejected 
      hnew = h/min(_facc1,fac11/_safe);
      reject = 1;
      if ( _naccpt >= 1 ) 
	_nrejct++;
      last = 0;
    }
    _h = h = hnew; 
  }
  return 1;

}

//=============================================================================
