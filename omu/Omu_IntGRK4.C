/*
 * Omu_IntGRK4.C --
 *   -- integrate a (stiff) ODE over a stage using a linear-implicit RK method
 *   -- class implementation
 *
 * Reference: FORTRAN code ros4.f:
 *
 *    NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
 *    SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y).
 *    THIS IS AN EMBEDDED ROSENBROCK METHOD OF ORDER (3)4  
 *    (WITH STEP SIZE CONTROL).
 *    C.F. SECTION IV.7
 *
 *    AUTHORS: E. HAIRER AND G. WANNER
 *             UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
 *             CH-1211 GENEVE 24, SWITZERLAND 
 *             E-MAIL:  HAIRER@CGEUGE51.BITNET,  WANNER@CGEUGE51.BITNET
 *    
 *    THIS CODE IS PART OF THE BOOK:
 *        E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
 *        EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
 *        SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
 *        SPRINGER-VERLAG (1990)               
 *     
 *    VERSION OF NOVEMBER 17, 1992
 *
 * E. Arnold   2000-05-25 C++ version
 *             2003-01-03
 *
 */

/*
    Copyright (C) 2000-2003  Eckhard Arnold

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
#include <If_Bool.h>
#include <If_Real.h>
#include <If_Class.h>

#include "Omu_IntGRK4.h"
extern "C" {
#include <meschach/matrix2.h>
}

IF_CLASS_DEFINE("GRK4", Omu_IntGRK4, Omu_Integrator);

//--------------------------------------------------------------------------

Omu_IntGRK4::Omu_IntGRK4()
{

    _y = v_get(1);
    _ynew = v_get(1);
    _fx = v_get(1);
    _dy = v_get(1);
    k1 = v_get(1);
    k1_ori = v_get(1);
    k2     = v_get(1);
    k3     = v_get(1);
    k4     = v_get(1);
    _u     = v_get(1);
    _yjac  = v_get(1);
    _yjacp = v_get(1);
    _yy    = m_get(1, 1);
    _yyn   = m_get(1, 1);
    _yu    = m_get(1, 1);
    _yun   = m_get(1, 1);
    _ppivot= px_get(1);

    _nmax = 100000;
    _hinit = 0.0;
    _hmax = 0.0;
    _uround = MACHEPS;
    _fac1 = 5.0;
    _fac2 = 1.0/6.0;

    _coeffs = 5;   // 1: SHAMPINE
                   // 2: GRK4A OF KAPS-RENTROP 
                   // 3: GRK4T OF KAPS-RENTROP
                   // 4: VAN VELDHUIZEN (GAMMA=1/2)
                   // 5: VAN VELDHUIZEN ("D-STABLE")
                   // 6: AN L-STABLE METHOD

    _sensrk4 = 0;  // sensitivities by 0: IMP, 1: RK4

    _ifList.append(new If_Real("prg_int_hinit", &_hinit));
    _ifList.append(new If_Real("prg_int_hmax", &_hmaxinit));
    _ifList.append(new If_Int("prg_int_coeffs", &_coeffs));
    _ifList.append(new If_Bool("prg_int_sensrk4", &_sensrk4));
    _res_evals = 0;
    _jac_evals = 0;
}

//--------------------------------------------------------------------------

Omu_IntGRK4::~Omu_IntGRK4()
{

    V_FREE(_y);
    V_FREE(_ynew);
    V_FREE(_fx);
    V_FREE(_dy);
    V_FREE(k1);
    V_FREE(k1_ori);
    V_FREE(k2);
    V_FREE(k3);
    V_FREE(k4);
    V_FREE(_u);
    V_FREE(_yjac);
    V_FREE(_yjacp);
    M_FREE(_yy);
    M_FREE(_yyn);
    M_FREE(_yu);
    M_FREE(_yun);
    PX_FREE(_ppivot);

}

//--------------------------------------------------------------------------

void Omu_IntGRK4::init_stage(int k,
			     const Omu_States &x, const Omu_Vector &u,
			     bool sa)
{
    int result = 0;

    Omu_IntODE::init_stage(k, x, u, sa);

    if ( _nv ) {
	fprintf(stderr,  "Omu_IntGRK4::init_stage");
	m_error(E_SIZES, "Omu_IntGRK4::init_stage");
    }
    
    _npar = _nd + _nu;

    resize();

    // check parameters
    if ( _rtol <= 10.0*_uround || _atol <= 0.0 ) {
	fprintf(stderr, "Tolerances are too small: rtol = %g, atol = %g\n", 
		_rtol, _atol);
	result = -1;
    }
    if ( _nmax < 0 ) {
	fprintf(stderr, "Wrong input for nmax = %ld\n", _nmax);
	result = -1;
    }
    if ( _uround <= 1e-35 || _uround >= 1.0 ) {
	fprintf(stderr, "Wrong input for uround = %g\n", _uround);
	result = -1;
    }

    if ( _fac2 < 0.0 || _fac1 < 1.0 || _fac2 > 1.0 ) {
	fprintf(stderr, 
		"Wrong input for step size parameters fac1 = %g, fac2 = %g\n",
		_fac1, _fac2);
	result = -1;
    } 

    // set coefficients according to _coeffs
    if ( !result )
	result = coeffs();

    if ( result < 0 ) {
	fprintf(stderr, "Omu_IntGRK4::init_stage");
	m_error(E_UNKNOWN, "Omu_IntGRK4::init_stage");
    }
}

//--------------------------------------------------------------------------

void Omu_IntGRK4::resize() 
{

    // resize _y: without sensitivity analysis _y->dim=_n
    if ( _sa ) 
	v_resize(_y, _n*(1+_n+_npar));
    else
	v_resize(_y, _n);
    v_resize(_ynew, _n);
    v_resize(_fx, _n);
    v_resize(_dy, _n);
    v_resize(_u, _npar);
    v_resize(_yjac, _n*(1+_n+_npar));
    v_resize(_yjacp, _n*(1+_n+_npar));
    m_resize(_yy, _n, _n);
    m_resize(_yyn, _n, _n);
    m_resize(_yu, _n, _npar);
    m_resize(_yun, _n, _npar);
    px_resize(_ppivot, _n);

    // ensure that k1, ..., k4, y1 are of the correct size 
    v_resize(k1, _n);
    v_resize(k1_ori, _n);
    v_resize(k2, _n);
    v_resize(k3, _n);
    v_resize(k4, _n);
}

//--------------------------------------------------------------------------

void Omu_IntGRK4::ode_solve(double tstart, VECP y, const VECP u, double tend)
{

    int i, result;

    result = 0;
    _x = tstart;
    _xend = tend;

    for ( i = 0; i < (int)_y->dim; i++ )
	_y[i] = y[i];
    _u = v_copy(u, _u);

    // start simulation
    _h = fabs(_hinit);
    if ( _h == 0.0 ) 
	_h = fabs(tend-tstart)/10.0;
    //  _hmax = fabs(_hmaxinit);
    result = simulation();
    if ( result < 0 ) {
	fprintf(stderr, "Omu_IntGRK4::ode_solve");
	m_error(E_UNKNOWN, "Omu_IntGRK4::ode_solve");
    }

    for ( i = 0; i < (int)_y->dim; i++ )
	y[i] = _y[i];
}

//-----------------------------------------------------------------------------
// calculate ODE rhs
void Omu_IntGRK4::sys(double t, const VECP x, VECP xp)
{
    bool sa_old;

    v_zero(xp);
    sa_old = _sa;
    _sa = 0;
    syseq(t, x, _u, xp);
    _sa = sa_old;
}

//-----------------------------------------------------------------------------
// calculate Jacobian of ODE rhs wrt x
// modifies _yjac, _yjacp
void Omu_IntGRK4::sys_jac(double t, const VECP x, VECP xp, MATP fx)
{
    int i, j;
    bool sa_old;

    v_zero(_yjac);
    v_zero(_yjacp);
    for ( i = 0; i < _n; i++ ) {
	_yjac[i] = x[i];
	_yjac[(1+_nd+i)*_n+i] = 1.0;
    }
    sa_old = _sa;
    _sa = 1;
    syseq(t, _yjac, _u, _yjacp);
    _sa = sa_old;
    _jac_evals++;

    for ( i = 0; i < _n; i++ ) {
	xp[i] = _yjacp[i];
	for ( j = 0; j < _n; j++ )
	    fx[i][j] = _yjacp[(1+_nd+j)*_n+i];
    }
}

//-----------------------------------------------------------------------------
// calculate Jacobian of ODE rhs wrt x and u
// modifies _yjac, _yjacp
void Omu_IntGRK4::sys_jac(double t, const VECP x, VECP xp, MATP fx, MATP fu)
{
    int i, j, ii;

    sys_jac( t, x, xp, fx);
    for ( i = 0, ii = 0; i < _n+_npar; i++ ) {
	if ( i < _nd || i >= _nd+_n ) {
	    for ( j = 0; j < _n; j++ ) 
		fu[j][ii] = _yjacp[_n*(1+i)+j];
	    ii++;
	}
    }
}

//-----------------------------------------------------------------------------
// LU factorize I*gamma-fx
// modifies _ppivot
int Omu_IntGRK4::lufac_jac(double gamma, MATP fx)
{
    int i, result = 0;

    sm_mlt(-1.0, fx, fx);
    for ( i = 0; i < _n; i++ ) 
	    fx[i][i] += gamma;
    m_catchall(LUfactor(fx, _ppivot);,
	       result = -1;
	       _nsing++; 
	       if ( _nsing == 5 ) {
		   fprintf(stderr, 
			   "Omu_IntGRK4::lufac_jac singular Jacobian");
		   m_error(E_SING, 
			   "Omu_IntGRK4::lufac_jac singular Jacobian");
	       }
	);
    return result;
}

//-----------------------------------------------------------------------------
// fx*(s+fac*sp)+fu
// modifies _fx, _dy
void Omu_IntGRK4::update_sens(const MATP fx, const VECP s, double fac, 
			      const VECP ds, const MATP fu, VECP sp)
{
    int i, ii, j;

    for ( i = 0, ii = 0; i < _n+_npar; i++ ) {
	for ( j = 0; j < _n; j++ )
	    _fx[j] = s[_n*(1+i)+j]+fac*ds[_n*(1+i)+j];
	mv_mlt(fx, _fx, _dy);
	if ( i < _nd || i >= _nd+_n ) {
	    for ( j = 0; j < _n; j++ )
		_dy[j] += fu[j][ii];
	    ii++;
	}
	for ( j = 0; j < _n; j++ )
	    sp[_n*(1+i)+j] = _dy[j];
    }
}

//-----------------------------------------------------------------------------
// LU solve Ax=b
void Omu_IntGRK4::lusolve_jac(MATP A, VECP b, VECP x)
{
    LUsolve(A, _ppivot, b, x);
}

//-----------------------------------------------------------------------------

int Omu_IntGRK4::simulation()
{

    int i, ii, j;
    int reject = 0, reject2 = 0;
    double facrej = 0.1;
    double h, hnew, xdelt, err;
    bool last = 0, jac_ok = 0;
    double hc21, hc31, hc32, hc41, hc42, hc43, fac, hd1, hd2, hd3, hd4;

    _nstep = _naccpt = _nrejct = _nsing = 0;
  
    h = _h;
    if ( _hmax == 0.0 ) 
	_hmax = fabs(_xend-_x);
    if ( (_xend > _x && h <= 0) || (_xend < _x && h > 0.0) ) 
	_h = h = -h;
  
    _posneg = (_xend-_x)/fabs(_xend-_x);
  
    _xold = _x; 
  
    // basic integration step 
    while ( ! last ) {

//	printf(" x=%g, h=%g, hnew=%g, xend=%g\n", _x, h, hnew, _xend);

	if ( _nstep > _nmax ) {
	    fprintf(stderr, 
		    "GRK4: more than nmax = %ld steps are needed.\n", _nmax); 
	    return -2;
	} 
	if ( _x+0.1*h == _x || fabs(h) <= _uround ) {
	    fprintf(stderr, "GRK4: stepsize to small, h = %g\n", h);
	    return -3;
	}     
	last = ((_x+1.01*h-_xend)*_posneg > 0.0); 
	if ( last ) 
	    _h = h = _xend-_x;
	_nstep++;

	// coefficients
	hc21 = C21/h;
	hc31 = C31/h;
	hc32 = C32/h;
	hc41 = C41/h;
	hc42 = C42/h;
	hc43 = C43/h;
	fac = 1.0/(h*GAMMA);
	hd1 = h*D1;
	hd2 = h*D2;
	hd3 = h*D3;
	hd4 = h*D4;

	if ( !reject && !jac_ok ) {
	    // calculate Jacobian
	    sys_jac(_x, _y, k1, _yy, _yu);
	    v_copy(k1, k1_ori);
	    m_copy(_yy, _yyn);
	    m_copy(_yu, _yun);
	} else {
	    // restore Jacobian
	    v_copy(k1_ori, k1);
	    m_copy(_yyn, _yy);
	    m_copy(_yun, _yu);
	    jac_ok = 0;
	}

	// matrix factorization I/(h*GAMMA)-_yy
	if ( !lufac_jac(fac, _yy) )
	    _nsing = 0;
	else {
	    h *= 0.5;
	    continue;
	}

	// derivative w.r.t. independent variable by finite differences
	xdelt = sqrt(_uround*max(1.0e-5, fabs(_x)));
	sys(_x+xdelt, _y, _fx);
	v_sub(_fx, k1, _fx);
	sv_mlt(1.0/xdelt, _fx, _fx);

	// k1
	v_mltadd(k1, _fx, hd1, k1);
  	lusolve_jac(_yy, k1, k1);

	// k2
	for ( i = 0; i < _n; i++ )
	    _ynew[i] = _y[i]+A21*k1[i];
	sys(_x+C2*h, _ynew, _dy);
	v_linlist(k2, _dy, 1.0, _fx, hd2, k1, hc21, NULL);
	lusolve_jac(_yy, k2, k2);

	// k3
	for ( i = 0; i < _n; i++ )
	    _ynew[i] = _y[i]+A31*k1[i]+A32*k2[i];
  	sys(_x+C3*h, _ynew, _dy);
	v_linlist(k3, _dy, 1.0, _fx, hd3, k1, hc31, k2, hc32, NULL);
	lusolve_jac(_yy, k3, k3);

	// k4
	v_linlist(k4, _dy, 1.0, _fx, hd4, k1, hc41, k2, hc42, k3, hc43, NULL);
	lusolve_jac(_yy, k4, k4);

	// _ynew
	for ( i = 0; i < _n; i++ )
	    _ynew[i] = _y[i]+B1*k1[i]+B2*k2[i]+B3*k3[i]+B4*k4[i];

	// error estimation
	for ( i = 0, err = 0.0; i < _n; i++ ) 
	    err += square((E1*k1[i]+E2*k2[i]+E3*k3[i]+E4*k4[i])/
			  (_atol+_rtol*max(fabs(_y[i]), fabs(_ynew[i]))));
	err = sqrt(err/(double) _n);
	hnew = h/max(_fac2, min(_fac1, pow(err, 0.25)/0.9));

	if ( err <= 1.0 ) {
	    // step accepted
	    for ( i = 0; i < _n; i++ ) {
		_dy[i] = _y[i]; // old _y for derivative
		_y[i] = _ynew[i];
	    }
	    _xold = _x;
	    _x += h;
	    hnew = _posneg*min(fabs(hnew), _hmax);
	    if ( reject )
		hnew = _posneg*min(fabs(hnew), fabs(h));
	    reject = reject2 = 0;
	    h = hnew; 
	} else {
	    // step rejected
	    if ( reject2 )
		hnew = h*facrej;
	    if ( reject )
		reject2 = 1;
	    reject = 1;
	    h = hnew;
	}
	
	// calculation of sensitivities
	if ( _sa && !reject ) {
	    if ( _sensrk4 ) {
		// calculation of sensitivities by RK4
		// Hermite interpolation 
		// y(t+h/2)=(y(t)+y(t+h))/2+h/8*(\dot y(t)-\dot y(t+h))
		sys_jac(_x, _y, k1, _yy, _yu);
		v_copy(k1, k1_ori);
		for ( i = 0; i < _n; i++ )
		    _ynew[i] = (_y[i]+_dy[i])/2.0
			+(_x-_xold)/8.0*(k1_ori[i]-k1[i]);
		v_resize(k1, _y->dim);
		v_resize(k2, _y->dim);
		v_resize(k3, _y->dim);
		v_resize(k4, _y->dim);

		update_sens(_yyn, _y, 0.0, _y, _yun, k1);
		sys_jac(0.5*(_x+_xold), _ynew, _dy, _yyn, _yun);
		update_sens(_yyn, _y, (_x-_xold)/2.0, k1, _yun, k2);
		update_sens(_yyn, _y, (_x-_xold)/2.0, k2, _yun, k3);
		update_sens(_yy, _y, (_x-_xold), k3, _yu, k4);
		for ( i = _n; i < (int) _y->dim; i++ )
		    _y[i] += (_x-_xold)/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);

		v_resize(k1, _n);
		v_resize(k2, _n);
		v_resize(k3, _n);
		v_resize(k4, _n);

		m_copy(_yy, _yyn);
		m_copy(_yu, _yun);
		jac_ok = 1;

	    } else {
		// calculation of sensitivities by IMP
		// calculate ODE r.h.s. and Jacobians _yy and _yu
		for ( i = 0; i < _n; i++ ) 
		    _ynew[i] = (_y[i]+_dy[i])/2.0;
		sys_jac(0.5*(_x+_xold), _ynew, k1, _yy, _yu);
		m_copy(_yy, _yyn);
		
		// factorize 2/dt*I-_yy
		lufac_jac(2.0/(_x-_xold), _yyn);
		
		// derivatives
		for ( i = 0, ii = 0; i < _n+_npar; i++ ) {
		    for ( j = 0; j < _n; j++ )
			k1[j] = _y[_n*(1+i)+j];
		    mv_mlt(_yy, k1, k2);
		    if ( i < _nd || i >= _nd+_n ) {
			for ( j = 0; j < _n; j++ )
			    k2[j] += _yu[j][ii];
			ii++;
		    }
		    lusolve_jac(_yyn, k2, k2);
		    for ( j = 0; j < _n; j++ )
			_y[_n*(1+i)+j] = _y[_n*(1+i)+j]+2.0*k2[j];
		}
	    }
	}
    }
    return 1;
}

//--------------------------------------------------------------------------
// coefficients for method GRK4 according to _coeffs
int Omu_IntGRK4::coeffs() 
{
    int result = 0;

    switch ( _coeffs ) {
	case 1: // METHOD OF SHAMPINE
	    A21 = 2.0E+00;
	    A31 = 48.0E+00/25.0E+00;
	    A32 = 6.0E+00/25.0E+00;
	    C21 = -8.0E+00;
	    C31 = 372.0E+00/25.0E+00;
	    C32 = 12.0E+00/5.0E+00;
	    C41 = -112.0E+00/125.0E+00;
	    C42 = -54.0E+00/125.0E+00;
	    C43 = -2.0E+00/5.0E+00;
	    B1 = 19.0E+00/9.0E+00;
	    B2 = 1.0E+00/2.0E+00;
	    B3 = 25.0E+00/108.0E+00;
	    B4 = 125.0E+00/108.0E+00;
	    E1 = 17.0E+00/54.0E+00;
	    E2 = 7.0E+00/36.0E+00;
	    E3 = 0.0E+00;
	    E4 = 125.0E+00/108.0E+00;
	    GAMMA = 0.5E+00;
	    C2 =  0.1000000000000000E+01;
	    C3 =  0.6000000000000000E+00;
	    D1 =  0.5000000000000000E+00;
	    D2 = -0.1500000000000000E+01;
	    D3 =  0.2420000000000000E+01;
	    D4 =  0.1160000000000000E+00;
	    break;
	case 2: // method GRK4A of KAPS-RENTROP
	    A21= 0.1108860759493671E+01;
	    A31= 0.2377085261983360E+01;
	    A32= 0.1850114988899692E+00;
	    C21=-0.4920188402397641E+01;
	    C31= 0.1055588686048583E+01;
	    C32= 0.3351817267668938E+01;
	    C41= 0.3846869007049313E+01;
	    C42= 0.3427109241268180E+01;
	    C43=-0.2162408848753263E+01;
	    B1= 0.1845683240405840E+01;
	    B2= 0.1369796894360503E+00;
	    B3= 0.7129097783291559E+00;
	    B4= 0.6329113924050632E+00;
	    E1= 0.4831870177201765E-01;
	    E2=-0.6471108651049505E+00;
	    E3= 0.2186876660500240E+00;
	    E4=-0.6329113924050632E+00;
	    GAMMA= 0.3950000000000000E+00;
	    C2= 0.4380000000000000E+00;
	    C3= 0.8700000000000000E+00;
	    D1= 0.3950000000000000E+00;
	    D2=-0.3726723954840920E+00;
	    D3= 0.6629196544571492E-01;
	    D4= 0.4340946962568634E+00;
	    break;
	case 3: // method GRK4T of KAPS-RENTROP
	    A21= 0.2000000000000000E+01;
	    A31= 0.4524708207373116E+01;
	    A32= 0.4163528788597648E+01;
	    C21=-0.5071675338776316E+01;
	    C31= 0.6020152728650786E+01;
	    C32= 0.1597506846727117E+00;
	    C41=-0.1856343618686113E+01;
	    C42=-0.8505380858179826E+01;
	    C43=-0.2084075136023187E+01;
	    B1= 0.3957503746640777E+01;
	    B2= 0.4624892388363313E+01;
	    B3= 0.6174772638750108E+00;
	    B4= 0.1282612945269037E+01;
	    E1= 0.2302155402932996E+01;
	    E2= 0.3073634485392623E+01;
	    E3=-0.8732808018045032E+00;
	    E4=-0.1282612945269037E+01;
	    GAMMA= 0.2310000000000000E+00;
	    C2= 0.4620000000000000E+00;
	    C3= 0.8802083333333334E+00;
	    D1= 0.2310000000000000E+00;
	    D2=-0.3962966775244303E-01;
	    D3= 0.5507789395789127E+00;
	    D4=-0.5535098457052764E-01;
	    break;
	case 4: // METHOD OF VAN VELDHUIZEN (GAMMA=1/2)
	    A21= 0.2000000000000000E+01;
	    A31= 0.1750000000000000E+01;
	    A32= 0.2500000000000000E+00;
	    C21=-0.8000000000000000E+01;
	    C31=-0.8000000000000000E+01;
	    C32=-0.1000000000000000E+01;
	    C41= 0.5000000000000000E+00;
	    C42=-0.5000000000000000E+00;
	    C43= 0.2000000000000000E+01;
	    B1= 0.1333333333333333E+01;
	    B2= 0.6666666666666667E+00;
	    B3=-0.1333333333333333E+01;
	    B4= 0.1333333333333333E+01;
	    E1=-0.3333333333333333E+00;
	    E2=-0.3333333333333333E+00;
	    E3=-0.0000000000000000E+00;
	    E4=-0.1333333333333333E+01;
	    GAMMA= 0.5000000000000000E+00;
	    C2= 0.1000000000000000E+01;
	    C3= 0.5000000000000000E+00;
	    D1= 0.5000000000000000E+00;
	    D2=-0.1500000000000000E+01;
	    D3=-0.7500000000000000E+00;
	    D4= 0.2500000000000000E+00;
	    break;
	case 5: // METHOD OF VAN VELDHUIZEN ("D-STABLE")
	    A21= 0.2000000000000000E+01;
	    A31= 0.4812234362695436E+01;
	    A32= 0.4578146956747842E+01;
	    C21=-0.5333333333333331E+01;
	    C31= 0.6100529678848254E+01;
	    C32= 0.1804736797378427E+01;
	    C41=-0.2540515456634749E+01;
	    C42=-0.9443746328915205E+01;
	    C43=-0.1988471753215993E+01;
	    B1= 0.4289339254654537E+01;
	    B2= 0.5036098482851414E+01;
	    B3= 0.6085736420673917E+00;
	    B4= 0.1355958941201148E+01;
	    E1= 0.2175672787531755E+01;
	    E2= 0.2950911222575741E+01;
	    E3=-0.7859744544887430E+00;
	    E4=-0.1355958941201148E+01;
	    GAMMA= 0.2257081148225682E+00;
	    C2= 0.4514162296451364E+00;
	    C3= 0.8755928946018455E+00;
	    D1= 0.2257081148225682E+00;
	    D2=-0.4599403502680582E-01;
	    D3= 0.5177590504944076E+00;
	    D4=-0.3805623938054428E-01;
	    break;
	case 6:  // AN L-STABLE METHOD
	    A21= 0.2000000000000000E+01;
	    A31= 0.1867943637803922E+01;
	    A32= 0.2344449711399156E+00;
	    C21=-0.7137615036412310E+01;
	    C31= 0.2580708087951457E+01;
	    C32= 0.6515950076447975E+00;
	    C41=-0.2137148994382534E+01;
	    C42=-0.3214669691237626E+00;
	    C43=-0.6949742501781779E+00;
	    B1= 0.2255570073418735E+01;
	    B2= 0.2870493262186792E+00;
	    B3= 0.4353179431840180E+00;
	    B4= 0.1093502252409163E+01;
	    E1=-0.2815431932141155E+00;
	    E2=-0.7276199124938920E-01;
	    E3=-0.1082196201495311E+00;
	    E4=-0.1093502252409163E+01;
	    GAMMA= 0.5728200000000000E+00;
	    C2= 0.1145640000000000E+01;
	    C3= 0.6552168638155900E+00;
	    D1= 0.5728200000000000E+00;
	    D2=-0.1769193891319233E+01;
	    D3= 0.7592633437920482E+00;
	    D4=-0.1049021087100450E+00;
	    break;
	default:
	    fprintf(stderr, "Wrong input for coeffs = %d\n", _coeffs);
	    result = -1;
    }
    return result;
}

//=============================================================================
