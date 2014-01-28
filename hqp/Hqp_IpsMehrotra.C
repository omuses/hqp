/*
 * Hqp_IpsMehrotra.C -- class definition
 *
 * rf, 5/28/94
 *
 * rf, 11/5/95: inserted _max_warm_iters
 * rf, 2/12/97: treat singular matrix errors
 * rf, 4/20/97: treat numerical overflow for _qp->x and _gap
 * ea, 5/9/97
 *   - Mehrotra's primal-dual predictor-corrector method for QP problems.
 *     Mehrotra, S.:  On the implementation of a primal-dual interior point
 *                    method. SIAM J. Optimization 2(1992)4, 575-601.
 *     Wright, S. J.: Primal-dual interior-point methods.
 *                    SIAM, Philadelphia, 1997.
 *     Czyzyk, J. and Mehrotra, S. and Wright, S. J.: PCx User Guide.
 *                    Technical Report OTC 96/01, 
 *                    Optimization Technology Center, 1997.  
 *     Gertz, E. M. and Wright, S. J.: 
 *                    Object-oriented Software for Quadratic Programming
 *     Salahi M., Peng, J. and Terlaky, T.:
 *                    On Mehrotra-Type Predictor-Corrector Algorithms.
 *                    McMaster University, Hamilton, Canada, 2005.
 *
 * Problem notation:
 *     min{ 0.5x'Qx + c'x | Ax+b = 0, Cx+d >= 0 }
 *
 * _matrix->factor, _matrix->solve solves:
 *     | -Q  A'  C'  0 | |dx|   |r1|
 *     |  A  0   0   0 | |dy|   |r2|
 *     |  C  0   0  -I | |dz| = |r3|
 *     |  0  0   W   Z | |dw|   |r4|
 *
 * 08/25/98 - _logging
 *            IpSolver -> IpsMehrotra
 * 2006-04-21 remove bug in corrector step calculation
 *            reformulate corrector step calculation
 *            mu -> mu/_m
 *            remove bug in Mehrotra's stepsize heuristic 
 *            alfa -> -alfa
 *            modify termination conditions
 *            Terlaky's modifications
 *
 */

/*
    Copyright (C) 1994--2014  Eckhard Arnold and Ruediger Franke

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

#include <If_Int.h>
#include <If_Real.h>
#include <If_Method.h>
#include <If_Module.h>

#include "Hqp_Program.h"
#include "Hqp_IpRedSpBKP.h"
#include "Hqp_IpsMehrotra.h"

IF_CLASS_DEFINE("Mehrotra", Hqp_IpsMehrotra, Hqp_Solver);

typedef If_Method<Hqp_IpsMehrotra> If_Cmd;

//--------------------------------------------------------------------------
Hqp_IpsMehrotra::Hqp_IpsMehrotra()
{
  _n = _me = _m = 0;

  _w = VNULL;
  _r1 = VNULL;
  _r2 = VNULL;
  _r3 = VNULL;
  _r4 = VNULL;
  _dx = VNULL;
  _dy = VNULL;
  _dz = VNULL;
  _dw = VNULL;
  _dxa = VNULL;
  _dya = VNULL;
  _dza = VNULL;
  _dwa = VNULL;
  _z_hot = VNULL;
  _w_hot = VNULL;
  _phimin = VNULL;
  _d1 = VNULL;
  _d2 = VNULL;

  _matrix = new Hqp_IpRedSpBKP;

  _alpha = 1.0;
  _gammaf = 0.01;
  //  _mu0 = 0.0;      // not used
  //  _beta = 0.995;   // not used
  _fail_iters = 0;
  _max_warm_iters = 25;
  _init_method = 0;
  _logging = 0;

  _ifList.append(new If_Real("qp_gap", &_gap));
  _ifList.append(new If_Real("qp_alpha", &_alpha));
  //  _ifList.append(new If_Real("qp_beta", &_beta));
  //  _ifList.append(new If_Real("qp_rhomin", &_rhomin));
  //  _ifList.append(new If_Real("qp_mu0", &_mu0));
  //  _ifList.append(new If_Real("qp_Ltilde", &_Ltilde));
  _ifList.append(new If_Int("qp_fail_iters", &_fail_iters));
  _ifList.append(new If_Int("qp_max_warm_iters", &_max_warm_iters));
  _ifList.append(new If_Int("qp_logging", &_logging));
  _ifList.append(new If_Int("qp_init_method", &_init_method));
  _ifList.append(new IF_MODULE("qp_mat_solver", &_matrix, Hqp_IpMatrix));

  _ifList.append(new If_Cmd("qp_init", &Hqp_IpsMehrotra::init, this));
  _ifList.append(new If_Cmd("qp_update", &Hqp_IpsMehrotra::update, this));
  _ifList.append(new If_Cmd("qp_cold_start", &Hqp_IpsMehrotra::cold_start, this));
  _ifList.append(new If_Cmd("qp_hot_start", &Hqp_IpsMehrotra::hot_start, this));
  _ifList.append(new If_Cmd("qp_step", &Hqp_IpsMehrotra::step, this));
  _ifList.append(new If_Cmd("qp_solve", &Hqp_IpsMehrotra::solve, this));
}

//--------------------------------------------------------------------------
Hqp_IpsMehrotra::~Hqp_IpsMehrotra()
{
  v_free(_w);
  v_free(_r1);
  v_free(_r2);
  v_free(_r3);
  v_free(_r4);
  v_free(_dx);
  v_free(_dy);
  v_free(_dz);
  v_free(_dw);
  v_free(_dxa);
  v_free(_dya);
  v_free(_dza);
  v_free(_dwa);
  v_free(_z_hot);
  v_free(_w_hot);
  v_free(_phimin);
  v_free(_d1);
  v_free(_d2);
  delete _matrix;
}

//--------------------------------------------------------------------------
void Hqp_IpsMehrotra::init()
{
  assert(_qp != NULL);

  // allocate matrices and vectors

  _n = _qp->Q->n;
  _me = _qp->A->m;
  _m = _qp->C->m;

  _y = v_resize(_y, _me);
  _z = v_resize(_z, _m);
  _w = v_resize(_w, _m);
  _r1 = v_resize(_r1, _n);
  _r2 = v_resize(_r2, _me);
  _r3 = v_resize(_r3, _m);
  _r4 = v_resize(_r4, _m);
  _dx = v_resize(_dx, _n);
  _dy = v_resize(_dy, _me);
  _dz = v_resize(_dz, _m);
  _dw = v_resize(_dw, _m);
  _dxa = v_resize(_dxa, _n);
  _dya = v_resize(_dya, _me);
  _dza = v_resize(_dza, _m);
  _dwa = v_resize(_dwa, _m);
  _z_hot = v_resize(_z_hot, _m);
  _w_hot = v_resize(_w_hot, _m);
  _d1 = v_resize(_d1, _m);
  _d2 = v_resize(_d2, _m);

  // fill up internal data

  _matrix->init(_qp);

  if ( _logging )
    printf("\nHqp_IpsMehrotra::init\n");
}

//--------------------------------------------------------------------------
void Hqp_IpsMehrotra::update()
{
  _matrix->update(_qp);
  _phimin = v_resize(_phimin, _max_iters+1);

  if ( _logging )
    printf("\nHqp_IpsMehrotra::update\n");
}

//--------------------------------------------------------------------------
void Hqp_IpsMehrotra::cold_start()
{
  int  i, izmin, iwmin, code;
  Real delz, delw, residuum;

  _iter = 0;
  v_zero(_phimin);
  _alpha = 1.0;

  _hot_started = 0;
  _result = Hqp_Infeasible;

  if ( _logging )
    printf("\nHqp_IpsMehrotra::cold_start\n");

  if ( _m > 0 ) {

    // initialization of _z and _w
    if ( _init_method == 1 ) {
      v_set(_z, 1.0);
      v_set(_w, max(v_norm_inf(_qp->d),1e-10)*sp_norm_inf(_qp->Q)/
	    sp_norm_inf(_qp->C));
    } else if ( _init_method == 2 ) {
      v_set(_z, 1.0);
      v_set(_w, sp_norm_inf(_qp->C)/max(v_norm_inf(_qp->d),1e-10)/
	    sp_norm_inf(_qp->Q));
    } else {
      v_ones(_z);
      v_ones(_w);
    }

    //  init right hand side

    v_copy(_qp->c, _r1);
    v_copy(_qp->b, _r2);
    sv_mlt(-1.0, _r2, _r2);
    v_copy(_qp->d, _r3);
    sv_mlt(-1.0, _r3, _r3);
    v_star(_z, _w, _r4);
    sv_mlt(-1.0, _r4, _r4);
    if ( !_init_method )
      v_zero(_r4);

    // factorize and solve
    
#ifdef m_catch
    m_catch(E_SING,
	    // try
	    _matrix->factor(_qp, _z, _w);
	    residuum =
	      _matrix->solve(_qp, _z, _w,
			     _r1, _r2, _r3, _r4, _dx, _dy, _dz, _dw),
	    // catch(E_SING)
	    if ( _logging )
	      printf("\nHqp_Degenerate: vmin(_z) = %g, vmin(_w) = %g\n", 
		     v_min(_z, &izmin), v_min(_w, &iwmin));
	    _result = Hqp_Degenerate;
	    return);
#else
    _matrix->factor(_qp, _z, _w);
    if ((code = setjmp(restart)) != 0) {
      set_err_flag(EF_EXIT);	// avoid recursive error calls
      if (code == E_SING) {
	if ( _logging )
	  printf("\nHqp_Degenerate: vmin(_z) = %g, vmin(_w) = %g\n", 
		 v_min(_z, &izmin), v_min(_w, &iwmin));
	_result = Hqp_Degenerate;
	return;
      }
      else
	error(code, "Hqp_IpsMehrotra::step");
    }
    else {
#   ifdef DEBUG
      set_err_flag(EF_JUMP);
#   else
      set_err_flag(EF_SILENT);
#   endif
    }
    residuum = _matrix->solve(_qp, _z, _w,
			      _r1, _r2, _r3, _r4, _dx, _dy, _dz, _dw);
#endif
    v_copy(_dx, _qp->x);
    v_copy(_dy, _y);

    if ( _init_method == 3 ) {
      v_add(_dz, _z, _dz);
      v_add(_dw, _w, _dw);
    }

    if ( v_norm_inf(_dz) == 0.0 )
      v_set(_dz, 1.0e-10);
    if ( v_norm_inf(_dw) == 0.0 )
      v_set(_dw, 1.0e-10);
    delz = max(-1.5*v_min(_dz, &i), 0.0);
    delw = max(-1.5*v_min(_dw, &i), 0.0);
    v_set(_d1, delz);
    v_add(_d1, _dz, _d1);
    v_set(_d2, delw);
    v_add(_d2, _dw, _d2);
    _gap = in_prod(_d1, _d2);
    delz +=0.5*_gap/(v_sum(_dw)+_m*delw);
    delw +=0.5*_gap/(v_sum(_dz)+_m*delz);
    v_set(_z, delz);
    v_add(_z, _dz, _z);
    v_set(_w, delw);
    v_add(_w, _dw, _w);
    if ( _logging )
      printf("delz = %g, delw = %g\n", delz, delw);
    v_set(_z_hot, 1.0);
    v_set(_w_hot, 1.0);
  }

  // initialize a program without inequality constraints
  else {
    v_zero(_qp->x);
    v_zero(_y);
  }
}

//--------------------------------------------------------------------------
void Hqp_IpsMehrotra::hot_start()
{
  _iter = 0;
  v_zero(_phimin);
  _alpha = 1.0;

  _hot_started = 1;
  _result = Hqp_Infeasible;

  if ( _logging )
    printf("\nHqp_IpsMehrotra::hot_start\n");

  if (_m > 0) {  
    v_copy(_z_hot, _z); // it is not necessary to correct _z or _w because of
    v_copy(_w_hot, _w); // Mehrotra's adaptive step size!
  }

  // initialize a program without inequality constraints
  else {
    v_zero(_qp->x);
    v_zero(_y);
  }
}

//--------------------------------------------------------------------------
void Hqp_IpsMehrotra::step()
{
  int	i, izmin, iwmin;
  Real 	residuum, smm, phi, pcost, norm_r, pm, pm30, t, gamma;
  Real  alpha_aff, alpha_corr, mu_aff, mu, sigma, zmin, wmin, fpd, mu_pl;
  int 	code, mod_terlaky = 1;

  if ( ( v_min(_z, &izmin) <= 0.0 ) || ( v_min(_w, &iwmin) <= 0.0 ) )
    printf("Should never occure: min(_z)=%g, min(_w)=%g\n", 
	   _z->ve[izmin], _w->ve[iwmin]);

  //   actual residuum of KKT conditions

  sp_mv_symmlt(_qp->Q, _qp->x, _r1);
  pcost = 0.5*in_prod(_qp->x, _r1)+in_prod(_qp->x, _qp->c);
  v_add(_r1, _qp->c, _r1);
  _gap = in_prod(_qp->x, _r1)+in_prod(_y, _qp->b)+in_prod(_z, _qp->d);
  sp_vm_mltadd(_r1, _y, _qp->A, -1.0, _r1);
  sp_vm_mltadd(_r1, _z, _qp->C, -1.0, _r1);
  //-  sv_mlt(-1.0, _r1, _r1);

  sp_mv_mlt(_qp->A, _qp->x, _r2);
  v_add(_r2, _qp->b, _r2);
  sv_mlt(-1.0, _r2, _r2); //-

  sp_mv_mlt(_qp->C, _qp->x, _r3);
  v_add(_r3, _qp->d, _r3);
  v_sub(_r3, _w, _r3);
  sv_mlt(-1.0, _r3, _r3); //-

  v_star(_z, _w, _r4);
  sv_mlt(-1.0, _r4, _r4); //-

  mu = in_prod(_z, _w)/_m;

  // convergence measure
  if ( _logging && 0 ) {
    printf("norm_r: %g %g %g \n", 
	   v_norm_inf(_r1), v_norm_inf(_r2), v_norm_inf(_r3));
    printf("z,w: %g %g %g %g %g %g\n", 
	   v_min(_z, &i), v_max(_z, &i), v_min(_w, &i), v_max(_w, &i), 
	   v_min(_r4, &i), v_max(_r4, &i));
  }
  norm_r = max(max(v_norm_inf(_r1), v_norm_inf(_r2)), v_norm_inf(_r3));
  if ( _iter == 0 ) {
    _mu0 = mu;
    _norm_r0 = norm_r;
    _norm_data = max(max(max(max(max(sp_norm_inf(_qp->Q), sp_norm_inf(_qp->A)),
				 sp_norm_inf(_qp->C)), v_norm_inf(_qp->c)), 
			 v_norm_inf(_qp->b)), v_norm_inf(_qp->d));
    if ( _logging )
      printf("eps=%g, norm_r0=%g, norm_data=%g\n", _eps, _norm_r0, _norm_data);
  }
  _test = phi = _phimin->ve[_iter] = (norm_r+fabs(_gap))/_norm_data;
  if ( _logging ) {
    if ( !( _iter % 15 ) )
      printf("\niter       _gap         mu        phi    sigma    alpha\n");
    printf("%4d %10.4g %10.4g %10.4g ", _iter, _gap, mu, phi);
  }

  //   prepare hot start

  //   5/5/98  0.5 ---> 0.3333
  if ( _test > pow(_eps, 0.3333) ) {
    v_copy(_z, _z_hot);
    v_copy(_w, _w_hot);
  }

  //   check for termination

  if ( ( mu <= _eps ) && ( norm_r <= _eps*_norm_data ) ) {
    _result = Hqp_Optimal;
    return;
  }

  //   check for infeasibility

  for ( i = 1, pm = _phimin->ve[0]; i <= _iter; i++ )
    pm = min(pm, _phimin->ve[i]);
  if ( ( phi > _eps ) && ( phi >= 1.0e4*pm ) ) {
    //    _result = Hqp_Infeasible;   // should be "infeasible"
    _result = Hqp_Suboptimal;      // should be "infeasible"
    if ( _logging )
      printf("\nHqp_Suboptimal: phi = %g, min(phi) = %g\n", phi, pm);
    return;
  }

  //   check for slow convergence

  if ( _iter >= 30 ) {
    for ( i = 2, pm30 = _phimin->ve[1]; i <= _iter-30; i++ )
      pm30 = min(pm30, _phimin->ve[i]);
    if ( pm >= 0.5*pm30 ) {
      _result = Hqp_Suboptimal;   // should be "unknown (slow progress)"
      if ( _logging )
	printf("\nHqp_Suboptimal: phi = %g, min(phi) = %g, min30(phi) = %g\n",
	       phi, pm, pm30);
      return;
    }
  }

  //   check for blowup in infeasibility-to-duality ratio

  if ( ( norm_r > _eps*_norm_data ) && ( norm_r/mu >= 1.0e8*_norm_r0/_mu0 ) ) {
    _result = Hqp_Suboptimal;   // should be "unknown (blow up)"
    if ( _logging )
      printf("\nHqp_Suboptimal: norm_r/mu = %g\n", norm_r/mu);
  }

  //   factorize the matrix
  //   and predictor (affine) step calculation

#ifdef m_catch
  m_catch(E_SING,
	  // try
	  _matrix->factor(_qp, _z, _w);
	  residuum =
	    _matrix->solve(_qp, _z, _w,
			   _r1, _r2, _r3, _r4, _dxa, _dya, _dza, _dwa),
	  // catch(E_SING)
	  if ( _logging )
	    printf("\nHqp_Degenerate: vmin(_z) = %g, vmin(_w) = %g\n", 
		   v_min(_z, &izmin), v_min(_w, &iwmin));
	  _result = Hqp_Degenerate;
	  return);

#else
  _matrix->factor(_qp, _z, _w);
  if ((code = setjmp(restart)) != 0) {
    set_err_flag(EF_EXIT);	// avoid recursive error calls
    if (code == E_SING) {
      _result = Hqp_Degenerate;

      if ( _logging )
	printf("\nHqp_Degenerate: vmin(_z) = %g, vmin(_w) = %g\n", 
	       v_min(_z, &izmin), v_min(_w, &iwmin));
      return;
    }
    else
      error(code, "Hqp_IpsMehrotra::step");
  }
  else {
#   ifdef DEBUG
      set_err_flag(EF_JUMP);
#   else
      set_err_flag(EF_SILENT);
#   endif
  }
  residuum = _matrix->solve(_qp, _z, _w,
			    _r1, _r2, _r3, _r4, _dxa, _dya, _dza, _dwa);
#endif

  //   predictor step size determination (find maximal feasible step)

  alpha_aff = 1.0;
  for ( i = 0; i < _m; i++ ) {
    if ( _dza->ve[i] < 0.0 ) //-
      alpha_aff = min(alpha_aff, -_z->ve[i]/_dza->ve[i]);
    if ( _dwa->ve[i] < 0.0 ) //-
      alpha_aff = min(alpha_aff, -_w->ve[i]/_dwa->ve[i]);
  }

  alpha_aff = max(0.0, min(alpha_aff, 1.0));

  //   centering parameter sigma

  if ( !mod_terlaky ) {
    // Mehrotra's original version
    v_mltadd(_z, _dza, alpha_aff, _d1); //-
    v_mltadd(_w, _dwa, alpha_aff, _d2); //-
    mu_aff = in_prod(_d1, _d2)/_m;
    sigma = pow(mu_aff/mu, 3.0);
  } else {
    // Terlaky's modification
    gamma = pow(1.0e-4, 0.25);
    for ( i = 0, t= 0.0; i < _m; i++ )
      if ( _dza->ve[i]*_dwa->ve[i] > 0.0 )
	t = max(t, _dza->ve[i]*_dwa->ve[i]/_z->ve[i]/_w->ve[i]);
    sigma = gamma*(t+1.0-alpha_aff)/(1.0-gamma);
  }
  smm = sigma*mu;

  //   centering and corrector step calculation

  if ( ( !mod_terlaky ) || ( alpha_aff >= 0.1 ) ) {
    v_star(_z, _w, _r4);
    for (i = 0; i < _m; i++ )
      _r4->ve[i] += _dza->ve[i]*_dwa->ve[i]-smm;
    sv_mlt(-1.0, _r4, _r4); //-
    residuum = _matrix->solve(_qp, _z, _w,
			      _r1, _r2, _r3, _r4, _dx, _dy, _dz, _dw);
  }
  if ( mod_terlaky ) {
    alpha_corr = 1.0;
    for ( i = 0; i < _m; i++ ) {
      if ( _dz->ve[i] < 0.0 ) //-
	alpha_corr = min(alpha_corr, -_z->ve[i]/_dz->ve[i]);
      if ( _dw->ve[i] < 0.0 ) //-
	alpha_corr = min(alpha_corr, -_w->ve[i]/_dw->ve[i]);
    }
    alpha_corr = max(0.0, min(alpha_corr, 1.0));
    if ( ( alpha_aff < 0.1 ) || ( alpha_corr < gamma*gamma/2.0/_m/_m ) ) {
      sigma = gamma/(1.0-gamma);
      smm = sigma*mu;
      v_star(_z, _w, _r4);
      sv_mlt(-1.0, _r4, _r4); //-
      for (i = 0; i < _m; i++ )
	_r4->ve[i] -= _dza->ve[i]*_dwa->ve[i]-smm; //-
      residuum = _matrix->solve(_qp, _z, _w,
				_r1, _r2, _r3, _r4, _dx, _dy, _dz, _dw);
    }
  }

  //   step size determination (Mehrotra´s adaptive algorithm)

  zmin = HUGE_VAL;
  izmin = -1;
  wmin = HUGE_VAL;
  iwmin = -1;
  for ( i = 0; i < _m; i++ ) {
    if ( ( _logging ) && ( _z->ve[i] <= 0.0 ) )
      printf("_z[%d] = %g\n", i, _z->ve[i]);
    if ( ( _logging ) && ( _w->ve[i] <= 0.0 ) )
      printf("_w[%d] = %g\n", i, _w->ve[i]);
    if ( _dz->ve[i] < 0.0 ) //-
      if ( -_z->ve[i]/_dz->ve[i] < zmin ) {
	izmin = i;
	zmin = -_z->ve[i]/_dz->ve[i];
      }
    if ( _dw->ve[i] < 0.0 )  //-
      if ( -_w->ve[i]/_dw->ve[i] < wmin ) {
	iwmin = i;
	wmin = -_w->ve[i]/_dw->ve[i];
      }
  }

  if ( ( izmin < 0 ) && ( iwmin < 0 ) )
    _alpha = 1.0;
  else {
    if ( izmin < 0 )
      _alpha = wmin;
    else if ( iwmin < 0 )
      _alpha = zmin;
    else 
      _alpha = min(zmin, wmin);
    v_mltadd(_z, _dz, _alpha, _d1); //-
    v_mltadd(_w, _dw, _alpha, _d2); //-
    mu_pl = in_prod(_d1, _d2)/_m;
    if ( ( _alpha == wmin ) && ( _z->ve[iwmin] > -_alpha*_dz->ve[iwmin] ) )
      fpd = (_gammaf*mu_pl/(_z->ve[iwmin]+_alpha*_dz->ve[iwmin])-_w->ve[iwmin])/
	(_alpha*_dw->ve[iwmin]); //-
    else if ( ( _alpha == zmin ) && ( _w->ve[izmin] > -_alpha*_dw->ve[izmin] ) )
      fpd = (_gammaf*mu_pl/(_w->ve[izmin]+_alpha*_dw->ve[izmin])-_z->ve[izmin])/
	(_alpha*_dz->ve[izmin]);
    else
      fpd = 0;
    _alpha = max(1-_gammaf, fpd)*_alpha;
  }
  if ( _logging ) 
    printf("%8.2g %8.2g %8.2g %8.2g %8.2g\n", 
	   sigma, _alpha, alpha_aff, alpha_corr, fpd);

  //   perform step

  v_mltadd(_qp->x, _dx, _alpha, _dx); //-
  v_mltadd(_y, _dy, _alpha, _y); //-
  v_mltadd(_z, _dz, _alpha, _z); //-
  v_mltadd(_w, _dw, _alpha, _w); //-

  mu = in_prod(_z, _w)/_m;

  if (!is_finite(mu) || !is_finite(v_norm_inf(_dx))) {
    _result = Hqp_Degenerate;

    if ( _logging )
      printf("\nHQP_Degenerate: mu = %g, _dx = %g\n", mu, v_norm_inf(_dx));
    return;
  }

  v_copy(_dx, _qp->x); 
  _iter++;
}

//--------------------------------------------------------------------------
void Hqp_IpsMehrotra::solve()
{
  Real test1 = 0.0;

  _fail_iters = 0;

  do {
    do {
      step();
      if (_hot_started) {
	if (_iter == 1) {
	  test1 = _test;
	} else
	  //   5/5/98   2.0 --> 1.2
	  if ( ( _test > test1/pow(1.2,_iter-1.0) ) || ( _alpha < 1.0e-5 ) ) {
	    //printf("Restarted cold after %d iters (%g)\n", _iter, _test);
	    _fail_iters += _iter;
	    cold_start();
	  }
      }

      if (_iter + _fail_iters >= _max_iters) break;
      if (_hot_started && _iter >= _max_warm_iters) break;

    } while (_result != Hqp_Optimal
	     && _result != Hqp_Suboptimal && _result != Hqp_Degenerate);

    if (_hot_started && _result != Hqp_Optimal) {
      //fprintf(stderr, "Bad hot-start, lost %d iters\n", _iter);
      _fail_iters += _iter;
      cold_start();
    } else
      break;

  } while (1);

  _iter += _fail_iters;
}

//==========================================================================


