/*
 * Hqp_IpMehrotra.C -- class definition
 *
 * rf, 5/28/94
 *
 * rf, 11/5/95: inserted _max_warm_iters
 * rf, 2/12/97: treat singular matrix errors
 * rf, 4/20/97: treat numerical overflow for _qp->x and _gap
 * ea, 5/9/97
 *   - Mehrotra´s primal-dual predictor-corrector method for QP problems.
 *     Mehrotra, S.:  On the implementation of a primal-dual interior point
 *                    method. SIAM J. Optimization 2(1992)4, 575-601.
 *     Wright, S. J.: Primal-dual interior-point methods.
 *                    SIAM, Philadelphia, 1997.
 *     Czyzyk, J. and Mehrotra, S. and Wright, S. J.: PCx User Guide.
 *                    Technical Report OTC 96/01, 
 *                    Optimization Technology Center, 1997.  
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
 */

/*
    Copyright (C) 1994--1998  Eckhard Arnold and Ruediger Franke

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
#include <If_Float.h>
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

  _matrix = new Hqp_IpRedSpBKP;

  _alpha = 1.0;
  _gammaf = 0.01;
  _mu0 = 0.0;      // not used
  _beta = 0.995;   // not used
  _fail_iters = 0;
  _max_warm_iters = 25;
  _logging = 0;

  _ifList.append(new If_Float("qp_gap", &_gap));
  _ifList.append(new If_Float("qp_alpha", &_alpha));
  _ifList.append(new If_Float("qp_beta", &_beta));
  //  _ifList.append(new If_Float("qp_rhomin", &_rhomin));
  _ifList.append(new If_Float("qp_mu0", &_mu0));
  //  _ifList.append(new If_Float("qp_Ltilde", &_Ltilde));
  _ifList.append(new If_Int("qp_fail_iters", &_fail_iters));
  _ifList.append(new If_Int("qp_max_warm_iters", &_max_warm_iters));
  _ifList.append(new If_Int("qp_logging", &_logging));
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
  delete _matrix;
}

//--------------------------------------------------------------------------
int Hqp_IpsMehrotra::init(IF_CMD_ARGS)
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

  // fill up internal data

  _matrix->init(_qp);

  if ( _logging )
    printf("\nHqp_IpsMehrotra::init\n");

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_IpsMehrotra::update(IF_CMD_ARGS)
{
  _matrix->update(_qp);
  _phimin = v_resize(_phimin, _max_iters+1);

  if ( _logging )
    printf("\nHqp_IpsMehrotra::update\n");

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_IpsMehrotra::cold_start(IF_CMD_ARGS)
{
  int  i, izmin, iwmin, code;
  Real delz, delw, residuum;
  int c_new = 0;

  if ( _m > 0 ) {

    //   ´very cold´ start
    //    v_zero(_qp->x);
    //    v_zero(_y);
    //    v_set(_z, 1.0);
    //    v_set(_w, 1.0);

    v_set(_z, 1.0);
    //    v_set(_w, max(v_norm_inf(_qp->d),1e-10)*sp_norm_inf(_qp->Q)/sp_norm_inf(_qp->C));
    v_set(_w, sp_norm_inf(_qp->C)/max(v_norm_inf(_qp->d),1e-10)/
	  sp_norm_inf(_qp->Q));
    if ( c_new )
      v_set(_w, 1e10);

#ifdef m_catch
    //  init right hand side

    v_copy(_qp->c, _r1);
    v_copy(_qp->b, _r2);
    sv_mlt(-1.0, _r2, _r2);
    v_copy(_qp->d, _r3);
    sv_mlt(-1.0, _r3, _r3);
    v_star(_z, _w, _r4);
    sv_mlt(-1.0, _r4, _r4);

    // factorize and solve
    
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
	    return IF_OK);
#else
    //   factorize the matrix

    _matrix->factor(_qp, _z, _w);
    if ((code = setjmp(restart)) != 0) {
      set_err_flag(EF_EXIT);	// avoid recursive error calls
      if (code == E_SING) {
	if ( _logging )
	  printf("\nHqp_Degenerate: vmin(_z) = %g, vmin(_w) = %g\n", 
		 v_min(_z, &izmin), v_min(_w, &iwmin));
	_result = Hqp_Degenerate;
	return IF_OK;
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

    //   right hand side and solve

    v_copy(_qp->c, _r1);
    v_copy(_qp->b, _r2);
    sv_mlt(-1.0, _r2, _r2);
    v_copy(_qp->d, _r3);
    sv_mlt(-1.0, _r3, _r3);
    v_star(_z, _w, _r4);
    sv_mlt(-1.0, _r4, _r4);
    
    residuum = _matrix->solve(_qp, _z, _w,
			      _r1, _r2, _r3, _r4, _dx, _dy, _dz, _dw);
#endif
    v_copy(_dx, _qp->x);
    v_copy(_dy, _y);

    if ( c_new )
      for ( i = 0; i < _m; i++ ) {
	if ( _dw->ve[i] < 0.0 ) {
	  _dz->ve[i] = -_dw->ve[i];
	  _dw->ve[i] = 0.0;
	} else
	  _dz->ve[i] = 0.0;
      }
    v_add(_dz, _z, _dz);
    v_add(_dw, _w, _dw);
    delz = v_min(_dz, &izmin);
    delw = v_min(_dw, &iwmin);
    if ( _logging )
      printf("_dz[%d] = %g, _dw[%d] = %g, ", izmin, delz, iwmin, delw);
    if ( v_norm_inf(_dz) == 0.0 )
      v_set(_dz, 1.0e-10);
    if ( v_norm_inf(_dw) == 0.0 )
      v_set(_dw, 1.0e-10);
    delz = max(-1.5*v_min(_dz, &i), 0.0);
    delw = max(-1.5*v_min(_dw, &i), 0.0);
    v_set(_r3, delz);
    v_add(_r3, _dz, _r3);  // 5/5/98
    v_set(_r4, delw);
    v_add(_r4, _dw, _r4);  // 5/5/98
    _gap = in_prod(_r3, _r4);
    delz +=0.5*_gap/(v_sum(_dw)+_m*delw);
    delw +=0.5*_gap/(v_sum(_dz)+_m*delz);
    for ( i = 0; i < _m; i++ ) {
      _z->ve[i] = _dz->ve[i]+delz;
      _w->ve[i] = _dw->ve[i]+delw;
    }
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

  _iter = 0;
  v_zero(_phimin);
  _alpha = 1.0;

  _hot_started = 0;
  _result = Hqp_Infeasible;

  if ( _logging )
    printf("\nHqp_IpsMehrotra::cold_start\n");

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_IpsMehrotra::hot_start(IF_CMD_ARGS)
{
  if (_m > 0) {  
    v_copy(_z_hot, _z); // it is not necessary to correct _z or _w because of
    v_copy(_w_hot, _w); // Mehrotra's adaptive step size!
  }

  // initialize a program without inequality constraints
  else {
    v_zero(_qp->x);
    v_zero(_y);
  }

  _iter = 0;
  v_zero(_phimin);
  _alpha = 1.0;

  _hot_started = 1;
  _result = Hqp_Infeasible;

  if ( _logging )
    printf("\nHqp_IpsMehrotra::hot_start\n");

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_IpsMehrotra::step(int, char *[], char **ret)
{
  int	i, izmin, iwmin;
  Real 	residuum, smm, phi, rmu0, pcost, rmu;
  Real  alfa_aff, mu_aff, mu, sigma, zmin, wmin, fpd, mu_pl;
  int 	code;

  //   actual residuum of KKT conditions

  sp_mv_symmlt(_qp->Q, _qp->x, _r1);
  pcost = 0.5*in_prod(_qp->x, _r1)+in_prod(_qp->x, _qp->c);
  v_add(_r1, _qp->c, _r1);
  _gap = in_prod(_qp->x, _r1)+in_prod(_y, _qp->b)+in_prod(_z, _qp->d);
  sp_vm_mltadd(_r1, _y, _qp->A, -1.0, _r1);
  sp_vm_mltadd(_r1, _z, _qp->C, -1.0, _r1);
  sv_mlt(-1.0, _r1, _r1);

  sp_mv_mlt(_qp->A, _qp->x, _r2);
  v_add(_r2, _qp->b, _r2);
  
  sp_mv_mlt(_qp->C, _qp->x, _r3);
  v_add(_r3, _qp->d, _r3);
  v_sub(_r3, _w, _r3);

  v_star(_z, _w, _r4);

  _test = phi = fabs(_gap)/
    max(1.0, max(v_norm_inf(_qp->c), 
		 max(v_norm_inf(_qp->b), v_norm_inf(_qp->d))))+
    v_norm_inf(_r1)/max(1.0, v_norm_inf(_qp->c))+
    v_norm_inf(_r2)/max(1.0, v_norm_inf(_qp->b))+
    v_norm_inf(_r3)/max(1.0, v_norm_inf(_qp->d));
  mu = in_prod(_z, _w);

  if ( _logging ) {
    if ( !( _iter % 15 ) )
      printf("\niter       _gap         mu        phi        rmu    sigma      alfa\n");
    printf("%4d %10.4g %10.4g %10.4g ", _iter, _gap, mu, phi);
  }

  //   prepare hot start

  //   5/5/98  0.5 ---> 0.3333
  if ( _test > pow(_eps, 0.3333) ) {
  //  if ( ( _test > pow(_eps, 0.5) ) && ( mu > pow(_eps, 0.5) ) || 
  //       ( _iter == 0 ) ) {
    v_copy(_z, _z_hot);
    v_copy(_w, _w_hot);
  }

  //   check for termination
  if ( ( v_norm_inf(_r1)/(1.0+v_norm_inf(_qp->c)) <= _eps ) &&
       ( v_norm_inf(_r2)/(1.0+v_norm_inf(_qp->b)) <= _eps ) &&
       ( v_norm_inf(_r3)/(1.0+v_norm_inf(_qp->d)) <= _eps ) &&
       ( fabs(_gap)/(1.0+fabs(pcost)) <= _eps ) ) {
    _result = Hqp_Optimal;
    return IF_OK;
  }

  //   _phimin

  if ( _iter == 0 )
    _phimin->ve[0] = phi;
  else
    _phimin->ve[_iter] = min(_phimin->ve[_iter-1], phi);

  //   check for infeasibility

  if ( phi >= max(_eps, 1.0e5*_phimin->ve[_iter]) ) {
    //    _result = Hqp_Infeasible;   // should be "infeasible"
    _result = Hqp_Suboptimal;  

    if ( _logging )
      printf("\nHqp_Suboptimal: _phi = %g, _phimin[_iter] = %g\n",
	     phi, _phimin->ve[_iter]);
    return IF_OK;
  }

  //   check for slow convergence

  if ( ( _iter >= 30 ) && 
       ( _phimin->ve[_iter] >= 0.5*_phimin->ve[_iter-30] ) ) {
    _result = Hqp_Suboptimal;   // should be "unknown"

    if ( _logging )
      printf("\nHqp_Suboptimal: _phimin[_iter] = %g, _phimin[_iter-30] = %g\n",
	     _phimin->ve[_iter], _phimin->ve[_iter-30]);
    return IF_OK;
  }

  //   check for blowup in infeasibility-to-duality ratio

  rmu = max(v_norm_inf(_r1), max(v_norm_inf(_r2), v_norm_inf(_r3)))/mu;
  if ( _logging )
    printf("%10.4g ", rmu);
  if ( _iter == 5 )
    rmu0 = rmu;
  else if ( ( _iter > 5 ) && 
	    ( ( v_norm_inf(_r1)/max(1.0, v_norm_inf(_qp->c)) > _eps ) ||
	      ( v_norm_inf(_r2)/max(1.0, v_norm_inf(_qp->b)) > _eps ) ||
	      ( v_norm_inf(_r3)/max(1.0, v_norm_inf(_qp->d)) > _eps ) ) &&
	    ( rmu/rmu0 >= 1.0e6 ) ) {
    //    printf("\nHqp_Suboptimal: _r1 = %g, _r2 = %g, _r3 = %g, mu = %g, rmu0 = %g\n",
    //	   v_norm_inf(_r1), v_norm_inf(_r2), v_norm_inf(_r3), mu, rmu0);
    //   5/5/98 removed check for 'blow up'
    //    _result = Hqp_Suboptimal;   // should be "unknown"
  }

#ifdef m_catch
  //   factorize the matrix
  //   and predictor (affine) step calculation
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
	  return IF_OK);

#else
  //   factorize the matrix

  _matrix->factor(_qp, _z, _w);
  if ((code = setjmp(restart)) != 0) {
    set_err_flag(EF_EXIT);	// avoid recursive error calls
    if (code == E_SING) {
      _result = Hqp_Degenerate;

      if ( _logging )
	printf("\nHqp_Degenerate: vmin(_z) = %g, vmin(_w) = %g\n", 
	       v_min(_z, &izmin), v_min(_w, &iwmin));
      return IF_OK;
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

  //   predictor (affine) step calculation

  residuum = _matrix->solve(_qp, _z, _w,
			    _r1, _r2, _r3, _r4, _dxa, _dya, _dza, _dwa);
#endif

  //   predictor step size determination (find maximal feasible step)

  alfa_aff = 1.0;
  for ( i = 0; i < _m; i++ ) {
    if ( _dza->ve[i] > 0.0 ) 
      alfa_aff = min(alfa_aff, _z->ve[i]/_dza->ve[i]);
    if ( _dwa->ve[i] > 0.0 ) 
      alfa_aff = min(alfa_aff, _w->ve[i]/_dwa->ve[i]);
  }

  alfa_aff = max(0.0, min(alfa_aff, 1.0));

  //   centering parameter sigma

  v_mltadd(_z, _dza, -alfa_aff, _r3);
  v_mltadd(_w, _dwa, -alfa_aff, _r4);
  mu_aff = in_prod(_r3, _r4);
  sigma = pow(mu_aff/mu, 3.0);

  v_star(_r3, _r4, _r4);
  smm = sigma*mu/_m;
  for (i = 0; i < _m; i++ )
    _r4->ve[i] -= smm;
  v_zero(_r1);  
  v_zero(_r2);  
  v_zero(_r3);  

  //   centering and corrector step calculation

  residuum = _matrix->solve(_qp, _z, _w,
			    _r1, _r2, _r3, _r4, _dx, _dy, _dz, _dw);
  v_add(_dx, _dxa, _dx);
  v_add(_dy, _dya, _dy);
  v_add(_dz, _dza, _dz);
  v_add(_dw, _dwa, _dw);

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
    if ( _dz->ve[i] > 0.0 ) 
      if ( _z->ve[i]/_dz->ve[i] < zmin ) {
	izmin = i;
	zmin = _z->ve[i]/_dz->ve[i];
      }
    if ( _dw->ve[i] > 0.0 ) 
      if ( _w->ve[i]/_dw->ve[i] < wmin ) {
	iwmin = i;
	wmin = _w->ve[i]/_dw->ve[i];
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
    v_mltadd(_z, _dz, -_alpha, _r3);
    v_mltadd(_w, _dw, -_alpha, _r4);
    mu_pl = in_prod(_r3, _r4)/_m;
    if ( ( _alpha == wmin ) && ( _z->ve[iwmin] > _alpha*_dz->ve[iwmin] ) )
      fpd = 1.0-_gammaf*mu_pl/_w->ve[iwmin]/
	(_z->ve[iwmin]-_alpha*_dz->ve[iwmin]);
    else if ( ( _alpha == zmin ) && ( _w->ve[izmin] > _alpha*_dw->ve[izmin] ) )
      fpd = 1.0-_gammaf*mu_pl/_z->ve[izmin]/
	(_w->ve[izmin]-_alpha*_dw->ve[izmin]);
    else
      fpd = 0;
    _alpha = min(max(1-_gammaf, fpd),1.0-MACHEPS)*_alpha;
  }

  _alpha = max(0.0, min(_alpha, 1.0));

  if ( _logging ) 
    printf("%8.2g  %8.2g\n", sigma, _alpha);

  //   perform step

  v_mltadd(_qp->x, _dx, -_alpha, _dx);
  v_mltadd(_y, _dy, -_alpha, _y);
  v_mltadd(_z, _dz, -_alpha, _z);
  v_mltadd(_w, _dw, -_alpha, _w);

  mu = in_prod(_z, _w);

  if (!is_finite(mu) || !is_finite(v_norm_inf(_dx))) {
    _result = Hqp_Degenerate;

    if ( _logging )
      printf("\nHQP_Degenerate: mu = %g, _dx = %g\n", mu, v_norm_inf(_dx));
    return IF_OK;
  }

  v_copy(_dx, _qp->x); 
  _iter++;

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_IpsMehrotra::solve(int, char *[], char **ret)
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

  return IF_OK;
}

//==========================================================================


