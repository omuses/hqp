/*
 * Hqp_SqpSchittkowski.C -- class definition
 *
 * rf, 6/8/94
 */

/*
    Copyright (C) 1994--2001  Ruediger Franke

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

#include <math.h>

#include <If_Real.h>
#include <If_Bool.h>
#include <If_Method.h>

#include "Hqp_SqpSchittkowski.h"
#include "Hqp_Solver.h"
#include "Hqp_SqpProgram.h"
#include "Hqp_Program.h"

typedef If_Method<Hqp_SqpSchittkowski> If_Cmd;

IF_CLASS_DEFINE("Schittkowski", Hqp_SqpSchittkowski, Hqp_SqpSolver);

//--------------------------------------------------------------------------
Hqp_SqpSchittkowski::Hqp_SqpSchittkowski()
{
  _x0 = v_resize(v_get(1), 0);
  _xk = v_resize(v_get(1), 0);
  _re = v_resize(v_get(1), 0);
  _r = v_resize(v_get(1), 0);
  _ve = v_resize(v_get(1), 0);
  _v = v_resize(v_get(1), 0);
  _v0 = v_resize(v_get(1), 0);
  _ve0 = v_resize(v_get(1), 0);
  _ue_ve = v_resize(v_get(1), 0);
  _u_v = v_resize(v_get(1), 0);
  _sgme = v_resize(v_get(1), 0);
  _sgm = v_resize(v_get(1), 0);

  _mu = 0.1;
  _beta = 0.1;
  _damped_multipliers = true;

  _ifList.append(new If_Real("sqp_mu", &_mu));
  _ifList.append(new If_Real("sqp_beta", &_beta));
  _ifList.append(new If_Bool("sqp_damped_multipliers",
			     &_damped_multipliers));
  _ifList.append(new If_Cmd("sqp_init", &Hqp_SqpSchittkowski::init, this));
}

//--------------------------------------------------------------------------
Hqp_SqpSchittkowski::~Hqp_SqpSchittkowski()
{
  v_free(_x0);
  v_free(_xk);
  v_free(_re);
  v_free(_r);
  v_free(_ve);
  v_free(_v);
  v_free(_v0);
  v_free(_ve0);
  v_free(_ue_ve);
  v_free(_u_v);
  v_free(_sgme);
  v_free(_sgm);
}

//--------------------------------------------------------------------------
int Hqp_SqpSchittkowski::init(IF_CMD_ARGS)
{
  int ret = Hqp_SqpSolver::init();

  Hqp_Program *qp = _prg->qp();

  int n = qp->c->dim;
  int me = qp->b->dim;
  int m = qp->d->dim;

  _x0 = v_resize(_x0, n);
  _xk = v_resize(_xk, n);
  _re = v_resize(_re, me);
  _r = v_resize(_r, m);
  _ve = v_resize(_ve, me);
  _v = v_resize(_v, m);
  _ve0 = v_resize(_ve0, me);
  _v0 = v_resize(_v0, m);
  _ue_ve = v_resize(_ue_ve, me);
  _u_v = v_resize(_u_v, m);
  _sgme = v_resize(_sgme, me);
  _sgm = v_resize(_sgm, m);

  v_ones(_re);
  v_ones(_r);
  v_zero(_ve);
  v_zero(_v);

  return ret;
}

//--------------------------------------------------------------------------
VEC *Hqp_SqpSchittkowski::update_sgm(const VEC *r, VEC *sgm)
{
  int i, m;
  const Real *r_ve;
  Real *sgm_ve;
  Real val;

  m = sgm->dim;
  r_ve = r->ve;
  sgm_ve = sgm->ve;
  for (i=0; i<m; i++, r_ve++, sgm_ve++) {
    val = (Real)_iter / sqrt(*r_ve);
    *sgm_ve = min(1.0, val);
  }

  return sgm;
}

//--------------------------------------------------------------------------
VEC *Hqp_SqpSchittkowski::update_r(const VEC *u, const VEC *v, const VEC *sgm,
				   Real dQd, VEC *r)
{
  int i, i_end;
  const Real *u_ve, *v_ve, *sgm_ve;
  Real *r_ve;
  Real m2, uv;
  Real val1, val2;

  m2 = 2.0 * (_re->dim + _r->dim);
  i_end = u->dim;
  u_ve = u->ve;
  v_ve = v->ve;
  sgm_ve = sgm->ve;
  r_ve = r->ve;
  for (i=0; i<i_end; i++, u_ve++, v_ve++, sgm_ve++, r_ve++) {
    val1 = *sgm_ve * *r_ve;
    uv = *u_ve - *v_ve;
    val2 = m2 * uv*uv / dQd;
    if (val2 > val1)	// says no for when val2 == NaN
      *r_ve = val2;
    else
      *r_ve = val1;
  }

  return r;
}

//--------------------------------------------------------------------------
Real Hqp_SqpSchittkowski::phi()
{
  int i, i_end;
  const Real *v_ve, *r_ve, *g_ve;
  Real g, v;
  Real ret;

  ret = _prg->f();

  v_ve = _ve->ve;
  r_ve = _re->ve;
  g_ve = _prg->qp()->b->ve;
  i_end = _re->dim;
  for (i=0; i<i_end; i++, v_ve++, r_ve++, g_ve++) {
    g = *g_ve;
    ret -= *v_ve * g - 0.5 * *r_ve * g*g;	// index set J
  }  

  v_ve = _v->ve;
  r_ve = _r->ve;
  g_ve = _prg->qp()->d->ve;
  i_end = _r->dim;
  for (i=0; i<i_end; i++, v_ve++, r_ve++, g_ve++) {
    g = *g_ve;
    v = *v_ve;
    if (g <= v / *r_ve)
      ret -= v * g - 0.5 * *r_ve * g*g;		// index set J
    else
      ret -= 0.5 * v*v / *r_ve;			// index set K
  }  

  return ret;
}

//--------------------------------------------------------------------------
//  dphi
//   Preconditions:
//    - _u_v, _ue_ve must be initialized
//
Real Hqp_SqpSchittkowski::dphi()
{
  Hqp_Program *qp = _prg->qp();
  int i, i_end;
  const Real *v_ve, *r_ve, *g_ve;
  Real val, *phiv_ve;
  VEC *phix, *phive, *phiv;
  VEC *v_rg = VNULL;
  Real ret;

  // (d Phi / d ve) belongs completely to intex set J

  phive = sv_mlt(-1.0, qp->b, VNULL);

  // initialize (d Phi / d x) with (d f / d x) and equality constraints

  phix = v_copy(qp->c, VNULL);
  v_rg = v_star(_re, qp->b, v_rg);
  v_sub(_ve, v_rg, v_rg);
  sp_vm_mltadd(phix, v_rg, qp->A, -1.0, phix);

  // do a loop to calculate (d Phi / d v) and to get index sets
  // use v_rg for calculating inequality-part of (d Phi / d x)

  v_rg = v_star(_r, qp->d, v_rg);
  v_sub(_v, v_rg, v_rg);

  phiv = v_get(_v->dim);
  phiv_ve = phiv->ve;
  g_ve = qp->d->ve;
  v_ve = _v->ve;
  r_ve = _r->ve;
  i_end = _v->dim;
  for (i=0; i<i_end; i++, g_ve++, v_ve++, r_ve++, phiv_ve++) {
    val = *v_ve / *r_ve;
    if (*g_ve <= val)		// index set J
      *phiv_ve = - *g_ve;
    else {			// index set K
      *phiv_ve = - val;
      v_rg->ve[i] = 0.0;
    }
  }    

  // actualize (d Phi / d x) with inequality constraints of index set J
  sp_vm_mltadd(phix, v_rg, qp->C, -1.0, phix);

  ret = in_prod(phix, qp->x);
  ret += in_prod(phive, _ue_ve);
  ret += in_prod(phiv, _u_v);
    
  v_free(phix);
  v_free(phive);
  v_free(phiv);
  v_free(v_rg);
  
  return ret;  
}

//--------------------------------------------------------------------------
void Hqp_SqpSchittkowski::update_vals()
{
  Hqp_Program *qp = _prg->qp();
  Real dphi0, phi0, phik;
  Real n_alpha;

  // update penalty coeffizients

  update_sgm(_re, _sgme);
  update_sgm(_r, _sgm);
  _y = v_copy(_solver->y(), _y);
  _z = v_copy(_solver->z(), _z);
  update_r(_solver->y(), _ve, _sgme, _sQs, _re);
  update_r(_solver->z(), _v, _sgm, _sQs, _r);

  // find step length and update _prg->x()

  _ue_ve = v_sub(_solver->y(), _ve, _ue_ve);
  _u_v = v_sub(_solver->z(), _v, _u_v);

  _x0 = v_copy(_prg->x(), _x0);
  _ve0 = v_copy(_ve, _ve0);
  _v0 = v_copy(_v, _v0);
  phi0 = phik = phi();
  dphi0 = dphi();

  if (dphi0 > 0.0) {
//    printf("No descending direction (%g).\n", dphi0);
    _alpha = _min_alpha;
  }
  else {
    _alpha = 1.0;
  }
  while (1) {
    _d = sv_mlt(_alpha, qp->x, _d);
    v_add(_x0, _d, _xk);
    v_mltadd(_ve0, _ue_ve, _alpha, _ve);
    v_mltadd(_v0, _u_v, _alpha, _v);
    if (_damped_multipliers && _alpha < 1.0) {
      v_copy(_ve, _y);
      v_copy(_v, _z);
    }
    _prg->x(_xk);
    _prg->update_fbd();
    if (!is_finite(_prg->f())) {
      _alpha *= 0.1;
      continue;
    }
    if (_alpha <= _min_alpha)
      break;
    phik = phi();
    if (phik <= (phi0 + _mu * _alpha * dphi0) || fabs(dphi0) <= _eps)
      break;
    n_alpha = 0.5 * dphi0 * _alpha*_alpha / (dphi0 * _alpha - (phik - phi0));
    /*
    if (fabs(_alpha - n_alpha) < _min_alpha)
      break;
    */
    if (!(n_alpha < _alpha))
      break;
    _alpha *= _beta;
    _alpha = max(_alpha, n_alpha);
  }
  _dphi = dphi0;
  _phi = phi0;
}


//==========================================================================
