/*
 * Hqp_IpBdLU.C -- class definition
 *
 * rf, 9/7/94
 */

/*
    Copyright (C) 1994--1998  Ruediger Franke

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

extern "C" {
#include <matrix2.h>
}

#include "Hqp_Program.h"
#include "Hqp_Structure.h"
#include "Hqp_IpBdLU.h"


//--------------------------------------------------------------------------
Hqp_IpBdLU::Hqp_IpBdLU()
{
  _n = _me = _m = 0;
  _CT = SMNULL;
  _J = BDNULL;
  _J_raw = BDNULL;
  _J_fct = BDNULL;
  _QP2J = PNULL;
  _J2QP = PNULL;
  _pivot = PNULL;
  _zw = VNULL;
  _scale = VNULL;
  _r12 = VNULL;
  _xy = VNULL;
  _test = VNULL;
  _CTC_structure = new Hqp_Structure;
}

//--------------------------------------------------------------------------
Hqp_IpBdLU::~Hqp_IpBdLU()
{
  sp_free(_CT);
  bd_free(_J);
  bd_free(_J_raw);
  bd_free(_J_fct);
  px_free(_QP2J);
  px_free(_J2QP);
  px_free(_pivot);
  v_free(_zw);
  v_free(_scale);
  v_free(_r12);
  v_free(_xy);
  v_free(_test);
  delete _CTC_structure;
}

//--------------------------------------------------------------------------
BAND *Hqp_IpBdLU::sub_CTC(const PERM *px, BAND *Q)
// return Q - _CT * diags(_zw) * _CT'
// read:  _CT, _zw
// write: _scale
{
  int i, j, j_idx, j_end;
  int qi, qj;
  SPROW *crow;
  Real sum, val;
  const IVEC *neigh;
  int lb;
  Real **Qmat;

  lb = Q->lb;
  Qmat = Q->mat->me;
  crow = _CT->row;
  for (i=0; i<_n; i++, crow++) {

    if (crow->len <= 0) {
      // just init scaling
      val = Qmat[lb][px->pe[i]];
      _scale->ve[i] = min(1.0, sqrt(-1.0 / val));
    }
    else {

      // calculate diagonal entry
      sum = sprow_inprod(crow, _zw, crow);
      val = Qmat[lb][px->pe[i]] -= sum;
      _scale->ve[i] = min(1.0, sqrt(-1.0 / val));

      // calculate resting entries
      neigh = _CTC_structure->neighbours(i);
      j_end = neigh->dim;
      for (j_idx=0; j_idx<j_end; j_idx++) {
	j = neigh->ive[j_idx];
	if (j < i) {
	  sum = sprow_inprod(crow, _zw, _CT->row + j);
	  qi = px->pe[i];
	  qj = px->pe[j];

	  // substract sum from Qij and Qji

	  Qmat[lb+qj-qi][qj] -= sum;
	  Qmat[lb+qi-qj][qi] -= sum;
	}
      }
    }
  }

  return Q;  
}

//--------------------------------------------------------------------------
void Hqp_IpBdLU::init(const Hqp_Program *qp)
{
  Hqp_Structure *sp = new Hqp_Structure;
  SPMAT *QCTC;
  SPROW *r1, *r2;
  int i, j;
  int dim;
  Real sum;
  const Hqp_BandSize *bds;

  _n = qp->c->dim;
  _me = qp->b->dim;
  _m = qp->d->dim;
  dim = _n + _me;

  // reallocations

  _pivot = px_resize(_pivot, dim);
  _zw = v_resize(_zw, _m);
  _scale = v_resize(_scale, _n);
  _r12 = v_resize(_r12, dim);
  _xy = v_resize(_xy, dim);
  _test = v_resize(_test, dim);

  // store C' for further computations
  // analyze structure of C'*C

  _CT = sp_transp(qp->C, _CT);
  sp_ones(_CT);
  v_ones(_zw);
  QCTC = sp_get(_n, _n, 10);
  r1 = _CT->row;
  for (i=0; i<_n; i++, r1++) {
    r2 = r1;
    for (j=i; j<_n; j++, r2++) {
      sum = sprow_inprod(r1, _zw, r2);
      if (sum != 0.0) {
	sp_set_val(QCTC, i, j, sum);
	if (i != j)
	  sp_set_val(QCTC, j, i, sum);
      }
    }
  }
  _CTC_structure->init_spmat(QCTC);

  // initialize structure of reduced qp

  QCTC = sp_add(qp->Q, QCTC, QCTC);
  sp->init_QA(QCTC, qp->A);
  sp->order_rcm();
  _QP2J = sp->px_get(_QP2J);
  _J2QP = px_inv(_QP2J, _J2QP);
  bds = sp->bd_size();
  _J = bd_resize(_J, bds->lb, bds->ub, dim);
  _J_raw = bd_resize(_J_raw, bds->lb, bds->ub, dim);
  _J_fct = bd_resize(_J_fct, bds->lb, bds->ub + bds->lb, dim);

  delete sp;
  sp_free(QCTC);

  // prepare iterations

  update(qp);
}

//--------------------------------------------------------------------------
void Hqp_IpBdLU::update(const Hqp_Program *qp)
{
  // prepare _J_raw
  m_zero(_J_raw->mat);
  sp_into_bd(qp->Q, -1.0, _J_raw, _QP2J, 0, 0);
  spT_into_bd(qp->A, 1.0, _J_raw, _QP2J, 0, _n);
  sp_into_bd(qp->A, 1.0, _J_raw, _QP2J, _n, 0);

  // update _CT
  _CT = sp_transp(qp->C, _CT);
}

//--------------------------------------------------------------------------
void Hqp_IpBdLU::factor(const VEC *z, const VEC *w)
{
  assert(z->dim == _m && w->dim == _m);

  int	i, i_end;
  int	j, j_end, k, l;
  int	lb, ub, bw;
  Real	scale;
  Real	**Jmat;
  int	dim = _n+_me;

  // copy _J_raw into _J

  bd_copy(_J_raw, _J);

  // augment _J
  v_slash(w, z, _zw);
  sub_CTC(_QP2J, _J);

  // diagonal scaling

  lb = _J->lb;
  ub = _J->ub;
  bw = _J->mat->m;
  Jmat = _J->mat->me;
  j_end = _n;
  for (j=0; j<j_end; j++) {
    scale = _scale->ve[j];
    k = _QP2J->pe[j];
    i_end = _J_raw->mat->m;	// scale k'th col of _J
    for (i=0; i<i_end; i++) {
      Jmat[i][k] *= scale;
    }
    i = max(0, lb-k);		// scale k'th row of _J
    l = max(0, k-lb);
    i_end = min(bw, bw + dim-k-1 - ub);
    for (; i<i_end; i++, l++) {
      Jmat[i][l] *= scale;
    }
  }  

  // copy _J to _J_fct for factorization
  // (_J is hold for calculating residuum)

  bd_resize(_J_fct, _J->lb, _J->ub, dim);	// no memory reallocation!
  bd_copy(_J, _J_fct);

  // factorization of _J_fct

  bdLUfactor(_J_fct, _pivot);
}

//--------------------------------------------------------------------------
Real Hqp_IpBdLU::solve(const VEC *r1, const VEC *r2, const VEC *r3,
		       VEC *dx, VEC *dy, VEC *dz)
{
  static VEC v;
  Real residuum;

  assert(r1->dim == _n && dx->dim == _n);
  assert(r2->dim == _me && dy->dim == _me);
  assert(r3->dim == _m && dz->dim == _m);

  // augment, copy, scale and permutate [r1;r2;r3] into _r12
  // calculate, permutate and scale x
  
  // augment r1
  v_part(_r12, 0, _n, &v);
  v_star(_zw, r3, dz);
  sp_mv_mlt(_CT, dz, &v);
  v_sub(r1, &v, &v);
  v_star(&v, _scale, &v);

  v_copy(r2, v_part(_r12, _n, _me, &v));

  px_vec(_J2QP, _r12, _r12);

  bdLUsolve(_J_fct, _pivot, _r12, _xy);

  bd_mv_mlt(_J, _xy, _test);
  v_sub(_r12, _test, _test);
  residuum = v_norm2(_test);

  px_vec(_QP2J, _xy, _xy);

  v_star(v_part(_xy, 0, _n, &v), _scale, dx);
  v_copy(v_part(_xy, _n, _me, &v), dy);

  sp_vm_mlt(_CT, dx, dz);
  v_sub(r3, dz, dz);
  v_star(_zw, dz, dz);

  return residuum;
}


//==========================================================================
