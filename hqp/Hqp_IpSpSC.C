/*
 * Hqp_IpSpSC.C -- class definition
 *
 * hl, 96/10/09
 */

/*
    Copyright (C) 1996--1998  Hartmut Linke

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
#include <meschach/sparse2.h>
}

extern "C" {
#include <meschach/matlab.h>
}

#include <If_Int.h>
#include <If_Real.h>

SPMAT *spLTMsolve(SPMAT *, SPMAT *,SPMAT *);
SPMAT *sp_mmt_mlt(const SPMAT *, SPMAT *);
SPMAT *sp_mmt2_mlt(const SPMAT *, const SPMAT *, SPMAT *);
SPMAT *sp_mmt2_mlt_u(const SPMAT *, const SPMAT *, SPMAT *);
SPMAT *spCHOLfac(SPMAT *);
SPMAT *spMODCHOLfac(SPMAT *, VEC *, double);
VEC *spCHOLsol(SPMAT *, VEC *, VEC *);

#include "sprcm.h"
#include "Hqp_Program.h"
#include "Hqp_IpSpSC.h"

IF_CLASS_DEFINE("SpSC", Hqp_IpSpSC, Hqp_IpMatrix);

/* sp_set_val3 -- sets the (i,j) entry of the sparse matrix A */
/* don't destroy flag_diag */
static double	sp_set_val3(SPMAT *A, int i, int j, double val)
{
   SPROW	*r;
   int	idx, idx2, new_len;
   
   if ( A == SMNULL )
     m_error(E_NULL,"sp_set_val3");
   if ( i < 0 || i >= A->m || j < 0 || j >= A->n )
     m_error(E_SIZES,"sp_set_val3");
   
   r = A->row+i;
   idx = sprow_idx(r,j);
   /* printf("sp_set_val: idx = %d\n",idx); */
   if ( idx >= 0 )
   {	r->elt[idx].val = val;	return val;	}
   /* else */ if ( idx < -1 )
   {
      /* Note: this destroys the column & diag access paths */
      /* rf: don't destroy flag_diag */
      A->flag_col = FALSE;
      /* shift & insert new value */
      idx = -(idx+2);	/* this is the intended insertion index */
      if ( r->len >= r->maxlen )
      {
	 r->len = r->maxlen;
	 new_len = max(2*r->maxlen+1,5);
	 if (mem_info_is_on()) {
	    mem_bytes(TYPE_SPMAT,A->row[i].maxlen*sizeof(row_elt),
			    new_len*sizeof(row_elt));
	 }

	 r->elt = RENEW(r->elt,new_len,row_elt);
	 if ( ! r->elt )	/* can't allocate */
	   m_error(E_MEM,"sp_set_val");
	 r->maxlen = new_len;
      }
      for ( idx2 = r->len-1; idx2 >= idx; idx2-- )
	MEM_COPY((char *)(&(r->elt[idx2])),
		 (char *)(&(r->elt[idx2+1])),sizeof(row_elt));
      /************************************************************
	if ( idx < r->len )
	MEM_COPY((char *)(&(r->elt[idx])),(char *)(&(r->elt[idx+1])),
	(r->len-idx)*sizeof(row_elt));
	************************************************************/
      r->len++;
      r->elt[idx].col = j;
      if (j == i)
	r->diag = j;
      else if (j < i)
	r->diag++;
      return r->elt[idx].val = val;
   }
   /* else -- idx == -1, error in index/matrix! */
   return 0.0;
}

//--------------------------------------------------------------------------
Hqp_IpSpSC::Hqp_IpSpSC()
{
  _n = _me = _m = 0;
  _piv_QCVC_flag = 0;
  _piv_AQCVCA_flag = 0;
  _sbw = -1;

  _macheps = MACHEPS;

  _Q = SMNULL;
  _CT = SMNULL;
  _AT = SMNULL;
  _QCVC = SMNULL;
  _QCVCA = SMNULL;
  _AQCVC = SMNULL;
  _AQCVCA = SMNULL;
  _AQCVCA_fac = SMNULL;

  _piv_QCVC = PNULL;
  _piv_AQCVCA = PNULL;

  _zw = VNULL;
  _v1 = VNULL;
  _v2 = VNULL;
  _scale = VNULL;

  _CTC_degree = IVNULL;
  _CTC_neigh_start = IVNULL;
  _CTC_neighs = IVNULL;

  //_ifList.append(new If_Int("mat_sbw", &_sbw)); // not supported
  _ifList.append(new If_Real("mat_macheps", &_macheps));
}

//--------------------------------------------------------------------------
Hqp_IpSpSC::~Hqp_IpSpSC()
{
  sp_free(_Q);
  sp_free(_CT);
  sp_free(_AT);
  sp_free(_QCVC);
  sp_free(_QCVCA);
  sp_free(_AQCVC);
  sp_free(_AQCVCA);
  sp_free(_AQCVCA_fac);

  px_free(_piv_QCVC);
  px_free(_piv_AQCVCA);

  v_free(_zw);
  v_free(_v1);
  v_free(_v2);
  v_free(_scale);

  iv_free(_CTC_degree);
  iv_free(_CTC_neigh_start);
  iv_free(_CTC_neighs);
}


//--------------------------------------------------------------------------
SPMAT *Hqp_IpSpSC::sub_CTC(const PERM *px, SPMAT *Q)
// return Q - _CT * diag(_zw) * _CT'
// read:  _CT, _zw
// write: _scale
{
  int i, j, j_idx, j_end;
  int qi, qj, qj_idx;
  SPROW *crow, *qrow;
  Real sum, val;
  IVEC neigh_header;
  IVEC *neigh = &neigh_header;

  assert(Q->n == Q->m);
  assert((int)px->size == Q->m && Q->m >= _n);

  if (!Q->flag_diag)
    sp_diag_access(Q);

  neigh_header.max_dim = 0;

  crow = _CT->row;
  for (i=0; i<_n; i++, crow++) {

    qrow = Q->row + px->pe[i];
    if (crow->len <= 0) {
      val = qrow->elt[qrow->diag].val;
      _scale->ve[i] = min(1.0, sqrt(-1.0 / val));
    }
    else {

      // calculate diagonal entry
      sum = sprow_inprod(crow, _zw, crow);
      j_idx = qrow->diag;
      val = qrow->elt[j_idx].val -= sum;
      _scale->ve[i] = min(1.0, sqrt(-1.0 / val));

      // calculate resting entries
      neigh->ive = _CTC_neighs->ive + _CTC_neigh_start->ive[i];
      neigh->dim = _CTC_neigh_start->ive[i + 1] - _CTC_neigh_start->ive[i];
      j_end = neigh->dim;
      for (j_idx = 0; j_idx < j_end; j_idx++) {
	j = neigh->ive[j_idx];
	if (j < i) {
	  sum = sprow_inprod(crow, _zw, _CT->row + j);
	  qi = px->pe[i];
	  qj = px->pe[j];

	  // substract sum from Qij or Qji (entry from upper part only)

	  if (qi < qj) {
	    qrow = Q->row + qi;
	    qj_idx = sprow_idx(qrow, qj);
	    if (qj_idx < 0) {
	      // the structure must already have been allocated in init()
	      m_error(E_INTERN, "Hqp_IpSpSC");
	    }
	    else {
	      qrow->elt[qj_idx].val -= sum;
	    }
	  }
	  else {
	    qrow = Q->row + qj;
	    qj_idx = sprow_idx(qrow, qi);
	    if (qj_idx < 0) {
	      // the structure must already have been allocated in init()
	      m_error(E_INTERN, "Hqp_IpSpSC");
	    }
	    else {
	      qrow->elt[qj_idx].val -= sum;
	    }
	  }
	}
      }
    }
  }

  return Q;  
}

//--------------------------------------------------------------------------
void Hqp_IpSpSC::init(const Hqp_Program *qp)
{

  SPROW *r1, *r2;
  int i, j;
  Real sum;

  _n = qp->c->dim;
  _me = qp->b->dim;
  _m = qp->d->dim;

  // reallocations

  _piv_QCVC = px_resize(_piv_QCVC, _n);
  _piv_AQCVCA = px_resize(_piv_AQCVCA, _me);

  px_ident(_piv_QCVC);

  _zw = v_resize(_zw, _m);
  _v1 = v_resize(_v1, _n);
  _v2 = v_resize(_v2, _me);
  _scale = v_resize(_scale, _n);

  // store C' for further computations
  // analyze structure of C'*C

  _CT = sp_get(_n, _m, 10);
  _CT = sp_transp(qp->C, _CT);

  sp_ones(_CT);
  v_ones(_zw);

  _Q = sp_get(_n, _n, 10);

  r1 = _CT->row;
  for (i=0; i<_n; i++, r1++) {
    r2 = r1;
    for (j=i; j<_n; j++, r2++) {
      sum = sprow_inprod(r1, _zw, r2);
      if (sum != 0.0) {
	sp_set_val(_Q, i, j, sum);
	if (i != j)
	  sp_set_val(_Q, j, i, sum);
      }
    }
  }


  _CTC_degree = iv_resize(_CTC_degree, _n);
  _CTC_neigh_start = iv_resize(_CTC_neigh_start, _n + 1);
  _CTC_neighs = sp_rcm_scan(_Q, SMNULL, SMNULL,
			    _CTC_degree, _CTC_neigh_start, _CTC_neighs);

  // initialize structure of reduced qp

  _QCVC = sp_resize(_QCVC,_n, _n);
  _AT = sp_resize(_AT,_n,_me);
  _QCVCA = sp_resize(_QCVCA,_n,_me);
  _AQCVC = sp_resize(_AQCVC,_me,_n);
  _AQCVCA = sp_resize(_AQCVCA,_me,_me);
  _AQCVCA_fac = sp_resize(_AQCVCA_fac,_me,_me);

  // prepare iterations

  update(qp);

}

//--------------------------------------------------------------------------
void Hqp_IpSpSC::update(const Hqp_Program *qp)
{
  // prepare QCVC
  sp_zero(_Q);
  symsp_into_symsp(qp->Q, -1.0, _Q, _piv_QCVC, 0);

  // update _CT
  _CT = sp_transp(qp->C, _CT);
  _AT = sp_transp(qp->A, _AT);
}

//--------------------------------------------------------------------------
void Hqp_IpSpSC::factor(const Hqp_Program *qp, const VEC *z, const VEC *w)
{

  assert((int)z->dim == _m && (int)w->dim == _m);

  _QCVC = sp_copy3(_Q,_QCVC);

  v_slash(w, z, _zw);
  sub_CTC(_piv_QCVC, _QCVC);
  sp_smlt(_QCVC, -1.0, _QCVC);

  sp_compact(_QCVC, 0.0);
}

//--------------------------------------------------------------------------
Real Hqp_IpSpSC::solve(const Hqp_Program *qp, const VEC *z, const VEC *w,
		       const VEC *r1, const VEC *r2, const VEC *r3,
		       const VEC *r4, VEC *dx, VEC *dy, VEC *dz, VEC *dw)
{
  assert((int)r1->dim == _n && (int)dx->dim == _n);
  assert((int)r2->dim == _me && (int)dy->dim == _me);
  assert((int)r3->dim == _m && (int)dz->dim == _m);
  assert((int)r4->dim == _m && (int)dw->dim == _m);

  _v1 = v_resize(_v1, r1->dim);

  // temporary store (W^{-1}r_4 + ZW^{-1}r_3) in dz
  v_slash(w, r4, dw);
  v_star(_zw, r3, dz);
  v_add(dw, dz, dz);

  // temporary store (r1-CT*W^(-1)*(Z*r3+r4))
  sp_mv_mlt(_CT, dz, _v1);
  v_sub(r1, _v1, _v1);
  v_copy(_v1, dx);

  //spCHOLfac(_QCVC);
  spMODCHOLfac(_QCVC, _v1, _macheps);
  spCHOLsol(_QCVC, _v1, _v1);
  sp_mv_mlt(qp->A, _v1, _v2);
  v_add(_v2, r2, _v2);

  _QCVCA = spLTMsolve(_QCVC, _AT, _QCVCA);
  _AQCVC = sp_transp(_QCVCA, _AQCVC);

  if(!_piv_AQCVCA_flag) {
    _AQCVCA = sp_mmt2_mlt(_AQCVC, _QCVCA, _AQCVCA);
    sp_symrcm(_AQCVCA, _piv_AQCVCA);
    _piv_AQCVCA_flag = 1;
    pxinv_sprows(_piv_AQCVCA, _AQCVCA, _AQCVCA);
    pxinv_spcols(_piv_AQCVCA, _AQCVCA, _AQCVCA);
  }
  else {
    _AQCVCA = sp_mmt2_mlt_u(_AQCVC, _QCVCA, _AQCVCA);
    pxinv_sprows(_piv_AQCVCA, _AQCVCA, _AQCVCA);
    pxinv_spcols(_piv_AQCVCA, _AQCVCA, _AQCVCA);
  }

  pxinv_vec(_piv_AQCVCA, _v2, _v2);

  sp_compact(_AQCVCA, 0.0);
  _AQCVCA_fac = sp_copy3(_AQCVCA, _AQCVCA_fac);
  spMODCHOLfac(_AQCVCA_fac, _v2, _macheps);
  // spCHOLfac(_AQCVCA_fac);
  spCHOLsol(_AQCVCA_fac, _v2, dy);

  px_vec(_piv_AQCVCA, dy, dy);

  // compute dx
  sp_vm_mlt(qp->A, dy, _v1);
  v_sub(_v1, dx, _v1);
  spCHOLsol(_QCVC, _v1, dx);

  // compute dz 
  sp_mv_mlt(qp->C, dx, dw);
  v_star(_zw, dw, dw);
  v_sub(dz, dw, dz);

  // compute dw
  sv_mlt(-1.0, r3, dw);
  sp_mv_mltadd(dw, dx, qp->C, 1.0, dw);

  return residuum(qp, z, w, r1, r2, r3, r4, dx, dy, dz, dw);
}


//==========================================================================
