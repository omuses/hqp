/*
 * Hqp_IpRedPardiso.C -- class definition
 *
 * hl, 2006/11/22
 */

/*
    Copyright (C) 2006    Hartmut Linke

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

#include "Hqp_IpRedPardiso.h"

#include <assert.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <omp.h>

// extern int omp_get_max_threads();
/* PARDISO prototype. */
extern "C" {
#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif
extern int PARDISO
	(void *, int *, int *, int *, int *, int *,
	double *, int *, int *, int *, int *, int *,
	int *, double *, double *, int *);
}

extern "C" {
#include <meschach/sparse2.h>
}

#include <If_Int.h>
#include <If_Real.h>

#include "sprcm.h"
#include "Hqp_Program.h"

IF_CLASS_DEFINE("RedPardiso", Hqp_IpRedPardiso, Hqp_IpMatrix);  

//--------------------------------------------------------------------------
Hqp_IpRedPardiso::Hqp_IpRedPardiso()
{
  _n = _me = _m = 0;
  _sbw = -1;
  _tol = 1.0;
  _CT = SMNULL;
  _pivot = PNULL;
  _zw = VNULL;
  _scale = VNULL;
  _r12 = VNULL;
  _xy = VNULL;

  _pivot_strategy = 2;

  _iv = iv_resize(iv_get(1),0);
  _jv = iv_resize(iv_get(1),0);
  _ipivot = iv_resize(iv_get(1),0);	
  _v = v_resize(v_get(1),0);
  _v_raw = v_resize(v_get(1),0);

  _ifList.append(new If_Int("mat_sbw", &_sbw));
  _ifList.append(new If_Int("mat_pivot_strategy", &_pivot_strategy));
  _ifList.append(new If_Real("mat_tol", &_tol));

}

//--------------------------------------------------------------------------
Hqp_IpRedPardiso::~Hqp_IpRedPardiso()
{

  //  Termination and release of memory of Pardiso solver
  double ddum;
  int idum;
  int dim = _n + _me;

  _phase = -1;        // Release internal memory

  PARDISO (_pardiso_pt, &_maxfct, &_mnum, &_mtype, &_phase,
	   &dim, &ddum, _iv->ive, _jv->ive, &idum, &_nrhs,
	   _pardiso_parm, &_msglvl, &ddum, &ddum, &_error);

  sp_free(_CT);
  px_free(_pivot);
  v_free(_zw);
  v_free(_scale);
  v_free(_r12);
  v_free(_xy);

  v_free(_v);
  iv_free(_iv);
  iv_free(_jv);
  iv_free(_ipivot);

}

//--------------------------------------------------------------------------
void Hqp_IpRedPardiso::sub_CTC( )
{

  int i, j, j_idx, j_end, col_idx;
  SPROW *crow;
  Real sum;
 
  crow = _CT->row;
  for ( i = 0; i < _n; i++, crow++) {

    // determine index of diagonal entry
	j_idx = _iv->ive[i];
	j_idx--;

	// starting column entry of each row must be diagonal entry
	if(_jv->ive[j_idx]-1 != i)
	  m_error(E_INTERN, "Hqp_IpRedPardiso");

	if (crow->len > 0) {

      // calculate diagonal entry
      sum = sprow_inprod(crow, _zw, crow);
      _v->ve[j_idx] -= sum;

      // calculate resting entries
      j_end = _iv->ive[i+1]-1;

      for (j = j_idx+1; j < j_end; j++) {
	     col_idx = _jv->ive[j]-1;

        if (col_idx < _n) {
	      sum = sprow_inprod(crow, _zw, _CT->row + col_idx);
	      _v->ve[j] -= sum;
	    }
      }
    }
  }

  return;  
}

//--------------------------------------------------------------------------
void Hqp_IpRedPardiso::init(const Hqp_Program *qp)
{

  int i;

  assert(_pivot_strategy >= 0 && _pivot_strategy <= 2);

  _n = qp->c->dim;
  _me = qp->b->dim;
  _m = qp->d->dim;
  _dim = _n + _me;

  // reallocations

  _pivot = px_resize(_pivot, _dim);
  _ipivot = iv_resize(_ipivot, _dim);
  _zw = v_resize(_zw, _m);
  _scale = v_resize(_scale, _n);
  _r12 = v_resize(_r12, _dim);
  _xy = v_resize(_xy, _dim);

  // store C' for further computations
  // analyze structure of C'*C

  _CT = sp_transp(qp->C, _CT);

  // Setup Pardiso control parameters
  for (i = 0; i < 64; i++) {
    _pardiso_parm[i] = 0;
  }
  _pardiso_parm[0] = 1;        // No solver default

  // 0 -> MMD reordering,  2 -> Fill-in reordering from METIS
  _pardiso_parm[1] = _pivot_strategy;

  // Numbers of processors, value of OMP_NUM_THREADS
  // _pardiso_parm[2] = omp_get_max_threads();
  std::cout << "max number of threads: " << omp_get_max_threads() << std::endl;
  _pardiso_parm[2] = omp_get_max_threads();
  _pardiso_parm[3] = 0;        // No iterative-direct algorithm

  if(_pivot_strategy == 1)
	_pardiso_parm[4] = 1;      // user fill-in reducing permutation
  else
	_pardiso_parm[4] = 0;      // no user fill-in reducing permutation

  _pardiso_parm[5] = 0;        // Write solution into x
  _pardiso_parm[6] = 0;        // Not in use
  _pardiso_parm[7] = 0;        // Max numbers of iterative refinement steps
  _pardiso_parm[8] = 0;        // Not in use
  _pardiso_parm[9] = 8;        // Perturb the pivot elements with 1E-8
  _pardiso_parm[10] = 1;       // Use nonsymmetric permutation and scaling MPS
  _pardiso_parm[11] = 0;       // Not in use
  _pardiso_parm[12] = 1;       // Not in use
  _pardiso_parm[13] = 0;       // Output: Number of perturbed pivots
  _pardiso_parm[14] = 0;       // Not in use
  _pardiso_parm[15] = 0;       // Not in use
  _pardiso_parm[16] = 0;       // Not in use 
  _pardiso_parm[17] = 0;       // -1 -> Output: Number of nonzeros in the factor LU
  _pardiso_parm[18] = 0;       // -1 -> Output: Mflops for LU factorization
  _pardiso_parm[19] = 0;       // Output: Numbers of CG Iterations
  _pardiso_parm[20] = 1;       // allows for 1x1 and 2x2 pivotizations

  _maxfct = 1;          // Maximum number of numerical factorizations
  _mnum = 1;            // Which factorization to use
  _msglvl = 0;          // Print statistical information in file
  _error = 0;           // Initialize error flag
  _mtype = -2;          // Real symmetric matrix
  _nrhs = 1;            // Number of right hand sides

  for (i = 0; i < 64; i++) {
    _pardiso_pt[i] = 0;
  }

  _reinit = -1;

  // prepare iterations

  update(qp);
}

//--------------------------------------------------------------------------
void Hqp_IpRedPardiso::update(const Hqp_Program *qp)
{
  // prepare coefficient matrix in Pardiso format

  SPMAT  *AT, *QCTC;
  SPROW *r1, *r2;
  IVEC *degree, *neigh_start, *neighs;
  int i, j, j_end, j_idx, k, nnz_row, col;
  double ddum;
  int idum;
  Real sum;
  const row_elt *elt;

  // store C' for further computations
  // analyze structure of C'*C

  _CT = sp_transp(qp->C, _CT);
  sp_ones(_CT);
  v_set(_zw,1.0e-5);
  QCTC = sp_get(_n, _n, 10);
  r1 = _CT->row;
  for (i=0; i<_n; i++, r1++) {
    r2 = r1;
    for (j=i; j<_n; j++, r2++) {
      sum = sprow_inprod(r1, _zw, r2);
    if (sum != 0.0) 
	  sp_set_val(QCTC, i, j, sum);
	if ((i == j) && sum == 0.0)
	  sp_set_val(QCTC, i, i, 0.0);
    }
  }

  // initialize structure of reduced qp

  QCTC = sp_add(QCTC, qp->Q, QCTC);

  AT = sp_get(_n, _me, 10);
  AT = sp_transp(qp->A, AT);

  _nnz = 0;
  for (i = 0; i < _n; i++) {
	nnz_row = 0;
    nnz_row += QCTC->row[i].len;
    nnz_row += AT->row[i].len;

    if((_reinit == 0) && (nnz_row != (_iv->ive[i+1] - _iv->ive[i])))
      _reinit = 1;

	_nnz += nnz_row;

  }

   // determine RCM ordering
 if(_pivot_strategy == 1 && _reinit != 0) {
    degree = iv_get(_dim);
	neigh_start = iv_get(_dim + 1);
	neighs = sp_rcm_scan(QCTC, qp->A, SMNULL, degree, neigh_start, IVNULL);
	_pivot = sp_rcm_order(degree, neigh_start, neighs, _pivot);
	_sbw = sp_rcm_sbw(neigh_start, neighs, _pivot);
	iv_resize(_ipivot,_pivot->size);
	// switch to Fortran numbering schema
	for(i = 0; i < (int)_pivot->size; i++)
	  _ipivot->ive[i] = _pivot->pe[i]+1;
 }

  // printf("\n nnz: %d    reinit: %d\n",_nnz,_reinit);


 if(_reinit != 0) {
	iv_resize(_iv,_dim+1);
	iv_zero(_iv);

	iv_resize(_jv,_nnz+_me);
	iv_zero(_jv);

	v_resize(_v,_nnz+_me);
	v_zero(_v);

	v_resize(_v_raw,_nnz+_me);
	v_zero(_v_raw);
 }
	
 _iv->ive[_dim] = _nnz+_me;

  // fill up data

  k = 0;

  for( i = 0; i < _n; i++) {

    _iv->ive[i] = k;
    elt = QCTC->row[i].elt;
    j_end = QCTC->row[i].len;
    for (j_idx=0; j_idx < j_end; j_idx++, elt++) {
      if (elt->col >= i) {
		col = _jv->ive[k] - 1;
	    // look for changes of the hessian structure
	    // -> reinit pardiso solver
	    if((_reinit == 0) && (elt->col != col))
	      _reinit = 1;
	_jv->ive[k] = elt->col;
	_v_raw->ve[k++] = -elt->val;
	  }
    }

    elt = AT->row[i].elt;
    j_end = AT->row[i].len;
    for (j_idx=0; j_idx < j_end; j_idx++, elt++) {
	  _jv->ive[k] = elt->col+_n;
	  _v_raw->ve[k++] = elt->val;
    }

  }

  for( i = _n; i < _n+_me; i++) {
	_iv->ive[i] = k;
	_jv->ive[k] = i;
	_v_raw->ve[k++] = 0.0;
  }

  v_copy(_v_raw,_v);

  // switch to Fortran numerbering schema
  for( i = 0; i < (int) _iv->dim; i++)
	_iv->ive[i] = _iv->ive[i]+1;

  for( i = 0; i < (int) _jv->dim; i++)
	_jv->ive[i] = _jv->ive[i]+1;

  _phase = 11;

  if(_reinit != 0) {
    printf("Reinitialization of RedPardiso solver!\n");
	PARDISO (_pardiso_pt, &_maxfct, &_mnum, &_mtype, &_phase,
			 &_dim, _v->ve, _iv->ive, _jv->ive, _ipivot->ive, &_nrhs,
			 _pardiso_parm, &_msglvl, &ddum, &ddum, &_error);
  }

  _reinit = 0;

  if (_error != 0) {
    printf("\nERROR during symbolic factorization: %d", _error);
    exit(1);
  }

  // fill up data again

  v_zero(_v_raw);
  QCTC = sp_zero(QCTC);
  QCTC = sp_add(QCTC, qp->Q, QCTC);

  k = 0;

  for( i = 0; i < _n; i++) {

    elt = QCTC->row[i].elt;
    j_end = QCTC->row[i].len;
    for (j_idx=0; j_idx < j_end; j_idx++, elt++) {
      if (elt->col >= i) {
	_v_raw->ve[k++] = -elt->val;
	  }
    }

    elt = AT->row[i].elt;
    j_end = AT->row[i].len;
    for (j_idx=0; j_idx < j_end; j_idx++, elt++) {
	  _v_raw->ve[k++] = elt->val;
    }

  }

  v_copy(_v_raw,_v);

  sp_free(QCTC);
  sp_free(AT);

  // update _CT
  _CT = sp_transp(qp->C, _CT);
}

//--------------------------------------------------------------------------
void Hqp_IpRedPardiso::factor(const Hqp_Program *qp,
			    const VEC *z, const VEC *w)
{
  assert((int)z->dim == _m && (int)w->dim == _m);

  double  ddum;
  int     idum;

  v_copy(_v_raw,_v);

  // augment _J
  v_slash(w, z, _zw);
  sub_CTC( );

  // PARDISO factor
  _phase = 22;

  clock_t start, finish;
  double  duration;

  start = clock();

  PARDISO (_pardiso_pt, &_maxfct, &_mnum, &_mtype, &_phase,
	   &_dim, _v->ve, _iv->ive, _jv->ive, _ipivot->ive, &_nrhs,
	   _pardiso_parm, &_msglvl, &ddum, &ddum, &_error);

  finish = clock();
  duration = (double)(finish - start) / CLOCKS_PER_SEC;
  // printf( "Hqp_IpRedPardiso::factor():  %2.8f seconds\n", duration );

  if (_error != 0) {
	printf("\nERROR during numerical factorization: %d", _error);
	exit(2);
  }

}

//--------------------------------------------------------------------------
void Hqp_IpRedPardiso::step(const Hqp_Program *qp, const VEC *z, const VEC *w,
			  const VEC *r1, const VEC *r2, const VEC *r3,
			  const VEC *r4, VEC *dx, VEC *dy, VEC *dz, VEC *dw)
{

  int     idum;
  VEC     v;

  assert((int)r1->dim == _n && (int)dx->dim == _n);
  assert((int)r2->dim == _me && (int)dy->dim == _me);
  assert((int)r3->dim == _m && (int)dz->dim == _m);
  assert((int)r4->dim == _m && (int)dw->dim == _m);

  // augment, copy, scale and permutate [r1;r2;r3] into _r12
  // calculate, permutate and scale x
  
  // augment r1
  // temporary store (W^{-1}r_4 + ZW^{-1}r_3) in dz
  v_part(_r12, 0, _n, &v);
  v_slash(w, r4, dw);
  v_star(_zw, r3, dz);
  v_add(dw, dz, dz);
  sp_mv_mlt(_CT, dz, &v);
  v_sub(r1, &v, &v);
  v_copy(r2, v_part(_r12, _n, _me, &v));

  // PARDISO solve

  _phase = 33;

  clock_t start, finish;
  double  duration;

  start = clock();


  PARDISO (_pardiso_pt, &_maxfct, &_mnum, &_mtype, &_phase,
	   &_dim, _v->ve, _iv->ive, _jv->ive, _ipivot->ive, &_nrhs,
	   _pardiso_parm, &_msglvl, _r12->ve, _xy->ve, &_error);

  finish = clock();
  duration = (double)(finish - start) / CLOCKS_PER_SEC;
  // printf( "Hqp_IpRedPardiso::solve():  %2.8f seconds\n", duration );

  if (_error != 0) {
	printf("\nERROR during solution: %d", _error);
	exit(3);
  }

  v_copy(v_part(_xy, 0, _n, &v), dx);
  v_copy(v_part(_xy, _n, _me, &v), dy);

  sp_vm_mlt(_CT, dx, dw);
  v_star(_zw, dw, dw);
  v_sub(dz, dw, dz);

  // calculate dw

  sv_mlt(-1.0, r3, dw);
  sp_mv_mltadd(dw, dx, qp->C, 1.0, dw);
  // usage of _CT is numerically worse!
  //sp_vm_mltadd(dw, dx, _CT, 1.0, dw);
}


//==========================================================================
