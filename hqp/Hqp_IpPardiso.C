/*
 * Hqp_IpPardiso.C -- class definition
 *
 * hl, 2006/11/22
 *
 * rf, 2010/08/14: add dynamic load of Pardiso solver
 *
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

#include "Hqp_IpPardiso.h"
#include "Hqp_Program.h"

#include <assert.h>
#include <math.h>
#include <time.h>

//#include <omp.h> // Note: need automatic configuration if available
//extern "C" int omp_get_max_threads();

extern "C" {
#include <meschach/sparse2.h>
}
#include "sprcm.h"

#include <If_Int.h>
#include <If_Real.h>
#include <If_String.h>

#define GET_SET_CB(vartype, prefix, name) \
  GET_CB(vartype, prefix, name), \
  IF_SET_CB(vartype, Hqp_IpPardiso, set_##name)

#define GET_CB(vartype, prefix, name) \
  prefix#name, \
  IF_GET_CB(vartype, Hqp_IpPardiso, name)

IF_CLASS_DEFINE("Pardiso", Hqp_IpPardiso, Hqp_IpMatrix);

//--------------------------------------------------------------------------
Hqp_IpPardiso::Hqp_IpPardiso()
{
  _n = _me = _m = 0;
  _sbw = -1;
  _tol = 1.0;
  _pivot = PNULL;
  _scale = VNULL;
  _r123 = VNULL;
  _xyz = VNULL;

  _pardiso_libname = strdup("mkl_core"); // MKL 10.2.5.035
  _pardiso_fname = strdup("mkl_pds_pardiso"); // MKL 10.2.5.035
  _pardiso_fp = NULL;
  _nparallel = 2; // Note: could use omp_get_max_threads();

  _reinit = -1;

  _iv = iv_resize(iv_get(1),0);
  _jv = iv_resize(iv_get(1),0);
  _v = v_resize(v_get(1),0);
  _v_raw = v_resize(v_get(1),0);

  _ifList.append(new If_Int("mat_sbw", &_sbw));
  _ifList.append(new If_Real("mat_tol", &_tol));
  _ifList.append(new If_Int(GET_SET_CB(int, "mat_", nparallel)));
  _ifList.append(new If_String(GET_SET_CB(const char *, "mat_", pardiso_libname)));
  _ifList.append(new If_String(GET_SET_CB(const char *, "mat_", pardiso_fname)));
}

//--------------------------------------------------------------------------
Hqp_IpPardiso::~Hqp_IpPardiso()
{

  free_pardiso();

  px_free(_pivot);
  v_free(_scale);
  v_free(_r123);
  v_free(_xyz);

  v_free(_v);
  v_free(_v_raw);
  iv_free(_iv);
  iv_free(_jv);

  free(_pardiso_libname);
  free(_pardiso_fname);
}

//--------------------------------------------------------------------------
void Hqp_IpPardiso::free_pardiso()
{

  //  Termination and release of memory of Pardiso solver
  double ddum;
  int idum;

  _phase = -1;        // Release internal memory

  if (_pardiso_fp) {
    (*_pardiso_fp) (_pardiso_pt, &_maxfct, &_mnum, &_mtype, &_phase,
                    &_dim, &ddum, _iv->ive, _jv->ive, &idum, &_nrhs,
                    _pardiso_parm, &_msglvl, &ddum, &ddum, &_error);
  }
}

//--------------------------------------------------------------------------
void Hqp_IpPardiso::reinit_pardiso()
{

  double ddum;
  int idum;

  if (!_pardiso_fp)
    m_error(E_NULL, "Hqp_IpPardiso::reinit_pardiso");

  if(_reinit == 1)
    free_pardiso();

  _phase = 11;

  (*_pardiso_fp)(_pardiso_pt, &_maxfct, &_mnum, &_mtype, &_phase,
                 &_dim, _v->ve, _iv->ive, _jv->ive, &idum, &_nrhs,
                 _pardiso_parm, &_msglvl, &ddum, &ddum, &_error);

  if (_error != 0) {
    fprintf(stderr, "\nERROR during symbolic factorization: %d", _error);
    m_error(E_CONV, "Pardiso symbolic factorization");
  }

  _reinit = 0;

  return;

}

//--------------------------------------------------------------------------
void Hqp_IpPardiso::set_pardiso_libname(const char *value)
{
  if (strcmp(value, _pardiso_libname) != 0) {
    free(_pardiso_libname);
    _pardiso_libname = strdup(value);
    _pardiso_fp = NULL; // needs to be loaded in Hqp_IpPardiso::init
    _dl.close();
  }
}

//--------------------------------------------------------------------------
void Hqp_IpPardiso::set_pardiso_fname(const char *value)
{
  if (strcmp(value, _pardiso_fname) != 0) {
    free(_pardiso_fname);
    _pardiso_fname = strdup(value);
    _pardiso_fp = NULL; // needs to be loaded in Hqp_IpPardiso::init
  }
}

//--------------------------------------------------------------------------
void Hqp_IpPardiso::init(const Hqp_Program *qp)
{
  int i;

  // load solver
  if (!_pardiso_fp) {
    _dl.open(_pardiso_libname);
    _pardiso_fp = (pardiso_ft*)_dl.symbol(_pardiso_fname);
  }
  if (!_pardiso_fp)
    m_error(E_NULL, "Hqp_IpPardiso::init");

  _n = qp->c->dim;
  _me = qp->b->dim;
  _m = qp->d->dim;
  _dim = _n + _me + _m;

  // reallocations

  _pivot = px_resize(_pivot, _dim);
  _scale = v_resize(_scale, _m);
  _r123 = v_resize(_r123, _dim);
  _xyz = v_resize(_xyz, _dim);

  // Setup Pardiso control parameters
  for (i = 0; i < 64; i++) {
    _pardiso_parm[i] = 0;
  }
  _pardiso_parm[0] = 1;        // No solver default
  _pardiso_parm[1] = 2;        // 1 -> MMD reordering algorithm; 2 -> Fill-in reordering from METIS
  _pardiso_parm[2] = _nparallel;// Number of processors
  _pardiso_parm[3] = 0;        // No iterative-direct algorithm
  _pardiso_parm[4] = 0;        // No user fill-in reducing permutation
  _pardiso_parm[5] = 0;        // Write solution into x
  _pardiso_parm[6] = 0;        // Not in use
  _pardiso_parm[7] = 0;        // Max numbers of iterative refinement steps
  _pardiso_parm[8] = 0;        // Not in use
  _pardiso_parm[9] = 8;        // Perturb the pivot elements with 1E-8
  _pardiso_parm[10] = 2;       // Use nonsymmetric permutation and scaling MPS
  _pardiso_parm[11] = 0;       // Not in use
  _pardiso_parm[12] = 1;       // Not in use
  _pardiso_parm[13] = 0;       // Output: Number of perturbed pivots
  _pardiso_parm[14] = 0;       // Not in use
  _pardiso_parm[15] = 0;       // Not in use
  _pardiso_parm[16] = 0;       // Not in use 
  _pardiso_parm[17] = 0;      // Output: Number of nonzeros in the factor LU
  _pardiso_parm[18] = 0;      // Output: Mflops for LU factorization
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

  // fill up data

  update(qp);
}

//--------------------------------------------------------------------------
void Hqp_IpPardiso::update(const Hqp_Program *qp)
{

  SPMAT  *AT, *CT;
  int i, j_end, j_idx, k, nnz_row, col;
  const row_elt *elt;
  double ddum;
  int idum;

  AT = sp_get(_n, _me, 10);
  AT = sp_transp(qp->A, AT);

  CT = sp_get(_n, _m, 10);
  CT = sp_transp(qp->C, CT);

  _nnz = 0;
  for (i = 0; i < _n; i++) {
    nnz_row = 0;
    nnz_row += qp->Q->row[i].len;
    nnz_row += AT->row[i].len;
    nnz_row += CT->row[i].len;

    if((_reinit == 0) && (nnz_row != (_iv->ive[i+1] - _iv->ive[i])))
      _reinit = 1;
    _nnz += nnz_row;

  }

  // printf("\n nnz: %d    reinit: %d\n",_nnz,_reinit);

  // number of non-zero entrys in qp-matrix has changed
  // this should only occur according to changes of the hessian
  // structure  -> reinit pardiso solver

  if(_reinit != 0) {
    iv_resize(_iv,_dim+1);
	iv_zero(_iv);

	iv_resize(_jv,_nnz+_me+_m);
	iv_zero(_jv);

	_iv->ive[_dim] = _nnz+_me+_m;

	v_resize(_v,_nnz+_me+_m);
	v_zero(_v);

	v_resize(_v_raw,_nnz+_me+_m);
	v_zero(_v_raw);
  }
  else 
	_iv->ive[_dim] = _nnz+_me+_m;

  // fill up data

  k = 0;

  for( i = 0; i < _n; i++) {

    _iv->ive[i] = k;
    elt = qp->Q->row[i].elt;
    j_end = qp->Q->row[i].len;
    for (j_idx=0; j_idx < j_end; j_idx++, elt++) {
      if (elt->col >= i) {
		  col = _jv->ive[k] - 1;
	// look for changes of the hessian structure
	// -> reinit pardiso solver
	if((_reinit == 0) && (elt->col != col)) {
	  _reinit = 1;
	  printf("row:  %d  col: %d   col_old: %d \n",i,elt->col, col);
	}
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

    elt = CT->row[i].elt;
    j_end = CT->row[i].len;
    for (j_idx=0; j_idx < j_end; j_idx++, elt++) {
	  _jv->ive[k] = elt->col+_n+_me;
	  _v_raw->ve[k++] = elt->val;
    }
  }

  for( i = _n; i < _n+_me; i++) {
	_iv->ive[i] = k;
	_jv->ive[k] = i;
	_v_raw->ve[k++] = 0.0;
  }

   for( i = _n+_me; i < _n+_me+_m; i++) {
	_iv->ive[i] = k;
	_jv->ive[k] = i;
	_v_raw->ve[k++] = 1.0;
  }

  v_copy(_v_raw,_v);

  // switch to Fortran numerbering schema
  for( i = 0; i < (int) _iv->dim; i++)
	_iv->ive[i] = _iv->ive[i]+1;

  for( i = 0; i < (int) _jv->dim; i++)
	_jv->ive[i] = _jv->ive[i]+1;

  _phase = 11;

  if(_reinit != 0) {
    printf("Reinitialization of Pardiso solver!\n");
    if (!_pardiso_fp)
      m_error(E_NULL, "Hqp_IpPardiso::init");
    (*_pardiso_fp) (_pardiso_pt, &_maxfct, &_mnum, &_mtype, &_phase,
                    &_dim, _v->ve, _iv->ive, _jv->ive, &idum, &_nrhs,
                    _pardiso_parm, &_msglvl, &ddum, &ddum, &_error);
  }

  _reinit = 0;

  if (_error != 0) {
    fprintf(stderr, "\nERROR during symbolic factorization: %d", _error);
    m_error(E_CONV, "Pardiso symbolic factorization");
  }

  sp_free(CT);
  sp_free(AT);

}

//--------------------------------------------------------------------------
void Hqp_IpPardiso::factor(const Hqp_Program *, const VEC *z, const VEC *w)
{
  assert((int)z->dim == _m && (int)w->dim == _m);

  if (!_pardiso_fp)
    m_error(E_NULL, "Hqp_IpPardiso::factor");

  double  ddum, wz;
  int     i, i0, idum;

  v_copy(_v_raw,_v);

  // insert slacks
  i0 = _iv->ive[_n+_me]-1;
  for (i = 0; i < _m; i++) {
    wz = w->ve[i] / z->ve[i];
    _v->ve[i+i0] = wz;
  }

  // PARDISO factor
  _phase = 22;

  clock_t start, finish;
  double  duration;

  start = clock();

  (*_pardiso_fp) (_pardiso_pt, &_maxfct, &_mnum, &_mtype, &_phase,
                  &_dim, _v->ve, _iv->ive, _jv->ive, &idum, &_nrhs,
                  _pardiso_parm, &_msglvl, &ddum, &ddum, &_error);

  finish = clock();
  duration = (double)(finish - start) / CLOCKS_PER_SEC;
  // printf( "Hqp_IpRedPardiso::factor():  %2.8f seconds\n", duration );

  if (_error != 0) {
    fprintf(stderr, "\nERROR during numerical factorization: %d", _error);
    m_error(E_CONV, "Pardiso numerical factorization");
  }

}

//--------------------------------------------------------------------------
void Hqp_IpPardiso::step(const Hqp_Program *qp, const VEC *z, const VEC *w,
		       const VEC *r1, const VEC *r2, const VEC *r3,
		       const VEC *r4, VEC *dx, VEC *dy, VEC *dz, VEC *dw)
{

  int idum;
  VEC v;

  assert((int)r1->dim == _n && (int)dx->dim == _n);
  assert((int)r2->dim == _me && (int)dy->dim == _me);
  assert((int)r3->dim == _m && (int)dz->dim == _m);
  assert((int)r4->dim == _m && (int)dw->dim == _m);

  if (!_pardiso_fp)
    m_error(E_NULL, "Hqp_IpPardiso::step");

  // copy, scale and permutate [r1;r2;r3] into _r123
  // calculate, permutate and scale x
  
  v_copy(r1, v_part(_r123, 0, _n, &v));
  v_copy(r2, v_part(_r123, _n, _me, &v));
  v_part(_r123, _n+_me, _m, &v);
  v_slash(z, r4, &v);
  v_add(&v, r3,	&v);

  _phase = 33;

  clock_t start, finish;
  double  duration;

  start = clock();

  (*_pardiso_fp) (_pardiso_pt, &_maxfct, &_mnum, &_mtype, &_phase,
                  &_dim, _v->ve, _iv->ive, _jv->ive, &idum, &_nrhs,
                  _pardiso_parm, &_msglvl, _r123->ve, _xyz->ve, &_error);

  finish = clock();
  duration = (double)(finish - start) / CLOCKS_PER_SEC;
  // printf( "Hqp_IpRedPardiso::solve():  %2.8f seconds\n", duration );

  if (_error != 0) {
    fprintf(stderr, "\nERROR during solution: %d", _error);
    m_error(E_CONV, "Pardiso solution");
  }

  v_copy(v_part(_xyz, 0, _n, &v), dx);
  v_copy(v_part(_xyz, _n, _me, &v), dy);
  v_copy(v_part(_xyz, _n+_me, _m, &v), dz);

  // calculate dw

  sv_mlt(-1.0, r3, dw);
  sp_mv_mltadd(dw, dx, qp->C, 1.0, dw);
}


//==========================================================================
