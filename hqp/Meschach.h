/*
 * Meschach.h -- Meschach declarations and some extensions
 *
 * rf, 8/18/94
 */

/*
    Copyright (C) 1994--2000  Ruediger Franke and Eckhard Arnold

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

#ifndef Meschach_H
#define Meschach_H

#include <assert.h>

#if defined(_MSC_VER)
// Include math.h prior to matrix.h
// as matrix.h would include it otherwise
// and as extern "C" specification does not work for MSC's math.h.
#include <math.h>
#endif

extern "C" {
#include <matrix.h>
#include <sparse.h>
}

extern IVEC *iv_set(IVEC *, int);
extern IVEC *iv_part(IVEC *iv, int offs, int dim, IVEC *header);
extern IVEC *iv_expand(IVEC *iv, int nel, int granul);

extern VEC *v_set(VEC *, Real);
extern VEC *v_part(VEC *v, int offs, int dim, VEC *header);
extern VEC *v_expand(VEC *v, int nel, int granul);
extern VEC *bd_mv_mlt(const BAND *, const VEC *, VEC *);

extern MAT *m_mltadd(const MAT *, const MAT *, Real, MAT *);

extern Real sp_norm_inf(SPMAT *);
extern SPMAT *sp_copy3(const SPMAT *, SPMAT *);
extern int sp_update_val(SPMAT *, int, int, Real);
extern void sp_insert_mat(SPMAT *dst, int i_offs, int j_offs, const MAT *src);
extern void symsp_insert_symmat(SPMAT *dst, int offs, const MAT *src);
extern void sp_update_mat(SPMAT *dst, int i_offs, int j_offs, const MAT *src);
extern void sp_extract_mat(const SPMAT *src, int i_offs, int j_offs, MAT *dst);
extern void symsp_extract_mat(const SPMAT *src, int offs, MAT *dst);
extern void sp_insert_mrow(SPMAT *dst, int i_offs, int j_offs,
			   const MAT *src, int i);
extern void sp_update_mrow(SPMAT *dst, int i_offs, int j_offs,
			   const MAT *src, int i);
extern void sp_extract_mrow(const SPMAT *src, int i_offs, int j_offs,
			    MAT *dst, int i);
extern SPMAT *sp_ident(SPMAT *);
extern SPMAT *sp_ones(SPMAT *);
extern VEC *sp_mv_mltadd(const VEC *v1, const VEC *v2, const SPMAT *A,
			 Real alpha, VEC *out);
extern VEC *sp_vm_mltadd(const VEC *v1, const VEC *v2, const SPMAT *A,
			 Real alpha, VEC *out);
extern VEC *sp_mv_symmlt(SPMAT *A, const VEC *v, VEC *out);

extern Real sprow_inprod(const SPROW *r1, const VEC *inner, const SPROW *r2);
extern void sprow_zero(SPROW *row);

extern SPMAT *spLUfactor2(SPMAT *A, PERM *px);
extern SPMAT *sp_transp(const SPMAT *, SPMAT *);

/*
 * BKP factor and solve routines
 */

extern MAT *matBKPfactor(MAT *A, PERM *pivot);
extern VEC *matBKPsolve(const MAT *A, const PERM *pivot,
			const VEC *b, VEC *x);

extern BAND *bdBKPfactor(BAND *A, PERM *pivot, PERM *relief);
extern VEC *bdBKPsolve(const BAND *A, const PERM *pivot, const PERM *relief,
		       const VEC *b, VEC *x);

extern SPMAT *spBKPfactor(SPMAT *, PERM *pivot, Real tol);
extern VEC *spBKPsolve(const SPMAT *, const PERM *pivot,
		       const VEC *b, VEC *x);

/*
 * routines for copying sparse matrices into sparse/band matrices
 *  -- "sym" means, that only the upper part of a symmetric matrix is filled
 */

extern void sp_into_sp(const SPMAT *src, Real s, SPMAT *dst,
		       const PERM *px, int i_offs, int j_offs);
extern void spT_into_sp(const SPMAT *src, Real s, SPMAT *dst,
			const PERM *px, int i_offs, int j_offs);
extern void sp_into_symsp(const SPMAT *src, Real s, SPMAT *dst,
			  const PERM *px, int i_offs, int j_offs);
extern void symsp_into_symsp(const SPMAT *src, Real s, SPMAT *dst,
			     const PERM *px, int offs);
extern void spT_into_symsp(const SPMAT *src, Real s, SPMAT *dst,
			   const PERM *px, int i_offs, int j_offs);

extern void sp_into_bd(const SPMAT *sp, Real s, BAND *bd,
		       const PERM *px, int i_offs, int j_offs);
extern void spT_into_bd(const SPMAT *sp, Real s, BAND *bd,
			const PERM *px, int i_offs, int j_offs);

// overload operators for Meschach data structures (esp. operator [])
// idea: 
//   - define for any STRUCT a class STRUCTP
//   - data representation of STRUCTP is compatible to STRUCT*
//   - overload operators as needed
// a C synonym can be thought as: 
//   typedef STRUCT* STRUCTP
// bugs:
//   - Meschach prototypes have no const's --> bad const handling
//     (VEC *v_add(VEC *a, VEC *b, VEC *out) instead of
//      VEC *v_add(const VEC *a, const VEC *b, VEC *out)

#ifdef DEBUG
#define MESCH_BOUNDS_CHECK(i, first_index, next_index) \
   assert((int)(first_index) <= (int)(i) && (int)(i) < (int)(next_index))
#else
#define MESCH_BOUNDS_CHECK(i, first_index, next_index)
#endif

#ifdef DEBUG
#define MESCH_NULL_CHECK(ptr) \
   assert (ptr != NULL)
#else
#define MESCH_NULL_CHECK(ptr)
#endif

/** Wrapper for VEC*. */
class VECP {

 protected:
  VEC *_v;

 public:
  // constructors and assignments
  VECP() {_v = VNULL;}
  VECP(VEC *cv) {_v = cv;}

  VEC *operator = (VEC *nv) {_v = nv; return _v;}

  // overloaded operators
  Real &operator [] (int j)
    {
      MESCH_NULL_CHECK(_v);
      MESCH_BOUNDS_CHECK(j, 0, _v->dim);
      return _v->ve[j];
    }
  VEC *operator -> () {return _v;}
  operator VEC*() {return _v;}
  // additional conversion to avoid warning of gcc2.95
  operator const VEC*() {return _v;} 

  // operators for const VECP
  const Real &operator [] (int j) const
    {
      MESCH_NULL_CHECK(_v);
      MESCH_BOUNDS_CHECK(j, 0, _v->dim);
      return _v->ve[j];
    }
  const VEC *operator -> () const {return _v;}
  operator const VEC*() const {return _v;}
};

/** Wrapper for IVEC*. */
class IVECP {

 protected:
  IVEC *_v;

 public:
  // constructors and assignment
  IVECP() {_v = IVNULL;}
  IVECP(IVEC *cv) {_v = cv;}

  IVEC *operator = (IVEC *nv) {_v = nv; return _v;}

  // operators
  int &operator [] (int j)
    {
      MESCH_NULL_CHECK(_v);
      MESCH_BOUNDS_CHECK(j, 0, _v->dim);
      return _v->ive[j];
    }
  IVEC *operator -> () {return _v;}
  operator IVEC*() {return _v;}
  // additional conversion to avoid warning of gcc2.95
  operator const IVEC*() {return _v;}

  // operators for const IVECP
  const int &operator [] (int j) const
    {
      MESCH_NULL_CHECK(_v);
      MESCH_BOUNDS_CHECK(j, 0, _v->dim);
      return _v->ive[j];
    }
  const IVEC *operator -> () const {return _v;}
  operator const IVEC*() const {return _v;}
};

/** Wrapper for PERM*. */
class PERMP {

 protected:
  PERM *_v;

 public:
  // constructors and assignment
  PERMP() {_v = PNULL;}
  PERMP(PERM *cv) {_v = cv;}

  PERM *operator = (PERM *nv) {_v = nv; return _v;}

  // operators
  u_int &operator [] (int j)
    {
      MESCH_NULL_CHECK(_v);
      MESCH_BOUNDS_CHECK(j, 0, _v->size);
      return _v->pe[j];
    }
  PERM *operator -> () {return _v;}
  operator PERM*() {return _v;}
  // additional conversion to avoid warning of gcc2.95
  operator const PERM*() {return _v;}

  // operators for const PERMP
  const u_int &operator [] (int j) const
    {
      MESCH_NULL_CHECK(_v);
      MESCH_BOUNDS_CHECK(j, 0, _v->size);
      return _v->pe[j];
    }
  const PERM *operator -> () const {return _v;}
  operator const PERM*() const {return _v;}
};

/** Wrapper for row in MAT*. */
class MATROWP {

 protected:
  Real *_row;
  int _dim;

 public:
  MATROWP(Real *row, int dim) {_row = row; _dim = dim;}

  Real &operator [] (int j)
    {
      MESCH_BOUNDS_CHECK(j, 0, _dim);
      return _row[j];
    }
  operator Real*() {return _row;}

  const Real &operator [] (int j) const
    {
      MESCH_BOUNDS_CHECK(j, 0, _dim);
      return _row[j];
    }
  operator const Real*() const {return _row;}
};

/** Wrapper for MAT*. */
class MATP {

 protected:
  MAT *_m;

 public:
  // constructors and assignment
  MATP() {_m = MNULL;}
  MATP(MAT *cm) {_m = cm;}

  MAT *operator = (MAT *nm) {_m = nm; return _m;}

  // operators
# ifdef DEBUG
  MATROWP operator [] (int i)
    {
      MESCH_NULL_CHECK(_m);
      MESCH_BOUNDS_CHECK(i, 0, _m->m);
      return MATROWP(_m->me[i], _m->n);
    }
# else
  Real *operator [] (int i)
    {
      return _m->me[i];
    }
# endif
  MAT *operator -> () {return _m;}
  operator MAT*() {return _m;}
  // additional conversion to avoid warning of gcc2.95
  operator const MAT*() {return _m;}

  // operators for const MATP
# ifdef DEBUG
  const MATROWP operator [] (int i) const
    {
      MESCH_NULL_CHECK(_m);
      MESCH_BOUNDS_CHECK(i, 0, _m->m);
      return MATROWP(_m->me[i], _m->n);
    }
# else
  const Real *operator [] (int i) const
    {
      return _m->me[i];
    }
# endif
  const MAT *operator -> () const {return _m;}
  operator const MAT*() const {return _m;}
};

typedef SPMAT* SPMATP;
typedef BAND* BANDP;

#endif
