/**
 * @file Meschach.h
 *   Meschach declarations and some extensions
 *
 * rf, 8/18/94
 */

/*
    Copyright (C) 1994--2002  Ruediger Franke and Eckhard Arnold

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

/** Avoid multiple inclusion */
#if !defined(Meschach_H)
#define Meschach_H

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
#include <limits> // for NaN
#if !defined(max)
#define	max(a,b)	((a) > (b) ? (a) : (b))
#endif
#if !defined(min)
#define	min(a,b)	((a) < (b) ? (a) : (b))
#endif
#if !defined(max)
#define	max(a,b)	((a) > (b) ? (a) : (b))
#endif
#include <string.h> // for sscanf

#include <assert.h>

#include <meschach/matrix.h>
#include <meschach/sparse.h>
#include <meschach/addon2_hqp.h>

/**
 * @name BKP factor and solve routines
 *   Reference: J.R.Bunch, L.Kaufman, and B.N.Parlett:
 *      Decomposition of a Symmetric Matrix, Numer Math 27, 95--109 (1976).
 */
//@{
extern MAT *matBKPfactor(MAT *A, PERM *pivot);
extern VEC *matBKPsolve(const MAT *A, const PERM *pivot,
			const VEC *b, VEC *x);

extern BAND *bdBKPfactor(BAND *A, PERM *pivot, PERM *relief);
extern VEC *bdBKPsolve(const BAND *A, const PERM *pivot, const PERM *relief,
		       const VEC *b, VEC *x);

extern SPMAT *spBKPfactor(SPMAT *, PERM *pivot, Real tol);
extern VEC *spBKPsolve(const SPMAT *, const PERM *pivot,
		       const VEC *b, VEC *x);
//@}

/**
   Wrappers for Meschach data structures allowing use of C++ operators.

   Idea: 
     - define for any STRUCT a class STRUCTP
     - data representation of STRUCTP is compatible to STRUCT*
     - overload operators as needed, esp. operator []

   C synonym: typedef STRUCT* STRUCTP

   Known Problem:
     - Meschach prototypes have no const's --> bad const handling
       (VEC *v_add(VEC *a, VEC *b, VEC *out) instead of
        VEC *v_add(const VEC *a, const VEC *b, VEC *out)
*/
namespace Mesch {

/** Bounds check for wrappers of Meschach types if compiled with DEBUG flag */
#ifdef DEBUG
#define MESCH_BOUNDS_CHECK(i, first_index, next_index) \
   assert((int)(first_index) <= (int)(i) && (int)(i) < (int)(next_index))
#else
#define MESCH_BOUNDS_CHECK(i, first_index, next_index)
#endif

/** NULL check for wrappers of Meschach types if compiled with DEBUG flag */
#ifdef DEBUG
#define MESCH_NULL_CHECK(ptr) \
   assert (ptr != NULL)
#else
#define MESCH_NULL_CHECK(ptr)
#endif

/** Wrapper for Meschach VEC* */
class VECP {

 protected:
  VEC *_v;	///< wrapped VEC*

 public:
  /// @name Constructors and assignments
  //@{
  VECP() {_v = VNULL;}
  VECP(VEC *cv) {_v = cv;}

  VEC *operator = (VEC *nv) {_v = nv; return _v;}
  //@}

  /// @name Operators for VECP
  //@{
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
  //@}

  /// @name Operators for const VECP
  //@{
  const Real &operator [] (int j) const
    {
      MESCH_NULL_CHECK(_v);
      MESCH_BOUNDS_CHECK(j, 0, _v->dim);
      return _v->ve[j];
    }
  const VEC *operator -> () const {return _v;}
  operator const VEC*() const {return _v;}
  //@}
};

/** Wrapper for Meschach IVEC* */
class IVECP {

 protected:
  IVEC *_v;	///< wrapped IVEC*

 public:
  /// @name Constructors and assignment
  //@{
  IVECP() {_v = IVNULL;}
  IVECP(IVEC *cv) {_v = cv;}

  IVEC *operator = (IVEC *nv) {_v = nv; return _v;}
  //@}

  /// @name Operators for IVECP
  //@{
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
  //@}

  /// @name Operators for const IVECP
  //@{
  const int &operator [] (int j) const
    {
      MESCH_NULL_CHECK(_v);
      MESCH_BOUNDS_CHECK(j, 0, _v->dim);
      return _v->ive[j];
    }
  const IVEC *operator -> () const {return _v;}
  operator const IVEC*() const {return _v;}
  //@}
};

/** Wrapper for Meschach PERM* */
class PERMP {

 protected:
  PERM *_v; 	///< wrapped PERMP

 public:
  /// @name Constructors and assignments
  //@{
  PERMP() {_v = PNULL;}
  PERMP(PERM *cv) {_v = cv;}

  PERM *operator = (PERM *nv) {_v = nv; return _v;}
  //@}

  /// @name Operators for PERMP
  //@{
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
  //@}

  /// @name Operators for const PERMP
  //@{
  const u_int &operator [] (int j) const
    {
      MESCH_NULL_CHECK(_v);
      MESCH_BOUNDS_CHECK(j, 0, _v->size);
      return _v->pe[j];
    }
  const PERM *operator -> () const {return _v;}
  operator const PERM*() const {return _v;}
  //@}
};

/** Wrapper for a row in Meschach MAT.
    It is used if compiled with DEBUG flag. */
class MATROWP {

 protected:
  Real *_row; 	///< pointer to data
  int _dim;	///< dimension of row

 public:
  /// @name Constructors and assignments
  //@{
  MATROWP(Real *row, int dim) {_row = row; _dim = dim;}
  //@}

  /// @name Operators for MATROWP
  //@{
  Real &operator [] (int j)
    {
      MESCH_BOUNDS_CHECK(j, 0, _dim);
      return _row[j];
    }
  operator Real*() {return _row;}
  //@}

  /// @name Operators for const MATROWP
  //@{
  const Real &operator [] (int j) const
    {
      MESCH_BOUNDS_CHECK(j, 0, _dim);
      return _row[j];
    }
  operator const Real*() const {return _row;}
  //@}
};

/** Wrapper for Meschach MAT* */
class MATP {

 protected:
  MAT *_m; 	///< wrapped MAT*

 public:
  /// @name Constructors and assignments
  //@{
  MATP() {_m = MNULL;}
  MATP(MAT *cm) {_m = cm;}

  MAT *operator = (MAT *nm) {_m = nm; return _m;}
  //@}

  /// @name Operators for MATP
  //@{
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
  //@}

  /// @name Operators for const MATP
  //@{
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
  //@}
};

/** Place holder for wrapper of Meschach SPMAT* */
typedef SPMAT* SPMATP;

/** Place holder for wrapper of Meschach BAND* */
typedef BAND* BANDP;

#undef Inf
/** Infinity for non existing constraints and numerical overflow */
const Real Inf = (Real)std::numeric_limits<double>::infinity();

/** check if a number is finite */
inline bool is_finite(Real x)
{
  return -Inf < x && x < Inf;
}

/** check for not a number */
inline bool is_nan(Real x)
{
  return (double)x == std::numeric_limits<double>::quiet_NaN();
}

/** Scan a string for a real number, including infinity (Inf, -Inf);
    Return the number or NaN in case of error. */
inline Real sscan_real(const char *str)
{
  Real val;
  float valf;
  if (sscanf(str, "%g", &valf)) {
    val = (Real)valf;
  } else {
    if (strncmp(str, "Inf", 3) == 0 || strncmp(str, "+Inf", 4) == 0)
      val = Inf;
    else if (strncmp(str, "-Inf", 4) == 0)
      val = -Inf;
    else
      val = (Real)std::numeric_limits<double>::quiet_NaN();
  }
  return val;
}

}; // namespace Mesch

using namespace Mesch;

#endif
