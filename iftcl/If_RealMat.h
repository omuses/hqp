/**
 *  @file If_RealMat.h
 *     Interface variables for matrices of reals (currently Mesch::MATP)
 *
 *  rf, 2/7/97
 */

/*
    Copyright (C) 1994--2002  Ruediger Franke

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

#ifndef If_RealMat_H
#define If_RealMat_H

// Include file with declaration of vector and matrix type.
#include <Meschach.h>

#include "If_Variable.h"

/** Interface real matrix type */
typedef const Mesch::MATP If_RealMat_t;

/** Interface variable of real matrix type. */
class IF_API If_RealMat: public If_Variable<If_RealMat_t> {

 protected:
  // conversion of internal data from and to a Tcl object
  int getTclObj(Tcl_Interp *);
  int setTclObj(Tcl_Interp *, Tcl_Obj *CONST objPtr);

 public:
  /** Constructor taking callback methods as arguments. */
  If_RealMat(const char *ifName, If_GetIf<If_RealMat_t> *getCb,
	 If_SetIf<If_RealMat_t> *setCb = NULL)
    :If_Variable<If_RealMat_t>(ifName, getCb, setCb) {}

  /** Alternative constructor for direct read/write access to a MATP. */
  If_RealMat(const char *ifName, MATP *varPtr)
    :If_Variable<If_RealMat_t>(ifName,
			      new If_GetPt<If_RealMat_t>(varPtr),
			      new SetMAT((MAT **)varPtr)) {}

  /** Alternative constructor for direct read/write access to a MAT*. */
  If_RealMat(const char *ifName, MAT **varPtr)
    :If_Variable<If_RealMat_t>(ifName,
			      new If_GetPt2<If_RealMat_t, MAT*>(varPtr),
			      new SetMAT(varPtr)) {}

 private:
  /** Write callback for direct access to a real matrix */
  class SetMAT: public If_SetIf<If_RealMat_t> {
   protected:
    MAT **_varPtr; 	///< pointer to accessed variable

   public:
    /** constructor */
    SetMAT(MAT **varPtr): _varPtr(varPtr) {}
    /** set method ensuring that size of matrix does not change */
    void set(If_RealMat_t value) {m_copy_elements(value, *_varPtr);}
  };
};

#endif
