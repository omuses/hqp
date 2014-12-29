/**
 *  @file If_RealVec.h
 *     Interface variables for vectors of reals (currently Mesch::VECP)
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

#ifndef If_RealVec_H
#define If_RealVec_H

// Include file with declaration of vector and matrix type.
#include <Meschach.h>

#include "If_Variable.h"

/** Interface real vector type */
typedef const Mesch::VECP If_RealVec_t;

/** Interface variable of real vector type. */
class IF_API If_RealVec: public If_Variable<If_RealVec_t> {

 protected:
  // conversion of internal data from and to a Tcl object
  int getTclObj(Tcl_Interp *);
  int setTclObj(Tcl_Interp *, Tcl_Obj *CONST objPtr);

 public:
  /** Constructor taking callback methods as arguments. */
  If_RealVec(const char *ifName, If_GetIf<If_RealVec_t> *getCb,
	 If_SetIf<If_RealVec_t> *setCb = NULL)
    :If_Variable<If_RealVec_t>(ifName, getCb, setCb) {}

  /** Alternative constructor for direct read/write access to a VECP. */
  If_RealVec(const char *ifName, VECP *varPtr)
    :If_Variable<If_RealVec_t>(ifName,
			      new If_GetPt<If_RealVec_t>(varPtr),
			      new SetVEC((VEC **)varPtr)) {}

  /** Alternative constructor for direct read/write access to a VEC*. */
  If_RealVec(const char *ifName, VEC **varPtr)
    :If_Variable<If_RealVec_t>(ifName,
			      new If_GetPt2<If_RealVec_t, VEC*>(varPtr),
			      new SetVEC(varPtr)) {}

 private:
  /** Write callback for direct access to a real vector */
  class SetVEC: public If_SetIf<If_RealVec_t> {
   protected:
    VEC **_varPtr; 	///< pointer to accessed variable

   public:
    /** constructor */
    SetVEC(VEC **varPtr): _varPtr(varPtr) {}
    /** set method ensuring that size of vector does not change */
    void set(If_RealVec_t value) {v_copy_elements(value, *_varPtr);}
  };
};

#endif
