/**
 *  @file If_IntVec.h
 *    Interface variables for integer vectors of the Meschach library.
 *
 *  rf, 1/31/97
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

#ifndef If_IntVec_H
#define If_IntVec_H

// Include file with declaration of vector and matrix type.
#include "Meschach.h"

#include "If_Variable.h"

/** Interface integer vector type */
typedef const Mesch::IVECP If_IntVec_t;

/** Interface variable of integer vector type. */
class IF_API If_IntVec: public If_Variable<If_IntVec_t> {

 protected:
  // conversion of internal data from and to a Tcl object
  int getTclObj(Tcl_Interp *);
  int setTclObj(Tcl_Interp *, Tcl_Obj *CONST objPtr);

 public:
  /** Constructor taking callback methods as arguments. */
  If_IntVec(const char *ifName, If_GetIf<If_IntVec_t> *getCb,
	 If_SetIf<If_IntVec_t> *setCb = NULL)
    :If_Variable<If_IntVec_t>(ifName, getCb, setCb) {}

  /** Alternative constructor for direct read/write access to an IVECP. */
  If_IntVec(const char *ifName, IVECP *varPtr)
    :If_Variable<If_IntVec_t>(ifName,
			      new If_GetPt<If_IntVec_t>(varPtr),
			      new SetIVEC((IVEC **)varPtr)) {}

  /** Alternative constructor for direct read/write access to an IVEC*. */
  If_IntVec(const char *ifName, IVEC **varPtr)
    :If_Variable<If_IntVec_t>(ifName,
			      new If_GetPt2<If_IntVec_t, IVEC*>(varPtr),
			      new SetIVEC(varPtr)) {}

 private:
  /** Write callback for direct access to an integer vector */
  class SetIVEC: public If_SetIf<If_IntVec_t> {
   protected:
    IVEC **_varPtr; 	///< pointer to accessed variable

   public:
    /** constructor */
    SetIVEC(IVEC **varPtr): _varPtr(varPtr) {}
    /** set method ensuring that size of vector does not change */
    void set(If_IntVec_t value) {iv_copy_elements(value, *_varPtr);}
  };
};

#endif
