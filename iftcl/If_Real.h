/**
 *  @file If_Real.h
 *     Interface variable for doubles.
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

#ifndef If_Real_H
#define If_Real_H

#include "If_Variable.h"

/** Interface real type */
typedef double If_Real_t;

/** Interface variable of real type. */
class IF_API If_Real: public If_Variable<If_Real_t> {

 protected:
  // conversion of internal data from and to a Tcl object
  int getTclObj(Tcl_Interp *);
  int setTclObj(Tcl_Interp *, Tcl_Obj *CONST objPtr);

 public:
  /** Constructor taking callback methods as arguments. */
  If_Real(const char *ifName, If_GetIf<If_Real_t> *getCb,
	  If_SetIf<If_Real_t> *setCb = NULL)
    :If_Variable<If_Real_t>(ifName, getCb, setCb) {}

  /** Alternative constructor for direct access to a variable pointer. */
  If_Real(const char *ifName, If_Real_t *varPtr)
    :If_Variable<If_Real_t>(ifName, new If_GetPt<If_Real_t>(varPtr),
			    new If_SetPt<If_Real_t>(varPtr)) {}
};

#endif
