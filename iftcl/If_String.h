/*
 *  If_String.h
 *     If_Variable for C strings (const char *)
 *
 *  rf, 3/17/01
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

#ifndef If_String_H
#define If_String_H

#include "If_Variable.h"

/** Interface type for C string */
typedef const char * If_String_t;

/** Interface variable of C string type. */
class IF_API If_String: public If_Variable<If_String_t> {

 protected:
  // conversion of internal data from and to a Tcl object
  int getTclObj(Tcl_Interp *);
  int setTclObj(Tcl_Interp *, Tcl_Obj *CONST objPtr);

 public:
  /** Constructor taking callback methods as arguments. */
  If_String(const char *ifName, If_GetIf<If_String_t> *getCb,
	    If_SetIf<If_String_t> *setCb = NULL)
    :If_Variable<If_String_t>(ifName, getCb, setCb) {}

  /** Alternative constructor for direct read access through a pointer. */
  If_String(const char *ifName, If_String_t *varPtr)
    :If_Variable<If_String_t>(ifName, new If_GetPt<If_String_t>(varPtr)) {}
};

#endif
