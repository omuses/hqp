/**
 *  @file If_StdString.h
 *     If_Variable for standard C++ string type.
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

#ifndef If_StdString_H
#define If_StdString_H

#include "If_Variable.h"

#include <string>

/** Interface type for standard string */
typedef const std::string & If_StdString_t;

/** Interface variable of standard string type. */
class If_StdString: public If_Variable<If_StdString_t> {

 protected:
  // conversion of internal data from and to a Tcl object
  int getTclObj(Tcl_Interp *);
  int setTclObj(Tcl_Interp *, Tcl_Obj *CONST objPtr);

 public:
  /** Constructor taking callback methods as arguments. */
  If_StdString(const char *ifName, If_GetIf<If_StdString_t> *getCb,
	       If_SetIf<If_StdString_t> *setCb = NULL)
    :If_Variable<If_StdString_t>(ifName, getCb, setCb) {}
};

/*
 *  If_StdString.C -- class definition
 *     (Everything is inlined so that only programs that use it
 *      have to link a Std C++ library.)
 */

//--------------------------------------------------------------------------
int If_StdString::setTclObj(Tcl_Interp *, Tcl_Obj *CONST objPtr)
{
  char *str_value;
  int str_length;

  // get the new value
  //------------------
  str_value = Tcl_GetStringFromObj(objPtr, &str_length);
  if (str_value == NULL)
    return TCL_ERROR;

  If_StdString_t newValue(str_value);

  // use the new value
  //------------------
  set(newValue);

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_StdString::getTclObj(Tcl_Interp *interp)
{
  const std::string &curValue = get();
  // (need to cast to CONST char* as Tcl may not use const qualifier)
  Tcl_Obj *objPtr = Tcl_NewStringObj((CONST char *)curValue.c_str(),
				     curValue.length());

  Tcl_SetObjResult(interp, objPtr);

  return TCL_OK;
}

//==========================================================================

#endif
