/*
 *  If_String.C -- class definition
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

#include "If_String.h"

#include <string.h>

//--------------------------------------------------------------------------
int If_String::setTclObj(Tcl_Interp *, Tcl_Obj *CONST objPtr)
{
  If_String_t newValue;
  int newLength = 0;

  // get the new value
  //------------------
  newValue = Tcl_GetStringFromObj(objPtr, &newLength);
  if (newValue == NULL)
    return TCL_ERROR;

  // use the new value
  //------------------
  set(newValue);

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_String::getTclObj(Tcl_Interp *interp)
{
  If_String_t curValue = get();
  // (need to cast to CONST char* as Tcl may not use const qualifier)
  Tcl_Obj *objPtr = Tcl_NewStringObj((CONST char *)curValue,
				     strlen(curValue));

  Tcl_SetObjResult(interp, objPtr);

  return TCL_OK;
}


//==========================================================================
