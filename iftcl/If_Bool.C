/*
 *  If_Bool.C -- class definiton
 *
 *  rf, 1/17/97
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

#include "If_Bool.h"

//--------------------------------------------------------------------------
int If_Bool::setTclObj(Tcl_Interp *interp, Tcl_Obj *CONST objPtr)
{
  int value;

  // parse the new value
  //--------------------
  if (Tcl_GetBooleanFromObj(interp, objPtr, &value) != TCL_OK)
    return TCL_ERROR;

  // use the new value
  //------------------
  set(value);

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_Bool::getTclObj(Tcl_Interp *interp)
{
  Tcl_Obj *objPtr = Tcl_NewBooleanObj(get());

  Tcl_SetObjResult(interp, objPtr);

  return TCL_OK;
}


//==========================================================================
