/*
 *  If_Float.C -- class definiton
 *
 *  rf, 1/25/97
 */

/*
    Copyright (C) 1994--1998  Ruediger Franke

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

#include <math.h>

#include "If_Float.h"

//--------------------------------------------------------------------------
If_Float::If_Float(char *ifName, If_Float_t *varPtr,
		   If_FloatWriteIf *callback)
:If_Variable(ifName)
{
  _varPtr = varPtr;
  _callback = callback;
}

//--------------------------------------------------------------------------
If_Float::~If_Float()
{
  delete _callback;
}

//--------------------------------------------------------------------------
int If_Float::put(Tcl_Obj *CONST objPtr)
{
  If_Float_t value;

  // parse the new value
  //--------------------
  if (Tcl_GetDoubleFromObj(theInterp, objPtr, &value) != TCL_OK)
    return TCL_ERROR;

  // use the new value
  //------------------
  if (_callback) {
    if (_callback->write(value) != IF_OK) {
      Tcl_AppendResult(theInterp, "writing of ", _ifName, " rejected", NULL);
      return TCL_ERROR;
    }
  }
  else
    *_varPtr = value;

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_Float::get()
{
  Tcl_Obj *objPtr = Tcl_NewDoubleObj(*_varPtr);

  Tcl_SetObjResult(theInterp, objPtr);

  return TCL_OK;
}


//==========================================================================
