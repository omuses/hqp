/*
 *  If_Bool.C -- class definiton
 *
 *  rf, 1/17/97
 */

/*
    Copyright (C) 1994--2001  Ruediger Franke

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

#include <stdio.h>
#include <string.h>

#include "If_Bool.h"


//--------------------------------------------------------------------------
If_Bool::If_Bool(const char *ifName, If_Bool_t *varPtr, const char *mode)
  :If_Variable(ifName, mode)
{
  _varPtr = varPtr;
  _callback = NULL;
}

//--------------------------------------------------------------------------
If_Bool::If_Bool(const char *ifName, If_Bool_t *varPtr,
		 If_BoolWriteIf *callback)
  :If_Variable(ifName)
{
  _varPtr = varPtr;
  _callback = callback;
}

//--------------------------------------------------------------------------
If_Bool::~If_Bool()
{
  delete _callback;
}

//--------------------------------------------------------------------------
int If_Bool::put(Tcl_Obj *CONST objPtr)
{
  int value;

  // parse the new value
  //--------------------
  if (Tcl_GetBooleanFromObj(theInterp, objPtr, &value) != TCL_OK)
    return TCL_ERROR;

  // use the new value
  //------------------
  if (_callback) {
    if (_callback->write((If_Bool_t)value) != IF_OK) {
      Tcl_AppendResult(theInterp, "writing of ", _ifName, " rejected", NULL);
      return TCL_ERROR;
    }
  }
  else
    *_varPtr = (If_Bool_t)value;

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_Bool::get()
{
  int value = (int)*_varPtr;
  Tcl_Obj *objPtr = Tcl_NewBooleanObj(*_varPtr);

  Tcl_SetObjResult(theInterp, objPtr);

  return TCL_OK;
}


//==========================================================================
