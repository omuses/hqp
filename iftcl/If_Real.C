/*
 *  If_Real.C -- class definiton
 *
 *  rf, 1/25/97
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

#include <math.h>

#include "If_Real.h"

//--------------------------------------------------------------------------
If_Real::If_Real(const char *ifName, If_Real_t *varPtr, const char *mode)
  :If_Variable(ifName, mode)
{
  _varPtr = varPtr;
  _callback = NULL;
}

//--------------------------------------------------------------------------
If_Real::If_Real(const char *ifName, If_Real_t *varPtr,
		 If_RealWriteIf *callback)
  :If_Variable(ifName)
{
  _varPtr = varPtr;
  _callback = callback;
}

//--------------------------------------------------------------------------
If_Real::~If_Real()
{
  delete _callback;
}

//--------------------------------------------------------------------------
int If_Real::put(Tcl_Obj *CONST objPtr)
{
  If_Real_t value;

  // parse the new value
  //--------------------
  if (Tcl_GetDoubleFromObj(theInterp, objPtr, &value) != TCL_OK) {
    // in case of error check for Inf, +Inf, -Inf
    int len;
    const char *str = Tcl_GetStringFromObj(objPtr, &len);
    value = 2.0;
    if (len == 4 && str[0] == '+')
      value = 1.0;
    else if (len == 4 && str[0] == '-')
      value = -1.0;
    if (value != 2.0) {
      len--;
      str++;
    }
    if (len == 3 && str[0] == 'I' && str[1] == 'n' && str[2] == 'f') {
      Tcl_ResetResult(theInterp);
      value *= Inf;
    }
    else
      return TCL_ERROR;
  }

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
int If_Real::get()
{
  Tcl_Obj *objPtr = Tcl_NewDoubleObj(*_varPtr);

  Tcl_SetObjResult(theInterp, objPtr);

  return TCL_OK;
}


//==========================================================================
