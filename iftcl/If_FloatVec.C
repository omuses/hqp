/*
 *  If_FloatVec.C -- class definition
 *
 *  rf, 6/22/94
 *
 *  rf, 2/2/97
 *   -- check NULL
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

#include <stdio.h>
#include <string.h>
#include <malloc.h>

#include "If_FloatVec.h"

//--------------------------------------------------------------------------
If_FloatVec::If_FloatVec(char *ifName, VEC **varPtr,
			 If_FloatVecWriteIf *callback)
:If_Variable(ifName)
{
  _varPtr = varPtr;
  _callback = callback;
}

//--------------------------------------------------------------------------
If_FloatVec::If_FloatVec(char *ifName, VECP *varPtr,
			 If_FloatVecWriteIf *callback)
:If_Variable(ifName)
{
  _varPtr = (VEC **)varPtr;
  _callback = callback;
}

//--------------------------------------------------------------------------
If_FloatVec::~If_FloatVec()
{
  delete _callback;
}

//--------------------------------------------------------------------------
int If_FloatVec::put(Tcl_Obj *CONST objPtr)
{
  int  j;
  int  ncols=0;
  Tcl_Obj **cols;
  double element;

  // split up objPtr into list elements
  //-----------------------------------
  if (Tcl_ListObjGetElements(theInterp, objPtr, &ncols, &cols) != TCL_OK) {
    return TCL_ERROR;
  }

  // check dimension
  //----------------
  if (*_varPtr == VNULL || ncols != (int)(*_varPtr)->dim) {
    Tcl_AppendResult(theInterp, "wrong dimension", NULL);
    return TCL_ERROR;
  }

  // create and fill a new vector
  //-----------------------------
  VEC *newVec = v_get(ncols);

  for (j = 0; j < ncols; j++) {
      
    // parse new value
    //----------------
    if(Tcl_GetDoubleFromObj(theInterp, cols[j], &element) != TCL_OK) {
      v_free(newVec);
      return TCL_ERROR;
    }
    newVec->ve[j] = element;
  }

  // use the new vector
  //-------------------
  if (_callback) {
    if (_callback->write(newVec) != IF_OK) {
      Tcl_AppendResult(theInterp, "writing of ", _ifName, " rejected", NULL);
      v_free(newVec);
      return TCL_ERROR;
    }
  }
  else
    *_varPtr = v_copy(newVec, *_varPtr);

  v_free(newVec);

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_FloatVec::get()
{
  VEC	*vec = *_varPtr;
  int   j, jend = vec? vec->dim: 0;
  Tcl_Obj *listPtr = Tcl_NewListObj(0, NULL);
  Tcl_Obj *objPtr;

  for (j = 0; j < jend; j++) {
    objPtr = Tcl_NewDoubleObj(vec->ve[j]);
    Tcl_ListObjAppendElement(theInterp, listPtr, objPtr);
  }

  Tcl_SetObjResult(theInterp, listPtr);

  return TCL_OK;
}


//==========================================================================
