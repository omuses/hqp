/*
 *  If_IntVEC.C -- class definition
 *
 *  rf, 1/31/97
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

#include "If_IntVec.h"

//--------------------------------------------------------------------------
If_IntVec::If_IntVec(char *ifName, IVEC **varPtr,
		     If_IntVecWriteIf *callback)
:If_Variable(ifName)
{
  _varPtr = varPtr;
  _callback = callback;
}

//--------------------------------------------------------------------------
If_IntVec::If_IntVec(char *ifName, IVECP *varPtr,
		     If_IntVecWriteIf *callback)
:If_Variable(ifName)
{
  _varPtr = (IVEC **)varPtr;
  _callback = callback;
}

//--------------------------------------------------------------------------
If_IntVec::If_IntVec(char *ifName, IVEC **varPtr, const char *mode)
:If_Variable(ifName, mode)
{
  _varPtr = varPtr;
  _callback = NULL;
}

//--------------------------------------------------------------------------
If_IntVec::If_IntVec(char *ifName, IVECP *varPtr, const char *mode)
:If_Variable(ifName, mode)
{
  _varPtr = (IVEC **)varPtr;
  _callback = NULL;
}

//--------------------------------------------------------------------------
If_IntVec::~If_IntVec()
{
  delete _callback;
}

//--------------------------------------------------------------------------
int If_IntVec::put(Tcl_Obj *CONST objPtr)
{
  int  j;
  int  ncols = 0;
  Tcl_Obj **cols;
  int element;

  // split up objPtr into list elements
  //-----------------------------------
  if (Tcl_ListObjGetElements(theInterp, objPtr, &ncols, &cols) != TCL_OK) {
    return TCL_ERROR;
  }

  // check dimension
  //----------------
  if (*_varPtr == IVNULL || ncols != (int)(*_varPtr)->dim) {
    Tcl_AppendResult(theInterp, "wrong dimension", NULL);
    return TCL_ERROR;
  }

  // create and fill a new vector
  //-----------------------------
  IVEC *newIVEC = iv_get(ncols);

  for (j = 0; j < ncols; j++) {
      
    // parse new value
    //----------------
    if(Tcl_GetIntFromObj(theInterp, cols[j], &element) != TCL_OK) {
      iv_free(newIVEC);
      return TCL_ERROR;
    }
    newIVEC->ive[j] = element;
  }

  // use the new vector
  //-------------------
  if (_callback) {
    if (_callback->write(newIVEC) != IF_OK) {
      Tcl_AppendResult(theInterp, "writing to ", _ifName, " rejected", NULL);
      iv_free(newIVEC);
      return TCL_ERROR;
    }
  }
  else
    *_varPtr = iv_copy(newIVEC, *_varPtr);

  iv_free(newIVEC);

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_IntVec::get()
{
  IVEC *ivec = *_varPtr;
  int j, jend = ivec? ivec->dim: 0;
  Tcl_Obj *listPtr = Tcl_NewListObj(0, NULL);
  Tcl_Obj *objPtr;

  for (j = 0; j < jend; j++) {
    objPtr = Tcl_NewIntObj(ivec->ive[j]);
    Tcl_ListObjAppendElement(theInterp, listPtr, objPtr);
  }

  Tcl_SetObjResult(theInterp, listPtr);

  return TCL_OK;
}


//==========================================================================
