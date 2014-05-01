/*
 *  If_IntVEC.C -- class definition
 *
 *  rf, 1/31/97
 *
 *  rf, 2/2/97
 *   -- check NULL
 */

/*
    Copyright (C) 1994--2014  Ruediger Franke

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

#include "If_IntVec.h"

//--------------------------------------------------------------------------
int If_IntVec::setTclObj(Tcl_Interp *interp, Tcl_Obj *CONST objPtr)
{
  const IVEC *curIVec = get();
  int  j;
  int  ncols = 0;
  Tcl_Obj **cols;
  int element;

  // split up objPtr into list elements
  //-----------------------------------
  if (Tcl_ListObjGetElements(interp, objPtr, &ncols, &cols) != TCL_OK) {
    return TCL_ERROR;
  }

  // create and fill a new vector
  //-----------------------------
  IVEC *newIVec = iv_get(ncols);

  for (j = 0; j < ncols; j++) {
      
    // parse new value
    //----------------
    if(Tcl_GetIntFromObj(interp, cols[j], &element) != TCL_OK) {
      iv_free(newIVec);
      return TCL_ERROR;
    }
    newIVec->ive[j] = element;
  }

  // use the new vector
  //-------------------
  set(newIVec);
  iv_free(newIVec);

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_IntVec::getTclObj(Tcl_Interp *interp)
{
  const IVEC *curIVec = get();
  int j, jend = curIVec? curIVec->dim: 0;
  Tcl_Obj *listPtr = Tcl_NewListObj(0, NULL);
  Tcl_Obj *objPtr;

  for (j = 0; j < jend; j++) {
    objPtr = Tcl_NewIntObj(curIVec->ive[j]);
    Tcl_ListObjAppendElement(interp, listPtr, objPtr);
  }

  Tcl_SetObjResult(interp, listPtr);

  return TCL_OK;
}


//==========================================================================
