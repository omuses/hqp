/*
 *  If_RealVec.C -- class definition
 *
 *  rf, 6/22/94
 *
 *  rf, 2/2/97
 *   -- check NULL
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

#include "If_RealVec.h"

//--------------------------------------------------------------------------
int If_RealVec::setTclObj(Tcl_Interp *interp, Tcl_Obj *CONST objPtr)
{
  const VEC *curVec = get();
  int  j;
  int  ncols=0;
  Tcl_Obj **cols;
  double element;

  // split up objPtr into list elements
  //-----------------------------------
  if (Tcl_ListObjGetElements(interp, objPtr, &ncols, &cols) != TCL_OK) {
    return TCL_ERROR;
  }

  // check dimension
  //----------------
  if (curVec == VNULL || ncols != (int)(curVec)->dim) {
    Tcl_AppendResult(interp, "wrong dimension for ", ifName(), NULL);
    return TCL_ERROR;
  }

  // create and fill a new vector
  //-----------------------------
  VEC *newVec = v_get(ncols);

  for (j = 0; j < ncols; j++) {
      
    // parse new value
    //----------------
    if(Tcl_GetDoubleFromObj(interp, cols[j], &element) != TCL_OK) {
      // in case of error check for Inf, +Inf, -Inf
      int len;
      const char *str = Tcl_GetStringFromObj(cols[j], &len);
      element = 2.0;
      if (len == 4 && str[0] == '+')
	element = 1.0;
      else if (len == 4 && str[0] == '-')
	element = -1.0;
      if (element != 2.0) {
	len--;
	str++;
      }
      if (len == 3 && str[0] == 'I' && str[1] == 'n' && str[2] == 'f') {
	Tcl_ResetResult(interp);
	element *= Inf;
      }
      else {
	v_free(newVec);
	return TCL_ERROR;
      }
    }
    newVec->ve[j] = element;
  }

  // use the new vector
  //-------------------
  set(newVec);
  v_free(newVec);

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_RealVec::getTclObj(Tcl_Interp *interp)
{
  const VEC *curVec = get();
  int   j, jend = curVec? curVec->dim: 0;
  Tcl_Obj *listPtr = Tcl_NewListObj(0, NULL);
  Tcl_Obj *objPtr;

  for (j = 0; j < jend; j++) {
    objPtr = Tcl_NewDoubleObj(curVec->ve[j]);
    Tcl_ListObjAppendElement(interp, listPtr, objPtr);
  }

  Tcl_SetObjResult(interp, listPtr);

  return TCL_OK;
}


//==========================================================================
