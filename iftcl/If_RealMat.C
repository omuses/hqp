/*
 *  If_RealMat.C -- class definition
 *
 *  rf, 8/19/94
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

#include "If_RealMat.h"

//--------------------------------------------------------------------------
int If_RealMat::setTclObj(Tcl_Interp *interp, Tcl_Obj *CONST objPtr)
{
  const MAT *curMat = get();
  int  i, j;
  int  nrows = 0, ncols = 0;
  int  dummy;
  Tcl_Obj **rows, **cols;
  double element;

  // split up the matrix into it's rows
  //-----------------------------------
  if (Tcl_ListObjGetElements(interp, objPtr, &nrows, &rows) != TCL_OK)
    return TCL_ERROR;

  // split up the first row -- to get the number of cols
  //----------------------------------------------------
  if (nrows > 0)
    if (Tcl_ListObjGetElements(interp, rows[0], &ncols, &cols) != TCL_OK) {
      return TCL_ERROR;
    }

  // check dimension
  //----------------
  if (curMat == MNULL
      || nrows != (int)(curMat)->m || ncols != (int)(curMat)->n) {
    Tcl_AppendResult(interp, "wrong dimension for ", _ifName, NULL);
    return TCL_ERROR;
  }

  // create and fill a new matrix
  //-----------------------------
  MAT *newMat = m_get(nrows, ncols);

  for (i = 0; i < nrows; i++) {
    if (i > 0) {
      if (Tcl_ListObjGetElements(interp, rows[i], &dummy, &cols) != TCL_OK) {
	m_free(newMat);
	return TCL_ERROR;
      }
      if (dummy != ncols) {
	Tcl_AppendResult(interp, "different sizes of rows in ",
			 _ifName, "!", NULL);
	m_free(newMat);
	return TCL_ERROR;
      }
    }
    for ( j = 0; j < ncols; j++) {
      
      // parse the new value
      //--------------------
      if (Tcl_GetDoubleFromObj(interp, cols[j], &element) != TCL_OK) {
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
	  m_free(newMat);
	  return TCL_ERROR;
	}
      }
      newMat->me[i][j] = element;
    }
  }

  // use the new matrix
  //-------------------
  set(newMat);
  m_free(newMat);

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_RealMat::getTclObj(Tcl_Interp *interp)
{
  const MAT *curMat = get();
  int i, iend = curMat? curMat->m: 0;
  int j, jend = curMat? curMat->n: 0;
  Tcl_Obj *listPtr = Tcl_NewListObj(0, NULL);
  Tcl_Obj *rowPtr, *objPtr;

  for (i = 0; i < iend; i++) {
    rowPtr = Tcl_NewListObj(0, NULL);
    for (j = 0; j < jend; j++) {
      objPtr = Tcl_NewDoubleObj(curMat->me[i][j]);
      Tcl_ListObjAppendElement(interp, rowPtr, objPtr);
    }
    Tcl_ListObjAppendElement(interp, listPtr, rowPtr);
  }

  Tcl_SetObjResult(interp, listPtr);

  return TCL_OK;
}


//==========================================================================
