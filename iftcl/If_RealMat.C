/*
 *  If_RealMat.C -- class definition
 *
 *  rf, 8/19/94
 *
 *  rf, 2/2/97
 *   -- check NULL
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
#include <malloc.h>

#include "If_RealMat.h"

//--------------------------------------------------------------------------
If_RealMat::If_RealMat(const char *ifName, MAT **varPtr, const char *mode)
  :If_Variable(ifName, mode)
{
  _varPtr = varPtr;
  _callback = NULL;
}

//--------------------------------------------------------------------------
If_RealMat::If_RealMat(const char *ifName, MATP *varPtr, const char *mode)
  :If_Variable(ifName, mode)
{
  _varPtr = (MAT **)varPtr;
  _callback = NULL;
}

//--------------------------------------------------------------------------
If_RealMat::If_RealMat(const char *ifName, MAT **varPtr,
		       If_RealMatWriteIf *callback)
  :If_Variable(ifName)
{
  _varPtr = varPtr;
  _callback = callback;
}

//--------------------------------------------------------------------------
If_RealMat::If_RealMat(const char *ifName, MATP *varPtr,
		       If_RealMatWriteIf *callback)
  :If_Variable(ifName)
{
  _varPtr = (MAT **)varPtr;
  _callback = callback;
}

//--------------------------------------------------------------------------
If_RealMat::~If_RealMat()
{
  delete _callback;
}

//--------------------------------------------------------------------------
int If_RealMat::put(Tcl_Obj *CONST objPtr)
{
  int  i, j;
  int  nrows = 0, ncols = 0;
  int  dummy;
  Tcl_Obj **rows, **cols;
  double element;

  // split up the matrix into it's rows
  //-----------------------------------
  if (Tcl_ListObjGetElements(theInterp, objPtr, &nrows, &rows) != TCL_OK)
    return TCL_ERROR;

  // split up the first row -- to get the number of cols
  //----------------------------------------------------
  if (nrows > 0)
    if (Tcl_ListObjGetElements(theInterp, rows[0], &ncols, &cols) != TCL_OK) {
      return TCL_ERROR;
    }

  // check dimension
  //----------------
  if (*_varPtr == MNULL
      || nrows != (int)(*_varPtr)->m || ncols != (int)(*_varPtr)->n) {
    Tcl_AppendResult(theInterp, "wrong dimension", NULL);
    return TCL_ERROR;
  }

  // create and fill a new matrix
  //-----------------------------
  MAT *newMat = m_get(nrows, ncols);

  for (i = 0; i < nrows; i++) {
    if (i > 0) {
      if (Tcl_ListObjGetElements(theInterp, rows[i], &dummy, &cols) != TCL_OK) {
	m_free(newMat);
	return TCL_ERROR;
      }
      if (dummy != ncols) {
	Tcl_AppendResult(theInterp, "Different sizes of rows!", NULL);
	m_free(newMat);
	return TCL_ERROR;
      }
    }
    for ( j = 0; j < ncols; j++) {
      
      // parse the new value
      //--------------------
      if (Tcl_GetDoubleFromObj(theInterp, cols[j], &element) != TCL_OK) {
	m_free(newMat);
	return TCL_ERROR;
      }
      newMat->me[i][j] = element;
    }
  }

  // use the new matrix
  //-------------------
  if (_callback) {
    if (_callback->write(newMat) != IF_OK) {
      Tcl_AppendResult(theInterp, "writing of ", _ifName, " rejected", NULL);
      m_free(newMat);
      return TCL_ERROR;
    }
  }
  else
    *_varPtr = m_copy(newMat, *_varPtr);

  m_free(newMat);

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_RealMat::get()
{
  MAT *mat = *_varPtr;
  int i, iend = mat? mat->m: 0;
  int j, jend = mat? mat->n: 0;
  Tcl_Obj *listPtr = Tcl_NewListObj(0, NULL);
  Tcl_Obj *rowPtr, *objPtr;

  for (i = 0; i < iend; i++) {
    rowPtr = Tcl_NewListObj(0, NULL);
    for (j = 0; j < jend; j++) {
      objPtr = Tcl_NewDoubleObj(mat->me[i][j]);
      Tcl_ListObjAppendElement(theInterp, rowPtr, objPtr);
    }
    Tcl_ListObjAppendElement(theInterp, listPtr, rowPtr);
  }

  Tcl_SetObjResult(theInterp, listPtr);

  return TCL_OK;
}


//==========================================================================
