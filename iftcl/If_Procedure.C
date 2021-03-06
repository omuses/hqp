/*
 *  If_Procedure.C -- class definition
 *
 *  rf, 7/20/94
 *
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

#include "If_Procedure.h"

//--------------------------------------------------------------------------
If_Procedure::If_Procedure(const char *ifName, If_Procedure_t *proc)
  :If_Element(ifName)
{
  _proc = proc;
}

//--------------------------------------------------------------------------
int If_Procedure::invoke(Tcl_Interp *interp, int objc, Tcl_Obj *CONST [])
{
  if (objc > 1) {
    Tcl_AppendResult(interp, "wrong # args, should be: ",
		     ifName(), NULL);
  }
  (*_proc)();
  return TCL_OK;
}


//==========================================================================
