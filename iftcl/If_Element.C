/*
 *  If_Element.C -- class definition
 *
 *  rf, 6/22/94
 */

/*
    Copyright (C) 1994--2017  Ruediger Franke

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

#include "If_ListElement.h"
#include "If_Element.h"

#include <meschach/err.h>

//--------------------------------------------------------------------------
Tcl_Interp *theInterp = NULL;

//--------------------------------------------------------------------------
If_Element::If_Element(const char *ifName)
{
  // (note: need to cast away const for Tcl < 8.4)
  _token = Tcl_CreateObjCommand(theInterp, (char *)ifName,
				&If_Element::tclCmd, (ClientData)this,
				&If_Element::tclCmdDeleted);
  _deleted = false;
}

//--------------------------------------------------------------------------
If_Element::~If_Element()
{
  if (!_deleted)
    Tcl_DeleteCommandFromToken(theInterp, _token);
}

//--------------------------------------------------------------------------
const char *If_Element::ifName()
{
  return Tcl_GetCommandName(theInterp, _token);
}

//--------------------------------------------------------------------------
void If_Element::tclCmdDeleted(ClientData cld)
{
  ((If_Element *)cld)->_deleted = true;
}

//--------------------------------------------------------------------------
int If_Element::tclCmd(ClientData cld, Tcl_Interp *interp,
		       int objc, Tcl_Obj *CONST objv[])
{
  If_Element *element = (If_Element *)cld;
  int ret = TCL_ERROR;

  m_catchall(// try
#ifdef DEBUG
	     // use tractcatch to get error message printed to stderr
	     m_tracecatch(// try
#endif
			  ret = element->invoke(interp, objc, objv),
#ifdef DEBUG
			  // catch and throw
			  m_error_description()),
#endif
	     // catch
	     Tcl_AppendResult(interp, "error invoking ",
			      element->ifName(), ": ", m_error_description(),
			      NULL));

  return ret;
}


//==========================================================================
