/*
 *  If_Variable.C
 *
 *
 *  rf, 6/22/94
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

#include <stdio.h>
#include <string.h>

#include "If_Variable.h"

#define IF_READ  1
#define IF_WRITE 2

//--------------------------------------------------------------------------
If_VariableCmd::If_VariableCmd(const char *ifName, const char *mode)
  :If_Element(ifName)
{
  Tcl_CreateObjCommand(theInterp, _ifName,
		       tclCmd, (ClientData)this, NULL);

  int len = strlen(mode);
  _mode = 0;
  if (len > 0) {
    if (mode[0] == 'r')
      _mode |= IF_READ;
    if (mode[0] == 'w' || len > 1 && mode[1] == 'w')
      _mode |= IF_WRITE;
  }
}

//--------------------------------------------------------------------------
If_VariableCmd::~If_VariableCmd()
{
  Tcl_DeleteCommand(theInterp, _ifName);
}

//--------------------------------------------------------------------------
int If_VariableCmd::tclCmd(ClientData cld, Tcl_Interp *interp,
			   int objc, Tcl_Obj *CONST objv[])
{
  If_VariableCmd *var = (If_VariableCmd *)cld;

  m_catchall(// try
	     // use tractcatch to get error message printed to stderr
	     m_tracecatch(// try
			  switch (objc) {

			  case 1:
			    if (var->_mode & IF_READ)
			      return var->getTclObj(interp);
			    else {
			      Tcl_AppendResult(interp,
					       "read permission denied for ",
					       var->_ifName, NULL);
			      return TCL_ERROR;
			    }

			  case 2:
			    if (var->_mode & IF_WRITE)
			      return var->setTclObj(interp, objv[1]);
			    else {
			      Tcl_AppendResult(interp,
					       "write permission denied for ",
					       var->_ifName, NULL);
			      return TCL_ERROR;
			    }

			  default:
			    Tcl_AppendResult(interp,
					     "wrong # args, should be: ",
					     var->_ifName, " [new value]",
					     NULL);
			    return TCL_ERROR;
			  },
			  // catch and throw
			  "If_Variable::tclCmd"),
	     // catch
	     Tcl_AppendResult(interp, "error accessing ",
			      var->_ifName, NULL);
	     return TCL_ERROR);

  return TCL_ERROR; // this shall never be reached
}


//==========================================================================
