/**
 *  @file If_Element.h
 *     interface element to Tcl
 *
 *  rf, 22/6/94
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

#ifndef If_Element_H
#define If_Element_H

#include <tcl.h>

#include "If.h"
#include "If_ListElement.h"

/** global Tcl interpreter that must be initialized by the application */
extern Tcl_Interp *theInterp;

/**
 *  Abstract base class for interface elements to Tcl, like variables and
 *  procedures. A Tcl command is created for each interface element.
 */
class IF_API If_Element: public If_ListElement {

 private:
  Tcl_Command 	_token;		///< token identifying the interface element
  bool 		_deleted; 	///< mark if command was deleted

  /** Callback for Tcl if command is invoked */
  static  int	tclCmd(ClientData, Tcl_Interp*,
		       int objc, Tcl_Obj *CONST objv[]);

  /** Callback for Tcl if command is being deleted */
  static  void	tclCmdDeleted(ClientData);

 protected:
  /** Interface for derived classes to process command invocation */
  virtual int invoke(Tcl_Interp *, int objc, Tcl_Obj *CONST objv[]) = 0;

 public:
  If_Element(const char *ifName); 	///< constructor
  virtual ~If_Element(); 		///< destructor
  const char *ifName(); 		///< get interface name
};


#endif
