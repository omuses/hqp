/*
 *  If_Variable.h
 *     - superclass: If_Element
 *     - abstract base class for interface-variables to Tcl
 *     - a derived class converts Tcl data to internal C++ representation
 *       for one C++ data type
 *       (see class If_Real for example)
 *     - has direct access to variables for reading
 *     - writing can be done with direct access or with a Callback procedure
 *
 *  rf, 6/22/94
 *
 *  rf, 8/13/98
 *   - use typed Tcl 8 objects instead of strings
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

#ifndef If_Variable_H
#define If_Variable_H

#include "If_Element.h"

class If_Variable: public If_Element {

 protected:
  If_Variable(const char *ifName, const char *mode = "rw");

  // interface to Tcl
  //-----------------
  static  int	tclCmd(ClientData, Tcl_Interp*,
		       int objc, Tcl_Obj *CONST objv[]);

  virtual int	put(Tcl_Obj *CONST objPtr) = 0;
  virtual int	get()=0;

  int		_mode;

 public:
  virtual	~If_Variable();
};


#endif
