/*
 * If_Element.h --
 *     - abstract base class for interface elements to Tcl
 *       (e.g. variables, commands)
 *     - provides:
 *        + a textual name for every interface element
 *        + a static Tcl_Interp for all elements
 *
 *  rf, 22/6/94, rev. 1.0
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

#ifndef If_Element_H
#define If_Element_H

extern "C" {
#ifdef VARARGS
#undef VARARGS
#endif
#include <tcl.h>
#undef VARARGS
}

#include "If.h"
#include "If_ListElement.h"


// global Tcl interpreter, that must be initialized by an application
//-------------------------------------------------------------------
extern Tcl_Interp *theInterp;

//-----------------------------------------------------
class If_Element: public If_ListElement {

 protected:
  char                  *_ifName;

 public:
                        If_Element(const char *ifName);
  virtual               ~If_Element();
};


#endif
