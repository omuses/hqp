/**
 *  @file If_Procedure.h
 *     binds a function pointer to an interface element
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

#ifndef If_Procedure_H
#define If_Procedure_H

#include "If_Element.h"

/** Interface procedure type */
typedef void If_Procedure_t();

/** Interface procedure */
class IF_API If_Procedure: public If_Element {

 protected:
  /** pointer to callback procedure */
  If_Procedure_t	*_proc;

  /** invoke callback procedure */
  int invoke(Tcl_Interp *, int objc, Tcl_Obj *CONST objv[]);

 public:
  /** constructor */
  If_Procedure(const char *ifName, If_Procedure_t *proc);
};


#endif
