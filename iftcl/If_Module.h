/*
 *  If_Module.h
 *   - binds a pointer to a method to an interface command
 *   - class template (inline code to be compatible with most compilers)
 *
 *  rf, 7/19/94
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

#ifndef If_Module_H
#define If_Module_H

#include "If_Element.h"
#include "If_Class.h"

// use this macro to create a new If_Module
//-----------------------------------------
#define IF_MODULE(id, ptr, base) \
  If_Module<base>(id, ptr, theIf_ClassList_##base)


template<class X>
class If_Module: public If_Element {

 protected:

  X	**_module;
  If_ClassList<X>*& _list;

  /** create requested module */
  int invoke(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
  {
    char *fallback;

    if (objc == 1) {
      if (*_module)
	Tcl_AppendResult(interp, (*_module)->name(), NULL);
      else
	Tcl_AppendResult(interp, "No module chosen.", NULL);
      return IF_OK;
    }
    else if (objc == 2) {
      if (!_list) {
	Tcl_AppendResult(interp, "No modules available.", NULL);
	return IF_ERROR;
      }

      // store the current module name
      if (*_module)
	fallback = (*_module)->name();
      else
	fallback = NULL;

      // delete the current module to free interface elements
      delete *_module;

      // create the requested module
      *_module = _list->createObject(Tcl_GetString(objv[1]));
      if (!*_module) {
	Tcl_AppendResult(interp, _list->cand_names(), NULL);
	if (fallback)
	  *_module = _list->createObject(fallback);
	return IF_ERROR;
      }
      return IF_OK;
    }

    Tcl_AppendResult(interp, "wrong # of args, should be: ",
		     ifName(), " [new module name]", NULL);
    return IF_ERROR;
  }	

 public:

  /** constructor */
  If_Module(const char *ifName, X **module, If_ClassList<X>*& list)
    :If_Element(ifName), _list(list) {
    _module = module;
  }
};


#endif
