/**
 *  @file If_Method.h
 *     binds a pointer to a method to an interface element
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

#ifndef If_Method_H
#define If_Method_H

#include "If_Element.h"

/** Enable old signatures for callback methods. */
#define IF_METHOD_WITH_DEPRECATED 1

/** Interface method */
template<class ClassType>
class If_Method: public If_Element {

 protected:
  ClassType	*_object; 		 ///< object the method belongs to
  void		(ClassType::*_method)(); ///< method pointer

#if defined(IF_METHOD_WITH_DEPRECATED)
  /** deprecated pointer to a method taking argc/argv arguments */
  int	(ClassType::*_dep_method)(int, char *[], char **);
#endif

  /** invoke callback method */
  int invoke(Tcl_Interp *interp, int objc, Tcl_Obj *CONST []) {
    if (objc > 1) {
      Tcl_AppendResult(interp, "wrong # args, should be: ",
		       ifName(), NULL);
      return TCL_ERROR;
    }
#if defined(IF_METHOD_WITH_DEPRECATED)
    if (_dep_method) {
      return (_object->*_dep_method)(0, NULL, &interp->result);
    }
#endif
    (_object->*_method)();
    return TCL_OK;
  }	

 public:
  /** constructor */
  If_Method(const char *ifName,
	    void (ClassType::*method)(), ClassType *object)
    :If_Element(ifName) {
    _object = object;
    _method = method;
#if defined(IF_METHOD_WITH_DEPRECATED)
    _dep_method = NULL;
#endif
  }

#if defined(IF_METHOD_WITH_DEPRECATED)
  /** Deprecated constructor. It temporarily provides backwards
      compatibility for methods using call arguments and returning a result.
      However, no values will be passed to the method, that is argc == 0,
      argv == NULL; result can be used as before. Please redefine methods
      to not having call arguments and to throw errors with m_error(). */
  If_Method(const char *ifName,
	    int (ClassType::*method)(int argc, char *argv[], char **result),
	    ClassType *object)
    :If_Element(ifName) {
    _object = object;
    _method = NULL;
    _dep_method = method;
  }
#endif
};


#endif
