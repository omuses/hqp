/*
 *  If_String.h
 *     - If_Variable for C++ string type
 *       (Everything is inlined so that only programs that use it
 *        have to link a Std C++ library.)
 *
 *  rf, 3/17/01
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

#ifndef If_String_H
#define If_String_H

#include "If_Variable.h"

#include <string>

typedef std::string If_String_t;

// callback for write-access
//--------------------------
class If_StringWriteIf {

 public:
  virtual ~If_StringWriteIf() {}
  virtual int write(If_String_t newVal)=0;
};

template <class X>
class If_StringWriteCB: public If_StringWriteIf {

 protected:
  X	*_object;
  int	(X::*_write)(If_String_t);

 public:
  If_StringWriteCB(int (X::*n_write)(If_String_t), X *n_object)
    {
      assert(n_write != NULL && n_object != NULL);
      _write = n_write;
      _object = n_object;
    }
  int write(If_String_t newVal)
    {
      return (_object->*_write)(newVal);
    }
};

// class declaration
//------------------
class If_String: public If_Variable {

 protected:
  If_String_t		*_varPtr;
  If_StringWriteIf	*_callback;

  // define abstract methods of If_Variable
  //---------------------------------------
  int                  put(Tcl_Obj *CONST objPtr);
  int                  get();

 public:
  If_String(const char *ifName, If_String_t *varPtr,
	    const char *mode = "rw");
  If_String(const char *ifName, If_String_t *varPtr,
	    If_StringWriteIf *callback);
  ~If_String();
};

/*
 *  If_String.C -- class definition
 */

//--------------------------------------------------------------------------
If_String::If_String(const char *ifName, If_String_t *varPtr,
		     const char *mode)
  :If_Variable(ifName, mode)
{
  _varPtr = varPtr;
  _callback = NULL;
}

//--------------------------------------------------------------------------
If_String::If_String(const char *ifName, If_String_t *varPtr,
		     If_StringWriteIf *callback)
  :If_Variable(ifName)
{
  _varPtr = varPtr;
  _callback = callback;
}

//--------------------------------------------------------------------------
If_String::~If_String()
{
  delete _callback;
}

//--------------------------------------------------------------------------
int If_String::put(Tcl_Obj *CONST objPtr)
{
  char *str_value;
  int str_length;

  // parse the new value
  //--------------------
  str_value = Tcl_GetStringFromObj(objPtr, &str_length);
  if (str_value == NULL)
    return TCL_ERROR;

  If_String_t value(str_value);

  // use the new value
  //------------------
  if (_callback) {
    if (_callback->write(value) != IF_OK) {
      Tcl_AppendResult(theInterp, "writing of ", _ifName, " rejected", NULL);
      return TCL_ERROR;
    }
  }
  else
    *_varPtr = value;

  return TCL_OK;
}

//--------------------------------------------------------------------------
int If_String::get()
{
  // (need to cast to char* as Tcl does not use const qualifier)
  Tcl_Obj *objPtr = Tcl_NewStringObj((char *)(*_varPtr).c_str(),
				     (*_varPtr).length());

  Tcl_SetObjResult(theInterp, objPtr);

  return TCL_OK;
}


//==========================================================================


#endif
