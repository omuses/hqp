/*
 * If_Float.h --
 *     - an If_Variable for floats (currently supported: double)
 *
 *  rf, 2/7/97
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

#ifndef If_Float_H
#define If_Float_H

#include "If_Variable.h"

typedef double If_Float_t;

// callback for write-access
//--------------------------
class If_FloatWriteIf {

 public:
  virtual ~If_FloatWriteIf() {}
  virtual int write(If_Float_t newVal)=0;
};

template <class X>
class If_FloatWriteCB: public If_FloatWriteIf {

 protected:
  X	*_object;
  int	(X::*_write)(If_Float_t);

 public:
  If_FloatWriteCB(int (X::*n_write)(If_Float_t), X *n_object)
    {
      assert(n_write != NULL && n_object != NULL);
      _write = n_write;
      _object = n_object;
    }
  int write(If_Float_t newVal)
    {
      return (_object->*_write)(newVal);
    }
};


// class declaration
//------------------
class If_Float: public If_Variable {

 protected:

  If_Float_t          *_varPtr;
  If_FloatWriteIf     *_callback;

  // define abstract methods of If_Variable
  //---------------------------------------
  int                put(Tcl_Obj *CONST objPtr);
  int                get();

 public:

  If_Float(char *ifName, If_Float_t *varPtr,
	    If_FloatWriteIf *callback=NULL);
  ~If_Float();
};


#endif
