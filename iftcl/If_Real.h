/*
 * If_Real.h --
 *     - an If_Variable for real numbers (currently supported: double)
 *
 *  rf, 2/7/97
 *
 *  rf, 8/13/98
 *   - use typed Tcl 8 objects instead of strings
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

#ifndef If_Real_H
#define If_Real_H

#include "If_Variable.h"

typedef double If_Real_t;

// callback for write-access
//--------------------------
class If_RealWriteIf {

 public:
  virtual ~If_RealWriteIf() {}
  virtual int write(If_Real_t newVal)=0;
};

template <class X>
class If_RealWriteCB: public If_RealWriteIf {

 protected:
  X	*_object;
  int	(X::*_write)(If_Real_t);

 public:
  If_RealWriteCB(int (X::*n_write)(If_Real_t), X *n_object)
    {
      assert(n_write != NULL && n_object != NULL);
      _write = n_write;
      _object = n_object;
    }
  int write(If_Real_t newVal)
    {
      return (_object->*_write)(newVal);
    }
};


// class declaration
//------------------
class If_Real: public If_Variable {

 protected:

  If_Real_t          *_varPtr;
  If_RealWriteIf     *_callback;

  // define abstract methods of If_Variable
  //---------------------------------------
  int                put(Tcl_Obj *CONST objPtr);
  int                get();

 public:

  If_Real(const char *ifName, If_Real_t *varPtr, const char *mode = "rw");
  If_Real(const char *ifName, If_Real_t *varPtr, If_RealWriteIf *callback);
  ~If_Real();
};


#endif
