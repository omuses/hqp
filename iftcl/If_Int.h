/*
 *  If_Int.h
 *     - an If_Variable for type Int
 *     - writing may be protected by an optinal write callback,
 *       a pointer to it is passed to the constructor and
 *       deleted by the destuctor
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

#ifndef If_Int_H
#define If_Int_H

#include "If_Variable.h"

typedef int If_Int_t;

// callback for write-access
//--------------------------
class If_IntWriteIf {

 public:
  virtual ~If_IntWriteIf() {}
  virtual int write(If_Int_t newVal)=0;
};

template <class X>
class If_IntWriteCB: public If_IntWriteIf {

 protected:
  X	*_object;
  int	(X::*_write)(If_Int_t);

 public:
  If_IntWriteCB(int (X::*n_write)(If_Int_t), X *n_object)
    {
      assert(n_write != NULL && n_object != NULL);
      _write = n_write;
      _object = n_object;
    }
  int write(If_Int_t newVal)
    {
      return (_object->*_write)(newVal);
    }
};

// class declaration
//------------------
class If_Int: public If_Variable {

 protected:
  If_Int_t	*_varPtr;
  If_IntWriteIf	*_callback;

  // define abstract methods of If_Variable
  //---------------------------------------
  int                  put(Tcl_Obj *CONST objPtr);
  int                  get();

 public:
  If_Int(char *ifName, If_Int_t *varPtr, If_IntWriteIf *callback=NULL);
  ~If_Int();
};


#endif
