/*
 *  If_IntVec.h
 *   - integer vectors of the Meschach library
 *
 *  rf, 1/31/97
 *
 *  rf, 8/13/98
 *   - use typed Tcl 8 objects instead of strings
 */

/*
    Copyright (C) 1994--2000  Ruediger Franke

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

#ifndef If_IntVec_H
#define If_IntVec_H

// Include file with declaration of vector and matrix type.
#include "Meschach.h"

#include "If_Variable.h"

// typedef for a callback for write-access
//----------------------------------------
class If_IntVecWriteIf {

 public:
  virtual ~If_IntVecWriteIf() {}
  virtual int write(IVEC *newIVEC)=0;
};

template <class X>
class If_IntVecWriteCB: public If_IntVecWriteIf {

 protected:
  X	*_object;
  int	(X::*_write)(IVEC *);

 public:
  If_IntVecWriteCB(int (X::*n_write)(IVEC *), X *n_object)
    {
      assert(n_write != NULL && n_object != NULL);
      _write = n_write;
      _object = n_object;
    }
  int write(IVEC *newIVEC)
    {
      return (_object->*_write)(newIVEC);
    }
};

// foreward declaration for a IVECP
class IVECP;

// class declaration
//------------------
class If_IntVec: public If_Variable {

 protected:

  IVEC			**_varPtr;
  If_IntVecWriteIf  	*_callback;

  // define abstract methods of If_Variable
  //---------------------------------------
  int                  put(Tcl_Obj *CONST objPtr);
  int                  get();

 public:

  If_IntVec(char *ifName, IVEC **varPtr, If_IntVecWriteIf *callback=NULL);
  If_IntVec(char *ifName, IVECP *varPtr, If_IntVecWriteIf *callback=NULL);
  If_IntVec(char *ifName, IVEC **varPtr, const char *mode);
  If_IntVec(char *ifName, IVECP *varPtr, const char *mode);
  ~If_IntVec();
};

#endif
