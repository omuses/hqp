/*
 *  If_FloatVec.h
 *   - vectors of floats (currently supported: Meschach VEC*)
 *
 *  rf, 2/7/97
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

#ifndef If_FloatVec_H
#define If_FloatVec_H

// Include file with declaration of vector and matrix type.
#include "Meschach.h"

#include "If_Variable.h"

// typedef for a callback for write-access
//----------------------------------------
class If_FloatVecWriteIf {

 public:
  virtual ~If_FloatVecWriteIf() {}
  virtual int write(VEC *newVEC)=0;
};

template <class X>
class If_FloatVecWriteCB: public If_FloatVecWriteIf {

 protected:
  X	*_object;
  int	(X::*_write)(VEC *);

 public:
  If_FloatVecWriteCB(int (X::*n_write)(VEC *), X *n_object)
    {
      assert(n_write != NULL && n_object != NULL);
      _write = n_write;
      _object = n_object;
    }
  int write(VEC *newVEC)
    {
      return (_object->*_write)(newVEC);
    }
};

// foreward declaration for a VECP
class VECP;

// class declaration
//------------------
class If_FloatVec: public If_Variable {

 protected:

  VEC			**_varPtr;
  If_FloatVecWriteIf  	*_callback;

  // define abstract methods of If_Variable
  //---------------------------------------
  int                  put(Tcl_Obj *CONST objPtr);
  int                  get();

 public:

  If_FloatVec(char *ifName, VEC **varPtr,
	      If_FloatVecWriteIf *callback=NULL);
  If_FloatVec(char *ifName, VECP *varPtr,
	      If_FloatVecWriteIf *callback=NULL);
  ~If_FloatVec();
};

#endif
