/*
 *  If_RealMat.h
 *   - matrices of reals (currently supported: Meschach MAT* and MATP)
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

#ifndef If_RealMat_H
#define If_RealMat_H

// Include file with declaration of vector and matrix type.
#include "Meschach.h"

#include "If_Variable.h"

// callback for write-access
//--------------------------
class If_RealMatWriteIf {

 public:
  virtual ~If_RealMatWriteIf() {}
  virtual int write(MAT *newMAT)=0;
};

template <class X>
class If_RealMatWriteCB: public If_RealMatWriteIf {

 protected:
  X	*_object;
  int	(X::*_write)(MAT *);

 public:
  If_RealMatWriteCB(int (X::*n_write)(MAT *), X *n_object)
    {
      assert(n_write != NULL && n_object != NULL);
      _write = n_write;
      _object = n_object;
    }
  int write(MAT *newMAT)
    {
      return (_object->*_write)(newMAT);
    }
};


// class declaration
//------------------
class If_RealMat: public If_Variable {

 protected:

  MAT		**_varPtr;
  If_RealMatWriteIf *_callback;

  // define abstract methods of If_Variable
  //---------------------------------------
  int           put(Tcl_Obj *CONST objPtr);
  int           get();

 public:

  If_RealMat(const char *ifName, MAT **varPtr, const char *mode = "rw");
  If_RealMat(const char *ifName, MATP *varPtr, const char *mode = "rw");
  If_RealMat(const char *ifName, MAT **varPtr, If_RealMatWriteIf *callback);
  If_RealMat(const char *ifName, MATP *varPtr, If_RealMatWriteIf *callback);
  ~If_RealMat();
};

#endif
