/*
 *  If_FloatMat.h
 *   - matrices of floats (currently supported: Meschach MAT*)
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

#ifndef If_FloatMat_H
#define If_FloatMat_H

// Include file with declaration of vector and matrix type.
#include "Meschach.h"

#include "If_Variable.h"

// callback for write-access
//--------------------------
class If_FloatMatWriteIf {

 public:
  virtual ~If_FloatMatWriteIf() {}
  virtual int write(MAT *newMAT)=0;
};

template <class X>
class If_FloatMatWriteCB: public If_FloatMatWriteIf {

 protected:
  X	*_object;
  int	(X::*_write)(MAT *);

 public:
  If_FloatMatWriteCB(int (X::*n_write)(MAT *), X *n_object)
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

// foreward declaration for a MATP
class MATP;

// class declaration
//------------------
class If_FloatMat: public If_Variable {

 protected:

  MAT		**_varPtr;
  If_FloatMatWriteIf *_callback;

  // define abstract methods of If_Variable
  //---------------------------------------
  int           put(Tcl_Obj *CONST objPtr);
  int           get();

 public:

  If_FloatMat(char *ifName, MAT **varPtr,
	      If_FloatMatWriteIf *callback=NULL);
  If_FloatMat(char *ifName, MATP *varPtr,
	      If_FloatMatWriteIf *callback=NULL);
  ~If_FloatMat();
};

#endif
