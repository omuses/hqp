/**
 * @file Hxi_mxArray.h
 *   Native %mxArray for compiling a Simulink(R) S-function
 *   with Hqp. mxArray is used for passing parameters to an S-function.
 *
 * (Simulink is a registered trademarks of The MathWorks, Inc.)
 *
 * rf, 05/06/2001
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

/** Avoid multiple inclusion */
#if !defined(Hxi_mxArray_H)
#define Hxi_mxArray_H

#include "Hxi_sfun_types.h"
#include <vector>

using namespace std;

/// @name Data types
//@{ 
#define mxREAL		0x0001
//@}

/// @name Access macros
//@{
#define mxCreateDoubleMatrix(m, n, ty) 	new mxArray(m, n, ty)
#define mxDestroyArray(a) 		delete (a)
#define mxGetPr(a) 			(a)->getPr()
#define mxSetNumberOfElements(a, num) 	(a)->setNumberOfElements(num)
#define mxGetNumberOfElements(a) 	(a)->getNumberOfElements()
#define mxSetM(a, m) 			(a)->setM(m)
#define mxGetM(a) 			(a)->getM()
#define mxSetN(a, n) 			(a)->setN(n)
#define mxGetN(a) 			(a)->getN()
#define mxIsEmpty(a) 			(a)->isEmpty()
#define mxIsSparse(a) 			(a)->isSparse()
#define mxIsComplex(a) 			(a)->isComplex()
#define mxIsNumeric(a) 			(a)->isNumeric()
//@}

namespace Hxi {

/** Native %mxArray for Hqp. */
class mxArray {
protected:
  real_T 	 	_dummy; ///< argument for ADOL-C memory management
  vector<real_T> 	_data; 	///< data vector
  int_T 		_m;	///< number of rows (more generally first dim)
  int_T 		_n;	///< number of cols (more generally other dims)

public:
  /** Constuctor. */
  mxArray(int_T m, int_T n, int_T) : _data(m*n, _dummy) {
    _m = m;
    _n = n;
    for (int i = 0; i < (int)_data.size(); i++)
      _data[i] = 0.0;
  }

  /** Get address of first element of real data. */
  real_T *getPr() {
    return _data.size() > 0? &_data[0]: NULL;
  }
  /** Get address of first element of const real data. */
  const real_T *getPr() const {
    return _data.size() > 0? &_data[0]: NULL;
  }

  /** Set number of data elements. */
  int_T setNumberOfElements(int_T num) {
    _data.resize(num, _dummy);
    return _data.size();
  }
  /** Get number of data elements. */
  int_T getNumberOfElements() const {
    return _data.size();
  }
  /** Set number of rows. */
  int_T setM(int_T m) {
    return _m = m;
  }
  /** Get number of rows. */
  int_T getM() const {
    return _m;
  }
  /** Set number of columns. */
  int_T setN(int_T n) {
    return _n = n;
  }
  /** Get number of columns. */
  int_T getN() const {
    return _n;
  }

  /** Check if data object is empty. */
  bool isEmpty() const {
    return (_data.size() == 0);
  }
  /** Check if data object is sparse. */
  bool isSparse() const {
    return false;
  }
  /** Check if data object is complex. */
  bool isComplex() const {
    return false;
  }
  /** Check if data object is of numeric type. */
  bool isNumeric() const {
    return true;
  }
};

}; // namespace Hxi

#endif
