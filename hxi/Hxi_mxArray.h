/*
 * Hxi_mxArray.h:
 *   mxArray for compiling a Simulink(R) S-function for HQP
 *
 * (Simulink is a registered trademark of The MathWorks, Inc.)
 *
 * rf, 05/06/2001
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

#if !defined(Hxi_mxArray_H)
#define Hxi_mxArray_H

#include "Hxi_sfun_types.h"

// Macros to access an mxArray (required for parameters)
#define mxSetPr(a, pr) 			(a)->setPr(pr)
#define mxGetPr(a) 			(a)->getPr()
#define mxSetNumberOfElements(a, num) 	(a)->setNumberOfElements(num)
#define mxGetNumberOfElements(a) 	(a)->getNumberOfElements()
#define mxIsEmpty(a) 			(a)->isEmpty()
#define mxIsSparse(a) 			(a)->isSparse()
#define mxIsComplex(a) 			(a)->isComplex()
#define mxIsNumeric(a) 			(a)->isNumeric()

/** mxArray for HQP. */
class mxArray {
protected:
  real_T 	*_data;
  int_T 	_size;

public:
  /** Constuctor. */
  mxArray() {
    _data = NULL;
    _size = 0;
  }

  /** Set address of first element of real data. */
  real_T *setPr(real_T *pr) {
    return _data = pr;
  }
  /** Get address of first element of real data. */
  real_T *getPr() const {
    return _data;
  }

  /** Set number of data elements. */
  int_T setNumberOfElements(int_T num) {
    return _size = num;
  }
  /** Get number of data elements. */
  int_T getNumberOfElements() const {
    return _size;
  }

  /** Check if data object is empty. */
  bool isEmpty() const {
    return (_data == NULL || _size == 0);
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

#endif
