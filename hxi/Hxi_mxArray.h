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
    Copyright (C) 1994--2005  Ruediger Franke

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
#include <string>
#include <stdlib.h>

using namespace std;

/// @name Data types
//@{ 
#define mxREAL		0x0001
#define mxCHAR		0x0002
//@}

/// @name Access macros
//@{
#define mxCreateDoubleMatrix(m, n, ty) 	new mxArray(m, n, ty)
#define mxCreateString(s) 		new mxArray(s)
#define mxDestroyArray(a) 		delete (a)
#define mxArrayToString(a) 		(a)->toString()
#define mxFree(p) 			free(p)
#define mxGetPr(a) 			(a)->getPr()
#define mxSetNumberOfElements(a, num) 	(a)->setNumberOfElements(num)
#define mxGetNumberOfElements(a) 	(a)->getNumberOfElements()
#define mxSetM(a, m) 			(a)->setM(m)
#define mxGetM(a) 			(a)->getM()
#define mxSetN(a, n) 			(a)->setN(n)
#define mxGetN(a) 			(a)->getN()
#define mxIsEmpty(a) 			(a)->isEmpty()
#define mxIsChar(a) 			(a)->isChar()
#define mxIsDouble(a) 			(a)->isDouble()
#define mxIsSparse(a) 			(a)->isSparse()
#define mxIsComplex(a) 			(a)->isComplex()
#define mxIsNumeric(a) 			(a)->isNumeric()
//@}

namespace Hxi {

/** Native %mxArray for Hqp. */
class mxArray {
protected:
  real_T 	 	_dummy; ///< argument for ADOL-C memory management
  int_T			_type; 	///< element data type
  vector<real_T> 	_realData; 	///< data vector for array of reals
  string 		_charData; 	///< data for char array
  int_T 		_m;	///< number of rows (more generally first dim)
  int_T 		_n;	///< number of cols (more generally other dims)

public:
  /** Constuctor for real array. */
  mxArray(int_T m, int_T n, int_T) : _realData(m*n, _dummy) {
    _type = mxREAL;
    _m = m;
    _n = n;
    for (int i = 0; i < (int)_realData.size(); i++)
      _realData[i] = 0.0;
  }

  /** Constuctor for char array. */
  mxArray(const char *s) : _charData(s) {
    _type = mxCHAR;
  }

  /** Get address of first element of real data. */
  char *toString() const {
    char *str;
    int i;
    switch (_type) {
    case mxCHAR:
      str = (char *)malloc(_charData.size() + 1);
      for (i = 0; i < _charData.size(); i++)
        str[i] = _charData[i];
      str[i] = '\0';
      break;
    default:
      str = NULL;
      break;
    }
    return str;
  }

  /** Get address of first element of real data. */
  real_T *getPr() {
    return _realData.size() > 0? &_realData[0]: NULL;
  }
  /** Get address of first element of const real data. */
  const real_T *getPr() const {
    return _realData.size() > 0? &_realData[0]: NULL;
  }

  /** Set number of data elements. */
  int_T setNumberOfElements(int_T num) {
    _realData.resize(num, _dummy);
    return _realData.size();
  }
  /** Get number of data elements. */
  int_T getNumberOfElements() const {
    return _realData.size();
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
    return (_realData.size() == 0 && _charData.size() == 0);
  }
  /** Check if data object is char array. */
  bool isChar() const {
    return (_type == mxCHAR);
  }
  /** Check if data object is real array. */
  bool isDouble() const {
    return (_type == mxREAL);
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
