/*
 * If.h: public declarations of the interface library
 *
 * rf, 10/4/95
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

#ifndef If_H
#define If_H

#include "Meschach.h"

/** define IF_API when compiling a Dynamic Link Library (DLL) */
#ifndef IF_API
#define IF_API 
#endif

/** Return code for successful calculation. */
#define IF_OK		0

/** Return code for failed calculation. */
#define IF_ERROR 	1


extern "C" {
  /** Return size of Int. */
  IF_API int If_SizeOfInt();

  /** Set an Int value. */
  IF_API int If_SetInt(const char *name, int val);

  /** Get an Int value. */
  IF_API int If_GetInt(const char *name, int &val);

  /** Return size of Real. */
  IF_API int If_SizeOfReal();

  /** Set a Real value. */
  IF_API int If_SetReal(const char *name, Real val);

  /** Get a Real value. */
  IF_API int If_GetReal(const char *name, Real &val);

  /** Evaluate a command. */
  IF_API int If_Eval(const char *command);

  /** Evaluate a command with string argument */
  IF_API int If_EvalStringArg(const char *command, const char *arg);

  /** Return the result string produced by the last If function call.
      After a failed calculation, the corresponding error message
      is returned. */
  IF_API const char *If_ResultString();
}

#endif
