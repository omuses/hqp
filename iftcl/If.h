/**
 * @file If.h
 *   public declarations of the interface library
 *
 * rf, 10/4/95
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

#ifndef If_H
#define If_H

#include <Meschach.h>
#include <tcl.h>

/** define IF_API when compiling a Dynamic Link Library (DLL) */
#ifndef IF_API
#define IF_API 
#endif

/** Return code for success. */
#define IF_OK		0

/** Return code for fail. */
#define IF_ERROR 	1

extern "C" {
  /** Initialize the command interface to use the specified Tcl_Interp.
      This function needs to be called by an extension being loaded
      to a Tcl application. It either returns TCL_OK or TCL_ERROR
      if the initialization fails. */
  IF_API int If_Init(Tcl_Interp *interp);

  /** Create a Tcl interpreter and initialize the path for the Tcl library
      based on argv[0]. This function needs to be called by a non Tcl
      application, i.e. an application that does not enter Tcl_Main.
      It is inlined as in this way the application determines,
      which Tcl version to use. */
  inline int If_CreateInterp(int argc, char *argv[])
  {
    char argv0 = '\0';
    Tcl_FindExecutable(argc > 0? argv[0]: &argv0);
    Tcl_Interp *interp = Tcl_CreateInterp();
    if (If_Init(interp) != TCL_OK || Tcl_Init(interp) != TCL_OK)
      return IF_ERROR;
    return IF_OK;
  }

  /** Return the Tcl interpreter used by the interface. */
  IF_API Tcl_Interp *If_Interp();

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

  /** Set a variables value from a string representation. */
  IF_API int If_SetString(const char *name, const char *val);

  /** Get the string representation of a variables value. */
  IF_API int If_GetString(const char *name, const char *&val);

  /** Evaluate a command. */
  IF_API int If_Eval(char *command);

  /** Return the result string produced by the last If function call.
      After a failed calculation, the corresponding error message
      is returned. */
  IF_API const char *If_ResultString();
}

#endif
