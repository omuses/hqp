/*
 * If.C: implementation of public interface functions
 *
 * rf, 11/09/00
 */

/*
    Copyright (C) 1994--2004  Ruediger Franke

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

#include "If.h"
#include "If_Element.h"

#include <string.h>
#include <stdarg.h>

//-----------------------------------------------------------------------
extern "C" void If_Log(const char *category, const char *message, ...)
{
  static char *logMessage = NULL;
  static int logMessageLength = 256;
  va_list args;
  int n;

  va_start(args, message);
  if (logMessage == NULL)
    logMessage = (char *)malloc(logMessageLength);
  n = vsnprintf(logMessage, logMessageLength, message, args);
  /* enlarge logMessage up to max 16 KB */
  if (n >= logMessageLength && n < 16*1024) {
    free(logMessage);
    logMessage = (char *)malloc(n+1);
    logMessageLength = n+1;
    vsnprintf(logMessage, logMessageLength, message, args);
  }
  va_end(args);

  Tcl_VarEval(theInterp, "puts {", category, ": ", logMessage, "}", NULL);
  Tcl_Eval(theInterp, "update idletasks");
}

//--------------------------------------------------------------------------
extern "C" int If_Init(Tcl_Interp *interp)
{
  theInterp = interp;
  return theInterp? TCL_OK: TCL_ERROR;
}

//--------------------------------------------------------------------------
extern "C" Tcl_Interp *If_Interp()
{
  return theInterp;
}

//-----------------------------------------------------------------------
extern "C" int If_SizeOfInt()
{
  return sizeof(int);
}

//-----------------------------------------------------------------------
extern "C" int If_SetInt(const char *name, int val)
{
  if (!theInterp)
    return IF_ERROR;

#if 0
  // unfortunately Tcl_EvalObjv was not available under Tcl 8.0
  Tcl_Obj *objv[2];

  objv[0] = Tcl_NewStringObj((char *)name, -1);
  objv[1] = Tcl_NewIntObj(val);

  int retcode;
  retcode = Tcl_EvalObjv(theInterp, 2, objv, 0);

  Tcl_DecrRefCount(objv[0]);
  Tcl_DecrRefCount(objv[1]);

  if (retcode != TCL_OK)
    return IF_ERROR;
#else
  char valstr[50];
  sprintf(valstr, "%d", val);
  if (Tcl_VarEval(theInterp, (char *)name, " ", valstr, NULL) != TCL_OK)
    return IF_ERROR;
#endif

  Tcl_ResetResult(theInterp); // reset result as val was accepted
  return IF_OK;
}

//-----------------------------------------------------------------------
extern "C" int If_GetInt(const char *name, int *val)
{
  if (!theInterp)
    return IF_ERROR;

  if (Tcl_Eval(theInterp, (char *)name) != TCL_OK)
    return IF_ERROR;

  if (Tcl_GetIntFromObj(theInterp, Tcl_GetObjResult(theInterp), val)
      != TCL_OK)
    return IF_ERROR;

  Tcl_ResetResult(theInterp); // reset result as val is already returned
  return IF_OK;
}

//-----------------------------------------------------------------------
extern "C" int If_SizeOfReal()
{
  return sizeof(Real);
}

//-----------------------------------------------------------------------
extern "C" int If_SetReal(const char *name, Real val)
{
  if (!theInterp)
    return IF_ERROR;

#if 0
  // unfortunately Tcl_EvalObjv was not available under Tcl 8.0
  Tcl_Obj *objv[2];

  objv[0] = Tcl_NewStringObj((char *)name, -1);
  objv[1] = Tcl_NewDoubleObj(val);

  int retcode;
  retcode = Tcl_EvalObjv(theInterp, 2, objv, 0);

  Tcl_DecrRefCount(objv[0]);
  Tcl_DecrRefCount(objv[1]);

  if (retcode != TCL_OK)
    return IF_ERROR;
#else
  char valstr[50];
  sprintf(valstr, "%g", val);
  if (Tcl_VarEval(theInterp, (char *)name, " ", valstr, NULL) != TCL_OK)
    return IF_ERROR;
#endif

  Tcl_ResetResult(theInterp); // reset result as val was accepted
  return IF_OK;
}

//-----------------------------------------------------------------------
extern "C" int If_GetReal(const char *name, Real *val)
{
  if (!theInterp)
    return IF_ERROR;

  if (Tcl_Eval(theInterp, (char *)name) != TCL_OK)
    return IF_ERROR;

  if (Tcl_GetDoubleFromObj(theInterp, Tcl_GetObjResult(theInterp), val)
      != TCL_OK)
    return IF_ERROR;

  Tcl_ResetResult(theInterp); // reset result as val is already returned
  return IF_OK;
}

//-----------------------------------------------------------------------
extern "C" int If_SetString(const char *name, const char *val)
{
  if (!theInterp)
    return IF_ERROR;

  if (Tcl_VarEval(theInterp, (char *)name, " {", (char *)val, "}", 
		  NULL) != TCL_OK)
    return IF_ERROR;
  
  return IF_OK;
}

//-----------------------------------------------------------------------
extern "C" int If_GetString(const char *name, const char **val)
{
  if (!theInterp)
    return IF_ERROR;

  if (Tcl_Eval(theInterp, (char *)name) != TCL_OK) {
    val = NULL;
    return IF_ERROR;
  }
  *val = Tcl_GetStringResult(theInterp);
  return IF_OK;
}

//-----------------------------------------------------------------------
extern "C" int If_Eval(const char *command)
{
  if (!theInterp)
    return IF_ERROR;

  if (Tcl_Eval(theInterp, command) != TCL_OK)
    return IF_ERROR;

  return IF_OK;
}


//=======================================================================
