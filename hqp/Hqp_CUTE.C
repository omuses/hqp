/*
 * Tcl_AppInit and main function for CUTE
 *
 * rf, 4/10/97
 */

/*
    Copyright (C) 1994--1998  Ruediger Franke

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

#include <string.h>
#include <tcl.h>

/*
 * initialization of a Tcl application that runs HQP
 */

extern "C" int Hqp_Init (Tcl_Interp *interp);
extern "C" int Tcl_AppInit(Tcl_Interp *interp)
{
  if (Tcl_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }

  if (Hqp_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }

  return TCL_OK;
}

/*
 * a Tcl main function to be called from FORTRAN
 */

extern "C" void tclmn_(char *arg1)
{
  char cmdline[256];
  char *argv[2];
  int idx;

  strcpy(cmdline, "hqpmin");
  cmdline[6] = ' ';
  idx = 0;
  while (arg1[idx] != ' ' && idx < 248) {
    cmdline[7 + idx] = arg1[idx];
    idx++;
  }
  cmdline[7 + idx] = '\0';
  
  argv[0] = cmdline;
  argv[1] = cmdline + 7;
  
  Tcl_Main(2, argv, Tcl_AppInit);
}

extern "C" void tclmn(char *arg1)
{
  tclmn_(arg1);
}
