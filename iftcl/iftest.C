/* 
 * Odc_Main.c --
 *
 *    main function for the Omuses demo collection
 *    (this file may be replaced by tclAppInit.c or tkAppInit.c of
 *     a Tcl/Tk distribution, provided that the modules Hqp, Omu, 
 *     and Odc are initialized)
 *
 *  rf, 2/6/97
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

#include <tcl.h>
#include "If_Int.h"
#include "If_Bool.h"
#include "If_Float.h"
#include "If_FloatVec.h"
#include "If_FloatMat.h"
#include "If_IntVec.h"

/*
 * the main function
 */

int
main(int argc, char **argv)
{
  Tcl_Main(argc, argv, Tcl_AppInit);
  return 0;
}

/*
 * Tcl_AppInit function
 */

Tcl_Interp *theInterp = NULL;
int if_int = 1234;
double if_float = 1.0/0.0;
bool if_bool = false;
VEC *if_floatVec = VNULL;
MAT *if_floatMat = MNULL;
IVEC *if_intVec = IVNULL;

int
Tcl_AppInit(Tcl_Interp *interp)
{
  if (Tcl_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }

  theInterp = interp;

  new If_Int("if_int", &if_int);
  new If_Bool("if_bool", &if_bool);
  new If_Float("if_float", &if_float);
  if_floatVec = v_get(5);
  new If_FloatVec("if_floatVec", &if_floatVec);
  if_floatMat = m_get(2,3);
  new If_FloatMat("if_floatMat", &if_floatMat);
  if_intVec = iv_get(3);
  new If_IntVec("if_intVec", &if_intVec);

  Tcl_SetVar(interp, "tcl_rcFileName", "~/.tclshrc", TCL_GLOBAL_ONLY);
  return TCL_OK;
}
