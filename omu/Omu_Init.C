/*
 * Omu_Init: initialize the Omuses package
 *
 * rf, 1/10/97
 *
 * rf, 03/22/00
 *  add support for Tcl package handling
 */

/*
    Copyright (C) 1997--2014  Ruediger Franke

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
#include <If_String.h>

#include "Omu.h"

#ifdef IF_CLASS_STATIC
//--------------------------------------------------------------------------
// ensure linkage of modules
//--------------------------------------------------------------------------
#include "Omu_IntEuler.h"
#include "Omu_IntRK4.h"
#ifdef OMU_WITH_ADOLC
#include "Omu_IntOdeTs.h"
#endif
#include "Omu_IntDopri5.h"
#include "Omu_IntIMP.h"
#include "Omu_IntGRK4.h"
#include "Omu_IntSDIRK.h"
#ifdef OMU_WITH_FORTRAN
#include "Omu_IntRKsuite.h"
#include "Omu_IntDASPK.h"
#endif
#include "Hqp_Omuses.h"
#include "Prg_DynamicOpt.h"
#include "Prg_DynamicEst.h"

static void Omu_ClassAlloc()
{
  IF_CLASS_ALLOC("Euler", Omu_IntEuler, Omu_Integrator);
  IF_CLASS_ALLOC("RK4", Omu_IntRK4, Omu_Integrator);
#ifdef OMU_WITH_ADOLC
  IF_CLASS_ALLOC("OdeTs", Omu_IntOdeTs, Omu_Integrator);
#endif
  IF_CLASS_ALLOC("Dopri5", Omu_IntDopri5, Omu_Integrator);
  IF_CLASS_ALLOC("IMP", Omu_IntIMP, Omu_Integrator);
  IF_CLASS_ALLOC("GRK4", Omu_IntGRK4, Omu_Integrator);
  IF_CLASS_ALLOC("SDIRK", Omu_IntSDIRK, Omu_Integrator);
#ifdef OMU_WITH_FORTRAN
  IF_CLASS_ALLOC("DASPK", Omu_IntDASPK, Omu_Integrator);
  IF_CLASS_ALLOC("RKsuite", Omu_IntRKsuite, Omu_Integrator);
#endif
  IF_CLASS_ALLOC("Omuses", Hqp_Omuses, Hqp_SqpProgram);
  IF_CLASS_ALLOC("DynamicOpt", Prg_DynamicOpt, Omu_Program);
  IF_CLASS_ALLOC("DynamicEst", Prg_DynamicEst, Omu_Program);
  IF_CLASS_ALLOC("SFunctionOpt", Prg_SFunctionOpt, Omu_Program);
  IF_CLASS_ALLOC("SFunctionEst", Prg_SFunctionEst, Omu_Program);
}
#endif

//--------------------------------------------------------------------------
const char *Omu_Version = VERSION;

//--------------------------------------------------------------------------
extern "C" int Hqp_Init(Tcl_Interp *interp);
extern "C" int OMU_API Omu_Init(Tcl_Interp *interp)
{
  // initialize stubs
  if (Tcl_InitStubs(interp, "8.1", 0) == NULL) {
    return TCL_ERROR;
  }

  // initialize Hqp if it is not already present
  // (for now do the initialization directly by calling Hqp_Init,
  //  instead of using Tcl_PkgRequire, to force linkage of Hqp)
  if (Tcl_PkgPresent(interp, "Hqp", VERSION, 1) == NULL) {
    if (Hqp_Init(interp) != TCL_OK)
      return TCL_ERROR;
  }
  // provide Tcl package Omuses
  if (/* Tcl_PkgRequire(interp, "Hqp", VERSION, 1) == NULL || */
      Tcl_PkgProvide(interp, "Omuses", (char *)Omu_Version) != TCL_OK) {
    return TCL_ERROR;
  }

  // initialize version command
  new If_String("omu_version", &Omu_Version);

  // allocate interface modules
# ifdef IF_CLASS_STATIC
  Omu_ClassAlloc();
# endif

  // use one command prg_name for Hqp and Omuses
  if (Tcl_Eval(interp, "rename prg_name _prg_name_hqp; \
                        proc prg_name {{name {}}} { \
                          if [catch {eval _prg_name_hqp $name} result] { \
                            _prg_name_hqp Omuses; \
                            _prg_name_omu $name \
                          } else { \
                            if {$result == {Omuses}} { \
                              return [_prg_name_omu] \
                            } else { \
                              return $result \
                            } \
                          } \
                        }") != TCL_OK) {
    return TCL_ERROR;
  }

  // evaluate ../hxi/fmi.tcl
  extern char *fmi;
  if (Tcl_Eval(interp, fmi) == TCL_ERROR) {
    fprintf(stderr,
	    "Evaluation of built-in code failed: %s\n",
            Tcl_GetStringResult(interp));
  }

  return TCL_OK;
}

//--------------------------------------------------------------------------
extern "C" int OMU_API Omu_Unload(Tcl_Interp *interp, int flags)
{
  // do package specific unloads

  return TCL_OK;
}
