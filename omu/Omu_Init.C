/*
 * Omu_Init: initialize the Omuses package
 *
 * rf, 1/10/97
 *
 * rf, 03/22/00
 *  add support for Tcl package handling
 */

/*
    Copyright (C) 1997--2002  Ruediger Franke

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
#include <If_Proc.h>

#include "Hqp_Omuses.h"

#ifdef IF_CLASS_STATIC
//--------------------------------------------------------------------------
// ensure linkage of modules
//--------------------------------------------------------------------------
#include "Omu_IntEuler.h"
#include "Omu_IntRK4.h"
#include "Omu_IntOdeTs.h"
#include "Omu_IntDopri5.h"
#include "Omu_IntIMP.h"
#include "Omu_IntSDIRK.h"
#ifdef OMU_WITH_FORTRAN
#include "Omu_IntRKsuite.h"
#include "Omu_IntDASPK.h"
#endif
#ifdef OMU_WITH_MEX
#include "Prg_SFunctionOpt.h"
#include "Prg_SFunctionEst.h"
#endif

static void Omu_ClassAlloc()
{
  IF_CLASS_ALLOC("Euler", Omu_IntEuler, Omu_Integrator);
  IF_CLASS_ALLOC("RK4", Omu_IntRK4, Omu_Integrator);
  IF_CLASS_ALLOC("OdeTs", Omu_IntOdeTs, Omu_Integrator);
  IF_CLASS_ALLOC("Dopri5", Omu_IntDopri5, Omu_Integrator);
  IF_CLASS_ALLOC("IMP", Omu_IntIMP, Omu_Integrator);
  IF_CLASS_ALLOC("SDIRK", Omu_IntSDIRK, Omu_Integrator);
#ifdef OMU_WITH_FORTRAN
  IF_CLASS_ALLOC("DASPK", Omu_IntDASPK, Omu_Integrator);
  IF_CLASS_ALLOC("RKsuite", Omu_IntRKsuite, Omu_Integrator);
#endif
#ifdef OMU_WITH_MEX
  IF_CLASS_ALLOC("SFunctionOpt", Prg_SFunctionOpt, Omu_Program);
  IF_CLASS_ALLOC("SFunctionEst", Prg_SFunctionEst, Omu_Program);
#endif
}
#endif

//--------------------------------------------------------------------------
const char *Omu_Version = VERSION;

static int Omu_VersionCmd(int, char *[], char **result)
{
  *result = (char *)Omu_Version;
  return IF_OK;
}

//--------------------------------------------------------------------------
extern "C" int Omu_Init(Tcl_Interp *interp)
{
  // provide Tcl package Omuses
  if (Tcl_PkgRequire(interp, "Tcl", "8.0", 0) == NULL ||
      Tcl_PkgProvide(interp, "Omuses", (char *)Omu_Version) != TCL_OK) {
    return TCL_ERROR;
  }

  // initialize global reference to Tcl interpreter
  theInterp = interp;

  new If_Proc("omu_version", &Omu_VersionCmd);

  // allocate interface modules
# ifdef IF_CLASS_STATIC
  Omu_ClassAlloc();
# endif

  // create instance of Hqp_Omuses
  new Hqp_Omuses;

  return TCL_OK;
}
