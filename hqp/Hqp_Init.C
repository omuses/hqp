/*
 * central command for setting up several tasks
 * 
 * rf, 6/21/94
 *
 * rf, 4/13/99
 *  rename hqp.tcl to hqp_solve.tcl for initial evaluation
 *
 * rf, 12/25/99
 *  add support for Tcl package handling
 *
 * rf, 00/11/01
 *  add support for MSC
 *
 * rf, 09/04/11
 *  add mixed integer solver (Hqp_MipSolver)
 *
 */

/*
    Copyright (C) 1994--2014  Ruediger Franke

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

#include <signal.h>
#include <tcl.h>

#include <If_List.h>
#include <If_Element.h>
#include <If_Procedure.h>
#include <If_String.h>
#include <If_Class.h>
#include <If_Module.h>

#include "Hqp.h"
#include "Hqp_SqpProgram.h"
#include "Hqp_SqpSolver.h"
#include "Hqp_MipSolver.h"

#include "Hqp_SqpPowell.h"

/** static list of general Hqp interface elements */
static If_List _ifList;

/** Hold a pointer to the Hqp_SqpProgram */
Hqp_SqpProgram *theSqpProgram = NULL;

/** Hold a pointer to the Hqp_SqpSolver */
Hqp_SqpSolver *theSqpSolver = NULL;

/** Hold a pointer to the Hqp_MipSolver */
Hqp_MipSolver *theMipSolver = NULL;

#ifdef IF_CLASS_STATIC
//--------------------------------------------------------------------------
// Ensure linkage of modules
// where automatic inclusion of all objects does not work.
//--------------------------------------------------------------------------
#if defined(HQP_WITH_LPSOLVE)
#include "Hqp_LPSolve.h"
#endif
#if defined(HQP_WITH_ASCEND)
#include "Prg_ASCEND.h"
#endif
#if defined(HQP_WITH_CUTE)
#include "Prg_CUTE.h"
#include "Prg_CUTE_ST.h"
#endif
#include "Hqp_SqpSchittkowski.h"
#include "Hqp_HL_BFGS.h"
#include "Hqp_HL_DScale.h"
#include "Hqp_HL_Gerschgorin.h"
#include "Hqp_IpRedSpBKP.h"
#include "Hqp_IpSpBKP.h"
#include "Hqp_IpLQDOCP.h"
#include "Hqp_IpSpSC.h"
#include "Hqp_IpPardiso.h"
#include "Hqp_IpsFranke.h"
#include "Hqp_IpsMehrotra.h"
#include "Hqp_Client.h"

static void Hqp_ClassAlloc()
{
  IF_CLASS_ALLOC("Client", Hqp_Client, Hqp_Solver);
  IF_CLASS_ALLOC("Mehrotra", Hqp_IpsMehrotra, Hqp_Solver);
  IF_CLASS_ALLOC("Franke", Hqp_IpsFranke, Hqp_Solver);
  IF_CLASS_ALLOC("Pardiso", Hqp_IpPardiso, Hqp_IpMatrix);
  IF_CLASS_ALLOC("SpSC", Hqp_IpSpSC, Hqp_IpMatrix);
  IF_CLASS_ALLOC("LQDOCP", Hqp_IpLQDOCP, Hqp_IpMatrix);
  IF_CLASS_ALLOC("SpBKP", Hqp_IpSpBKP, Hqp_IpMatrix);
  IF_CLASS_ALLOC("RedSpBKP", Hqp_IpRedSpBKP, Hqp_IpMatrix);
  IF_CLASS_ALLOC("Gerschgorin", Hqp_HL_Gerschgorin, Hqp_HL);
  IF_CLASS_ALLOC("DScale", Hqp_HL_DScale, Hqp_HL);
  IF_CLASS_ALLOC("BFGS", Hqp_HL_BFGS, Hqp_HL);
  IF_CLASS_ALLOC("Schittkowski", Hqp_SqpSchittkowski, Hqp_SqpSolver);
  IF_CLASS_ALLOC("Powell", Hqp_SqpPowell, Hqp_SqpSolver);
#if defined(HQP_WITH_CUTE)
  IF_CLASS_ALLOC("CUTE_ST", Prg_CUTE_ST, Hqp_SqpProgram);
  IF_CLASS_ALLOC("CUTE", Prg_CUTE, Hqp_SqpProgram);
#endif
#if defined(HQP_WITH_ASCEND)
  IF_CLASS_ALLOC("ASCEND", Prg_ASCEND, Hqp_SqpProgram);
#endif
#if defined(HQP_WITH_LPSOLVE)
  IF_CLASS_ALLOC("LPSolve", Hqp_LPSolve, Hqp_MipSolver);
#endif
}
#endif

//--------------------------------------------------------------------------
// return HQP's version
//--------------------------------------------------------------------------
const char *Hqp_Version = VERSION;

//--------------------------------------------------------------------------
// fulfill the Meschach copyright
//--------------------------------------------------------------------------
static void m_version_cmd()
{
  m_version();
}

//--------------------------------------------------------------------------
static void signal_handler(int code)
{
  switch (code) {

  case SIGINT:
    Tcl_Eval(theInterp, "hqp_exit {signal interrupt}");
    fprintf(stderr, "HQP %s: %s\n", Hqp_Version, "signal interrupt");
    break;
#if !defined(_MSC_VER) && !defined(__MINGW32__)
  case SIGXCPU:
    signal(SIGXCPU, SIG_IGN);
    Tcl_Eval(theInterp, "hqp_exit {signal cputime}");
    fprintf(stderr, "HQP %s: %s\n", Hqp_Version, "signal cputime");
    break;
#endif
  default:
    break;
  }

  fprintf(stderr, "%s\n", Tcl_GetStringResult(theInterp));
  exit(0);
}

//--------------------------------------------------------------------------
extern "C" HQP_API int Hqp_InitSignalHandler() 
{
  // install a handler for signal interrupt
  signal(SIGINT, &signal_handler);
#if !defined(_MSC_VER) && !defined(__MINGW32__)
  signal(SIGXCPU, &signal_handler);
#endif

  // ignore signals from floating point arithmetics
  // (needed for ADOL-C 1.7 on Alpha with OSF3.2, OSF4.0)
  signal(SIGFPE, SIG_IGN);
}

//--------------------------------------------------------------------------
extern "C" HQP_API int Hqp_Init(Tcl_Interp *interp) 
{
  // provide Tcl package Hqp
  if (Tcl_InitStubs(interp, "8.1", 0) == NULL ||
      Tcl_PkgProvide(interp, "Hqp", (char *)Hqp_Version) != TCL_OK) {
    return TCL_ERROR;
  }

  // initialize global reference to Tcl interpreter
  theInterp = interp;

  // disable Meschach's error counting
  // (otherwise program would exit if counter reaches 100)
  count_errs(0);

  // allocate interface modules
# ifdef IF_CLASS_STATIC
  Hqp_ClassAlloc();
# endif

  // create initial state
  theSqpProgram = NULL;
  theSqpSolver = new Hqp_SqpPowell;

  _ifList.append(new IF_MODULE("prg_name", &theSqpProgram, Hqp_SqpProgram));
  _ifList.append(new IF_MODULE("sqp_solver", &theSqpSolver, Hqp_SqpSolver));
  _ifList.append(new IF_MODULE("mip_solver", &theMipSolver, Hqp_MipSolver));
  _ifList.append(new If_String("hqp_version", &Hqp_Version));
  _ifList.append(new If_Procedure("m_version", &m_version_cmd));

  // evaluate hqp_solve.tcl
  extern char *hqp_solve;
  if (Tcl_Eval(interp, hqp_solve) == TCL_ERROR) {
    fprintf(stderr,
	    "Evaluation of built-in code failed: %s\n", Tcl_GetStringResult(interp));
  }

  // evaluate file ~/.hqprc if it exists
  if (Tcl_Eval(interp, "if {[file exists ~/.hqprc]} {source ~/.hqprc}")
      != TCL_OK) {
    fprintf(stderr,
	    "Evaluation of ~/.hqprc failed: %s\n", Tcl_GetStringResult(interp));
  }

  return TCL_OK;
}
