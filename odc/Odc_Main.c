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

#include <tcl.h>

/*
 * prototypes for module initializations
 */

int Hqp_Init _ANSI_ARGS_((Tcl_Interp *interp));
int Omu_Init _ANSI_ARGS_((Tcl_Interp *interp));
int Odc_Init _ANSI_ARGS_((Tcl_Interp *interp));

/*
 * the main function
 */

int
main(argc, argv)
     int argc;
     char **argv;
{
  Tcl_Main(argc, argv, Tcl_AppInit);
  return 0;
}

/*
 * Tcl_AppInit function
 */

int
Tcl_AppInit(interp)
     Tcl_Interp *interp;
{
  if (Tcl_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }

  if (Hqp_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }
  Tcl_StaticPackage(interp, "Hqp", Hqp_Init, NULL);

  if (Omu_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }
  Tcl_StaticPackage(interp, "Omu", Omu_Init, NULL);

  if (Odc_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }
  Tcl_StaticPackage(interp, "Odc", Odc_Init, NULL);

  Tcl_SetVar(interp, "tcl_rcFileName", "~/.tclshrc", TCL_GLOBAL_ONLY);
  return TCL_OK;
}
