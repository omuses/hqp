/*
 * Odc_Init: initialize the package Odc (Omuses demo collection)
 *
 * rf, 1/14/97
 */

#include <tcl.h>
#include <If_String.h>

//--------------------------------------------------------------------------
const char *Odc_Version = "1.4";

//--------------------------------------------------------------------------
extern "C" int Odc_Init(Tcl_Interp *interp)
{
  // do package specific initializations

  new If_String("odc_version", &Odc_Version);

  return TCL_OK;
}
