/*
 * Odc_Init: initialize the package Odc (Omuses demo collection)
 *
 * rf, 1/14/97
 */

#include <tcl.h>
#include <If_Proc.h>

//--------------------------------------------------------------------------
const char *Odc_Version = "1.2";

static int Odc_VersionCmd(int, char *[], char **result)
{
  *result = (char *)Odc_Version;
  return IF_OK;
}

//--------------------------------------------------------------------------
extern "C" int Odc_Init(Tcl_Interp *interp)
{
  // do package specific initializations

  new If_Proc("odc_version", &Odc_VersionCmd);

  return TCL_OK;
}
