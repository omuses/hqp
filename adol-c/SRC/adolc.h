#ifndef _ADOLC_H_
#define _ADOLC_H_
/*
   ------------------------------------------------------------- 
   File adutils.h of ADOL-C version 1.8.0        as of Dec/01/98
   -------------------------------------------------------------
   Provides all C/C++ interfaces of ADOL-C.
   NOTICE: ALL C/C++ headers will be included DEPENDING ON 
           whether the source code is plain C or C/C++ code.
 
   Last changes: 
     981201 olvo     this new version

  ------------------------------------------------------------- 
*/

/****************************************************************************/
/*                                                            JUST INCLUDES */


#ifdef __cplusplus
/****************************************************************************/
/****************************************************************************/
/*                                                  Now the pure C++ THINGS */

/*--------------------------------------------------------------------------*/
/* Operator overloading things (active doubles & vectors) */
#include "adouble.h"
#include "avector.h"


#endif


/****************************************************************************/
/****************************************************************************/
/*                                                     Now the C/C++ THINGS */

/*--------------------------------------------------------------------------*/
/* interfaces to basic forward/reverse routines */
#include "interfaces.h"

/*--------------------------------------------------------------------------*/
/* interfaces to "Easy To Use" driver routines for ... */
#include "DRIVERS/drivers.h"    /* ... optimization & nonlinear equations */
#include "DRIVERS/taylor.h"     /* ... higher order tensors & inverse/implicit 
                                       functions */
#include "DRIVERS/odedrivers.h" /* ... ordinary differential equations */

/*--------------------------------------------------------------------------*/
/* interfaces to TAPEDOC package */
#include "TAPEDOC/tapedoc.h"

/*--------------------------------------------------------------------------*/
/* interfaces to SPARSE package */
#include "SPARSE/sparse.h"
#include "SPARSE/jacutils.h"

/*--------------------------------------------------------------------------*/
/* tape utilities */
#include "taputil.h"

/*--------------------------------------------------------------------------*/
/* allocation utilities */
#include "adalloc.h"


/****************************************************************************/
/*                                                               THAT'S ALL */
#endif



