/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     adolc.h
 Revision: $Id: adolc.h,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: Provides all C/C++ interfaces of ADOL-C.
           NOTICE: ALL C/C++ headers will be included DEPENDING ON 
           whether the source code is plain C or C/C++ code. 

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20040423 kowarz: adapted to configure - make - make install
          19981201 olvo:   this new version

----------------------------------------------------------------------------*/

#if !defined(ADOLC_ADOLC_H)
#define ADOLC_ADOLC_H 1

#include "common.h"

/****************************************************************************/
/*                                                  Now the pure C++ THINGS */
#if defined(__cplusplus)
/*--------------------------------------------------------------------------*/
/* Operator overloading things (active doubles & vectors) */
#  include "adouble.h"
#  include "avector.h"
#endif

/****************************************************************************/
/*                                                     Now the C/C++ THINGS */

/*--------------------------------------------------------------------------*/
/* interfaces to basic forward/reverse routines */
#include "interfaces.h"

/*--------------------------------------------------------------------------*/
/* interfaces to "Easy To Use" driver routines for ... */
#include "drivers/drivers.h"    /* optimization & nonlinear equations */
#include "drivers/taylor.h"     /* higher order tensors &
                                         inverse/implicit functions */
#include "drivers/odedrivers.h" /* ordinary differential equations */

/*--------------------------------------------------------------------------*/
/* interfaces to TAPEDOC package */
#include "tapedoc/tapedoc.h"

/*--------------------------------------------------------------------------*/
/* interfaces to SPARSE package */
#include "sparse/sparse.h"
#include "sparse/jacutils.h"

/*--------------------------------------------------------------------------*/
/* tape and value stack utilities */
#include "taputil.h"
#include "tayutil.h"

/*--------------------------------------------------------------------------*/
/* allocation utilities */
#include "adalloc.h"

/****************************************************************************/
#endif
