/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     common.h
 Revision: $Id: common.h,v 1.1 2004/10/13 14:18:11 e_arnold Exp $
 Contents: Common (global) ADOL-C header  

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.
----------------------------------------------------------------------------*/

#if !defined(ADOLC_COMMON_H)
#define ADOLC_COMMON_H 1

/*--------------------------------------------------------------------------*/
/* system dependend configuration */
#if HAVE_CONFIG_H
#  include "config.h"
#endif

/*--------------------------------------------------------------------------*/
/* developer and user parameters */
#include "dvlparms.h"
#include "usrparms.h"

/*--------------------------------------------------------------------------*/
/* standard includes */
#include <stdlib.h>
#include <stdio.h>

/*--------------------------------------------------------------------------*/
/* further helpful macros */
#if defined(__cplusplus)
#  define BEGIN_C_DECLS extern "C" {
#  define END_C_DECLS   }
#else
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif

#define maxinc(a,b) if ((a) < (b)) (a) = (b)
#define mindec(a,b) if ((a) > (b)) (a) = (b)

#if !defined(max)
#define max(a,b) ( (a)<(b)? (b):(a) )
#endif
#if !defined(min)
#define min(a,b) ( (a)>(b)? (b):(a) )
#endif

/*--------------------------------------------------------------------------*/
#endif
