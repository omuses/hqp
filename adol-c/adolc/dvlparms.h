/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     dvlparms.h
 Revision: $Id: dvlparms.h,v 1.1 2004/10/13 14:18:11 e_arnold Exp $
 Contents: Developer parameters:
           These parameters are intended for use by the developers and 
           maintainers of ADOL-C to specify library wide definitions.

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          19991122 olvo:  version 1.8.5
          19990816 olvo:  version 1.8.4
          19990310 olvo:  version 1.8.2
          19981130 olvo:  last check

----------------------------------------------------------------------------*/

#if !defined(ADOLC_DVLPARMS_H)
#define ADOLC_DVLPARMS_H 1

/*--------------------------------------------------------------------------*/
/* File names for the tapes */
#define FNAME3     "vs_tape"
#define FNAME2     "rl_tape."
#define FNAME1     "in_tape."
#define FNAME      "op_tape."

/*--------------------------------------------------------------------------*/
/* Standard output used for diagnostics by ADOL-C, e.g. stdout or stderr or */
/* whatever file identifier */
#define DIAG_OUT stdout

/*--------------------------------------------------------------------------*/
/* TAPE IDENTIFICATION (ADOLC & version check) */
#define statSpace   22
#define statSize    11
#define adolcIDSize  5
/* NOTE: adolcIDSize + statSize <= statSpace required! */

/*--------------------------------------------------------------------------*/
/* ADOL-C configuration (never change this) */
#define overwrite 1
#define compsize >

/*--------------------------------------------------------------------------*/
#endif
