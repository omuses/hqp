/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     tayutil.h
 Revision: $Id: tayutil.h,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: Management for the value stack tape (Taylors)

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20030306 olvo:   extracted from tayutil.h of ADOL-C 1.8.7
          20030305 andrea: clean up for vs_data
          19991122 olvo:   new op_codes eq_plus_prod eq_min_prod
                           for  y += x1 * x2
                           and  y -= x1 * x2
                           --> new: delete_scaylor(..)  
          19981130 olvo:   automatic cleanup from utils.C moved here
          19980921 olvo:   new interface of void overwrite_scaylor(..) to
                           allow correction of old overwrite in store
          19980708 olvo:   new:  void overwrite_scaylor(..)

----------------------------------------------------------------------------*/

#if !defined(ADOLC_TAYUTIL_H)
#define ADOLC_TAYUTIL_H 1

#include "common.h"

BEGIN_C_DECLS

/****************************************************************************/
/*                                                       DEALLOCATE VS_DATA */
extern void clean_vs_data(int);

/****************************************************************************/
END_C_DECLS

#endif
