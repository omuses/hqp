/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     tapedoc/tapedoc.h
 Revision: $Id: tapedoc.h,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
 Contents: Contains declaration of tapedoc driver.

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
          19981130 olvo:   newly created (taken from adutilsc.h)

----------------------------------------------------------------------------*/

#if !defined(ADOLC_TAPEDOC_TAPEDOC_H)
#define ADOLC_TAPEDOC_TAPEDOC_H 1

#include "../common.h"

BEGIN_C_DECLS

/****************************************************************************/
/*                                                                 tape_doc */ 
/* tape_doc(tag, m, n, x[n], y[m])                                          */ 

void tape_doc(short, int, int, double*, double*);


/****************************************************************************/
/*                                                               THAT'S ALL */

END_C_DECLS

#endif
