/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     hov_wk_forward.c
 Revision: $Id: hov_wk_forward.c,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
 Contents: hov_wk_forward (higher-order-vector forward mode with keep)

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.  
----------------------------------------------------------------------------*/
#define _HOV_WK_ 1
#define _KEEP_   1
#include "uni5_for.c"
#undef _KEEP_
#undef _HOV_WK_
