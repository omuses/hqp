/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     drivers/drivers.h
 Revision: $Id: hos_ov_reverse.c,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: hos_ov_reverse (higher-order-scalar reverse mode on vectors)

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20040417 kowarz: adapted to configure - make - make install
          20030304 andrea: first version

----------------------------------------------------------------------------*/

#define _HOS_OV_ 1
#include "ho_rev.c"
#undef _HOS_OV_
