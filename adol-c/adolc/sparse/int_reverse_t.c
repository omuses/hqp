/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     int_reverse_t.c
 Revision: $Id: int_reverse_t.c,v 1.2 2004/10/14 13:29:48 e_arnold Exp $
 Contents: int_reverse_tight
           ( first-order-vector reverse mode for bit patterns,
             checks all dependences on taylors and real values,
             more precize )

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.  
----------------------------------------------------------------------------*/

#define _INT_REV_TIGHT_ 1
#include "int_rev.c"
#undef _INT_REV_TIGHT_
