/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     tapeutilc.cpp
 Revision: $Id: taputilc.cpp,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: C++ interface for initialization and stopage of the taping 
           process

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
        990713 olvo: trace_on/off: default values for arguments 
        981130 olvo: newly created from utils.C
        
----------------------------------------------------------------------------*/

#include "taputil.h"
#include "adouble.h"

/****************************************************************************/
/*                                                                 TRACE_ON */
/* Trace_on:                                                             
   Initialization for the taping process.  Sets up the arrays op_tape,   
   int_tape, val_tape, and stats.  Op_tape, int_tape, val_tape are arrays
   of pointers to individual buffers for operations, integers (locints), 
   and values (doubles).  Also initializes buffers for this tape, sets   
   files names, and calls appropriate setup routines */
void trace_on( short tnum, int revals )
{ start_trace(tnum,revals);
  take_stock();   /* record all existing adoubles on the tape */
}

/****************************************************************************/
/*                                                                TRACE_OFF */
/* Stop Tracing.  Clean up, and turn off trace_flag */
void trace_off( int flag )
{ int locations;
  locations = keep_stock();     /* copy remaining live variables and turns */
                                /* off trace_flag  */
  stop_trace(locations,flag);   
  std::cout.flush();
}

/****************************************************************************/
/*                                                               THAT'S ALL */
