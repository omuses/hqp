#define _TAPUTILC_CPP_
#define _ADOLC_SRC_
/*
   --------------------------------------------------------------
   File taputilC.C of ADOL-C version 1.8.3        as of Jul/13/99
   --------------------------------------------------------------
   C++ interface for initialization and stopage of the taping 
   process

   Last changes:
        990713 olvo: trace_on/off: default values for arguments 
        981130 olvo: newly created from utils.C
        
   --------------------------------------------------------------
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "dvlparms.h"
#include "usrparms.h"
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
  cout.flush();
}


/****************************************************************************/
/*                                                               THAT'S ALL */
#undef _ADOLC_SRC_
#undef _TAPUTILC_CPP_
