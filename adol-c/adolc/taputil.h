/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     taputil.h
 Revision: $Id: taputil.h,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: Initialization & stopage of the taping process, as well as
           statistics gathering functions.

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20030304 andrea: new variable for value stack name
          20030306 olvo:   extracted from taputil.h of ADOL-C 1.8.7
          19991122 olvo:   new op_codes eq_plus_prod eq_min_prod
                           for  y += x1 * x2
                           and  y -= x1 * x2
                           --> new: upd_resloc_inc_prod(..)  
          19990713 olvo:   trace_on/off: default values for arguments 
          19981130 olvo:   newly created by unification of taputil?.h
                           and all tape stuff

 History of taputil1.h:
          19980914 olvo:   adolcIDSize 5 (check size of locints ..)
          19980825 olvo:   #defines instead of const (C-Code!)
          19980820 olvo:   Version check
          19980723 olvo:   taputil3.* moved here
          19980713 olvo:   (1) no write_... routines anymore!
                           (2) statistic stuff kept here only
          19980709 olvo:   void write_pos_sign_a(..)
                           void write_neg_sign_a(..)
          19980708 olvo:   void write_upd(..)
          19980707 olvo:   void write_dot_av_av(..)
          19980706 olvo:   void write_incr_decr_a(..)
          19980623 olvo:   new operation code: take_stock_op

 History of taputil2.h:
          19980517 olvo:   griewank's idea:
                           int upd_resloc(locint, locint);
                     
----------------------------------------------------------------------------*/

#if !defined(ADOLC_TAPUTIL_H)
#define ADOLC_TAPUTIL_H 1

#include "common.h"

BEGIN_C_DECLS

/****************************************************************************/
/*                                                       TRACING ON/OFF (C) */
void start_trace(short,int);
void stop_trace(int,int);

/****************************************************************************/
/*                                                          TAPE STATISTICS */
void tapestats(short,int *);

END_C_DECLS

/****************************************************************************/
/****************************************************************************/
/*                                                       Now the C++ THINGS */
#if defined(__cplusplus)

/****************************************************************************/
/*                                                     TRACING ON/OFF (C++) */
void trace_on( short, int = 0 ); 
void trace_off( int = 0 );

#endif

/****************************************************************************/
#endif
