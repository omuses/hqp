/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     taputil_p.h
 Revision: $Id: taputil_p.h,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: Preparation & gets & puts & closing of the taping process
           (ADOL-C internal routines)

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20030306 olvo  extracted from taputil.h of ADOL-C 1.8.7
          
----------------------------------------------------------------------------*/

#if !defined(ADOLC_TAPUTIL_P_H)
#define ADOLC_TAPUTIL_P_H 1

#include "common.h"

BEGIN_C_DECLS

/****************************************************************************/
/*                                                         GLOBAL VARIABLES */

/*--------------------------------------------------------------------------*/
/* Statistic stuff. */
extern int ind_ptr;
extern int dep_ptr;
extern int vs_ptr;
extern int revalso;

/*--------------------------------------------------------------------------*/
/* File Name */
extern char vs_file_name[20];

#ifndef HARDDEBUG
/*--------------------------------------------------------------------------*/
/* Buffers for the operation tape, location tape, real tape. */
extern unsigned char *op_codes;
extern locint        *loc_tape;
extern double        *real_tape;

/*--------------------------------------------------------------------------*/
/* Pointers into the operation tape, location tape, real tape */
extern unsigned char *g_op_ptr;
extern locint        *g_loc_ptr;
extern double        *g_real_ptr;

extern int op_ptr;
extern int loc_ptr;
extern int real_ptr;

/*--------------------------------------------------------------------------*/
/*                                                        MACRO or FUNCTION */
#define get_op_f() *g_op_ptr++ 
#define get_op_r() *(--g_op_ptr)

#define get_locint_f() *g_loc_ptr++
#define get_locint_r() *(--g_loc_ptr)

#define get_val_f() *g_real_ptr++
#define get_val_r() *(--g_real_ptr)

#else /* HARDDEBUG */
extern unsigned char get_op_f(void);
extern unsigned char get_op_r(void);

extern locint get_locint_f(void);
extern locint get_locint_r(void);

extern double get_val_f(void);
extern double get_val_r(void);
#endif 

/****************************************************************************/
/*                                        CONTROL STUFF (inits, ends, etc.) */
extern void init_for_sweep(int);
extern void init_rev_sweep(int);
extern void set_buf_size(int);
extern void set_buffers(char*,unsigned char*,char*, locint*,char*, double *);
extern void close_tape(int*, int);
extern void end_sweep(void);
extern void get_fstr( char*,short,char*);

/****************************************************************************/
/*                                                                     PUTs */

/*--------------------------------------------------------------------------*/
/* Operations */
extern void put_op(unsigned char);

/*--------------------------------------------------------------------------*/
/* Locations */
extern void put_locint(locint);

/*--------------------------------------------------------------------------*/
/* Values */
extern void put_val(double);
extern void put_vals_p(double *,int);
extern void put_vals_r(double *,int);

/*--------------------------------------------------------------------------*/
/* Update/correction of values or locations */
extern void reset_val_r(void);
extern int upd_resloc(locint, locint);
extern int upd_resloc_inc_prod(locint, locint, unsigned char);

/****************************************************************************/
/*                                                                     GETs */

/*--------------------------------------------------------------------------*/
/* Operations */
extern void get_op_block_f(void);
extern void get_op_block_r(void);

/*--------------------------------------------------------------------------*/
/* Locations */
extern void get_loc_block_f(void);
extern void get_loc_block_r(void);

/*--------------------------------------------------------------------------*/
/* Values */
extern int get_val_space(void);
extern void get_val_block_f(void);
extern void get_val_block_r(void);
extern double * get_val_v_f(locint);
extern double * get_val_v_r(locint);

/****************************************************************************/
/*                                                                    UTILs */
extern double make_nan(void);
extern double make_inf(void);

END_C_DECLS

/****************************************************************************/
#endif
