#ifndef _TAPUTIL_H_
#define _TAPUTIL_H_
/*
   --------------------------------------------------------------
   File taputil.h of ADOL-C version 1.8.5         as of Nov/22/99
   --------------------------------------------------------------
   Initialization, stopage, and gets & puts of the taping 
   process, as well as statistics gathering functions.

   Last changes:
        991122 olvo  new op_codes eq_plus_prod eq_min_prod
                     for  y += x1 * x2
                     and  y -= x1 * x2
                     --> new: upd_resloc_inc_prod(..)  
        990713 olvo: trace_on/off: default values for arguments 
        981130 olvo: newly created by unification of taputil?.h
                     and all tape stuff

   History of taputil1.h:
        980914 olvo: adolcIDSize 5 (check size of locints ..)
        980825 olvo: #defines instead of const (C-Code!)
        980820 olvo: Version check
        980723 olvo: taputil3.* moved here
        980713 olvo: (1) no write_... routines anymore!
                     (2) statistic stuff kept here only
        980709 olvo: void write_pos_sign_a(..)
                     void write_neg_sign_a(..)
        980708 olvo: void write_upd(..)
        980707 olvo: void write_dot_av_av(..)
        980706 olvo: void write_incr_decr_a(..)
        980623 olvo: new operation code: take_stock_op
   History of taputil2.h:
        980517 olvo: griewank's idea:
                     int upd_resloc(locint, locint);

   --------------------------------------------------------------
*/


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                           PUBLIC EXPORTS */


#ifdef __cplusplus
/****************************************************************************/
/****************************************************************************/
/*                                                       Now the C++ THINGS */

/****************************************************************************/
/*                                                     TRACING ON/OFF (C++) */
void trace_on( short, int = 0 ); 
void trace_off( int = 0 );


/****************************************************************************/
/****************************************************************************/
/*                                                         Now the C THINGS */
extern "C" {
#endif

/****************************************************************************/
/*                                                          TAPE STATISTICS */
extern void tapestats(short,int *);


/****************************************************************************/
/*                                                       TRACING ON/OFF (C) */
extern void start_trace(short,int);
extern void stop_trace(int,int);


#ifdef __cplusplus
}
#endif


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                  ADOL-C INTERNAL EXPORTS */
#ifdef _ADOLC_SRC_


#ifdef __cplusplus
/****************************************************************************/
/****************************************************************************/
/*                                                            No C++ THINGS */


/****************************************************************************/
/****************************************************************************/
/*                                                         Now the C THINGS */
extern "C" {
#endif


/****************************************************************************/
/*                                                         GLOBAL VARIABLES */

/*--------------------------------------------------------------------------*/
/* Statistic stuff. */
extern int ind_ptr;
extern int dep_ptr;
extern int vs_ptr;
extern int revalso;


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


/****************************************************************************/
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

#endif /* _ADOLC_SRC_ */

#endif

