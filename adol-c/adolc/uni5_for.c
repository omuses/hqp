/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     uni5_for.mc
 Revision: $Id: uni5_for.c,v 1.2 2004/10/14 13:29:48 e_arnold Exp $
 Contents: Contains the routines :
           zos_forward (zero-order-scalar forward mode):      define _ZOS_   
           fos_forward (first-order-scalar forward mode):     define _FOS_
           hos_forward (higher-order-scalar forward mode):    define _HOS_
           fov_forward (first-order-vector forward mode):     define _FOV_
           hov_forward (higher-order-vector forward mode):    define _HOV_
           hov_wk_forward (higher-order-vector forward mode): define _HOV_WK_

           Uses the preprocessor to compile the 6 different object files
           with/without "keep" parameter:                     define _KEEP_

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
          20030317 andrea: change pow_op
          20030305 andrea: change taylor_begin
          20010719 andrea: start of hov_wk_forward
          20000310 olvo:   (1) better error messages for discontinuity
                           (2) removed trigraphs stuff
          19991122 olvo:   new op_codes eq_plus_prod eq_min_prod
                           for  y += x1 * x2
                           and  y -= x1 * x2  
          19990816 olvo:   computation of size of value stack buffer
                           (=taylor buffer) changed in order to avoid
                           overflows allowing larger values of vs_size
                           (unsigned longs or even doubles for
                           counters considered for future developments)
          19980929 olvo:   (1) new locint checkSize
                           (2) changed strategy to allow reflexive
                               - eq_mult_av_a, mult_a_av, div_av_a
          19980924 olvo:   (1) allow reflexive operations for
                               - vector operations: eq_mult_av_a, 
                               dot_av_av, mult_a_av, div_av_a
                               - cond_assign, cond_assign_s
                           (2) new macros TQO*
                           (3) deleted all int_* opcodes
          19980923 olvo:   allow reflexive operations for
                               - min_op, abs_op
          19980922 olvo:   allow reflexive operations for
                               - gen_quad
                               - div_d_a, div_a_a
          19980921 olvo:   (1) changed save-order in sin/cos
                           (2) allow reflexive operations for
                               - sin/cos
                               - asin/acos/atan
                               - asinh/acosh/atanh
                               - erf, log, pow                          
          19980915 olvo.   (1) modified mult_a_a to allow x = x*x
                           (2) correction of c&p-error in macros
                               (TARG_*, TARG1_*, TARG2_*) 
                           (3) modified exp_op to allow x = exp(x)
          19980820 olvo:   new comparison strategy 
          19980818 olvo:   ec BREAK_ZOS 
          19980721 olvo:   (1) ec in min_op, abs_op, ceil_op, floor_op
                           (2) write of taylors in subscript and
                               subscript_m 
          19980714 olvo:   some error elimination
                           op-code mult_av_a removed  
          19980710 olvo:   sin/cos writes 2 taylors
          19980709 olvo:   new operation code: neg_sign_a
                                               pos_sign_a
          19980706 olvo:   new operation code: int_adb_d_one
                                               int_adb_d_zero
          19980703 olvo:   new operation code: assign_d_one
                                               assign_d_zero
          19980625 olvo:   vector operations & subscripts
          19980624 olvo:   rearrangements
          19980623 mitev/olvo: revision + stock stuff 
          19980615 olvo:   griewank's idea
          19980609 mitev  

----------------------------------------------------------------------------*/

#include "interfaces.h"
#include "adalloc.h"
#include "taputil.h"
#include "taputil_p.h"
#include "tayutil_p.h"
#include "oplate.h"

#include <stdio.h>
#include <math.h>
#include <malloc.h>

/****************************************************************************/
/*                                                                   MACROS */
#undef _ADOLC_VECTOR_
#undef _HIGHER_ORDER_

/*--------------------------------------------------------------------------*/
#if defined(_ZOS_)
#  define GENERATED_FILENAME "zos_forward" 

/*--------------------------------------------------------------------------*/
#else
#if defined(_FOS_)
#define GENERATED_FILENAME "fos_forward"        

#define ARGUMENT(indexi,l,i) argument[indexi]
#define TAYLORS(indexd,l,i)   taylors[indexd]

/*--------------------------------------------------------------------------*/
#else
#if defined(_FOV_)
#define GENERATED_FILENAME "fov_forward"         

#define _ADOLC_VECTOR_

#define ARGUMENT(indexi,l,i) argument[indexi][l]
#define TAYLORS(indexd,l,i)   taylors[indexd][l]

/*--------------------------------------------------------------------------*/
#else
#if defined(_HOS_)
#define GENERATED_FILENAME "hos_forward"         

#define _HIGHER_ORDER_

#define ARGUMENT(indexi,l,i) argument[indexi][i]
#define TAYLORS(indexd,l,i)   taylors[indexd][i]

/*--------------------------------------------------------------------------*/
#else
#if defined(_HOV_)
#define GENERATED_FILENAME "hov_forward"         

#define _ADOLC_VECTOR_
#define _HIGHER_ORDER_

#define ARGUMENT(indexi,l,i) argument[indexi][l][i]
#define TAYLORS(indexd,l,i)   taylors[indexd][l][i]

/*--------------------------------------------------------------------------*/
#else
#if defined(_HOV_WK_)
#define GENERATED_FILENAME "hov_wk_forward"         

#define _ADOLC_VECTOR_
#define _HIGHER_ORDER_

#define ARGUMENT(indexi,l,i) argument[indexi][l][i]
#define TAYLORS(indexd,l,i)   taylors[indexd][l][i]

#else
#error Error ! Define [_ZOS_ | _FOS_ | _HOS_ | _FOV_ | _HOV_ | _HOV_WK_ ] [{_KEEP_}] 
#endif
#endif
#endif
#endif
#endif
#endif

/*--------------------------------------------------------------------------*/
/*                                                               KEEP stuff */
#if defined(_KEEP_)

#if defined(_HOV_WK_) /* keep in this vector mode */
#define IF_KEEP_TAYLOR_CLOSE \
if (keep){\
  fprintf(DIAG_OUT,"Succeeding reverse sweep will fail!\n");\
  taylor_close(-1,depcheck,indcheck);\
} 
#define IF_KEEP_WRITE_TAYLOR(res,keep,k,p) if (keep){ write_scaylor(T0[res]); \
if (keep > 1) write_taylors(res,(keep-1),k,p);} 
#else
#if defined(_ADOLC_VECTOR_) /* otherwise no keep */
#define IF_KEEP_TAYLOR_CLOSE  
#define IF_KEEP_WRITE_TAYLOR(res,keep,k,p) 
#else /* _ZOS_, _FOS_, _HOS_ */
#define IF_KEEP_TAYLOR_CLOSE \
if (keep){\
  fprintf(DIAG_OUT,"Otherwise succeeding reverse sweep will fail!\n");\
  taylor_close(-1,depcheck,indcheck);\
}
#if defined(_ZOS_)
#define IF_KEEP_WRITE_TAYLOR(res,keep,k,p) if (keep) write_scaylor(T0[res]);
#else
#if defined(_FOS_)
#define IF_KEEP_WRITE_TAYLOR(res,keep,k,p)  if (keep){ write_scaylor(T0[res]); \
if (keep > 1) write_scaylor(T[res]);}
#else
#if defined(_HOS_) 
#define IF_KEEP_WRITE_TAYLOR(res,keep,k,p) if (keep){ write_scaylor(T0[res]); \
if (keep > 1) write_taylor(res,keep-1);}
#endif
#endif
#endif
#endif
#endif

#else  /* no _KEEP_ */
#define IF_KEEP_TAYLOR_CLOSE  
#define IF_KEEP_WRITE_TAYLOR(res,keep,k,p) 
#endif 

/*--------------------------------------------------------------------------*/
/*                                                      access to variables */
#if !defined(_ZOS_)
#if defined(_FOS_)
#define TRES         *Tres
#define TARG         *Targ
#define TARG1        *Targ1
#define TARG2        *Targ2
#define TQO          *Tqo

#define TRES_INC     *Tres
#define TARG_INC     *Targ
#define TARG1_INC    *Targ1
#define TARG2_INC    *Targ2
#define TQO_INC      *Tqo

#define TRES_DEC     *Tres
#define TARG_DEC     *Targ
#define TARG1_DEC    *Targ1
#define TARG2_DEC    *Targ2
#define TQO_DEC      *Tqo

#define TRES_FOINC   *Tres
#define TARG_FOINC   *Targ
#define TARG1_FOINC  *Targ1
#define TARG2_FOINC  *Targ2
#define TQO_FOINC    *Tqo

#define TRES_FODEC   *Tres
#define DEC_TRES_FO
#define TARG_FODEC   *Targ
#define TARG1_FODEC  *Targ1
#define TARG2_FODEC  *Targ2
#define TQO_FODEC    *Tqo

#define ASSIGN_T(a,b)  a = &b;
 
#else  /* _HOS_, _FOV_, _HOV_, _HOV_WK */
#define TRES         *Tres
#define TARG         *Targ
#define TARG1        *Targ1
#define TARG2        *Targ2
#define TQO          *Tqo

#define TRES_INC     *Tres++
#define TARG_INC     *Targ++
#define TARG1_INC    *Targ1++
#define TARG2_INC    *Targ2++
#define TQO_INC      *Tqo++

#define TRES_DEC     *Tres--
#define TARG_DEC     *Targ--
#define TARG1_DEC    *Targ1--
#define TARG2_DEC    *Targ2--
#define TQO_DEC      *Tqo--

#if defined(_FOV_)
#define TRES_FOINC   *Tres++
#define TARG_FOINC   *Targ++
#define TARG1_FOINC  *Targ1++
#define TARG2_FOINC  *Targ2++
#define TQO_FOINC    *Tqo++

#define TRES_FODEC   *Tres
#define DEC_TRES_FO  Tres--;
#define TARG_FODEC   *Targ--
#define TARG1_FODEC  *Targ1--
#define TARG2_FODEC  *Targ2--
#define TQO_FODEC    *Tqo--
#else /* _HOS_, _HOV_, _HOV_WK */
#define TRES_FOINC   *Tres
#define TARG_FOINC   *Targ
#define TARG1_FOINC  *Targ1
#define TARG2_FOINC  *Targ2
#define TQO_FOINC    *Tqo

#define TRES_FODEC   *Tres
#define DEC_TRES_FO
#define TARG_FODEC   *Targ
#define TARG1_FODEC  *Targ1
#define TARG2_FODEC  *Targ2
#define TQO_FODEC    *Tqo
#endif

#define ASSIGN_T(a,b)  a = b;
#endif
#endif


/*--------------------------------------------------------------------------*/
/*                                                               loop stuff */
#if defined(_ADOLC_VECTOR_)
#define FOR_0_LE_l_LT_p for (l=0; l<p; l++)  
#define FOR_p_GT_l_GE_0 for (l=p-1; l>=0; l--)  
#else
#define FOR_0_LE_l_LT_p 
#define FOR_p_GT_l_GE_0  
#endif

#if defined(_HIGHER_ORDER_)
#define FOR_0_LE_i_LT_k for (i=0; i<k; i++) 
#define FOR_k_GT_i_GE_0 for (i=k-1; i>=0; i--) 
#else
#define FOR_0_LE_i_LT_k  
#define FOR_k_GT_i_GE_0  
#endif

#if defined(_HOV_)
#define FOR_0_LE_l_LT_pk for (l=0; l<pk; l++)  
#define INC_pk_1(T)      T += pk-1;
#define VEC_INC(T,inc)   T += inc; 
#define HOV_INC(T,inc)   T += inc; 
#else
#if defined(_HOV_WK_)
#define FOR_0_LE_l_LT_pk for (l=0; l<pk; l++)  
#define INC_pk_1(T)      T += pk-1;
#define VEC_INC(T,inc)   T += inc; 
#define HOV_INC(T,inc)   T += inc; 
#else
#if defined(_FOV_)
#define FOR_0_LE_l_LT_pk for (l=0; l<p; l++)  
#define INC_pk_1(T)      T += p-1;
#define VEC_INC(T,inc)   T++; 
#define HOV_INC(T,inc) 
#else
#if defined(_HOS_)
#define FOR_0_LE_l_LT_pk for (l=0; l<k; l++)  
#define INC_pk_1(T)      T += k-1;
#define VEC_INC(T,inc) 
#define HOV_INC(T,inc) 
#else
#define FOR_0_LE_l_LT_pk
#define INC_pk_1(T) 
#define VEC_INC(T,inc)
#define HOV_INC(T,inc) 
#endif
#endif
#endif
#endif

/*--------------------------------------------------------------------------*/
/*                                                        higher order case */
#if defined(_HIGHER_ORDER_)
#define BREAK_FOR_I break; 
#else
#define BREAK_FOR_I ;  
#endif

/*--------------------------------------------------------------------------*/
/*                                                             other macros */
#define FMIN(x,y)  ((y<x)?y:x)
#define MINDEC(a,b) if ((a) > (b)) (a) = (b)

/* END Macros */

/****************************************************************************/
/*                                                             NOW THE CODE */
BEGIN_C_DECLS

/*--------------------------------------------------------------------------*/
/*                                                   Local Static Variables */
static short tag;

static int for_location_cnt;
static int dep_cnt;
static int ind_cnt;

#if !defined(_ADOLC_VECTOR_)
static int valstack;
#endif
#if defined(_HOV_WK_)
static int valstack;
#endif

#if defined(_ZOS_)
/****************************************************************************/
/* Zero Order Scalar version of the forward mode.                           */
/****************************************************************************/
#if defined(_KEEP_)
int  zos_forward(
#else
int  zos_forward_nk(
#endif
		 short  tnum,        /* tape id */
		 int    depcheck,    /* consistency chk on # of deps */
		 int    indcheck,    /* consistency chk on # of indeps */
#if defined(_KEEP_)
		 int    keep,        /* flag for reverse sweep */
#endif
		 double *basepoint,  /* independant variable values */
		 double *valuepoint) /* dependent variable values */

#else
#if defined(_FOS_)
/****************************************************************************/
/* First Order Scalar version of the forward mode.                          */
/****************************************************************************/
#if defined(_KEEP_)
int  fos_forward(
#else
int  fos_forward_nk(
#endif
		 short  tnum,        /* tape id */
		 int    depcheck,    /* consistency chk on # of deps */
		 int    indcheck,    /* consistency chk on # of indeps */
#if defined(_KEEP_)
		 int    keep,        /* flag for reverse sweep */
#endif
                 double *basepoint,  /* independent variable values */
		 double *argument,   /* Taylor coefficients (input) */
		 double *valuepoint, /* Taylor coefficients (output) */
                 double *taylors)    /* matrix of coifficient vectors */
/* the order of the indices in argument and taylors is [var][taylor] */

#else
#if defined(_FOV_)
/****************************************************************************/
/* First Order Vector version of the forward mode.                          */
/****************************************************************************/
int  fov_forward(
		 short  tnum,        /* tape id */
		 int    depcheck,    /* consistency chk on # of deps */
		 int    indcheck,    /* consistency chk on # of indeps */
		 int    p,           /* # of taylor series */
                 double *basepoint,  /* independent variable values */
		 double **argument,  /* Taylor coefficients (input) */
                 double *valuepoint, /* Taylor coefficients (output) */
		 double **taylors)   /* matrix of coifficient vectors */
/* the order of the indices in argument and taylors is [var][taylor] */

#else
#if defined(_HOS_)
/****************************************************************************/
/* Higher Order Scalar version of the forward mode.                         */
/****************************************************************************/
#if defined(_KEEP_)
int  hos_forward(
#else
int  hos_forward_nk(
#endif
		 short  tnum,        /* tape id */
		 int    depcheck,    /* consistency chk on # of dependents */
		 int    indcheck,    /* consistency chk on # of independents */
		 int    gdegree,     /* highest derivative degree */
#if defined(_KEEP_)
		 int    keep,        /* flag for reverse sweep */
#endif
                 double *basepoint,  /* independent variable values */
		 double **argument,  /* independant variable values */
                 double *valuepoint, /* Taylor coefficients (output) */
		 double **taylors)   /* matrix of coifficient vectors */


#else
/****************************************************************************/
/* Higher Order Vector version of the forward mode.                         */
/****************************************************************************/
#if defined(_KEEP_)
int  hov_wk_forward(
#else
int  hov_forward(
#endif
		 short  tnum,        /* tape id */
		 int    depcheck,    /* consistency chk on # of deps */
		 int    indcheck,    /* consistency chk on # of indeps */
		 int    gdegree,     /* highest derivative degree */
#if defined(_KEEP_)
		 int    keep,        /* flag for reverse sweep */
#endif
		 int    p,           /* # of taylor series */
                 double *basepoint,  /* independent variable values */
		 double ***argument, /* Taylor coefficients (input) */
                 double *valuepoint, /* Taylor coefficients (output) */
		 double ***taylors)  /* matrix of coifficient vectors */
/* the order of the indices in argument and taylors is [var][taylor][deriv] */

#endif
#endif
#endif
#endif

{
/****************************************************************************/
/*                                                            ALL VARIABLES */
  unsigned char operation; /* operation code */
  int tape_stats[11];      /* tape stats */
  int ret_c =3;            /* return value */

  locint size = 0;
  locint res  = 0;
  locint arg  = 0;
  locint arg1 = 0;
  locint arg2 = 0;
  locint checkSize;

  double coval = 0, *d = 0;

  int indexi = 0,  indexd = 0;

  /* loop indices */
#if !defined (_ZOS_)
  int i;
#endif
#if defined (_HIGHER_ORDER_)
  int j;
#endif
  int l=0, ls;

  /* other necessary variables */
#if !defined (_ZOS_)
  double r0=0.0, x, y, divs;
  int even, flag;
#endif
  int buffer;
#if defined (_HIGHER_ORDER_)
  static int kax;
#endif
#if defined (_ADOLC_VECTOR_)
  static int pax;
#endif
  static int fax;
  
  /* Taylor stuff */
  static double  *T0;
#if !defined(_ZOS_)
#if  defined(_FOS_)
  static double  *T;
  double  Ttemp;
#else
  static double  **T;
  static double  *Ttemp;
#endif
  double         *Tres, *Targ, *Targ1, *Targ2, *Tqo;

#if defined (_HIGHER_ORDER_)
  double         *TresOP2, *zOP;
#endif
  double         *TresOP, *TargOP, *Targ1OP, *Targ2OP;
#endif
  double         T0temp;
#define T0res  T0temp
#define T0arg  T0temp

#if defined(_HIGHER_ORDER_)
  static double *z;
  int k = gdegree;
#endif
  
#if !defined(_ADOLC_VECTOR_)
  int numoperations;
#endif

#if defined(_KEEP_)  
  int taylbuf=0;
#endif

#if defined(_HOV_)
  int pk = k*p;
#else
#if defined(_HOV_WK_)
  int pk = k*p;
  int numoperations;
#endif
#endif

#if defined(ADOLC_DEBUG)
/****************************************************************************/
/*                                                           DEBUG MESSAGES */
  fprintf(DIAG_OUT,"Call of %s(..) with tag: %d, n: %d, m %d,\n",
                   GENERATED_FILENAME, tnum, indcheck, depcheck);      
#if defined(_KEEP_)
  fprintf(DIAG_OUT,"                    keep: %d\n", keep);
#endif
#if defined(_HIGHER_ORDER_)
  fprintf(DIAG_OUT,"                    degree: %d\n",gdegree);
#endif
#if defined(_ADOLC_VECTOR_)
  fprintf(DIAG_OUT,"                    p: %d\n\n",p);
#endif
  
#endif

/****************************************************************************/
/*                                                                    INITs */

/* Set up stuff for the tape */
  tag = tnum;         /*tag is global which specifies which tape to look at */

  tapestats(tag,tape_stats);
  ind_cnt          = tape_stats[0];
  dep_cnt          = tape_stats[1];
  for_location_cnt = tape_stats[2];
  buffer           = tape_stats[4];
#if !defined(_ADOLC_VECTOR_)
  /* olvo 980615 deleted nl 
  valstack = tape_stats[3]; 
                 new n2l */
  valstack = tape_stats[3]; 
  numoperations    = tape_stats[5];
  /*valstack = numoperations */ /* Hopefully an upper bound - says griewank */ 
  /*        + for_location_cnt;*/ /* olvo 980618 ec think about that!!! */
#endif
#if defined(_HOV_WK_)
  valstack = tape_stats[3]; 
  numoperations    = tape_stats[5];
#endif
  set_buf_size(buffer);

  if ((depcheck != dep_cnt)||(indcheck != ind_cnt))
  { fprintf(DIAG_OUT,"ADOL-C error: forward sweep on tape %d  aborted!\n",tag);
    fprintf(DIAG_OUT,"Number of dependent and/or independent variables passed"
                   " to forward is\ninconsistant with number"
                   " recorded on tape %d \n",tag);
    exit (-1);
  }


/****************************************************************************/
/*                                                        MEMORY ALLOCATION */
/* olvo 980626 has to be revised for common blocks */ 

/*--------------------------------------------------------------------------*/
#if defined(_ZOS_)                                                   /* ZOS */
#if defined(_KEEP_) 
  if (keep>1){
    fprintf(DIAG_OUT,"\n ADOL-C error: zero order scalar forward cannot save"
                   " more\nthan zero order taylor coefficients!\n");
    exit (-1);
  }
#endif
  if (for_location_cnt compsize fax)
  { if (fax)
      free((char *) T0);
    T0 = myalloc1(for_location_cnt);
    fax = for_location_cnt;
  }
#if defined(_KEEP_) 
  /* olvo 990816 new: */
  if (keep)
  { taylbuf = TBUFSIZE/(keep*sizeof(revreal));
    taylbuf *= keep*sizeof(revreal); 
    taylor_begin(tnum,taylbuf,&T0,keep-1);
  }
#endif

/*--------------------------------------------------------------------------*/
#else                                                                /* FOS */
#if defined(_FOS_)
#if defined(_KEEP_)
  if (keep>2){
    fprintf(DIAG_OUT,"\n ADOL-C error: first order scalar forward cannot save"
                   " more  \nthan first order taylor coefficients!\n");
    exit (-1);
  }
#endif
  if (for_location_cnt compsize fax)
  { if (fax)
    { free((char*) T0);
      free((char*) T); 
    }
    T0 = myalloc1(for_location_cnt);
    T  = myalloc1(for_location_cnt);
    fax = for_location_cnt;
  }
#if defined(_KEEP_) 
  /* olvo 990816 new: */
  if (keep)
  { taylbuf = TBUFSIZE/(keep*sizeof(revreal));
    taylbuf *= keep*sizeof(revreal); 
    taylor_begin(tnum,taylbuf,&T,keep-1);
  }
#endif

/*--------------------------------------------------------------------------*/
#else                                                                /* FOV */
#if defined(_FOV_)
  if (p compsize pax || for_location_cnt compsize fax)
  { if (pax || fax) 
    { free((char*) T0);
      free((char*) *T); free((char*) T);
      free((char*) Ttemp);
    }
    T0    = myalloc1(for_location_cnt);
    T     = myalloc2(for_location_cnt,p);
    Ttemp = myalloc1(p);
    pax = p;
    fax = for_location_cnt;
  }

/*--------------------------------------------------------------------------*/
#else                                                                /* HOS */ 
#if defined(_HOS_)
  if (k compsize kax || for_location_cnt compsize fax)
  { if (kax || fax)
    { free((char*) T0); 
      free((char*) *T); free((char*) T);
      free((char*) z);
      free((char*) Ttemp);
    }
    T0 = myalloc1(for_location_cnt);
    T  = myalloc2(for_location_cnt,k);
    z  = myalloc1(k);
    Ttemp = myalloc1(k);
    kax = k;
    fax = for_location_cnt;
  }
#if defined(_KEEP_)
  /* olvo 990816 new: */
  if (keep)
  { taylbuf = TBUFSIZE/(keep*sizeof(revreal));
    taylbuf *= keep*sizeof(revreal); 
    taylor_begin(tnum,taylbuf,T,keep-1);
  }
#endif

/*--------------------------------------------------------------------------*/
#else                                                     /* HOV and HOV_WK */
  if (k compsize kax || for_location_cnt compsize fax || p compsize pax)
  { if (kax || pax || fax)
    { free((char*) T0); 
      free((char*) *T); free((char*) T);
      free((char*) z);
      free((char*) Ttemp);
    } 
    T0    = myalloc1(for_location_cnt);
    T     = myalloc2(for_location_cnt,p*k);
    z     = myalloc1(k);
    Ttemp = myalloc1(p*k);
    kax = k;
    fax = for_location_cnt;
    pax = p ;
  }

#if defined(_KEEP_)
  if (keep)
  { taylbuf = TBUFSIZE/(keep*sizeof(revreal));
    taylbuf *= p*keep*sizeof(revreal); 
    taylor_begin(tnum,taylbuf,T,keep-1);
  }
#endif

#endif
#endif
#endif
#endif

/****************************************************************************/
/*                                                            FORWARD SWEEP */

  /* Initialize the Forward Sweep */
  init_for_sweep(tag);

  operation=get_op_f();
  while (operation !=end_of_tape)
  { switch (operation){

      
/****************************************************************************/
/*                                                                  MARKERS */

/*--------------------------------------------------------------------------*/
      case end_of_op:                                          /* end_of_op */
        get_op_block_f();
        operation=get_op_f(); 
        /* Skip next operation, it's another end_of_op */
        break;

/*--------------------------------------------------------------------------*/
      case end_of_int:                                        /* end_of_int */
        get_loc_block_f();
        break;

/*--------------------------------------------------------------------------*/
      case end_of_val:                                        /* end_of_val */
        get_val_block_f();
        break;
/*--------------------------------------------------------------------------*/
      case start_of_tape:                                  /* start_of_tape */
      case end_of_tape:                                      /* end_of_tape */
	break;


/****************************************************************************/
/*                                                               COMPARISON */

/*--------------------------------------------------------------------------*/
      case eq_zero:                                              /* eq_zero */
        arg = get_locint_f();

        if (T0[arg] != 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator eq_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          IF_KEEP_TAYLOR_CLOSE
          end_sweep();
          return (-1);
        }
        ret_c = 0;
        break;

/*--------------------------------------------------------------------------*/
      case neq_zero:                                            /* neq_zero */
        arg = get_locint_f();

        if (T0[arg] == 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator neq_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          IF_KEEP_TAYLOR_CLOSE
          end_sweep();
          return (-1);
        }
        break;

/*--------------------------------------------------------------------------*/
      case le_zero:                                              /* le_zero */
        arg = get_locint_f();

        if (T0[arg] > 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator le_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          IF_KEEP_TAYLOR_CLOSE
          end_sweep();
          return (-1);
        }
        if (T0[arg] == 0)
          ret_c = 0;
        break;

/*--------------------------------------------------------------------------*/
      case gt_zero:                                              /* gt_zero */
        arg = get_locint_f();

        if (T0[arg] <= 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator gt_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          IF_KEEP_TAYLOR_CLOSE
          end_sweep();
          return (-1);
        }
        break;

/*--------------------------------------------------------------------------*/
      case ge_zero:                                              /* ge_zero */
        arg = get_locint_f();

        if (T0[arg] < 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator ge_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          IF_KEEP_TAYLOR_CLOSE
          end_sweep();
          return (-1);
        }
        if (T0[arg] == 0)
          ret_c = 0;
        break;

/*--------------------------------------------------------------------------*/
      case lt_zero:                                              /* lt_zero */
        arg = get_locint_f();

        if (T0[arg] >= 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator lt_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          IF_KEEP_TAYLOR_CLOSE
          end_sweep();
          return (-1);
        }
        break;


/****************************************************************************/
/*                                                              ASSIGNMENTS */

/*--------------------------------------------------------------------------*/
      case assign_a:           /* assign an adouble variable an    assign_a */
	                       /* adouble value. (=) */
        arg = get_locint_f();
        res = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Targ,T[arg])
        ASSIGN_T(Tres,T[res])        

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case assign_d:            /* assign an adouble variable a    assign_d */
	                        /* double value. (=) */
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	T0[res] = coval;
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres, T[res])

        FOR_0_LE_l_LT_pk
	  TRES_INC = 0;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case assign_d_zero:  /* assign an adouble variable a    assign_d_zero */
	                   /* double value. (0) (=) */
        res   = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	T0[res] = 0.0;
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres, T[res])

        FOR_0_LE_l_LT_pk
	  TRES_INC = 0;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case assign_d_one:    /* assign an adouble variable a    assign_d_one */
	                    /* double value. (1) (=) */
        res   = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	T0[res] = 1.0;
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres, T[res])

        FOR_0_LE_l_LT_pk
	  TRES_INC = 0;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case assign_ind:       /* assign an adouble variable an    assign_ind */
	                     /* independent double value (<<=) */
        res = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	T0[res] = basepoint[indexi];
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres, T[res])

        FOR_0_LE_l_LT_p
          FOR_0_LE_i_LT_k
  	    TRES_INC = ARGUMENT(indexi,l,i);

#endif /* ALL_TOGETHER_AGAIN */
	++indexi;
	break;

/*--------------------------------------------------------------------------*/
      case assign_dep:           /* assign a float variable a    assign_dep */
	                         /* dependent adouble value. (>>=) */
	res = get_locint_f();

        valuepoint[indexd] = T0[res];
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres, T[res])

        if (taylors != 0 )  /* ??? question: why here? */
          FOR_0_LE_l_LT_p
            FOR_0_LE_i_LT_k
	      TAYLORS(indexd,l,i) = TRES_INC;

#endif /* ALL_TOGETHER_AGAIN */
	indexd++;
	break;


/****************************************************************************/
/*                                                   OPERATION + ASSIGNMENT */

/*--------------------------------------------------------------------------*/
      case eq_plus_d:            /* Add a floating point to an    eq_plus_d */
	                         /* adouble. (+=) */
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] += coval;
	break;

/*--------------------------------------------------------------------------*/
      case eq_plus_a:             /* Add an adouble to another    eq_plus_a */
	                          /* adouble. (+=) */
        arg = get_locint_f();
        res = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] += T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC += TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case eq_min_d:       /* Subtract a floating point from an    eq_min_d */
                           /* adouble. (-=) */
        res = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	T0[res] -= coval;
	break;

/*--------------------------------------------------------------------------*/
      case eq_min_a:        /* Subtract an adouble from another    eq_min_a */
	                    /* adouble. (-=) */
        arg = get_locint_f();
        res = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	T0[res] -= T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])
        
        FOR_0_LE_l_LT_pk
	  TRES_INC -= TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_d:              /* Multiply an adouble by a    eq_mult_d */
	                           /* flaoting point. (*=) */
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	T0[res] *= coval;
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres, T[res])

        FOR_0_LE_l_LT_pk
	    TRES_INC *= coval;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_a:       /* Multiply one adouble by another    eq_mult_a */
	                    /* (*=) */
        arg = get_locint_f();
        res = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        INC_pk_1(Tres)
        INC_pk_1(Targ)

        FOR_p_GT_l_GE_0
	  FOR_k_GT_i_GE_0
     { TRES_FODEC = T0[res]*TARG_DEC + TRES*T0[arg];
       DEC_TRES_FO
#ifdef _HIGHER_ORDER_
            TresOP = Tres-i;
            TargOP = Targ;

            for (j=0;j<i;j++)
              *Tres += (*TresOP++) * (*TargOP--);
            Tres--;
#endif /* _HIGHER_ORDER_ */
	  }
	
#endif /* ALL_TOGETHER_AGAIN */
        T0[res] *= T0[arg];
	break;

/*--------------------------------------------------------------------------*/
      case incr_a:                        /* Increment an adouble    incr_a */
        res   = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res]++;
	break;

/*--------------------------------------------------------------------------*/
      case decr_a:                        /* Increment an adouble    decr_a */
        res   = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res]--;
	break;


/****************************************************************************/
/*                                                        BINARY OPERATIONS */

/*--------------------------------------------------------------------------*/
      case plus_a_a:                 /* : Add two adoubles. (+)    plus a_a */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = T0[arg1] + T0[arg2];
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres,  T[res])
	ASSIGN_T(Targ1, T[arg1])
	ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG1_INC + TARG2_INC;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case plus_d_a:             /* Add an adouble and a double    plus_d_a */
	                         /* (+) */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = T0[arg] + coval;
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case min_a_a:              /* Subtraction of two adoubles     min_a_a */
	                         /* (-) */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = T0[arg1] - T0[arg2];
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
	ASSIGN_T(Targ1, T[arg1])
	ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG1_INC - TARG2_INC;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case min_d_a:                /* Subtract an adouble from a    min_d_a */
	                           /* double (-) */
        arg =get_locint_f();
        res = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = coval - T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC = -TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case mult_a_a:               /* Multiply two adoubles (*)    mult_a_a */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        /* olvo 980915 now in reverse order to allow x = x*x etc. */
        INC_pk_1(Tres)
        INC_pk_1(Targ1)
        INC_pk_1(Targ2)

        FOR_p_GT_l_GE_0
          FOR_k_GT_i_GE_0
	  { TRES_FODEC = T0[arg1]*TARG2_DEC + TARG1_DEC*T0[arg2]; 
       DEC_TRES_FO
#if defined(_HIGHER_ORDER_) 
            Targ1OP = Targ1-i+1;
            Targ2OP = Targ2;
            
            for (j=0;j<i;j++)
              *Tres += (*Targ1OP++) * (*Targ2OP--);
            Tres--;
#endif /* _HIGHER_ORDER_ */
          } 
                  
#endif /* ALL_TOGETHER_AGAIN */
        T0[res] = T0[arg1] * T0[arg2];
	break;

/*--------------------------------------------------------------------------*/
      /* olvo 991122: new op_code with recomputation */
      case eq_plus_prod:   /* increment a product of           eq_plus_prod */
                           /* two adoubles (*) */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        /* olvo 980915 now in reverse order to allow x = x*x etc. */
        INC_pk_1(Tres)
        INC_pk_1(Targ1)
        INC_pk_1(Targ2)

        FOR_p_GT_l_GE_0
          FOR_k_GT_i_GE_0
	  { TRES_FODEC += T0[arg1]*TARG2_DEC + TARG1_DEC*T0[arg2]; 
       DEC_TRES_FO
#if defined(_HIGHER_ORDER_) 
            Targ1OP = Targ1-i+1;
            Targ2OP = Targ2;
            
            for (j=0;j<i;j++)
              *Tres += (*Targ1OP++) * (*Targ2OP--);
            Tres--;
#endif /* _HIGHER_ORDER_ */
          } 
                  
#endif /* ALL_TOGETHER_AGAIN */
        T0[res] += T0[arg1] * T0[arg2];
	break;

/*--------------------------------------------------------------------------*/
      /* olvo 991122: new op_code with recomputation */
      case eq_min_prod:    /* decrement a product of            eq_min_prod */
                           /* two adoubles (*) */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        /* olvo 980915 now in reverse order to allow x = x*x etc. */
        INC_pk_1(Tres)
        INC_pk_1(Targ1)
        INC_pk_1(Targ2)

        FOR_p_GT_l_GE_0
          FOR_k_GT_i_GE_0
	  { TRES_FODEC -= T0[arg1]*TARG2_DEC + TARG1_DEC*T0[arg2]; 
       DEC_TRES_FO
#if defined(_HIGHER_ORDER_) 
            Targ1OP = Targ1-i+1;
            Targ2OP = Targ2;
            
            for (j=0;j<i;j++)
              *Tres -= (*Targ1OP++) * (*Targ2OP--);
            Tres--;
#endif /* _HIGHER_ORDER_ */
          } 
                  
#endif /* ALL_TOGETHER_AGAIN */
        T0[res] -= T0[arg1] * T0[arg2];
	break;

/*--------------------------------------------------------------------------*/
      case mult_d_a:         /* Multiply an adouble by a double    mult_d_a */
                             /* (*) */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = T0[arg] * coval;
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG_INC * coval;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case div_a_a:           /* Divide an adouble by an adouble    div_a_a */
                              /* (/) */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

#if !defined(_ZOS_) /* BREAK_ZOS */
	divs = 1.0 / T0[arg2];	
#endif /* ALL_TOGETHER_AGAIN */

        T0[res] = T0[arg1] / T0[arg2];
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_p
          FOR_0_LE_i_LT_k
	  { /* olvo 980922 changed order to allow x = y/x */ 
#if defined(_HIGHER_ORDER_)
            zOP      = z+i;
            (*zOP--) = -(*Targ2) * divs;
#endif /* _HIGHER_ORDER_ */

            TRES_FOINC = TARG1_INC * divs + T0[res] * (-TARG2_INC * divs);

#if defined(_HIGHER_ORDER_)
            TresOP = Tres-i;

	    for (j=0;j<i;j++)
	      *Tres += (*TresOP++) * (*zOP--);
            Tres++;
#endif /* _HIGHER_ORDER_ */
	  }         

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case div_d_a:             /* Division double - adouble (/)    div_d_a */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	/* olvo 980922 necessary for reverse */
        if (arg == res)
        { IF_KEEP_WRITE_TAYLOR(arg,keep,k,p)
	}

#if !defined(_ZOS_) /* BREAK_ZOS */
	divs = 1.0 / T0[arg];	
#endif /* ALL_TOGETHER_AGAIN */

        T0[res] = coval / T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_p
          FOR_0_LE_i_LT_k
	  { /* olvo 980922 changed order to allow x = d/x */ 
#if defined(_HIGHER_ORDER_)
            zOP      = z+i;
            (*zOP--) = -(*Targ) * divs;
#endif /* _HIGHER_ORDER_ */

            TRES_FOINC = T0[res] * (-TARG_INC * divs);

#if defined(_HIGHER_ORDER_)
            TresOP = Tres-i;

	    for (j=0;j<i;j++)
	      *Tres += (*TresOP++) * (*zOP--);
            Tres++;
#endif /* _HIGHER_ORDER_ */
	  }         

#endif /* ALL_TOGETHER_AGAIN */
	break;


/****************************************************************************/
/*                                                         SIGN  OPERATIONS */

/*--------------------------------------------------------------------------*/
      case pos_sign_a:                                        /* pos_sign_a */
        arg   = get_locint_f();
        res   = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case neg_sign_a:                                        /* neg_sign_a */
        arg   = get_locint_f();
        res   = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = -T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC = -TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
	break;


/****************************************************************************/
/*                                                         UNARY OPERATIONS */

/*--------------------------------------------------------------------------*/
      case exp_op:                          /* exponent operation    exp_op */
        arg = get_locint_f();
        res = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = exp(T0[arg]);
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_p
          FOR_0_LE_i_LT_k
	  { /* olvo 980915 changed order to allow x = exp(x) */
#if defined(_HIGHER_ORDER_)            
            zOP      = z+i;
	    (*zOP--) = (i+1) * (*Targ);
#endif /* _HIGHER_ORDER_ */

            TRES_FOINC = T0[res] * TARG_INC;

#if defined(_HIGHER_ORDER_)            
            TresOP = Tres-i;

            *Tres *= (i+1);
	    for (j=0;j<i;j++)
	      *Tres += (*TresOP++) * (*zOP--);
	    *Tres++ /= (i+1); /* important only for i>0 */
#endif /* _HIGHER_ORDER_ */
          }
        

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case sin_op:                              /* sine operation    sin_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(arg2,keep,k,p) /* olvo 980710 covalue */
        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	T0[arg2] = cos(T0[arg1]); /* Note: always arg2 != arg1 */
	T0[res]  = sin(T0[arg1]);
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres,  T[res]) 
	ASSIGN_T(Targ1, T[arg1])	
        ASSIGN_T(Targ2, T[arg2])
        
        FOR_0_LE_l_LT_p
          FOR_0_LE_i_LT_k
	  { /* olvo 980921 changed order to allow x = sin(x) */
#if defined(_HIGHER_ORDER_)            
            zOP      = z+i;
	    (*zOP--) = (i+1) * (*Targ1);
#endif /* _HIGHER_ORDER_ */

            TARG2_FOINC = -T0[res]  * TARG1; /* Note: always arg2 != arg1 */ 
            TRES_FOINC  =  T0[arg2] * TARG1_INC; 

#if defined(_HIGHER_ORDER_)
            TresOP  = Tres-i;
            Targ2OP = Targ2-i;          

            *Tres  *= (i+1);
            *Targ2 *= (i+1);
	    for (j=0;j<i;j++)
	    { *Tres  += (*Targ2OP++) * (*zOP); 
              *Targ2 -= (*TresOP++)  * (*zOP--);
	    } 
	    *Targ2++ /= (i+1);
	    *Tres++  /= (i+1);
#endif /* _HIGHER_ORDER_ */
	  }

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case cos_op:                            /* cosine operation    cos_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(arg2,keep,k,p) /* olvo 980710 covalue */
        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[arg2] = sin(T0[arg1]); /* Note: always arg2 != arg1 */
        T0[res]  = cos(T0[arg1]);
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres,  T[res])
	ASSIGN_T(Targ1, T[arg1])
	ASSIGN_T(Targ2, T[arg2])
 
        FOR_0_LE_l_LT_p
          FOR_0_LE_i_LT_k
	  { /* olvo 980921 changed order to allow x = cos(x) */
#if defined(_HIGHER_ORDER_)            
            zOP      = z+i;
	    (*zOP--) = (i+1) * (*Targ1);
#endif /* _HIGHER_ORDER_ */

            TARG2_FOINC =  T0[res]  * TARG1; /* Note: always arg2 != arg1 */
            TRES_FOINC  = -T0[arg2] * TARG1_INC;

#if defined(_HIGHER_ORDER_)   
            TresOP  = Tres-i;
            Targ2OP = Targ2-i;

            *Tres  *= (i+1);
            *Targ2 *= (i+1);  
	    for (j=0;j<i;j++)
	    { *Tres  -= (*Targ2OP++) * (*zOP); 
              *Targ2 += (*TresOP++)  * (*zOP--);
	    }
	    *Targ2++ /= (i+1);
	    *Tres++  /= (i+1);
#endif /* _HIGHER_ORDER_ */
	  }

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case atan_op:                                              /* atan_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res]=atan(T0[arg1]);
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_p
        { FOR_0_LE_i_LT_k
	  { /* olvo 980921 changed order to allow x = atan(x) */
#if defined(_HIGHER_ORDER_)            
            zOP      = z+i;
	    (*zOP--) = (i+1) * (*Targ1);
#endif /* _HIGHER_ORDER_ */

            TRES_FOINC = T0[arg2] * TARG1_INC;

#if defined(_HIGHER_ORDER_)
            Targ2OP = Targ2;          

            *Tres *= (i+1); 
	    for (j=0;j<i;j++)
	      *Tres  += (*Targ2OP++) * (*zOP--); 
	    *Tres++ /= (i+1);
#endif /* _HIGHER_ORDER_ */
	  }
          HOV_INC(Targ2, k)
        }

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case asin_op:                                              /* asin_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = asin(T0[arg1]);
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        if (T0[arg1] == 1.0)
          FOR_0_LE_l_LT_p
          { FOR_0_LE_i_LT_k
              if (TARG1 > 0.0)
              { r0 = make_nan();
                VEC_INC(Targ1, k-i)
                BREAK_FOR_I
              }
              else 
                if (TARG1 < 0.0)
                { r0 = make_inf();
                  VEC_INC(Targ1, k-i)
                  BREAK_FOR_I
                }
                else 
                { r0 = 0.0;
                  Targ1++;
                }
            TRES = r0;
            VEC_INC(Tres, k)
          }
        else 
          if (T0[arg1] == -1.0)
            FOR_0_LE_l_LT_p
            { FOR_0_LE_i_LT_k
                if (TARG1 > 0.0)
                { r0 = make_inf();
                  VEC_INC(Targ1, k-i)
                  BREAK_FOR_I
                }
                else 
                  if (TARG1 < 0.0)
                  { r0 = make_nan();
                    VEC_INC(Targ1, k-i)
                    BREAK_FOR_I
                  }
                  else 
                  { r0 = 0.0;
                    Targ1++;
                  }
              TRES = r0;
              VEC_INC(Tres, k)
            }
          else 
            FOR_0_LE_l_LT_p
            { FOR_0_LE_i_LT_k
              { /* olvo 980921 changed order to allow x = asin(x) */
#if defined(_HIGHER_ORDER_)
                zOP      = z+i;            
                (*zOP--) = (i+1) * (*Targ1);
#endif /* _HIGHER_ORDER_ */

                TRES_FOINC = T0[arg2] * TARG1_INC;

#if defined(_HIGHER_ORDER_)
                Targ2OP = Targ2;

                *Tres *= (i+1);  
  	        for (j=0;j<i;j++)
	          *Tres += (*Targ2OP++) * (*zOP--); 
	        *Tres++ /= (i+1);
#endif /* _HIGHER_ORDER_ */
	      }
              HOV_INC(Targ2, k)
            }

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case acos_op:                                              /* acos_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = acos(T0[arg1]);
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        if (T0[arg1] == 1.0)
          FOR_0_LE_l_LT_p
          { FOR_0_LE_i_LT_k
              if (TARG1 > 0.0)
              { r0 = make_nan();
                VEC_INC(Targ1, k-i)
                BREAK_FOR_I
              } 
              else 
                if (TARG1 < 0.0)
                { r0 = -make_inf();
                  VEC_INC(Targ1, k-i)
                  BREAK_FOR_I
                }
                else 
                { r0 = 0.0;
                  Targ1++;
                }
            TRES = r0;  
            VEC_INC(Tres, k)
          }
        else 
          if (T0[arg1] == -1.0)
            FOR_0_LE_l_LT_p
            { FOR_0_LE_i_LT_k
                if (TARG1 > 0.0)
                { r0 = -make_inf();
                  VEC_INC(Targ1, k-i)
                  BREAK_FOR_I
                }
                else 
                  if (TARG1 < 0.0)
                  { r0 = make_nan();
                    VEC_INC(Targ1, k-i)
                    BREAK_FOR_I
                  }
                  else 
                  { r0 = 0.0;
                    Targ1++;
                  }
              TRES = r0;
              VEC_INC(Tres, k)
            }
          else 
            FOR_0_LE_l_LT_p
            { FOR_0_LE_i_LT_k
	      { /* olvo 980921 changed order to allow x = acos(x) */
#if defined(_HIGHER_ORDER_)
                zOP      = z+i;            
                (*zOP--) = (i+1) * (*Targ1);
#endif /* _HIGHER_ORDER_ */

                TRES_FOINC = T0[arg2] * TARG1_INC;

#if defined(_HIGHER_ORDER_)
                Targ2OP = Targ2;         

                *Tres *= (i+1);  
  	        for (j=0;j<i;j++)
	          *Tres += (*Targ2OP++) * (*zOP--); 
	        *Tres++ /= (i+1);
#endif /* _HIGHER_ORDER_ */
	      }
              HOV_INC(Targ2, k)
            }

#endif /* ALL_TOGETHER_AGAIN */
	break;

#ifdef ATRIG_ERF

/*--------------------------------------------------------------------------*/
      case asinh_op:                                            /* asinh_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = asinh(T0[arg1]);
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_p
        { FOR_0_LE_i_LT_k
          { /* olvo 980921 changed order to allow x = asinh(x) */
#if defined(_HIGHER_ORDER_)
            zOP      = z+i;            
            (*zOP--) = (i+1) * (*Targ1);
#endif /* _HIGHER_ORDER_ */

            TRES_FOINC = T0[arg2] * TARG1_INC;

#if defined(_HIGHER_ORDER_)
            Targ2OP = Targ2;
            
            *Tres *= (i+1);	    
	    for (j=0;j<i;j++)
	      *Tres += (*Targ2OP++) * (*zOP--); 
	    *Tres++ /= (i+1);
#endif /* _HIGHER_ORDER_ */
	  }
          HOV_INC(Targ2, k)
        }

#endif /* ALL_TOGETHER_AGAIN */
        break;

/*--------------------------------------------------------------------------*/
       case acosh_op:                                           /* acosh_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = acosh(T0[arg1]);
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])
 
        if (T0[arg1] == 1.0)
          FOR_0_LE_l_LT_p
          { FOR_0_LE_i_LT_k 
              if (TARG1 > 0.0)
              { r0 = make_inf();
                VEC_INC(Targ1, k-i)
                BREAK_FOR_I
              }
              else 
                if (TARG1 < 0.0) 
                { r0 = make_nan();
                  VEC_INC(Targ1, k-i)
                  BREAK_FOR_I
                }
                else 
                { r0 = 0.0;
                  Targ1++;
                }
            TRES_INC = r0;  
#if defined(_HIGHER_ORDER_)
            for (i=1;i<k;i++)
              *Tres++ = make_nan();
#endif /* _HIGHER_ORDER_ */
          }
        else 
          FOR_0_LE_l_LT_p
          { FOR_0_LE_i_LT_k
            { /* olvo 980921 changed order to allow x = acosh(x) */
#if defined(_HIGHER_ORDER_)
              zOP      = z+i;            
              (*zOP--) = (i+1) * (*Targ1);
#endif /* _HIGHER_ORDER_ */

              TRES_FOINC = T0[arg2] * TARG1_INC;

#if defined(_HIGHER_ORDER_)
              Targ2OP = Targ2;            

              *Tres *= (i+1);
	      for (j=0;j<i;j++)
	        *Tres += (*Targ2OP++) * (*zOP--); 
	      *Tres++ /= (i+1);
#endif /* _HIGHER_ORDER_ */
	    }
            HOV_INC(Targ2, k)
          }

#endif /* ALL_TOGETHER_AGAIN */
        break;

/*--------------------------------------------------------------------------*/
      case atanh_op:                                            /* atanh_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = atanh(T0[arg1]);
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        if (T0[arg1] == 1.0)
          FOR_0_LE_l_LT_p
          { FOR_0_LE_i_LT_k  
              if (TARG1 > 0.0)
              { r0 = make_nan();
                VEC_INC(Targ1, k-i)
                BREAK_FOR_I
              }
              else 
                if (TARG1 < 0.0)
                { r0 = make_inf();
                  VEC_INC(Targ1, k-i)
                  BREAK_FOR_I
                } 
                else 
                { r0 = 0.0;
                  Targ1++;
                }
            TRES_INC = r0;
#if defined(_HIGHER_ORDER_)
            for (i=1;i<k;i++) 
              *Tres++ = make_nan();
#endif /* _HIGHER_ORDER_ */
          }
        else 
          if (T0[arg1] == -1.0)
            FOR_0_LE_l_LT_p
            { FOR_0_LE_i_LT_k
                if (TARG1 > 0.0)
                { r0 = make_inf();
                  VEC_INC(Targ1, k-i)
                  BREAK_FOR_I
                }
                else 
                  if (TARG1 < 0.0)
                  { r0 = make_nan();
                    VEC_INC(Targ1, k-i)
                    BREAK_FOR_I
                  } 
                  else 
                  { r0 = 0.0;
                    Targ1++;
                  }
              TRES_INC = r0;
#if defined(_HIGHER_ORDER_)
              for (i=1;i<k;i++) 
                *Tres++ = make_nan();
#endif /* _HIGHER_ORDER_ */
            }
          else 
            FOR_0_LE_l_LT_p
            { FOR_0_LE_i_LT_k
              { /* olvo 980921 changed order to allow x = atanh(x) */
#if defined(_HIGHER_ORDER_)
                zOP      = z+i;            
                (*zOP--) = (i+1) * (*Targ1);
#endif /* _HIGHER_ORDER_ */

                TRES_FOINC = T0[arg2] * TARG1_INC;

#if defined(_HIGHER_ORDER_)
                Targ2OP = Targ2;            

                *Tres *= (i+1);
	        for (j=0;j<i;j++)
	          *Tres += (*Targ2OP++) * (*zOP--); 
	        *Tres++ /= (i+1);
#endif /* _HIGHER_ORDER_ */
              }
              HOV_INC(Targ2, k)
            }

#endif /* ALL_TOGETHER_AGAIN */
        break;

/*--------------------------------------------------------------------------*/
      case erf_op:                                                /* erf_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = erf(T0[arg1]);
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ1,T[arg1])
        ASSIGN_T(Targ2,T[arg2])
        
        FOR_0_LE_l_LT_p
        { FOR_0_LE_i_LT_k
          { /* olvo 980921 changed order to allow x = erf(x) */
#if defined(_HIGHER_ORDER_)
            zOP      = z+i;            
            (*zOP--) = (i+1) * (*Targ1);
#endif /* _HIGHER_ORDER_ */

            TRES_FOINC = T0[arg2] * TARG1_INC;

#if defined(_HIGHER_ORDER_)
            Targ2OP = Targ2;

            *Tres *= (i+1);
	    for (j=0;j<i;j++)
	      *Tres += (*Targ2OP++) * (*zOP--); 
	    *Tres++ /= (i+1);
#endif /* _HIGHER_ORDER_ */
	  }
          HOV_INC(Targ2, k)
        }

#endif /* ALL_TOGETHER_AGAIN */
        break;
           
#endif

/*--------------------------------------------------------------------------*/
      case log_op:                                                /* log_op */
        arg = get_locint_f();
        res = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        divs = 1.0 / T0[arg];
        FOR_0_LE_l_LT_p
        { if (T0[arg] == 0.0)
          { TargOP = Targ;
            FOR_0_LE_i_LT_k
            { if (*TargOP++ < 0.0)
              { divs = make_nan();
                BREAK_FOR_I
              }
	    }
          }

          /* olvo 980921 changed order to allow x = log(x) */
          FOR_0_LE_i_LT_k
	  { TRES_FOINC = TARG_INC * divs;
#if defined(_HIGHER_ORDER_)	    
            TresOP = Tres - i;
            zOP    = z+i;
   
            (*zOP--) = *Tres;
	    (*Tres) *= i+1;
 	    for (j=0;j<i;j++)
	      (*Tres) -= (*zOP--) * (*TresOP++) * (j+1); 
	    *Tres++ /= i+1;
#endif /* _HIGHER_ORDER_ */
          }
        }

#endif /* ALL_TOGETHER_AGAIN */
        T0[res] = log(T0[arg]);
	break;

/*--------------------------------------------------------------------------*/
      case pow_op:                                                /* pow_op */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	/* olvo 980921 necessary for reverse */
        if (arg == res)
        { IF_KEEP_WRITE_TAYLOR(arg,keep,k,p)
	}

#if !defined(_ZOS_) /* BREAK_ZOS */
        T0arg   = T0[arg];
#endif /* ALL_TOGETHER_AGAIN */ 

        T0[res] = pow(T0[arg], coval);
#if !defined(_ZOS_) /* BREAK_ZOS */
	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        if (T0arg == 0.0)
          FOR_0_LE_l_LT_p
          { x = TARG;
	    y = TARG;            
	    i = 0;
            r0 = 0;
            FOR_0_LE_i_LT_k
	      { /* current derivative is zero */
              if (coval - i > 1)
	       {   
                Targ++;
                TRES_INC = r0;
                x *= y;
               }
              else
		{ /* would result in 1/0 */
      	        if (coval <= 0.0)
                  r0 = make_nan();
                else
		  /* would result in 1/0 */
                 if (coval - i < 1)
                   r0 = make_inf();
		  /* current derivative is one times x_0^coval*/ 
                 else
                   r0 = x;   
                TRES_INC = r0;
                VEC_INC(Targ, k-i)
                BREAK_FOR_I
               }
            }
#if defined(_HIGHER_ORDER_)
            for (j=i+1;j<k;j++)
            { if ((coval <= 0.0) || (coval - i < 1))
                *Tres++ = make_nan();
              else
                *Tres++ = 0;
            }
#endif /* _HIGHER_ORDER_ */
          }
        else
  	{ r0 = 1.0 / T0arg;
          FOR_0_LE_l_LT_p
            FOR_0_LE_i_LT_k
            { /* olvo 980921 changed order to allow x = pow(x,n) */
#if defined(_HIGHER_ORDER_)
              /* printf(" hier Tres = %e \n",*Tres); */
              zOP      = z+i;
              (*zOP--) = (*Targ) * r0;
#endif /* _HIGHER_ORDER_ */

              TRES_FOINC = T0[res] * TARG_INC * coval * r0;

#if defined(_HIGHER_ORDER_)
              TresOP = Tres-i;            

 	      (*Tres) *= i+1;           
	      y = coval*i -1;            
              for (j=0;j<i;j++)
              { *Tres += (*TresOP++) * (*zOP--) * y;
		y -= coval + 1;
	      }
              /* printf(" hier Tres = %e \n",*Tres); */
	      *Tres++ /= (i+1);
#endif /* _HIGHER_ORDER_ */
	    }
        }

#endif /* ALL_TOGETHER_AGAIN */ 
	break;

/*--------------------------------------------------------------------------*/
      case sqrt_op:                                              /* sqrt_op */
        arg = get_locint_f();
        res = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = sqrt(T0[arg]);
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Targ, T[arg])
	ASSIGN_T(Tres, T[res])	

        FOR_0_LE_l_LT_p
        { TargOP = Targ;
          if (T0[arg] == 0.0)   /* Note: <=> T0[res] == 0.0 */
          { r0 = 0.0;
            FOR_0_LE_i_LT_k
	    { if (TARG>0.0)
              { r0 = make_inf();
                VEC_INC(Targ, k-i)
                BREAK_FOR_I
              }
              else 
               if (TARG<0.0)
               { r0 = make_nan();
                 VEC_INC(Targ, k-i)
                 BREAK_FOR_I
               }
               else
                 Targ++;
            }
          }
          else 
          { r0 = 0.5/T0[res];
          }
          Targ = TargOP;

	  even = 1;
          FOR_0_LE_i_LT_k
	  { TRES_FOINC = r0 * TARG_INC;
#if defined(_HIGHER_ORDER_)
            TresOP  = Tres-i;            
            TresOP2 = Tres-1;            

            x = 0;
            for (j=1;2*j-1<i;j++)
              x += (*TresOP++) * (*TresOP2--);
	    x *= 2;
	    if (!even) 
              x += (*TresOP) * (*TresOP2); /* !!! */
	    even = !even;
	    *Tres++ -= r0*x;
#endif /* _HIGHER_ORDER_ */
	  }
        }

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case gen_quad:                                            /* gen_quad */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

        if (get_val_f()!=T0[arg1])
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: forward sweep aborted; tape invalid!\n");
          IF_KEEP_TAYLOR_CLOSE
          end_sweep();
          return -2;
        }
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = coval;
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_p
        { FOR_0_LE_i_LT_k
	  { /* olvo 980922 changed order to allow x = gen_quad(x) */
#if defined(_HIGHER_ORDER_)
            zOP      = z+i;
	    (*zOP--) = (i+1) * (*Targ1);
#endif /* _HIGHER_ORDER_ */

            TRES_FOINC = T0[arg2] * TARG1_INC;

#if defined(_HIGHER_ORDER_)
            Targ2OP = Targ2;
            
            *Tres *= (i+1);
	    for (j=0;j<i;j++)
	      *Tres += (*Targ2OP++) * (*zOP--);  
	    *Tres++ /= (i+1);
#endif /* _HIGHER_ORDER_ */
	  }
          HOV_INC(Targ2, k)
        } 

#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case min_op:                                                /* min_op */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	/* olvo 980923 changed order to allow x = min(x,y) etc. */ 

        /* olvo/mitev 980721 return value (taken from below) */
        if (T0[arg1] > T0[arg2])
        { if (coval)
            MINDEC(ret_c,2);
        }
        else
          if (T0[arg1] < T0[arg2])
          { if (!coval)
              MINDEC(ret_c,2);
          }
          else
            if (arg1 != arg2)
              MINDEC(ret_c,1);

#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])
        ASSIGN_T(Tres,  T[res])

        Tqo = NULL;
        if (T0[arg1] > T0[arg2])
          Tqo = Targ2;
        else 
          if (T0[arg1] < T0[arg2])
            Tqo = Targ1;

        FOR_0_LE_l_LT_p
        { Targ = Tqo;
          if (Targ == NULL) /* e.g. T0[arg1] == T0[arg2] */
          { Targ1OP = Targ1;
            Targ2OP = Targ2;
            FOR_0_LE_i_LT_k
            { if (TARG1 > TARG2)
              { Targ = Targ2OP;
                VEC_INC(Targ1, k-i)
                VEC_INC(Targ2, k-i)
                BREAK_FOR_I
              }             
              else 
                if (TARG1 < TARG2)
                { Targ = Targ1OP;
                  VEC_INC(Targ1, k-i)
                  VEC_INC(Targ2, k-i)
                  BREAK_FOR_I
                }
              Targ1++; Targ2++;
            }          
            if (Targ == NULL) /* e.g. both are equal */
              Targ = Targ1OP;
	  }

          FOR_0_LE_i_LT_k
            TRES_INC = TARG_INC;

          if (Tqo)
          { VEC_INC(Tqo, k)
	  }
        }

#endif /* ALL_TOGETHER_AGAIN */
        T0[res] = FMIN(T0[arg1], T0[arg2]);
        break;

/*--------------------------------------------------------------------------*/
      case abs_val:                                              /* abs_val */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	/* olvo 980923 changed order to allow x = min(x,y) etc. */ 

        /* olvo/mitev 980721 ec n3l (taken from below) */
        if (T0[arg] < 0.0) 
        { if (coval)
            MINDEC(ret_c,2);
        }
        else
          if (T0[arg] > 0.0)
          { if (!coval)
              MINDEC(ret_c,2);
          }

#if !defined(_ZOS_) /* BREAK_ZOS */        
        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])
  
        y = 0.0;
        if (T0[arg] != 0.0) {
          if (T0[arg] < 0.0)
            y = -1.0;
          else
            y = 1.0;
        }

        FOR_0_LE_l_LT_p
        { x = y;
          FOR_0_LE_i_LT_k
          { if ((x == 0.0) && (TARG != 0.0))
            { MINDEC(ret_c,1);
              if (TARG < 0.0)
                x = -1.0;
              else
                x = 1.0;
            }
            TRES_INC = x * TARG_INC;
          }
        } 

#endif /* ALL_TOGETHER_AGAIN */
        T0[res] = fabs(T0[arg]);
        break;

/*--------------------------------------------------------------------------*/
      case ceil_op:                                              /* ceil_op */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res]=ceil(T0[arg]);
        /* olvo/mitev 980721 ec n2l (taken from below) */
        if (coval != T0[res])
          MINDEC(ret_c,2);
#if !defined(_ZOS_) /* BREAK_ZOS */        
        ASSIGN_T(Tres, T[res])

        FOR_0_LE_l_LT_pk
          TRES_INC = 0.0;

#endif /* ALL_TOGETHER_AGAIN */
        break;

/*--------------------------------------------------------------------------*/
      case floor_op:                 /* Compute ceil of adouble    floor_op */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = floor(T0[arg]);
        /* olvo/mitev 980721 ec n2l (taken from below) */
        if (coval != T0[res])
          MINDEC(ret_c,2);
#if !defined(_ZOS_) /* BREAK_ZOS */        
        ASSIGN_T(Tres, T[res])

        FOR_0_LE_l_LT_pk
          TRES_INC = 0.0;

#endif /* ALL_TOGETHER_AGAIN */
        break;


/****************************************************************************/
/*                                                             CONDITIONALS */

/*--------------------------------------------------------------------------*/
      case cond_assign:                                      /* cond_assign */
        arg   = get_locint_f();
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        /* olvo 980924 changed order to allow reflexive ops */
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        if (T0[arg] > 0)
          FOR_0_LE_l_LT_pk
            TRES_INC = TARG1_INC;
        else
          FOR_0_LE_l_LT_pk
            TRES_INC = TARG2_INC;
#endif /* ALL_TOGETHER_AGAIN */

        if (T0[arg] > 0)
        { if (coval <= 0.0)
            MINDEC(ret_c,2);
          T0[res] = T0[arg1];
        }
        else
        { if (coval > 0.0)
            MINDEC(ret_c,2);
          if (T0[arg] == 0)
            MINDEC(ret_c,0);
          T0[res] = T0[arg2];
        }
        break;

/*--------------------------------------------------------------------------*/
      case cond_assign_s:                                  /* cond_assign_s */
        arg   = get_locint_f();
        arg1  = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        /* olvo 980924 changed order to allow reflexive ops */
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
 
        if (T0[arg] > 0)
          FOR_0_LE_l_LT_pk
            TRES_INC = TARG1_INC;
#endif /* ALL_TOGETHER_AGAIN */

        if (T0[arg] > 0)
	{ if (coval <= 0.0)
            MINDEC(ret_c,2);
          T0[res] = T0[arg1];
        } 
        else
          if (T0[arg] == 0)
            MINDEC(ret_c,0); 
        break;


/****************************************************************************/
/*                                                       VECTOR ASSIGNMENTS */

/*--------------------------------------------------------------------------*/
      case assign_av:                                          /* assign_av */
        arg  = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();

        for (ls=0; ls<size; ls++)
        { /* code for assign_a */
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] = T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Targ, T[arg])
          ASSIGN_T(Tres, T[res])
            
          FOR_0_LE_l_LT_pk
	    TRES_INC = TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
          res++; arg++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case assign_dv:                                          /* assign_dv */
        size  = get_locint_f();
        res   = get_locint_f();
        d     = get_val_v_f(size);

        for (ls=0; ls<size; ls++)
        { /* code for assign_d */
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] = *d++;
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres, T[res])
   
          FOR_0_LE_l_LT_pk
	    TRES_INC = 0;

#endif /* ALL_TOGETHER_AGAIN */
          res++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case assign_indvec:                                  /* assign_indvec */
        size = get_locint_f();
        res  = get_locint_f();

        for (ls=0; ls<size; ls++)
        { /* code for assign_ind */
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] = basepoint[indexi];
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres, T[res])
   
          FOR_0_LE_l_LT_p
            FOR_0_LE_i_LT_k      
              TRES_INC = ARGUMENT(indexi,l,i);

#endif /* ALL_TOGETHER_AGAIN */
          ++indexi;
          res++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case assign_depvec:                                  /* assign_depvec */
        size = get_locint_f();
        res  = get_locint_f();

        for (ls=0; ls<size; ls++)
        { /* code for assign_dep */

          valuepoint[indexd] = T0[res];
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres, T[res])

          if (taylors != 0 )  
            FOR_0_LE_l_LT_p
              FOR_0_LE_i_LT_k      
	        TAYLORS(indexd,l,i) = TRES_INC;

#endif /* ALL_TOGETHER_AGAIN */
          indexd++;
          res++;
        } 
        break;


/****************************************************************************/
/*                                            VECTOR OPERATION + ASSIGNMENT */

/*--------------------------------------------------------------------------*/
      case eq_plus_av:                                        /* eq_plus_av */
        arg  = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();

        for (ls=0; ls<size; ls++)
        { /* code for eq_plus_a */
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] += T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */            
          ASSIGN_T(Tres, T[res])
          ASSIGN_T(Targ, T[arg])
            
          FOR_0_LE_l_LT_pk
	    TRES_INC += TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
          res++; arg++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case eq_min_av:                                          /* eq_min_av */
        arg  = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();

        for (ls=0; ls<size; ls++)
        { /* code for eq_min_a */
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] -= T0[arg]; 
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres, T[res])
          ASSIGN_T(Targ, T[arg])
            
          FOR_0_LE_l_LT_pk
	    TRES_INC -= TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
          res++; arg++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case eq_mult_av_d:                                    /* eq_mult_av_d */
        size  = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        for (ls=0; ls<size; ls++)
        { /* code for eq_mult_d*/
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] *= coval;
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres, T[res])
            
          FOR_0_LE_l_LT_pk
	      TRES_INC *= coval;

#endif /* ALL_TOGETHER_AGAIN */
          res++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case eq_mult_av_a:                                    /* eq_mult_av_a */
        arg  = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();

        /* olvo 980929 new strategy to check for overwrites 
           (changes computation order) */
        checkSize = res+size;  

        for (ls=0; ls<size; ls++)
        { if (res == arg) /* skip res==arg first */
            res++;
          if (res == checkSize) /* checks if arg==res was skipped */
            res = arg; 

          /* code for eq_mult_a*/
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres, T[res])
          ASSIGN_T(Targ, T[arg])
            
          INC_pk_1(Tres)
          INC_pk_1(Targ)

          FOR_p_GT_l_GE_0
	    FOR_k_GT_i_GE_0
	    { TRES_FODEC = T0[res]*TARG_DEC + TRES*T0[arg]; 
         DEC_TRES_FO
#if defined(_HIGHER_ORDER_)
              TresOP = Tres-i;
              TargOP = Targ;

              for (j=0;j<i;j++)
                *Tres += (*TresOP++) * (*TargOP--);
              Tres--;
#endif /* _HIGHER_ORDER_ */
	    }

#endif /* ALL_TOGETHER_AGAIN */
          T0[res] *= T0[arg];
          res++;
        } 
        break;


/****************************************************************************/
/*                                                 BINARY VECTOR OPERATIONS */

/*--------------------------------------------------------------------------*/
      case plus_av_av:                                        /* plus_av_av */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();

        for (ls=0; ls<size; ls++)
        { /* code for plus_a_a */
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] = T0[arg1] + T0[arg2];
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres,  T[res])
          ASSIGN_T(Targ1, T[arg1])
          ASSIGN_T(Targ2, T[arg2])
            
          FOR_0_LE_l_LT_pk
	    TRES_INC = TARG1_INC + TARG2_INC;

#endif /* ALL_TOGETHER_AGAIN */
          res++; arg1++; arg2++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case sub_av_av:                                          /* sub_av_av */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();

        for (ls=0; ls<size; ls++)
        { /* code for min_a_a */
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] = T0[arg1] - T0[arg2];
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres,  T[res])
          ASSIGN_T(Targ1, T[arg1])
          ASSIGN_T(Targ2, T[arg2])
            
          FOR_0_LE_l_LT_pk
	    TRES_INC = TARG1_INC - TARG2_INC;

#endif /* ALL_TOGETHER_AGAIN */
          res++; arg1++; arg2++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case dot_av_av:                                          /* dot_av_av */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

	T0res = 0.0;
        /* olvo 980924 check for overwrites -- if necessary use 
           tempories for res-stuff to allow reflexive ops */
#if !defined(_ZOS_) /* BREAK_ZOS */
        if (   ((res >= arg1) && (res < arg1+size))
            || ((res >= arg2) && (res < arg2+size)) )
	{ ASSIGN_T(TresOP, Ttemp)
          flag = 1;
	}
        else 
        { ASSIGN_T(TresOP, T[res])
          flag = 0;
	}

        Tres = TresOP;
        FOR_0_LE_l_LT_pk
          TRES_INC = 0.0;
#endif /* ALL_TOGETHER_AGAIN */

        for (ls=0; ls<size; ls++)
        { /* code for mult_a_a  */
#if !defined(_ZOS_) /* BREAK_ZOS */
          Tres = TresOP;
          ASSIGN_T(Targ1, T[arg1])
          ASSIGN_T(Targ2, T[arg2])
           
          /* olvo 980915 now in reverse order to allow x = x*x etc. */
          INC_pk_1(Tres)
          INC_pk_1(Targ1)
          INC_pk_1(Targ2)

          FOR_p_GT_l_GE_0
            FOR_k_GT_i_GE_0
	    { TRES_FODEC += T0[arg1]*TARG2_DEC + TARG1_DEC*T0[arg2]; 
         DEC_TRES_FO
#if defined(_HIGHER_ORDER_) 
              Targ1OP = Targ1-i+1;
              Targ2OP = Targ2;
            
              for (j=0;j<i;j++)
                *Tres += (*Targ1OP++) * (*Targ2OP--);
              Tres--;
#endif /* _HIGHER_ORDER_ */
            } 
#endif /* ALL_TOGETHER_AGAIN */
          T0res += T0[arg1] * T0[arg2];
          arg1++; arg2++;
        }

        /* copy results if necessary */
        T0[res] = T0res;
#if !defined(_ZOS_) /* BREAK_ZOS */
        if (flag)
	{ ASSIGN_T(Tres,T[res]) 
          Tqo = TresOP;
  
          FOR_0_LE_l_LT_pk
            TRES_INC = TQO_INC;
	}
#endif /* ALL_TOGETHER_AGAIN */
	break;

/*--------------------------------------------------------------------------*/
      case mult_a_av:                                          /* mult_a_av */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();
         
        /* olvo 980929 new strategy to check for overwrites 
           (changes computation order) */
        checkSize = res+size;  

        for (ls=0; ls<size; ls++)
        { if (res == arg2) /* skip res==arg2 first */
          { res++; arg1++; 
	  }
          if (res == checkSize) /* checks if arg2==res was skipped */
          { arg1 -= res-arg2;
            res = arg2; 
	  }

          /* code for mult_a_a */
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres,  T[res])
          ASSIGN_T(Targ1, T[arg1])
          ASSIGN_T(Targ2, T[arg2])

          /* olvo 980915 now in reverse order to allow x = x*x etc. */
          INC_pk_1(Tres)
          INC_pk_1(Targ1)
          INC_pk_1(Targ2)

          FOR_p_GT_l_GE_0
            FOR_k_GT_i_GE_0
	    { TRES_FODEC = T0[arg1]*TARG2_DEC + TARG1_DEC*T0[arg2]; 
         DEC_TRES_FO
#if defined(_HIGHER_ORDER_) 
              Targ1OP = Targ1-i+1;
              Targ2OP = Targ2;
            
              for (j=0;j<i;j++)
                *Tres += (*Targ1OP++) * (*Targ2OP--);
              Tres--;
#endif /* _HIGHER_ORDER_ */
            } 
                  
#endif /* ALL_TOGETHER_AGAIN */
          T0[res] = T0[arg1] * T0[arg2];
          res++; arg1++;
        }
        break;

/*--------------------------------------------------------------------------*/
      case mult_d_av:                                          /* mult_d_av */
        arg   = get_locint_f();
        size  = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();

        for (ls=0; ls<size; ls++)
        { /* code for mult_d_a */
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] = T0[arg] * coval;
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres, T[res])
          ASSIGN_T(Targ, T[arg])
            
          FOR_0_LE_l_LT_pk
	    TRES_INC = TARG_INC * coval;

#endif /* ALL_TOGETHER_AGAIN */
          res++; arg++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case div_av_a:                                            /* div_av_a */
        arg1   = get_locint_f();
        arg2   = get_locint_f();
        size   = get_locint_f();
        res    = get_locint_f();

#if !defined(_ZOS_) /* BREAK_ZOS */
        divs   = 1.0 / T0[arg2];
#endif /* ALL_TOGETHER_AGAIN */

        /* olvo 980929 new strategy to check for overwrites 
           (changes computation order) */
        checkSize = res+size;  

        for (ls=0; ls<size; ls++)
        { if (res == arg2) /* skip res==arg2 first */
          { res++; arg1++; 
	  }
          if (res == checkSize) /* checks if arg2==res was skipped */
          { arg1 -= res-arg2;
            res = arg2; 
	  }

          /* code for div_a_a */
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] = T0[arg1] / T0[arg2];
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres,  T[res])
          ASSIGN_T(Targ1, T[arg1])
          ASSIGN_T(Targ2, T[arg2])

          FOR_0_LE_l_LT_p
            FOR_0_LE_i_LT_k
	    { /* olvo 980922 changed order to allow x = y/x */ 
#if defined(_HIGHER_ORDER_)
              zOP      = z+i;
              (*zOP--) = -(*Targ2) * divs;
#endif /* _HIGHER_ORDER_ */
 
 
              TRES_FOINC = TARG1_INC * divs + T0[res] * (-TARG2_INC * divs);

#if defined(_HIGHER_ORDER_)
              TresOP = Tres-i;

	      for (j=0;j<i;j++)
	        *Tres += (*TresOP++) * (*zOP--);
              Tres++;
#endif /* _HIGHER_ORDER_ */
	    }         

#endif /* ALL_TOGETHER_AGAIN */
          res++; arg1++;
        } 
        break;


/****************************************************************************/
/*                                                               SUBSCRIPTS */

/*--------------------------------------------------------------------------*/
      case subscript:                                          /* subscript */
        arg2 = get_locint_f(); /* Base */
        arg1 = get_locint_f(); /* pointer to variable containing offset */
        res  = get_locint_f();

        arg  = arg2 + (int)(T0[arg1]);

        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          MINDEC(ret_c,2);

        /* olvo 980721 new nl */
        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])
        
        FOR_0_LE_l_LT_pk
          TRES_INC = TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
        break;

/*--------------------------------------------------------------------------*/
      case subscript_l:                                      /* subscript_l */
        arg2 = get_locint_f();  /* Base */
	arg1 = get_locint_f();
        arg  = get_locint_f();         

        res  = arg2 + (int)(T0[arg1]);

        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          MINDEC(ret_c,2);

        IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

        T0[res] = T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])
        
        FOR_0_LE_l_LT_pk
          TRES_INC = TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
        break;

/*--------------------------------------------------------------------------*/
      case subscript_ld:                                    /* subscript_ld */
        arg2 = get_locint_f(); /* Base */
        arg1 = get_locint_f(); /* pointer to variable containing offset */

        arg  = arg2+(int)(T0[arg1]);

        IF_KEEP_WRITE_TAYLOR(arg,keep,k,p)

        T0[arg] = get_val_f();

        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          MINDEC(ret_c,2);
#if !defined(_ZOS_) /* BREAK_ZOS */
        ASSIGN_T(Targ, T[arg])
        
        FOR_0_LE_l_LT_pk
          TARG_INC = 0;

#endif /* ALL_TOGETHER_AGAIN */
        break;

/*--------------------------------------------------------------------------*/
      case m_subscript:                                      /* m_subscript */
        arg2 = get_locint_f(); /* Base */
        arg1 = get_locint_f(); /* pointer to variable containing offset */
        size = get_locint_f();
        res  = get_locint_f();

        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          MINDEC(ret_c,2);

        arg = arg2 + (int)(T0[arg1])*size;
        for (ls=0; ls<size; ls++)
        { /* olvo 980721 new nl */
          IF_KEEP_WRITE_TAYLOR(res,keep,k,p)
          
          T0[res] = T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres, T[res])   
          ASSIGN_T(Targ, T[arg])
            
          FOR_0_LE_l_LT_pk
            TRES_INC = TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
          res++; arg++;
        }
        break;

/*--------------------------------------------------------------------------*/
      case m_subscript_l:                                  /* m_subscript_l */
        arg2 = get_locint_f();  /* Base LHS */
	arg1 = get_locint_f();  /* Offset LHS */
        size = get_locint_f();
        arg  = get_locint_f();  /* RHS */

        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          MINDEC(ret_c,2);

        res = arg2 + (int)(T0[arg1])*size;
        for (ls=0; ls<size; ls++)
        { IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] = T0[arg];
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres, T[res])
          ASSIGN_T(Targ, T[arg])
            
          FOR_0_LE_l_LT_pk
            TRES_INC = TARG_INC;

#endif /* ALL_TOGETHER_AGAIN */
          res++; arg++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case m_subscript_ld:                                /* m_subscript_ld */
        arg2 = get_locint_f(); /* Base */
        arg1 = get_locint_f(); /* pointer to variable containing offset */
        arg  = get_locint_f(); /* offset in the vector itself */
        size = get_locint_f();

        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          MINDEC(ret_c,2);

        d      = get_val_v_f(size);

        res = arg2 + (int)(T0[arg1])*size + arg;
        for (ls=0; ls<size; ls++)
        { IF_KEEP_WRITE_TAYLOR(res,keep,k,p)

          T0[res] = d[l];
#if !defined(_ZOS_) /* BREAK_ZOS */
          ASSIGN_T(Tres, T[res])
            
          FOR_0_LE_l_LT_pk
            TRES_INC = 0.0;

#endif /* ALL_TOGETHER_AGAIN */
          res++;
        } 
        break;


/****************************************************************************/
/*                                                          REMAINING STUFF */

/*--------------------------------------------------------------------------*/
      case take_stock_op:                                  /* take_stock_op */
        size = get_locint_f();
        res  = get_locint_f();
        d    = get_val_v_f(size);

        for (ls=0;ls<size;ls++)
        { T0[res]=*d;
#if !defined(_ZOS_) /* BREAK_ZOS */           
          ASSIGN_T(Tres,T[res])
          
          FOR_0_LE_l_LT_pk
            TRES_INC = 0;

#endif /* ALL_TOGETHER_AGAIN */
          res++; d++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case death_not:                                          /* death_not */
        arg1=get_locint_f();
        arg2=get_locint_f();

#ifdef _KEEP_
        if (keep) 
        { do 
            IF_KEEP_WRITE_TAYLOR(arg2,keep,k,p)
          while(arg1 < arg2-- ); 
        }
#endif
  	break;

/*--------------------------------------------------------------------------*/
      default:                                                   /* default */
	/* Die here, we screwed up */

        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " ("
                       __FILE__
                       ") : no such operation %d\n", operation);
	exit(-1);
	break;
	
      } /* endswitch */
      
      /* Read the next operation */
      operation=get_op_f();
    }  /* endwhile */

#if defined(_KEEP_)
  if (keep) taylor_close(taylbuf,depcheck,indcheck);
#endif

  end_sweep();
  return ret_c;
}

/****************************************************************************/

END_C_DECLS
