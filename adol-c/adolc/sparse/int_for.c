/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     sparse/int_for.mc
 Revision: $Id: int_for.c,v 1.2 2004/10/14 13:29:48 e_arnold Exp $
 Contents: int_forward_tight,
               (first-order-vector forward mode for bit patterns, 
                tight version = basepoint check,  more precize)
           int_forward_safe
               (first-order-vector forward mode for bit patterns, 
                safe version, no basepoint check, return always 3 )

           Attention : no proper subscript handling in the safe version !
 
 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History: 20040414 kowarz:  adapted to configure - make - make install
          20000310 olvo:    (1) better error messages for discontinuity
                            (2) removed trigraphs stuff
          20000214 olvo:    new op_codes 'eq_plus_prod' and 'eq_min_prod'
                            for  y += x1 * x2
                            and  y -= x1 * x2  
          19990308 christo: bit patterns : 
                            unsigned int -> unsigned long int  
          19981130 olvo:    last check (includes etc.)
          19981113 christo: valuepoint may be NULL, check
          19980930 olvo:    (1) new: locint checkSize; int flag;
                                     unsigned int *TresOP 
                            (2) new: static unsigned int Ttemp;
                            (3) new macros TQO*
                            (4) changed strategy to allow reflexive
                                - eq_mult_av_a, dot_av_av, mult_a_av, 
                                div_av_a
          19980925 olvo:    allow reflexive operations for
                            - cond_assign, cond_assign_s
          19980924 olvo:    deleted all int_* opcodes
          19980923 olvo:    updated all changes from uni5_for.c
                            from 980915  til 980923
                            (sin_op, cos_op, min_op, abs_op)
          19980917 olvo/christo: macros  
          19980915 christo: new comments
          19980820 olvo:    new comparison strategy 

----------------------------------------------------------------------------*/

/****************************************************************************/
/*                                                                   MACROS */

/*--------------------------------------------------------------------------*/
#ifdef _INT_FOR_TIGHT_
#define GENERATED_FILENAME "int_forward_t"
#define _TIGHT_
#undef _SAFE_

/*--------------------------------------------------------------------------*/
#elif _INT_FOR_SAFE_
#define GENERATED_FILENAME "int_forward_s"  
#undef _TIGHT_
#define _SAFE_

#else
#error Error !  Define [_INT_FOR_SAFE_ | _INT_FOR_TIGHT_] 
#endif
/*--------------------------------------------------------------------------*/

#define ARGUMENT(indexi,l,i) argument[indexi][l]
#define TAYLORS(indexd,l,i)   taylors[indexd][l]

/*--------------------------------------------------------------------------*/
/*                                                      access to variables */

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

#define TRES_FOINC   *Tres++
#define TARG_FOINC   *Targ++
#define TARG1_FOINC  *Targ1++
#define TARG2_FOINC  *Targ2++
#define TQO_FOINC    *Tqo++

#define TRES_FODEC   *Tres--
#define TARG_FODEC   *Targ--
#define TARG1_FODEC  *Targ1--
#define TARG2_FODEC  *Targ2--
#define TQO_FODEC    *Tqo--


#define ASSIGN_T(a,b)  a = b;

/*--------------------------------------------------------------------------*/
/*                                                               loop stuff */
#define FOR_0_LE_l_LT_p for (l=0; l<p; l++)  
#define FOR_p_GT_l_GE_0 for (l=p-1; l>=0; l--)  


#define FOR_0_LE_l_LT_pk for (l=0; l<p; l++)  
#define INC_pk_1(T)      T += p-1;
#define VEC_INC(T,inc)   T++;  


/*--------------------------------------------------------------------------*/
/*                                                             other macros */
#define FMIN(x,y)  ((y<x)?y:x)

/* END Macros */


/****************************************************************************/
/*                                                       NECESSARY INCLUDES */
#include "../sparse/jacutils.h"
#include "../oplate.h"
#include "../adalloc.h"
#include "../taputil.h"
#include "../tayutil.h"
#include "../taputil_p.h"

#include <malloc.h>
#include <math.h>

BEGIN_C_DECLS

/****************************************************************************/
/*                                                             NOW THE CODE */

/*--------------------------------------------------------------------------*/
/*                                                   Local Static Variables */
static short tag;

static int for_location_cnt;
static int dep_cnt;
static int ind_cnt;


#ifdef _INT_FOR_TIGHT_
/****************************************************************************/
/* First Order Vector version of the forward mode for bit patterns, tight   */
/****************************************************************************/
int int_forward_tight(
	short             tnum,     /* tape id                              */
        int               depcheck, /* consistency chk on # of dependents   */
        int               indcheck, /* consistency chk on # of independents */
        int               p,        /* # of taylor series, bit pattern      */
        double            *basepoint,  /* independent variable values   (in)*/
        unsigned long int **argument,  /* Taylor coeff.                 (in)*/
        double            *valuepoint, /* dependent variable values    (out)*/
        unsigned long int **taylors)   /* matrix of coefficient vectors(out)*/

/* int_forward_tight( tag, m, n, p, x[n], X[n][p], y[m], Y[m][p]),           
   
     nBV = number of Boolean Vectors to be packed
                      (see Chapter Dependence Analysis, ADOL-C Documentation)
     bits_per_long = 8*sizeof(unsigned long int)
     p = nBV / bits_per_long + ( (nBV % bits_per_long) != 0 )

     The order of the indices in argument and taylors is [var][taylor]
 
     For the full Jacobian matrix set 
     p = indep / bits_per_long + ((indep % bits_per_long) != 0)
     and pass a bit pattern version of the identity matrix as an argument   */


#elif _INT_FOR_SAFE_
/****************************************************************************/
/* First Order Vector version of the forward mode, bit pattern, safe        */
/****************************************************************************/
int int_forward_safe(
	short             tnum,     /* tape id                              */
        int               depcheck, /* consistency chk on # of dependents   */
        int               indcheck, /* consistency chk on # of independents */
        int               p,        /* # of taylor series, bit pattern      */
        unsigned long int **argument, /* Taylor coeff.                  (in)*/
        unsigned long int **taylors)  /* matrix of coefficient vectors (out)*/

/* int_forward_tight( tag, m, n, p, X[n][p], Y[m][p]),           
   
     nBV = number of Boolean Vectors to be packed
                      (see Chapter Dependence Analysis, ADOL-C Documentation)
     bits_per_long = 8*sizeof(unsigned long int)
     p = nBV / bits_per_long + ( (nBV % bits_per_long) != 0 )

     The order of the indices in argument and taylors is [var][taylor]
 
     For the full Jacobian matrix set 
     p = indep / bits_per_long + ((indep % bits_per_long) != 0)
     and pass a bit pattern version of the identity matrix as an argument    */

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

#ifdef _TIGHT_
  double coval = 0, *d = 0;
#endif /* _TIGHT_ */

  int indexi = 0,  indexd = 0;

  /* loop indices */ 
  int l=0, ls;

  /* other necessary variables */
#ifdef _TIGHT_
  double y; 
#endif
  int buffer, flag;
  static int fax, pax;
  
  /* Taylor stuff */
#ifdef _TIGHT_
  static double  *T0;
#endif /* _TIGHT_ */

  static unsigned long int  **T;
  static unsigned long int  *Ttemp;

  unsigned long int         *Tres, *Targ, *Targ1, *Targ2, *Tqo;
  unsigned long int         *TresOP;
#ifdef _TIGHT_
  unsigned long int         *Targ1OP, *Targ2OP;
  unsigned long int         T0temp;
#endif
#define T0res  T0temp


#ifdef DEBUG
/****************************************************************************/
/*                                                           DEBUG MESSAGES */
  fprintf(DIAG_OUT,"Call of %s(..) with tag: %d, n: %d, m %d,\n",
                   GENERATED_FILENAME, tnum, indcheck, depcheck);      
  fprintf(DIAG_OUT,"                    p: %d\n\n",p);
  
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

/*--------------------------------------------------------------------------*/

  if (p compsize pax || for_location_cnt compsize fax)
  { if ( pax || fax )
    { 
#ifdef _TIGHT_
      free((char*) T0);
#endif /* _TIGHT_ */
      free((char*) *T); free((char*) T);
      free((char*) Ttemp);
    }
#ifdef _TIGHT_
    T0= myalloc1(for_location_cnt);
#endif /* _TIGHT_ */
    T     = myalloc2_ulong(for_location_cnt,p);
    Ttemp = myalloc1_ulong(p);
    pax = p;
    fax = for_location_cnt;
  }


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
#ifdef _TIGHT_
        get_val_block_f();
#endif /* _TIGHT_ */
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

#ifdef _TIGHT_
        if (T0[arg] != 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator eq_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          end_sweep();
          return (-1);
        }
        ret_c = 0;
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case neq_zero:                                            /* neq_zero */
        arg = get_locint_f();

#ifdef _TIGHT_
        if (T0[arg] == 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator neq_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          end_sweep();
          return (-1);
        }
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case le_zero:                                              /* le_zero */
        arg = get_locint_f();

#ifdef _TIGHT_
        if (T0[arg] > 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator le_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          end_sweep();
          return (-1);
        }
        if (T0[arg] == 0)
          ret_c = 0;
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case gt_zero:                                              /* gt_zero */
        arg = get_locint_f();

#ifdef _TIGHT_
        if (T0[arg] <= 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator gt_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          end_sweep();
          return (-1);
        }
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case ge_zero:                                              /* ge_zero */
        arg = get_locint_f();

#ifdef _TIGHT_
        if (T0[arg] < 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator ge_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          end_sweep();
          return (-1);
        }
        if (T0[arg] == 0)
          ret_c = 0;
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case lt_zero:                                              /* lt_zero */
        arg = get_locint_f();

#ifdef _TIGHT_
        if (T0[arg] >= 0)
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: Branch switch detected in comparison "
                  "(operator lt_zero).\n"
                  "Forward sweep aborted! Retaping recommended!\n");
          end_sweep();
          return (-1);
        }
#endif /* _TIGHT_ */
        break;


/****************************************************************************/
/*                                                              ASSIGNMENTS */

/*--------------------------------------------------------------------------*/
      case assign_a:           /* assign an adouble variable an    assign_a */
	                       /* adouble value. (=) */
        arg = get_locint_f();
        res = get_locint_f();

#ifdef _TIGHT_
        T0[res] = T0[arg];
#endif /* _TIGHT_ */

	ASSIGN_T(Targ,T[arg])
        ASSIGN_T(Tres,T[res])        

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG_INC;

	break;

/*--------------------------------------------------------------------------*/
      case assign_d:            /* assign an adouble variable a    assign_d */
	                        /* double value. (=) */
        res   = get_locint_f();

#ifdef _TIGHT_
        coval = get_val_f();

	T0[res] = coval;
#endif /* _TIGHT_ */

        ASSIGN_T(Tres, T[res])

        FOR_0_LE_l_LT_pk
	  TRES_INC = 0;


	break;

/*--------------------------------------------------------------------------*/
      case assign_d_zero:  /* assign an adouble variable a    assign_d_zero */
	                   /* double value. (0) (=) */
        res   = get_locint_f();

#ifdef _TIGHT_
	T0[res] = 0.0;
#endif /* _TIGHT_ */

        ASSIGN_T(Tres, T[res])

        FOR_0_LE_l_LT_pk
	  TRES_INC = 0;


	break;

/*--------------------------------------------------------------------------*/
      case assign_d_one:    /* assign an adouble variable a    assign_d_one */
	                    /* double value. (1) (=) */
        res   = get_locint_f();

#ifdef _TIGHT_
	T0[res] = 1.0;
#endif /* _TIGHT_ */

        ASSIGN_T(Tres, T[res])

        FOR_0_LE_l_LT_pk
	  TRES_INC = 0;


	break;

/*--------------------------------------------------------------------------*/
      case assign_ind:       /* assign an adouble variable an    assign_ind */
	                     /* independent double value (<<=) */
        res = get_locint_f();

#ifdef _TIGHT_
	T0[res] = basepoint[indexi];
#endif /* _TIGHT_ */

	ASSIGN_T(Tres, T[res])

        FOR_0_LE_l_LT_p
  	    TRES_INC = ARGUMENT(indexi,l,i);

	++indexi;
	break;

/*--------------------------------------------------------------------------*/
      case assign_dep:           /* assign a float variable a    assign_dep */
	                         /* dependent adouble value. (>>=) */
	res = get_locint_f();

#ifdef _TIGHT_
        if ( valuepoint != NULL )
          valuepoint[indexd] = T0[res];
#endif /* _TIGHT_ */

	ASSIGN_T(Tres, T[res])

        if (taylors != 0 )  /* ??? question: why here? */
          FOR_0_LE_l_LT_p
	      TAYLORS(indexd,l,i) = TRES_INC;

	indexd++;
	break;


/****************************************************************************/
/*                                                   OPERATION + ASSIGNMENT */

/*--------------------------------------------------------------------------*/
      case eq_plus_d:            /* Add a floating point to an    eq_plus_d */
	                         /* adouble. (+=) */
        res   = get_locint_f();
#ifdef _TIGHT_
        coval = get_val_f();

        T0[res] += coval;
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case eq_plus_a:             /* Add an adouble to another    eq_plus_a */
	                          /* adouble. (+=) */
        arg = get_locint_f();
        res = get_locint_f();

#ifdef _TIGHT_
        T0[res] += T0[arg];
#endif /* _TIGHT_ */

        ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC |= TARG_INC;

	break;

/*--------------------------------------------------------------------------*/
      case eq_min_d:       /* Subtract a floating point from an    eq_min_d */
                           /* adouble. (-=) */
        res = get_locint_f();
#ifdef _TIGHT_
        coval = get_val_f();

	T0[res] -= coval;
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case eq_min_a:        /* Subtract an adouble from another    eq_min_a */
	                    /* adouble. (-=) */
        arg = get_locint_f();
        res = get_locint_f();

#ifdef _TIGHT_
	T0[res] -= T0[arg];
#endif /* _TIGHT_ */

	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])
        
        FOR_0_LE_l_LT_pk
	  TRES_INC |= TARG_INC;


	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_d:              /* Multiply an adouble by a    eq_mult_d */
	                           /* flaoting point. (*=) */
        res   = get_locint_f();
#ifdef _TIGHT_
        coval = get_val_f();

	T0[res] *= coval;
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_a:       /* Multiply one adouble by another    eq_mult_a */
	                    /* (*=) */
        arg = get_locint_f();
        res = get_locint_f();

	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        INC_pk_1(Tres)
        INC_pk_1(Targ)

        FOR_p_GT_l_GE_0
	  TRES_FODEC |= TARG_DEC; 

#ifdef _TIGHT_	
        T0[res] *= T0[arg];
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case incr_a:                        /* Increment an adouble    incr_a */
        res   = get_locint_f();

#ifdef _TIGHT_
        T0[res]++;
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case decr_a:                        /* Increment an adouble    decr_a */
        res   = get_locint_f();

#ifdef _TIGHT_
        T0[res]--;
#endif /* _TIGHT_ */
	break;


/****************************************************************************/
/*                                                        BINARY OPERATIONS */

/*--------------------------------------------------------------------------*/
      case plus_a_a:                 /* : Add two adoubles. (+)    plus a_a */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] = T0[arg1] + T0[arg2];
#endif /* _TIGHT_ */

	ASSIGN_T(Tres,  T[res])
	ASSIGN_T(Targ1, T[arg1])
	ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG1_INC | TARG2_INC;

	break;

/*--------------------------------------------------------------------------*/
      case plus_d_a:             /* Add an adouble and a double    plus_d_a */
	                         /* (+) */
        arg   = get_locint_f();
        res   = get_locint_f();
#ifdef _TIGHT_
        coval = get_val_f();

        T0[res] = T0[arg] + coval;
#endif /* _TIGHT_ */

	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG_INC;

	break;

/*--------------------------------------------------------------------------*/
      case min_a_a:              /* Subtraction of two adoubles     min_a_a */
	                         /* (-) */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] = T0[arg1] - T0[arg2];
#endif /* _TIGHT_ */

        ASSIGN_T(Tres,  T[res])
	ASSIGN_T(Targ1, T[arg1])
	ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG1_INC | TARG2_INC;

	break;

/*--------------------------------------------------------------------------*/
      case min_d_a:                /* Subtract an adouble from a    min_d_a */
	                           /* double (-) */
        arg =get_locint_f();
        res = get_locint_f();

#ifdef _TIGHT_
        coval = get_val_f();

        T0[res] = coval - T0[arg];
#endif /* _TIGHT_ */

        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG_INC;

	break;

/*--------------------------------------------------------------------------*/
      case mult_a_a:               /* Multiply two adoubles (*)    mult_a_a */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] = T0[arg1] * T0[arg2];
#endif /* _TIGHT_ */

        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_p
          TRES_FOINC = TARG2_INC | TARG1_INC; 
         
	break;

/*--------------------------------------------------------------------------*/
      /* olvo 20000214: new op_code with recomputation */
      case eq_plus_prod:   /* increment a product of           eq_plus_prod */
                           /* two adoubles (*) */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] += T0[arg1] * T0[arg2];
#endif /* _TIGHT_ */

        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_p
          TRES_FOINC |= TARG2_INC | TARG1_INC; 
         
	break;

/*--------------------------------------------------------------------------*/
      /* olvo 20000214: new op_code with recomputation */
      case eq_min_prod:     /* decrement a product of           eq_min_prod */
                            /* two adoubles (*) */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] -= T0[arg1] * T0[arg2];
#endif /* _TIGHT_ */

        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_p
          TRES_FOINC |= TARG2_INC | TARG1_INC; 
         
	break;

/*--------------------------------------------------------------------------*/
      case mult_d_a:         /* Multiply an adouble by a double    mult_d_a */
                             /* (*) */
        arg   = get_locint_f();
        res   = get_locint_f();

#ifdef _TIGHT_
        coval = get_val_f();

        T0[res] = T0[arg] * coval;
#endif /* _TIGHT_ */

        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG_INC;

	break;

/*--------------------------------------------------------------------------*/
      case div_a_a:           /* Divide an adouble by an adouble    div_a_a */
                              /* (/) */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] = T0[arg1] / T0[arg2];
#endif /* _TIGHT_ */

        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

        FOR_0_LE_l_LT_p
	  TRES_FOINC = TARG1_INC | TARG2_FOINC;        

	break;

/*--------------------------------------------------------------------------*/
      case div_d_a:             /* Division double - adouble (/)    div_d_a */
        arg   = get_locint_f();
        res   = get_locint_f();

#ifdef _TIGHT_
        coval = get_val_f();

        T0[res] = coval / T0[arg];
#endif /* _TIGHT_ */

        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_p
	  TRES_FOINC = TARG_FOINC;

	break;


/****************************************************************************/
/*                                                         SIGN  OPERATIONS */

/*--------------------------------------------------------------------------*/
      case pos_sign_a:                                        /* pos_sign_a */
        arg   = get_locint_f();
        res   = get_locint_f();

#ifdef _TIGHT_
        T0[res] = T0[arg];
#endif /* _TIGHT_ */

	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG_INC;

	break;

/*--------------------------------------------------------------------------*/
      case neg_sign_a:                                        /* neg_sign_a */
        arg   = get_locint_f();
        res   = get_locint_f();

#ifdef _TIGHT_
        T0[res] = -T0[arg];
#endif /* _TIGHT_ */

	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_pk
	  TRES_INC = TARG_INC;

	break;


/****************************************************************************/
/*                                                         UNARY OPERATIONS */

/*--------------------------------------------------------------------------*/
      case exp_op:                          /* exponent operation    exp_op */
        arg = get_locint_f();
        res = get_locint_f();

#ifdef _TIGHT_
        T0[res] = exp(T0[arg]);
#endif /* _TIGHT_ */

	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_p
          TRES_FOINC = TARG_FOINC;

	break;

/*--------------------------------------------------------------------------*/
      case sin_op:                              /* sine operation    sin_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
	T0[arg2] = cos(T0[arg1]);
	T0[res]  = sin(T0[arg1]);
#endif /* _TIGHT_ */

	ASSIGN_T(Tres,  T[res]) 
	ASSIGN_T(Targ1, T[arg1])	
        ASSIGN_T(Targ2, T[arg2])
        
        FOR_0_LE_l_LT_p
        { /* olvo 980923 changed order to allow x = sin(x) */
          TARG2_FOINC =  TARG1; 
          TRES_FOINC  =  TARG1_FOINC; 
        }


	break;

/*--------------------------------------------------------------------------*/
      case cos_op:                            /* cosine operation    cos_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[arg2] = sin(T0[arg1]);
        T0[res]  = cos(T0[arg1]);
#endif /* _TIGHT_ */

	ASSIGN_T(Tres,  T[res])
	ASSIGN_T(Targ1, T[arg1])
	ASSIGN_T(Targ2, T[arg2])
 
        FOR_0_LE_l_LT_p
        { /* olvo 980923 changed order to allow x = cos(x) */
          TARG2_FOINC = TARG1;
          TRES_FOINC  = TARG1_FOINC;
        }

	break;

/*--------------------------------------------------------------------------*/
      case atan_op:                                              /* atan_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] = atan(T0[arg1]);
#endif /* _TIGHT_ */

        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])

        FOR_0_LE_l_LT_p
	  TRES_FOINC = TARG1_FOINC;

	break;

/*--------------------------------------------------------------------------*/
      case asin_op:                                              /* asin_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] = asin(T0[arg1]);
#endif /* _TIGHT_ */

        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])

        FOR_0_LE_l_LT_p
          TRES_FOINC = TARG1_FOINC;
	    
	break;

/*--------------------------------------------------------------------------*/
      case acos_op:                                              /* acos_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] = acos(T0[arg1]);
#endif /* _TIGHT_ */

        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])

        FOR_0_LE_l_LT_p
 	  TRES_FOINC = TARG1_FOINC;

	break;

#ifdef ATRIG_ERF

/*--------------------------------------------------------------------------*/
      case asinh_op:                                            /* asinh_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] = asinh(T0[arg1]);
#endif /* _TIGHT_ */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])

        FOR_0_LE_l_LT_p
          TRES_FOINC = TARG1_FOINC;

        break;

/*--------------------------------------------------------------------------*/
       case acosh_op:                                           /* acosh_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] = acosh(T0[arg1]);
#endif /* _TIGHT_ */

        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
  
        FOR_0_LE_l_LT_p
          TRES_FOINC = TARG1_FOINC;

        break;

/*--------------------------------------------------------------------------*/
      case atanh_op:                                            /* atanh_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] = atanh(T0[arg1]);
#endif /* _TIGHT_ */

        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])

        FOR_0_LE_l_LT_p
          TRES_FOINC = TARG1_FOINC;

        break;

/*--------------------------------------------------------------------------*/
      case erf_op:                                                /* erf_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0[res] = erf(T0[arg1]);
#endif /* _TIGHT_ */

        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ1,T[arg1])
        
        FOR_0_LE_l_LT_p
          TRES_FOINC = TARG1_FOINC;

        break;
           
#endif

/*--------------------------------------------------------------------------*/
      case log_op:                                                /* log_op */
        arg = get_locint_f();
        res = get_locint_f();

#ifdef _TIGHT_
        T0[res] = log(T0[arg]);
#endif /* _TIGHT_ */

	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_p
	  TRES_FOINC = TARG_INC;

	break;

/*--------------------------------------------------------------------------*/
      case pow_op:                                                /* pow_op */
        arg   = get_locint_f();
        res   = get_locint_f();

#ifdef _TIGHT_
        coval = get_val_f();

        T0[res] = pow(T0[arg], coval);
#endif /* _TIGHT_ */

	ASSIGN_T(Tres, T[res])
	ASSIGN_T(Targ, T[arg])

        FOR_0_LE_l_LT_p
          TRES_FOINC = TARG_INC;

	break;



/*--------------------------------------------------------------------------*/
      case sqrt_op:                                              /* sqrt_op */
        arg = get_locint_f();
        res = get_locint_f();

#ifdef _TIGHT_
        T0[res] = sqrt(T0[arg]);
#endif /* _TIGHT_ */
        ASSIGN_T(Targ, T[arg])
	ASSIGN_T(Tres, T[res])	

        FOR_0_LE_l_LT_p
          TRES_FOINC = TARG_INC;

	break;

/*--------------------------------------------------------------------------*/
      case gen_quad:                                            /* gen_quad */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        if (get_val_f()!=T0[arg1])
        { fprintf(DIAG_OUT,
                  "ADOL-C Warning: forward sweep aborted; tape invalid!\n");
          end_sweep();
          return -2;
        }

        T0[res] = get_val_f();
#endif /* _TIGHT_ */

        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])

        FOR_0_LE_l_LT_p
	  TRES_FOINC = TARG1_FOINC;

	break;

/*--------------------------------------------------------------------------*/
      case min_op:                                                /* min_op */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();

#ifdef _TIGHT_
        coval = get_val_f();

        /* olvo 980923 changed order to allow reflexive ops */
        if (T0[arg1] > T0[arg2])
        { if (coval)
            mindec(ret_c,2);
        }
        else 
          if (T0[arg1] < T0[arg2])
          { if (!coval)
              mindec(ret_c,2);
          }
          else
            if (arg1 != arg2)
              mindec(ret_c,1);
#endif /* _TIGHT_ */

        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])
        ASSIGN_T(Tres,  T[res])

#ifdef _TIGHT_
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
            if (TARG1 > TARG2)
              Targ = Targ2OP;
            else 
              if (TARG1 < TARG2)
                Targ = Targ1OP;
            Targ1++; Targ2++;
            if (Targ == NULL) /* e.g. both are equal */
              Targ = Targ1OP;
          }
        
          TRES_INC = TARG_INC;
          
          if (Tqo)
            Tqo++;	   
        }

        T0[res] = FMIN(T0[arg1], T0[arg2]);
#endif /* _TIGHT_ */
#ifdef _SAFE_
          TRES_INC = TARG1_INC | TARG2_INC;
#endif /* _SAFE_ */

        break;

/*--------------------------------------------------------------------------*/
      case abs_val:                                              /* abs_val */
        arg   = get_locint_f();
        res   = get_locint_f();

#ifdef _TIGHT_
        coval = get_val_f();

        /* olvo 980923 changed order to allow reflexive ops */

        y = 0.0;
        if (T0[arg] != 0.0) {
          if (T0[arg] < 0.0)
          { if (coval)
              mindec(ret_c,2);
            y = -1.0;
          }
          else
          { if (!coval)
              mindec(ret_c,2);
            y = 1.0;
          }
        }
        
#endif /* _TIGHT_ */
        
        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])

#ifdef _TIGHT_  
        FOR_0_LE_l_LT_p
        { if ((y == 0.0) && (TARG != 0.0))
            mindec(ret_c,1);
           
          TRES_INC = TARG_INC;
        }

        T0[res] = fabs(T0[arg]);
#endif /* _TIGHT_ */
#ifdef _SAFE_
        FOR_0_LE_l_LT_p
            TRES_INC = TARG_INC;
#endif /* _SAFE_ */
        break;

/*--------------------------------------------------------------------------*/
      case ceil_op:                                              /* ceil_op */
        arg   = get_locint_f();
        res   = get_locint_f();
#ifdef _TIGHT_
        coval = get_val_f();

        T0[res]=ceil(T0[arg]);

        if ( coval != T0[res] )
          mindec(ret_c,2);
#endif /* _TIGHT_ */        
        ASSIGN_T(Tres, T[res])
   
        FOR_0_LE_l_LT_pk
          TRES_INC = 0;

        break;

/*--------------------------------------------------------------------------*/
      case floor_op:                 /* Compute ceil of adouble    floor_op */
        arg   = get_locint_f();
        res   = get_locint_f();
#ifdef _TIGHT_
        coval = get_val_f();

        T0[res] = floor(T0[arg]);

        if ( coval != T0[res] )
          mindec(ret_c,2);
#endif /* _TIGHT_ */
        
        ASSIGN_T(Tres, T[res])
          
        FOR_0_LE_l_LT_pk
          TRES_INC = 0;

        break;


/****************************************************************************/
/*                                                             CONDITIONALS */

/*--------------------------------------------------------------------------*/
      case cond_assign:                                      /* cond_assign */
        arg   = get_locint_f();
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();

        /* olvo 980925 changed order to allow reflexive ops */
        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])
        ASSIGN_T(Targ2, T[arg2])

#ifdef _TIGHT_
        coval = get_val_f();

        if (T0[arg] > 0)
          FOR_0_LE_l_LT_pk
            TRES_INC = TARG1_INC;
        else
          FOR_0_LE_l_LT_pk
            TRES_INC = TARG2_INC;

        if (T0[arg] > 0)
        { if (coval <= 0.0)
            mindec(ret_c,2);
          T0[res] = T0[arg1];
        }
        else
        { if (coval > 0.0)
            mindec(ret_c,2);
          if (T0[arg] == 0)
            mindec(ret_c,0);
          T0[res] = T0[arg2];
        }
#endif /* _TIGHT_ */

#ifdef _SAFE_
        FOR_0_LE_l_LT_pk
          TRES_INC = TARG1_INC | TARG2_INC;
#endif /* _SAFE_ */

        break;

/*--------------------------------------------------------------------------*/
      case cond_assign_s:                                  /* cond_assign_s */
        arg   = get_locint_f();
        arg1  = get_locint_f();
        res   = get_locint_f();

        ASSIGN_T(Tres,  T[res])
        ASSIGN_T(Targ1, T[arg1])

        /* olvo 980925 changed order to allow reflexive ops */
#ifdef _TIGHT_ 
        coval = get_val_f();

        if (T0[arg] > 0)
#endif /* _TIGHT_ */
          FOR_0_LE_l_LT_pk
            TRES_INC = TARG1_INC;

#ifdef _TIGHT_
        if (T0[arg] > 0)
	{ if (coval <= 0.0)
            mindec(ret_c,2);
          T0[res] = T0[arg1];
        } 
        else  
          if (T0[arg] == 0)
            mindec(ret_c,0); 
#endif /* _TIGHT_ */

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
#ifdef _TIGHT_
          T0[res] = T0[arg];
#endif /* _TIGHT_ */

          ASSIGN_T(Targ, T[arg])
          ASSIGN_T(Tres, T[res])
            
          FOR_0_LE_l_LT_pk
	    TRES_INC = TARG_INC;

          res++; arg++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case assign_dv:                                          /* assign_dv */
        size  = get_locint_f();
        res   = get_locint_f();
#ifdef _TIGHT_
        d     = get_val_v_f(size);
#endif /* _TIGHT_ */

        for (ls=0; ls<size; ls++)
        { /* code for assign_d */
#ifdef _TIGHT_
          T0[res] = *d++;
#endif /* _TIGHT_ */

          ASSIGN_T(Tres, T[res])
   
          FOR_0_LE_l_LT_pk
	    TRES_INC = 0;

          res++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case assign_indvec:                                  /* assign_indvec */
        size = get_locint_f();
        res  = get_locint_f();

        for (ls=0; ls<size; ls++)
        { /* code for assign_ind */
#ifdef _TIGHT_
          T0[res] = basepoint[indexi];
#endif /* _TIGHT_ */

          ASSIGN_T(Tres, T[res])
   
          FOR_0_LE_l_LT_p
            TRES_INC = ARGUMENT(indexi,l,i);

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
#ifdef _TIGHT_
          if ( valuepoint != NULL )
            valuepoint[indexd] = T0[res];
#endif /* _TIGHT_ */

          ASSIGN_T(Tres, T[res])

          if (taylors != 0 )  
            FOR_0_LE_l_LT_p
              TAYLORS(indexd,l,i) = TRES_INC;

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

#ifdef _TIGHT_
          T0[res] += T0[arg];
#endif /* _TIGHT_ */

          ASSIGN_T(Tres, T[res])
          ASSIGN_T(Targ, T[arg])
            
          FOR_0_LE_l_LT_pk
	    TRES_INC |= TARG_INC;

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

#ifdef _TIGHT_
          T0[res] -= T0[arg]; 
#endif /* _TIGHT_ */

          ASSIGN_T(Tres, T[res])
          ASSIGN_T(Targ, T[arg])
            
          FOR_0_LE_l_LT_pk
	    TRES_INC |= TARG_INC;

          res++; arg++;
        } 
        break;

/*--------------------------------------------------------------------------*/
      case eq_mult_av_d:                                    /* eq_mult_av_d */
        size  = get_locint_f();
        res   = get_locint_f();

#ifdef _TIGHT_
        coval = get_val_f();

        for (ls=0; ls<size; ls++)
        { /* code for eq_mult_d*/

          T0[res] *= coval;

          res++;
        } 
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case eq_mult_av_a:                                    /* eq_mult_av_a */
        arg  = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();

        /* olvo 980930 new strategy to check for overwrites 
           (changes computation order) */
        checkSize = res+size;  

        for (ls=0; ls<size; ls++)
        { if (res == arg) /* skip res==arg first */
            res++;
          if (res == checkSize) /* checks if arg==res was skipped */
            res = arg; 

          /* code for eq_mult_a*/

          ASSIGN_T(Tres, T[res])
          ASSIGN_T(Targ, T[arg])
            
          INC_pk_1(Tres)
          INC_pk_1(Targ)

          FOR_p_GT_l_GE_0
            TRES_FODEC |= TARG_DEC; 

#ifdef _TIGHT_
          T0[res] *= T0[arg];
#endif /* _TIGHT_ */

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

#ifdef _TIGHT_
          T0[res] = T0[arg1] + T0[arg2];
#endif /* _TIGHT_ */

          ASSIGN_T(Tres,  T[res])
          ASSIGN_T(Targ1, T[arg1])
          ASSIGN_T(Targ2, T[arg2])
            
          FOR_0_LE_l_LT_pk
	    TRES_INC = TARG1_INC | TARG2_INC;

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

#ifdef _TIGHT_
          T0[res] = T0[arg1] - T0[arg2];
#endif /* _TIGHT_ */

          ASSIGN_T(Tres,  T[res])
          ASSIGN_T(Targ1, T[arg1])
          ASSIGN_T(Targ2, T[arg2])
            
          FOR_0_LE_l_LT_pk
	    TRES_INC = TARG1_INC | TARG2_INC;

          res++; arg1++; arg2++;
        } 
        break;


/*--------------------------------------------------------------------------*/
      case dot_av_av:                                          /* dot_av_av */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        T0res = 0;
#endif /* _TIGHT_ */

        /* olvo 980930 check for overwrites -- if necessary use 
           tempories for res-stuff to allow reflexive ops */
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

        for (ls=0; ls<size; ls++)
        { /* code for mult_a_a  */

          Tres = TresOP;
          ASSIGN_T(Targ1, T[arg1])
          ASSIGN_T(Targ2, T[arg2])
           
          /* olvo 980915 now in reverse order to allow x = x*x etc. */
          INC_pk_1(Tres)
          INC_pk_1(Targ1)
          INC_pk_1(Targ2)

          FOR_p_GT_l_GE_0
	    TRES_FODEC |= TARG2_DEC | TARG1_DEC; 

#ifdef _TIGHT_
          T0res += T0[arg1] * T0[arg2];
#endif /* _TIGHT_ */
          arg1++; arg2++;
        }

        /* copy results if necessary */
#ifdef _TIGHT_
        T0[res] = T0res;
#endif /* _TIGHT_ */
        if (flag)
	{ ASSIGN_T(Tres,T[res]) 
          Tqo = TresOP;
  
          FOR_0_LE_l_LT_pk
            TRES_INC = TQO_INC;
	}
        break;

/*--------------------------------------------------------------------------*/
      case mult_a_av:                                          /* mult_a_av */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();

        /* olvo 980930 new strategy to check for overwrites 
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

          ASSIGN_T(Tres,  T[res])
          ASSIGN_T(Targ1, T[arg1])
          ASSIGN_T(Targ2, T[arg2])
            
          /* olvo 980915 now in reverse order to allow x = x*x etc. */
          INC_pk_1(Tres)
          INC_pk_1(Targ1)
          INC_pk_1(Targ2)

          FOR_p_GT_l_GE_0
	    TRES_FOINC = TARG2_INC | TARG1_INC; 

#ifdef _TIGHT_
          T0[res] = T0[arg1] * T0[arg2];
#endif /* _TIGHT_ */

          res++; arg1++;
        }
        break;

/*--------------------------------------------------------------------------*/
      case mult_d_av:                                          /* mult_d_av */
        arg   = get_locint_f();
        size  = get_locint_f();
        res   = get_locint_f();
#ifdef _TIGHT_
        coval = get_val_f();

        for (ls=0; ls<size; ls++)
        { /* code for mult_d_a */
    
          T0[res] = T0[arg] * coval;

          res++; arg++;
        } 
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case div_av_a:                                            /* div_av_a */
        arg1   = get_locint_f();
        arg2   = get_locint_f();
        size   = get_locint_f();
        res    = get_locint_f();

        /* olvo 980930 new strategy to check for overwrites 
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

#ifdef _TIGHT_
          T0[res] = T0[arg1] / T0[arg2];
#endif /* _TIGHT_ */

          ASSIGN_T(Tres,  T[res])
          ASSIGN_T(Targ1, T[arg1])
          ASSIGN_T(Targ2, T[arg2])

          FOR_0_LE_l_LT_p
	    TRES_FOINC = TARG1_INC | TARG2_FOINC;

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

#ifdef _TIGHT_
        arg  = arg2 + (int)(T0[arg1]);

        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          mindec(ret_c,2);

        T0[res] = T0[arg];

        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])
        
        FOR_0_LE_l_LT_pk
          TRES_INC = TARG_INC;
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern forward mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
        break;

/*--------------------------------------------------------------------------*/
      case subscript_l:                                      /* subscript_l */
        arg2 = get_locint_f();  /* Base */
	arg1 = get_locint_f();
        arg  = get_locint_f();         

#ifdef _TIGHT_
        res  = arg2 + (int)(T0[arg1]);

        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          mindec(ret_c,2);


        T0[res] = T0[arg];

        ASSIGN_T(Tres, T[res])
        ASSIGN_T(Targ, T[arg])
        
        FOR_0_LE_l_LT_pk
          TRES_INC = TARG_INC;
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern forward mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
        break;

/*--------------------------------------------------------------------------*/
      case subscript_ld:                                    /* subscript_ld */
        arg2 = get_locint_f(); /* Base */
        arg1 = get_locint_f(); /* pointer to variable containing offset */

#ifdef _TIGHT_
        arg  = arg2+(int)(T0[arg1]);

        T0[arg] = get_val_f();

        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          mindec(ret_c,2);

        ASSIGN_T(Targ, T[arg])
        
        FOR_0_LE_l_LT_pk
          TARG_INC = 0;
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern forward mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
        break;

/*--------------------------------------------------------------------------*/
      case m_subscript:                                      /* m_subscript */
        arg2 = get_locint_f(); /* Base */
        arg1 = get_locint_f(); /* pointer to variable containing offset */
        size = get_locint_f();
        res  = get_locint_f();

#ifdef _TIGHT_
        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          mindec(ret_c,2);

        arg = arg2 + (int)(T0[arg1])*size;
        for (ls=0; ls<size; ls++)
        { T0[res] = T0[arg];

          ASSIGN_T(Tres, T[res])   
          ASSIGN_T(Targ, T[arg])
            
          FOR_0_LE_l_LT_pk
            TRES_INC = TARG_INC;

          res++; arg++;
        }
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern forward mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
        break;

/*--------------------------------------------------------------------------*/
      case m_subscript_l:                                  /* m_subscript_l */
        arg2 = get_locint_f();  /* Base LHS */
	arg1 = get_locint_f();  /* Offset LHS */
        size = get_locint_f();
        arg  = get_locint_f();  /* RHS */

#ifdef _TIGHT_
        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          mindec(ret_c,2);

        res = arg2 + (int)(T0[arg1])*size;
        for (ls=0; ls<size; ls++)
        { 

          T0[res] = T0[arg];

          ASSIGN_T(Tres, T[res])
          ASSIGN_T(Targ, T[arg])
            
          FOR_0_LE_l_LT_pk
            TRES_INC = TARG_INC;

          res++; arg++;
        } 
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern forward mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
        break;

/*--------------------------------------------------------------------------*/
      case m_subscript_ld:                                /* m_subscript_ld */
        arg2 = get_locint_f(); /* Base */
        arg1 = get_locint_f(); /* pointer to variable containing offset */
        arg  = get_locint_f(); /* offset in the vector itself */
        size = get_locint_f();

#ifdef _TIGHT_
        if ( (int)(T0[arg1]) != (int)(get_val_f()) )
          mindec(ret_c,2);

        d      = get_val_v_f(size);

        res = arg2 + (int)(T0[arg1])*size + arg;
        for (ls=0; ls<size; ls++)
        { 

          T0[res] = d[l];

          ASSIGN_T(Tres, T[res])
            
          FOR_0_LE_l_LT_pk
            TRES_INC = 0;

          res++;
        } 
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern forward mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
        break;


/****************************************************************************/
/*                                                          REMAINING STUFF */

/*--------------------------------------------------------------------------*/
      case take_stock_op:                                  /* take_stock_op */
        size = get_locint_f();
        res  = get_locint_f();
#ifdef _TIGHT_
        d    = get_val_v_f(size);
#endif /* _TIGHT_ */

        for (ls=0;ls<size;ls++)
        { 
#ifdef _TIGHT_
          T0[res] = *d;

#endif /* _TIGHT_ */
           
          ASSIGN_T(Tres,T[res])
          
          FOR_0_LE_l_LT_pk
            TRES_INC = 0;

          res++; 
#ifdef _TIGHT_
          d++;
#endif /* _TIGHT_ */
        } 
        break;

/*--------------------------------------------------------------------------*/
      case death_not:                                          /* death_not */
        arg1=get_locint_f();
        arg2=get_locint_f();

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

  end_sweep();
  return ret_c;
}


/****************************************************************************/

END_C_DECLS
