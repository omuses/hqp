/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     sparse/int_rev.mc
 Revision: $Id: int_rev.c,v 1.2 2004/10/14 13:29:48 e_arnold Exp $
 Contents: int_reverse_tight,
                ( first-order-vector reverse mode for bit patterns,
                  checks all dependences on taylors and real values,
                  more precize) 
           int_reverse_safe,
                ( first-order-vector reverse mode for bit patterns, 
                  return always 3, 
                  no dependences on taylors and real values,
                  faster than tight) 
           Attention : no subscript handling !

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History: 20040414 kowarz:  adaption to configure - make - make install
          20000310 olvo:    removed trigraphs stuff
          20000214 olvo:    new op_codes 'eq_plus_prod' and 'eq_min_prod'
                            for  y += x1 * x2
                            and  y -= x1 * x2  
          19990308 christo: bit patterns : 
                            unsigned int -> unsigned long int 
          19981203 olvo:    untransposing reverse
          19981130 olvo:    last check (includes etc.)
          19980930 olvo:    allow reflexive
                            - eq_mult_av_a, mult_a_av, div_av_a
          19980925 olvo:    (1) allow reflexive operations for
                            - cond_assign, cond_assign_s
                            (2) new macros AQO*
          19980924 olvo:    deleted all int_* opcodes
          19980923 olvo:    updated all changes from fo_rev.c
                            from 980915  til 980923
                            (div_a_a, div_d_a, pow_op, sin_op, cos_op)  
          19980820 olvo:    new comparison strategy
          
----------------------------------------------------------------------------*/

/****************************************************************************/
/*                                                                   MACROS */

/*--------------------------------------------------------------------------*/
#ifdef _INT_REV_TIGHT_
#define GENERATED_FILENAME "int_reverse_t" 
#define _TIGHT_
#undef _SAFE_

/*--------------------------------------------------------------------------*/
#elif _INT_REV_SAFE_
#define GENERATED_FILENAME "int_reverse_s"        
#undef _TIGHT_
#define _SAFE_

#else
#error Error !  Define [_INT_REV_SAFE_ | _INT_REV_TIGHT_] 
#endif
/*--------------------------------------------------------------------------*/

#define RESULTS(indexi,l)  results[l][indexi]
#define LAGRANGE(indexd,l) lagrange[l][indexd] 


/*--------------------------------------------------------------------------*/
/*                                                     access to variables  */

#define ARES       *Ares
#define AARG       *Aarg
#define AARG1      *Aarg1
#define AARG2      *Aarg2
#define AQO        *Aqo

#define ARES_INC   *Ares++
#define AARG_INC   *Aarg++
#define AARG1_INC  *Aarg1++
#define AARG2_INC  *Aarg2++
#define AQO_INC    *Aqo++

#define ARES_INC_O  Ares++
#define AARG_INC_O  Aarg++
#define AARG1_INC_O Aarg1++
#define AARG2_INC_O Aarg2++
#define AQO_INC_O   Aqo++

#define ASSIGN_A(a,b)  a = b;

#define TRES       T[res]
#define TARG       T[arg]
#define TARG1      T[arg1]
#define TARG2      T[arg2]

#define ASSIGN_T(a,b)

/*--------------------------------------------------------------------------*/
/*                                                              loop stuff  */

#define FOR_0_LE_l_LT_p for (l=0; l<p; l++)  
#define FOR_p_GT_l_GE_0 for (l=p-1; l>=0; l--)  

#define FOR_0_LE_l_LT_pk1 for (l=0; l<p; l++)  

/* END Macros */


/****************************************************************************/
/*                                                       NECESSARY INCLUDES */
#include "../sparse/jacutils.h"
#include "../oplate.h"
#include "../adalloc.h"
#include "../taputil.h"
#include "../taputil_p.h"
#include "../tayutil.h"
#include "../tayutil_p.h"

#include <malloc.h>
#include <math.h>

BEGIN_C_DECLS

/****************************************************************************/
/*                                                             NOW THE CODE */

/*--------------------------------------------------------------------------*/
/*                                                   Local Static Variables */
static short tag;

static int rev_location_cnt;
static int dep_cnt;
static int ind_cnt;

/****************************************************************************/
/* First-Order Vector Reverse Pass for bit patterns.                        */
/****************************************************************************/

#ifdef _TIGHT_
int int_reverse_tight(
#endif /* _TIGHT_ */
#ifdef _SAFE_
int int_reverse_safe(
#endif /* _SAFE_ */
          short             tnum,  /* tape id                               */
          int               depen, /* consistency chk on # of deps          */
          int               indep, /* consistency chk on # of indeps        */
          int               nrows, /* # of Jacobian rows being calculated q */
          unsigned long int **lagrange,/* domain weight vector[var][row](in)*/
          unsigned long int **results) /* matrix of coeff. vectors[var][row]*/

/* int_reverse_...( tag, m, n, q, U[q][m], Z[q][n])                         

     nBV = number of Boolean Vectors to be packed
                      (see Chapter Dependence Analysis, ADOL-C Documentation)
     bits_per_long = 8*sizeof(unsigned long int)
     q = nBV / bits_per_long + ( (nBV % bits_per_long) != 0 )

     For the full Jacobian matrix set 
     q = depen / bits_per_long + ((depen % bits_per_long) != 0)
     and pass a bit pattern version of the identity matrix as an argument    */

{
/****************************************************************************/
/*                                                           ALL VARIABLES  */
  unsigned char operation; /* operation code */
  int tape_stats[11];      /* tape stats */
  int ret_c=3;             /* return value */

  locint size = 0;
  locint res  = 0;
  locint arg  = 0;
  locint arg1 = 0;
  locint arg2 = 0;

#ifdef _TIGHT_
  double coval = 0, *d = 0;
#endif /* _TIGHT_ */

  int indexi = 0,  indexd = 0;

  /* loop indices */ 
  int j, l, ls;

  /* other necessary variables */
  int buffer;
  static int rax, pax;
#ifdef _TIGHT_
  int taycheck;
  int numdep,numind;
#endif
  unsigned long int aTmp;

/*--------------------------------------------------------------------------*/
  /* Taylor stuff */

#ifdef _TIGHT_
  static revreal *T;
#endif /* _TIGHT_ */
 
/*--------------------------------------------------------------------------*/
  /* Adjoint stuff */

  static unsigned long int **A;

  static unsigned long int *Atemp;
  unsigned long int        *Ares, *Aarg, *Aarg1, *Aarg2, *Aqo;


/*--------------------------------------------------------------------------*/

  /* Notice: The "q" in the interfaces is reperesented internally by "p" */
  int p = nrows;


#ifdef DEBUG
/****************************************************************************/
/*                                                           DEBUG MESSAGES */
  fprintf(DIAG_OUT,"Call of %s(..) with tag: %d, n: %d, m %d,\n",
                   GENERATED_FILENAME, tnum, indep, depen);      
  fprintf(DIAG_OUT,"                    p: %d\n\n",nrows);
#endif


/****************************************************************************/
/*                                                                    INITs */

  /*------------------------------------------------------------------------*/
  /* Set up stuff for the tape */

  tag = tnum;   /*tag is global which indicates which tape to look at */
        
  tapestats(tag,tape_stats);
  ind_cnt          = tape_stats[0];
  dep_cnt          = tape_stats[1];
  rev_location_cnt = tape_stats[2];
  buffer           =  tape_stats[4];

  set_buf_size(buffer);

  if ((depen != dep_cnt)||(indep != ind_cnt))
  { fprintf(DIAG_OUT,"ADOL-C error: Reverse sweep on tape %d  aborted!\n",tag);
    fprintf(DIAG_OUT,"Number of dependent and/or independent variables "
             "passed to reverse is\ninconsistent with number "
             "recorded on tape %d \n",tag);
    exit (-1);
  }
  
  indexi = ind_cnt - 1;
  indexd = dep_cnt - 1;


/****************************************************************************/
/*                                                  MEMORY ALLOCATION STUFF */

  if (rev_location_cnt compsize  rax || p compsize pax)
  { if (rax || pax)
    { free((char*) Atemp);
#ifdef _TIGHT_
      free((char*) T);
#endif /* _TIGHT_ */
      free((char*) *A); free((char*) A);
    }
    Atemp = myalloc1_ulong(p);
#ifdef _TIGHT_
    T = (revreal *)malloc(sizeof(revreal)*rev_location_cnt);
    if (T == NULL)
    { fprintf(DIAG_OUT,"ADOL-C error: cannot allocate %i bytes!\n",
              sizeof(revreal)*rev_location_cnt);
      exit (-1);
    }
#endif /* _TIGHT_ */
    A = myalloc2_ulong(rev_location_cnt,p);
    rax = rev_location_cnt;
    pax = p;
  }


/****************************************************************************/
/*                                                    TAYLOR INITIALIZATION */
#ifdef _TIGHT_

  taylor_back(tag,T,&numdep,&numind,&taycheck);

  if (taycheck < 0)   
  { fprintf(DIAG_OUT,"\n ADOL-C error: reverse fails because it was not"
                   " preceeded\nby a forward sweep with degree>0, keep=1!\n");
    exit(-2);
  };

  if((numdep != depen)||(numind != indep))
  { fprintf(DIAG_OUT, "\n ADOL-C error: reverse fails on tape %d because the"
                    " number of\nindependent and/or dependent variables"
                    " given to reverse are\ninconsistent with that of the"
                    " internal taylor array.\n",tag);
    exit(-2);
  }

#endif /* _TIGHT_ */  
/****************************************************************************/
/*                                                            REVERSE SWEEP */

  /* Initialize the Reverse Sweep */
  init_rev_sweep(tag);

  operation=get_op_r();
  while (operation != start_of_tape) 
  { /* Switch statement to execute the operations in Reverse */
    switch (operation) {


/****************************************************************************/
/*                                                                  MARKERS */

/*--------------------------------------------------------------------------*/
      case end_of_op:                                          /* end_of_op */
        get_op_block_r();
        operation = get_op_r(); 
        /* Skip next operation, it's another end_of_op */
        break;

/*--------------------------------------------------------------------------*/
      case end_of_int:                                        /* end_of_int */
        get_loc_block_r(); /* Get the next int block */
        break;

/*--------------------------------------------------------------------------*/
      case end_of_val:                                        /* end_of_val */
#ifdef _TIGHT_
        get_val_block_r(); /* Get the next val block */
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case start_of_tape:                                  /* start_of_tape */
      case end_of_tape:                                      /* end_of_tape */
	break;


/****************************************************************************/
/*                                                               COMPARISON */

/*--------------------------------------------------------------------------*/
      case eq_zero  :                                            /* eq_zero */
        arg   = get_locint_r();

#ifdef _TIGHT_
        ret_c = 0;
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case neq_zero :                                           /* neq_zero */
      case gt_zero  :                                            /* gt_zero */
      case lt_zero :                                             /* lt_zero */
        arg   = get_locint_r();
        break;

/*--------------------------------------------------------------------------*/
      case ge_zero :                                             /* ge_zero */
      case le_zero :                                             /* le_zero */
        arg   = get_locint_r();

#ifdef _TIGHT_
        ASSIGN_T( Targ, T[arg])

        if (TARG == 0)
          ret_c = 0;
#endif /* _TIGHT_ */
        break;


/****************************************************************************/
/*                                                              ASSIGNMENTS */

/*--------------------------------------------------------------------------*/
      case assign_a:           /* assign an adouble variable an    assign_a */
	                       /* adouble value. (=) */
        res = get_locint_r();
        arg = get_locint_r();

	ASSIGN_A( Aarg, A[arg])
	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
        { AARG_INC |= ARES;
          ARES_INC = 0;
        }

#ifdef _TIGHT_
	get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case assign_d:            /* assign an adouble variable a    assign_d */
	                        /* double value. (=) */
        res   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();
#endif /* _TIGHT_ */

	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
          ARES_INC = 0;
#ifdef _TIGHT_        
        get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case assign_d_zero:  /* assign an adouble variable a    assign_d_zero */
      case assign_d_one:   /* double value (0 or 1). (=)       assign_d_one */
        res   = get_locint_r();

	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
          ARES_INC = 0;
    
#ifdef _TIGHT_    
        get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case assign_ind:       /* assign an adouble variable an    assign_ind */
	                     /* independent double value (<<=) */
        res = get_locint_r();

	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
          RESULTS(indexi,l) = ARES_INC;
    
#ifdef _TIGHT_    
        get_taylor(res);
#endif /* _TIGHT_ */
	indexi--;
	break;

/*--------------------------------------------------------------------------*/
      case assign_dep:           /* assign a float variable a    assign_dep */
	                         /* dependent adouble value. (>>=) */
        res = get_locint_r();

	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
          ARES_INC = LAGRANGE(indexd,l);

	indexd--;
	break;


/****************************************************************************/
/*                                                   OPERATION + ASSIGNMENT */

/*--------------------------------------------------------------------------*/
      case eq_plus_d:            /* Add a floating point to an    eq_plus_d */
	                         /* adouble. (+=) */
        res   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();

	get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case eq_plus_a:             /* Add an adouble to another    eq_plus_a */
	                          /* adouble. (+=) */
        res = get_locint_r();
        arg = get_locint_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg]);

        FOR_0_LE_l_LT_p
          AARG_INC |= ARES_INC;

#ifdef _TIGHT_
	get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case eq_min_d:       /* Subtract a floating point from an    eq_min_d */
                           /* adouble. (-=) */
        res   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();

	get_taylor(res);
#endif /* _TIGHT_ */
	break;
	
/*--------------------------------------------------------------------------*/
      case eq_min_a:        /* Subtract an adouble from another    eq_min_a */
	                    /* adouble. (-=) */
        res = get_locint_r();
        arg = get_locint_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
          AARG_INC |= ARES_INC;

#ifdef _TIGHT_
	get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_d:              /* Multiply an adouble by a    eq_mult_d */
	                           /* flaoting point. (*=)  */
        res   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();

	get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_a:       /* Multiply one adouble by another    eq_mult_a */
	                    /* (*=)                         */
        res = get_locint_r();
        arg = get_locint_r();
#ifdef _TIGHT_
 	get_taylor(res);
#endif /* _TIGHT_ */
	ASSIGN_A( Ares, A[res])
        ASSIGN_A( Aarg, A[arg])
    
        FOR_0_LE_l_LT_p
	  AARG_INC |= ARES_INC;

        break;
        
/*--------------------------------------------------------------------------*/
      case incr_a:                        /* Increment an adouble    incr_a */
      case decr_a:                        /* Increment an adouble    decr_a */
        res   = get_locint_r();

#ifdef _TIGHT_
	get_taylor(res);
#endif /* _TIGHT_ */
	break;


/****************************************************************************/
/*                                                        BINARY OPERATIONS */

/*--------------------------------------------------------------------------*/
      case plus_a_a:                 /* : Add two adoubles. (+)    plus a_a */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg1, A[arg1])
	ASSIGN_A( Aarg2, A[arg2])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG1_INC |= aTmp;
          AARG2_INC |= aTmp;          
        }
#ifdef _TIGHT_
      	get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case plus_d_a:             /* Add an adouble and a double    plus_d_a */
	                         /* (+) */
        res   = get_locint_r();
        arg   = get_locint_r();

#ifdef _TIGHT_
        coval = get_val_r();
#endif /* _TIGHT_ */

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG_INC |= aTmp;
        }
#ifdef _TIGHT_
      	get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case min_a_a:              /* Subtraction of two adoubles    min_a_a */
	                         /* (-) */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg1, A[arg1])
	ASSIGN_A( Aarg2, A[arg2])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG1_INC |= aTmp;
          AARG2_INC |= aTmp;
        }
#ifdef _TIGHT_
        get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case min_d_a:                /* Subtract an adouble from a    min_d_a */
	                           /* double (-) */
        res   = get_locint_r();
        arg   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();
#endif /* _TIGHT_ */

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG_INC |= aTmp;
        }
#ifdef _TIGHT_
        get_taylor(res);
#endif /* _TIGHT_ */
	break;
	
/*--------------------------------------------------------------------------*/
      case mult_a_a:               /* Multiply two adoubles (*)    mult_a_a */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();
#ifdef _TIGHT_
        get_taylor(res);
#endif /* _TIGHT_ */

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg2, A[arg2])
	ASSIGN_A( Aarg1, A[arg1])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG2_INC |= aTmp;
	  AARG1_INC |= aTmp;
        }
        break;

/*--------------------------------------------------------------------------*/
      /* olvo 20000214: new op_code with recomputation */
      case eq_plus_prod:   /* increment a product of           eq_plus_prod */
                           /* two adoubles (*) */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg2, A[arg2])
	ASSIGN_A( Aarg1, A[arg1])

#ifdef _TIGHT_
        /* Recomputation */
        ASSIGN_T( Targ1, T[arg1])
	ASSIGN_T( Targ2, T[arg2])
        ASSIGN_T( Tres,  T[res])

        TRES -= TARG1*TARG2;
#endif /* _TIGHT_ */

        FOR_0_LE_l_LT_p
        { AARG2_INC |= ARES;
	  AARG1_INC |= ARES_INC;
        }
        break;

/*--------------------------------------------------------------------------*/
      /* olvo 20000214: new op_code with recomputation */
      case eq_min_prod:    /* decrement a product of            eq_min_prod */
                           /* two adoubles (*) */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg2, A[arg2])
	ASSIGN_A( Aarg1, A[arg1])

#ifdef _TIGHT_
        /* Recomputation */
        ASSIGN_T( Targ1, T[arg1])
	ASSIGN_T( Targ2, T[arg2])
        ASSIGN_T( Tres,  T[res])

        TRES += TARG1*TARG2;
#endif /* _TIGHT_ */

        FOR_0_LE_l_LT_p
        { AARG2_INC |= ARES;
	  AARG1_INC |= ARES_INC;
        }
        break;

/*--------------------------------------------------------------------------*/
      case mult_d_a:         /* Multiply an adouble by a double    mult_d_a */
                             /* (*) */
        res   = get_locint_r();
        arg   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();
#endif /* _TIGHT_ */
	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG_INC |= aTmp;
        }

#ifdef _TIGHT_
        get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case div_a_a:           /* Divide an adouble by an adouble    div_a_a */
                              /* (/) */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg2, A[arg2])
	ASSIGN_A( Aarg1, A[arg1])

	/* olvo 980923 to allow x =y/x */
#ifdef _TIGHT_
        get_taylor(res);
#endif /* _TIGHT_ */

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG1_INC |= aTmp;
	  AARG2_INC |= aTmp;
        }
        break;

/*--------------------------------------------------------------------------*/
      case div_d_a:             /* Division double - adouble (/)    div_d_a */
        res   = get_locint_r();
        arg   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();
#endif /* _TIGHT_ */

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

	/* olvo 980923 to allow x =d/x */
#ifdef _TIGHT_
        if (arg == res)
          get_taylor(arg);
#endif /* _TIGHT_ */

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG_INC |= aTmp;
        }
#ifdef _TIGHT_
        get_taylor(res);
#endif /* _TIGHT_ */
        break;


/****************************************************************************/
/*                                                         SIGN  OPERATIONS */

/*--------------------------------------------------------------------------*/
      case pos_sign_a:                                        /* pos_sign_a */
        res   = get_locint_r();
        arg   = get_locint_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG_INC |= aTmp;
        }

#ifdef _TIGHT_
      	get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case neg_sign_a:                                        /* neg_sign_a */
        res   = get_locint_r();
        arg   = get_locint_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG_INC |= aTmp;
        }
#ifdef _TIGHT_
      	get_taylor(res);
#endif /* _TIGHT_ */
	break;


/****************************************************************************/
/*                                                         UNARY OPERATIONS */

/*--------------------------------------------------------------------------*/
      case exp_op:                          /* exponent operation    exp_op */
        res = get_locint_r();
        arg = get_locint_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG_INC |= aTmp;
        }
#ifdef _TIGHT_
 	get_taylor(res);
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case sin_op:                              /* sine operation    sin_op */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg1, A[arg1])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG1_INC |= aTmp;
        }

#ifdef _TIGHT_
        /* olvo 980923 changed order to allow x = sin(x) */
       	get_taylor(res);
       	get_taylor(arg2); /* olvo 980710 covalue */
	                  /* NOTE: A[arg2] should be 0 already */
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case cos_op:                            /* cosine operation    cos_op */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg1, A[arg1])
    
        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG1_INC |= aTmp;          
        }

#ifdef _TIGHT_
        /* olvo 980923 changed order to allow x = cos(x) */
       	get_taylor(res);
       	get_taylor(arg2); /* olvo 980710 covalue */
	                  /* NOTE A[arg2] should be 0 already */
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case atan_op:                                             /* atan_op  */
      case asin_op:                                             /* asin_op  */
      case acos_op:                                             /* acos_op  */
      case asinh_op:                                            /* asinh_op */
      case acosh_op:                                            /* acosh_op */
      case atanh_op:                                            /* atanh_op */
      case erf_op:                                              /* erf_op   */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();
#ifdef _TIGHT_
        get_taylor(res);
#endif /* _TIGHT_ */
        ASSIGN_A( Ares,  A[res])
        ASSIGN_A( Aarg1, A[arg1])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG1_INC |= aTmp;
        }
	break;

/*--------------------------------------------------------------------------*/
      case log_op:                                                /* log_op */
        res = get_locint_r();
        arg = get_locint_r();
#ifdef _TIGHT_
        get_taylor(res);
#endif /* _TIGHT_ */

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG_INC |= aTmp;
        }
	break;

/*--------------------------------------------------------------------------*/
      case pow_op:                                                /* pow_op */
        res   = get_locint_r();
        arg   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();
#endif /* _TIGHT_ */
	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

	/* olvo 980923 to allow x = pow(x,d) */
#ifdef _TIGHT_
        if (arg == res)
          get_taylor(arg);
#endif /* _TIGHT_ */

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG_INC |= aTmp;
        }
#ifdef _TIGHT_ 
        get_taylor(res);
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case sqrt_op:                                              /* sqrt_op */
        res = get_locint_r();
        arg = get_locint_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG_INC |= aTmp;
        }
#ifdef _TIGHT_
        get_taylor(res);
#endif /* _TIGHT_ */
        break;

/*--------------------------------------------------------------------------*/
      case gen_quad:                                            /* gen_quad */
        res   = get_locint_r();
        arg2  = get_locint_r();
        arg1  = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();
        coval = get_val_r();
#endif /* _TIGHT_ */

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg1, A[arg1])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0; 
          AARG1_INC |= aTmp;
        }
#ifdef _TIGHT_  
      	get_taylor(res);
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case min_op:                                                /* min_op */
        res   = get_locint_r();
        arg2  = get_locint_r();
        arg1  = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();
  
        get_taylor(res);
#endif /* _TIGHT_ */

        ASSIGN_A( Aarg1, A[arg1])
        ASSIGN_A( Aarg2, A[arg2])
        ASSIGN_A( Ares,  A[res])
#ifdef _TIGHT_
        ASSIGN_T( Targ1, T[arg1])
        ASSIGN_T( Targ2, T[arg2])

        if (TARG1 > TARG2)
          FOR_0_LE_l_LT_p
          { aTmp = ARES;
            ARES_INC = 0; 
            if ((coval) && (aTmp))
              mindec(ret_c,2);
            AARG2_INC |= aTmp;
          }
        else 
          if (TARG1 < TARG2)
            FOR_0_LE_l_LT_p
            { aTmp = ARES;
              ARES_INC = 0; 
              if ((!coval) && (aTmp))
                mindec(ret_c,2);
              AARG1_INC |= aTmp;
            }
          else
          { /* both are equal */
            FOR_0_LE_l_LT_p
            { aTmp = ARES;
              ARES_INC = 0; 
              AARG2_INC |= aTmp;
              AARG1_INC |= aTmp;
            }
            if (arg1 != arg2)
              mindec(ret_c,1);
          }
#endif /* _TIGHT_ */
#ifdef _SAFE_
        FOR_0_LE_l_LT_p
          { aTmp = ARES;
            ARES_INC = 0; 
            AARG1_INC |= aTmp;
            AARG2_INC |= aTmp;
          }
#endif /* _SAFE_ */
        break;

/*--------------------------------------------------------------------------*/
      case abs_val:                                              /* abs_val */
        res   = get_locint_r();
        arg   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();
  
        get_taylor(res);
#endif /* _TIGHT_ */

        ASSIGN_A( Ares, A[res])
        ASSIGN_A( Aarg, A[arg])
#ifdef _TIGHT_
        ASSIGN_T( Targ, T[arg])

        if (TARG < 0.0)
          FOR_0_LE_l_LT_p
          { aTmp = ARES;
            ARES_INC = 0; 
            if ((coval) && (aTmp))
              mindec(ret_c,2);
            AARG_INC |= aTmp;
          } 
        else 
          if (TARG > 0.0)
            FOR_0_LE_l_LT_p
            { aTmp = ARES;
              ARES_INC = 0; 
              if ((!coval) && (aTmp))
                mindec(ret_c,2);
              AARG_INC |= aTmp;
            }
          else 
            FOR_0_LE_l_LT_p 
            { aTmp = ARES;
              ARES_INC = 0; 
              if (aTmp)        
                mindec(ret_c,1);
            }
#endif /* _TIGHT_ */
#ifdef _SAFE_
        FOR_0_LE_l_LT_p
          { aTmp = ARES;
            ARES_INC = 0; 
            AARG_INC |= aTmp;
          } 
#endif /* _SAFE_ */
        break;

/*--------------------------------------------------------------------------*/
      case ceil_op:                                              /* ceil_op */
        res   = get_locint_r();
        arg   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();

        get_taylor(res);
#endif /* _TIGHT_ */
        ASSIGN_A( Ares, A[res])
#ifdef _TIGHT_
        ASSIGN_T( Targ, T[arg])

        coval = (coval != ceil(TARG) );
#endif /* _TIGHT_ */

        FOR_0_LE_l_LT_p
        { 
#ifdef _TIGHT_
          if ((coval) && (ARES))
            mindec(ret_c,2);
#endif /* _TIGHT_ */
          ARES_INC = 0;
        }
        break;

/*--------------------------------------------------------------------------*/
      case floor_op:                                            /* floor_op */
        res   = get_locint_r();
        arg   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();

        get_taylor(res);
#endif /* _TIGHT_ */

        ASSIGN_A( Ares, A[res])
        ASSIGN_A( Aarg, A[arg])
#ifdef _TIGHT_
        ASSIGN_T( Targ, T[arg])

        coval = ( coval != floor(TARG1) );
#endif /* _TIGHT_ */
        FOR_0_LE_l_LT_p
        { 
#ifdef _TIGHT_
          if ( (coval) && (ARES) )
            mindec(ret_c,2);
#endif /* _TIGHT_ */
          ARES_INC = 0;
	}
        break;    


/****************************************************************************/
/*                                                             CONDITIONALS */

/*--------------------------------------------------------------------------*/
      case cond_assign:                                      /* cond_assign */
        res    = get_locint_r();
        arg2   = get_locint_r();
        arg1   = get_locint_r();
        arg    = get_locint_r();
#ifdef _TIGHT_
        coval  = get_val_r(); 

	get_taylor(res);
#endif /* _TIGHT_ */

        ASSIGN_A( Aarg1, A[arg1])
	ASSIGN_A( Ares,  A[res])
        ASSIGN_A( Aarg2, A[arg2])
#ifdef _TIGHT_
        ASSIGN_T( Targ,  T[arg])

	/* olvo 980925 changed code a little bit */
        if (TARG > 0.0)
        { if (res != arg1)
            FOR_0_LE_l_LT_p
            { if ((coval <= 0.0) && (ARES))
                mindec(ret_c,2);
              AARG1_INC |= ARES;
              ARES_INC = 0.0;
	    }
          else
            FOR_0_LE_l_LT_p
              if ((coval <= 0.0) && (ARES_INC))
                mindec(ret_c,2);
	}
        else
        { if (res != arg2)
            FOR_0_LE_l_LT_p
            { if ((coval <= 0.0) && (ARES))
                mindec(ret_c,2);
              AARG2_INC |= ARES;
              ARES_INC = 0.0;
	    }
          else
            FOR_0_LE_l_LT_p
              if ((coval <= 0.0) && (ARES_INC))
                mindec(ret_c,2);
	}
#endif /* _TIGHT_ */

#ifdef _SAFE_
         if (res != arg1)
         { FOR_0_LE_l_LT_p
             AARG1_INC |= ARES_INC;
           ASSIGN_A( Ares,  A[res])
	 }
         if (res != arg2)
         { FOR_0_LE_l_LT_p
             AARG2_INC |= ARES_INC;
           ASSIGN_A( Ares,  A[res])
	 }
         if ((res != arg1) && (res != arg2))  
           FOR_0_LE_l_LT_p
             ARES_INC = 0;
#endif /* _SAFE_ */
        break;

/*--------------------------------------------------------------------------*/
      case cond_assign_s:                                  /* cond_assign_s */
        res   = get_locint_r();
        arg1  = get_locint_r();
        arg   = get_locint_r(); 
#ifdef _TIGHT_
        coval = get_val_r();

        get_taylor(res);
#endif /* _TIGHT_ */
        ASSIGN_A( Aarg1, A[arg1])
	ASSIGN_A( Ares,  A[res])
#ifdef _TIGHT_
        ASSIGN_T( Targ,  T[arg])
        
	/* olvo 980924 changed code a little bit */
        if (TARG > 0.0)
        { if (res != arg1)
            FOR_0_LE_l_LT_p
            { if ((coval <= 0.0) && (ARES))
                mindec(ret_c,2);
              AARG1_INC |= ARES;
              ARES_INC = 0.0;
	    }
          else
            FOR_0_LE_l_LT_p
              if ((coval <= 0.0) && (ARES_INC))
                mindec(ret_c,2);
	}
        else
          if (TARG == 0.0) /* we are at the tie */
            FOR_0_LE_l_LT_p
              if (ARES_INC)
                mindec(ret_c,0);
#endif /* _TIGHT_ */

#ifdef _SAFE_ 
        if (res != arg1)
          FOR_0_LE_l_LT_p
          { AARG1 |= ARES;
            ARES_INC = 0;
          }
#endif /* _SAFE_ */
        break;


/****************************************************************************/
/*                                                       VECTOR ASSIGNMENTS */

/*--------------------------------------------------------------------------*/
      case assign_av:                                          /* assign_av */
        res  = get_locint_r();
        size = get_locint_r();
        arg  = get_locint_r();

	res += size;
	arg += size;
        for (ls=size; ls>0; ls--)
	{ res--;             /* Location of left-hand-side  */
	  arg--;             /* Location of right-hand-side */

          /* code for assign_a */
          ASSIGN_A( Aarg, A[arg])
	  ASSIGN_A( Ares, A[res])

          FOR_0_LE_l_LT_p
	  { AARG_INC |= ARES;
	    ARES_INC = 0;
	  }
#ifdef _TIGHT_
	  get_taylor(res);
#endif /* _TIGHT_ */
	}
	break;

/*--------------------------------------------------------------------------*/
      case assign_dv:                                          /* assign_dv */
        res  = get_locint_r();
        size = get_locint_r();
#ifdef _TIGHT_
        d    = get_val_v_r(size);
#endif /* _TIGHT_ */
        res += size;

	for (ls=size; ls>0; ls--)
	{ res--;                /* Location of left-hand-side */     

          /* code for assign_d */
	  ASSIGN_A( Ares, A[res])

	  FOR_0_LE_l_LT_p
	    ARES_INC = 0;
#ifdef _TIGHT_
	  get_taylor(res);
#endif /* _TIGHT_ */
	}
	break;

/*--------------------------------------------------------------------------*/
      case assign_indvec:                                  /* assign_indvec */
        res  = get_locint_r();
        size = get_locint_r();

        res += size;
	for (ls=size; ls>0; ls--)
	{ res--;             /* Location of the left-hand-side */

          /* code for assign_ind */
	  ASSIGN_A( Ares, A[res])

	  FOR_0_LE_l_LT_p
	    RESULTS(indexi,l) = ARES_INC;
	  indexi--;
#ifdef _TIGHT_
	  get_taylor(res);
#endif /* _TIGHT_ */
        }
#ifdef _TIGHT_
	reset_val_r();
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case assign_depvec:                                  /* assign_depvec */
        res  = get_locint_r();
        size = get_locint_r();

        res += size;
	for (ls=size; ls>0; ls--)
	{ res--;             /* Location of the left-hand-side */

	  /* code for assign_dep */
	  ASSIGN_A( Ares, A[res])
            
          FOR_0_LE_l_LT_p
            ARES_INC = LAGRANGE(indexd,l);
	  indexd--;
	}
	break;


/****************************************************************************/
/*                                            VECTOR OPERATION + ASSIGNMENT */

/*--------------------------------------------------------------------------*/
      case eq_plus_av:                                        /* eq_plus_av */
        res  = get_locint_r();
        size = get_locint_r();
        arg  = get_locint_r();

        res += size;
        arg += size;
	for (ls=size; ls>0; ls--)
	{ res--;            /* Location of left-hand-side  */
	  arg--;            /* Location on right-hand-side */

          /* code for eq_plus_a */
          ASSIGN_A( Ares, A[res])
	  ASSIGN_A (Aarg, A[arg])

	  FOR_0_LE_l_LT_p
	    AARG_INC |= ARES_INC;
#ifdef _TIGHT_
          get_taylor(res);
#endif /* _TIGHT_ */
	}
	break;

/*--------------------------------------------------------------------------*/
      case eq_min_av:                                          /* eq_min_av */
        res  = get_locint_r();
        size = get_locint_r();
        arg  = get_locint_r();

        res += size;
        arg += size;
	for (ls=size; ls>0; ls--)
	{ res--;              /* Location of left-hand-side  */
	  arg--;              /* Location on right-hand-side */

	  /* code for eq_min_a */ 
	  ASSIGN_A( Ares, A[res])
	  ASSIGN_A( Aarg, A[arg])

	  FOR_0_LE_l_LT_p
	    AARG_INC |= ARES_INC;
#ifdef _TIGHT_
          get_taylor(res);
#endif /* _TIGHT_ */
	}
	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_av_d:                                    /* eq_mult_av_d */
        res   = get_locint_r();
        size  = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();

        res += size;
	for (ls=size; ls>0; ls--)
        { res--;            /* Location of the left-hand-side  */

	  /* code for eq_mult_d*/

	  get_taylor(res);
        }
#endif /* _TIGHT_ */
	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_av_a:                                    /* eq_mult_av_a */
        res  = get_locint_r();
        size = get_locint_r();
        arg  = get_locint_r();

        /* olvo 980930 new strategy to check for overwrites 
           (changes computation order) */
#ifdef _TIGHT_
        if ((arg >= res) && (arg < res+size)) 
	{ /* FIRST compute the case: res==arg */
          /* simplified code for eq_mult_a*/
	  get_taylor(arg);
     	}
#endif /* _TIGHT_ */

	res += size;
        for (ls=size; ls>0; ls--)
	{ res--;                 /* Location of the left-hand-side  */
	  /* arg    = fixed;        Location on the right-hand-side */
     
          if (res == arg) /* NOW skip this case */
            continue;

          /* code for eq_mult_a*/
#ifdef _TIGHT_
	  get_taylor(res);
#endif /* _TIGHT_ */
	
	  ASSIGN_A( Aarg, A[arg])
	  ASSIGN_A( Ares, A[res]) 

          FOR_0_LE_l_LT_p
            AARG_INC |= ARES_INC;

	}
	break;


/****************************************************************************/
/*                                                 BINARY VECTOR OPERATIONS */

/*--------------------------------------------------------------------------*/
      case plus_av_av:                                        /* plus_av_av */
        res  = get_locint_r();
        size = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	res  += size;
	arg1 += size;
	arg2 += size;
        for (ls=size; ls>0; ls--)
	{ arg2--;       /* Location of var 2  */
	  arg1--;       /* Location of var 1  */
	  res--;        /* Location of result */

          /* code for plus_a_a */
	  ASSIGN_A( Ares,  A[res])
	  ASSIGN_A( Aarg1, A[arg1])
	  ASSIGN_A( Aarg2, A[arg2])

	  FOR_0_LE_l_LT_p
	  { aTmp = ARES;
            ARES_INC = 0; 
            AARG1_INC |= aTmp;
            AARG2_INC |= aTmp;          
          }
#ifdef _TIGHT_
      	  get_taylor(res);
#endif /* _TIGHT_ */
	}
	break;

/*--------------------------------------------------------------------------*/
      case sub_av_av:                                          /* sub_av_av */
        res  = get_locint_r();
        size = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	res  += size;
	arg1 += size;
	arg2 += size;
        for (ls=size; ls>0; ls--)
	{ arg2--;       /* Location of var 2  */
	  arg1--;       /* Location of var 1  */
	  res--;        /* Location of result */

	  /* code for min_a_a */
	  ASSIGN_A( Ares,  A[res])
	  ASSIGN_A( Aarg1, A[arg1])
	  ASSIGN_A( Aarg2, A[arg2])

	  FOR_0_LE_l_LT_p
	  { aTmp = ARES;
            ARES_INC = 0; 
            AARG1_INC |= aTmp;
            AARG2_INC |= aTmp;          
          }
#ifdef _TIGHT_
          get_taylor(res);
#endif /* _TIGHT_ */
	}
	break;

/*--------------------------------------------------------------------------*/
      case dot_av_av:                                          /* dot_av_av */
        res  = get_locint_r();
        size = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();
#ifdef _TIGHT_
        get_taylor(res);
#endif /* _TIGHT_ */

        /* save Ares to Atemp */
        ASSIGN_A( Aqo,  Atemp)
	ASSIGN_A( Ares, A[res])
        FOR_0_LE_l_LT_p
	{ AQO_INC = ARES;
          ARES_INC = 0; 
	}

	for (ls=0; ls<size; ls++)
	{ /* code for mult_a_a  */
	  ASSIGN_A( Aarg2, A[arg2])
	  ASSIGN_A( Aarg1, A[arg1])
          ASSIGN_A( Aqo,  Atemp)

	  FOR_0_LE_l_LT_p
          { AARG2_INC |= AQO;
	    AARG1_INC |= AQO_INC;
          }
          
          arg1++; arg2++;
	}
	break;

/*--------------------------------------------------------------------------*/
      case mult_a_av:                                          /* mult_a_av */
        res  = get_locint_r();
        size = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

        /* olvo 980930 new strategy to check for overwrites 
           (changes computation order) */
        if ((arg2 >= res) && (arg2 < res+size)) 
	{ /* FIRST compute the case: res==arg2 */
	  /* simplified code for mult_a_a */
#ifdef _TIGHT_
    	  get_taylor(res);
#endif /* _TIGHT_ */

	  ASSIGN_A( Aarg2, A[arg2])
	  ASSIGN_A( Aarg1, A[arg1])

          FOR_0_LE_l_LT_p
            AARG1_INC |= AARG2_INC;
	}

        res  += size;
        arg1 += size;
	for (ls=size; ls>0; ls--)
	{ arg1--;    /* Location of rght hnd side vectore[l]  */
	  res--;     /* Location of the result */
	    
          if (res == arg2) /* NOW skip this case */
            continue;

	  /* code for mult_a_a */
#ifdef _TIGHT_
    	  get_taylor(res);
#endif /* _TIGHT_ */

	  ASSIGN_A( Ares,  A[res])
	  ASSIGN_A( Aarg2, A[arg2])
	  ASSIGN_A( Aarg1, A[arg1])

          FOR_0_LE_l_LT_p
          { aTmp = ARES;
            ARES_INC = 0; 
            AARG2_INC |= aTmp;
	    AARG1_INC |= aTmp;
          }
	}
	break;

/*--------------------------------------------------------------------------*/
      case mult_d_av:                                          /* mult_d_av */
        res   = get_locint_r();
        size  = get_locint_r();
        arg   = get_locint_r();
#ifdef _TIGHT_
        coval = get_val_r();
#endif /* _TIGHT_ */

        res += size;
        arg += size;
	for (ls=size; ls>0; ls--)
	{ arg--;     /* Location on the right-hand-side */
	  res--;     /* location of the result */
	  /* coval = Fixed double value */
	    
	  /* code for mult_d_a */
	  ASSIGN_A( Ares, A[res])
	  ASSIGN_A( Aarg, A[arg])

          FOR_0_LE_l_LT_p
          { aTmp = ARES;
            ARES_INC = 0; 
            AARG_INC |= aTmp;
          }
#ifdef _TIGHT_
          get_taylor(res);
#endif /* _TIGHT_ */
        }
	break;

/*--------------------------------------------------------------------------*/
      case div_av_a:                                            /* div_av_a */
        res  = get_locint_r();
        size = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

        /* olvo 980930 new strategy to check for overwrites 
           (changes computation order) */
        if ((arg2 >= res) && (arg2 < res+size)) 
	{ /* FIRST compute the case: res==arg2 */
	  /* simplified code for div_a_a */
	  ASSIGN_A( Aarg2, A[arg2])
	  ASSIGN_A( Aarg1, A[arg1])

          FOR_0_LE_l_LT_p
            AARG1_INC |= AARG2_INC;

#ifdef _TIGHT_
          get_taylor(res);
#endif /* _TIGHT_ */
	}

        res  += size;
        arg1 += size;
	for (ls=size; ls>0; ls--)
	{ arg1--;    /* Location of right-hand-side vector[l] */
	  res--;     /* Location of the result */
	    
          if (res == arg2) /* NOW skip this case */
            continue;

	  /* code for div_a_a */
	  ASSIGN_A( Ares,  A[res])
	  ASSIGN_A( Aarg2, A[arg2])
	  ASSIGN_A( Aarg1, A[arg1])

          FOR_0_LE_l_LT_p
          { aTmp = ARES;
            ARES_INC = 0; 
            AARG1_INC |= aTmp;
	    AARG2_INC |= aTmp;
          }
#ifdef _TIGHT_
          get_taylor(res);
#endif /* _TIGHT_ */
	}
	break;


/****************************************************************************/
/*                                                               SUBSCRIPTS */


/*--------------------------------------------------------------------------*/
      case subscript:                                          /* subscript */
        res   = get_locint_r();
        arg1  = get_locint_r();
        arg2  = get_locint_r(); 
#ifdef _TIGHT_
        coval = get_val_r();
         
        ASSIGN_T( Targ1, T[arg1])

        arg  = arg2 + (int)(TARG1);

	ASSIGN_A( Aarg, A[arg])
	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
        { 
          if (((int)(coval) != (int)(TARG1)) && (ARES))
            mindec(ret_c,2);

          AARG_INC |= ARES;
          if (arg != res)
            ARES_INC = 0;
          else 
            ARES_INC;
        }
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern reverse mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
	break;

/*--------------------------------------------------------------------------*/
      case subscript_l:                                      /* subscript_l */
        arg   = get_locint_r();
        arg1  = get_locint_r();
        arg2  = get_locint_r(); 
#ifdef _TIGHT_
        coval = get_val_r();

        ASSIGN_T( Targ1, T[arg1])

        res   = arg2 + (int)(TARG1);

	get_taylor(res);

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

	FOR_0_LE_l_LT_p
        { 
          if (((int)(coval) != (int)(TARG1)) && (ARES))
            mindec(ret_c,2);

          AARG_INC |= ARES;
          if(arg != res)
            ARES_INC = 0;
          else
            ARES_INC;
        }
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern reverse mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
	break;
      
/*--------------------------------------------------------------------------*/
      case subscript_ld:                                    /* subscript_ld */
        arg1  = get_locint_r();
        arg2  = get_locint_r(); 
#ifdef _TIGHT_
        coval = get_val_r();
        coval = get_val_r();

        ASSIGN_T( Targ1, T[arg1])

        arg = arg2 + (int)(TARG1);

        get_taylor(arg);

        if((int)(coval)!=(int)(TARG1))
          mindec(ret_c,2);

	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
	  AARG_INC = 0;
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern reverse mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
	break;
       
/*--------------------------------------------------------------------------*/
      case m_subscript:                                      /* m_subscript */
        res   = get_locint_r();
        size  = get_locint_r();
        arg1  = get_locint_r();
        arg2  = get_locint_r(); 
#ifdef _TIGHT_
        coval = get_val_r();

        ASSIGN_T( Targ1, T[arg1])

        arg = arg2 + ((int)(TARG1) + 1)*size;

        res += size;
        for (ls=size; ls>0; ls--)
        { res--; arg--;

          ASSIGN_A( Aarg, A[arg])
	  ASSIGN_A( Ares, A[res])

          FOR_0_LE_l_LT_p
          { 
            if (((int)(coval)!=(int)(TARG1)) && (ARES))
              mindec(ret_c,2);

            AARG_INC |= ARES;
            if (arg != res)
              ARES_INC = 0;
            else
              ARES_INC;
          }
        }
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern reverse mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
	break;

/*--------------------------------------------------------------------------*/
      case m_subscript_l:                                  /* m_subscript_l */
        arg   = get_locint_r();
        size  = get_locint_r();
        arg1  = get_locint_r();
        arg2  = get_locint_r();
#ifdef _TIGHT_ 
        coval = get_val_r();

        ASSIGN_T( Targ1, T[arg1])
  
        res = arg2 + ((int)(TARG1) + 1)*size;
        arg += size;
        for (ls=size; ls>0; ls--)
        { arg--; res--;

          get_taylor(res);

          ASSIGN_A( Aarg, A[arg])
          ASSIGN_A( Ares, A[res])

          FOR_0_LE_l_LT_p
          { 
            if (((int)(coval) != (int)(TARG1)) && (ARES))
              mindec(ret_c,2);
            
            AARG_INC |= ARES;
            if (arg != res)
              ARES_INC = 0;
            else
              ARES_INC;
          } 
        } 
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern reverse mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
	break;

/*--------------------------------------------------------------------------*/
      case m_subscript_ld:                                /* m_subscript_ld */
        size   = get_locint_r();
        arg    = get_locint_r();
        arg1   = get_locint_r();
        arg2   = get_locint_r(); 
#ifdef _TIGHT_
        d      = get_val_v_r(size);
        coval  = get_val_r();

	ASSIGN_T( Targ1, T[arg1])

        if ((int)(coval) != (int)(TARG1))
          mindec(ret_c,2);

        res = arg2 + ((int)(TARG1) + 1)*size + arg;
        for (ls=size; ls>0; ls--)
        { res--;
          
          get_taylor(res);
        
          ASSIGN_A( Ares, A[res])
      
          FOR_0_LE_l_LT_p
            ARES_INC = 0;
        }
#endif /* _TIGHT_ */
#ifdef _SAFE_
        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__ ") : no such operation %d\n\n"
                "active subscripts are not implemented in the safe version"
                " of bit pattern reverse mode\n"
                , operation);
        exit(-1);
#endif /* _SAFE_ */
	break;



/****************************************************************************/
/*                                                          REMAINING STUFF */

/*--------------------------------------------------------------------------*/
      case take_stock_op:                                  /* take_stock_op */
        res  = get_locint_r();
        size = get_locint_r();
#ifdef _TIGHT_
        d    = get_val_v_r(size);
#endif /* _TIGHT_ */

        res += size;
	for (ls=size; ls>0; ls--)
	{ res--;

          ASSIGN_A( Ares, A[res])

          FOR_0_LE_l_LT_p
            ARES_INC = 0;
        }
        break;

/*--------------------------------------------------------------------------*/
      case death_not:                                          /* death_not */
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	for (j=arg1;j<=arg2;j++)
	{ ASSIGN_A(Aarg1, A[j])

          FOR_0_LE_l_LT_p
            AARG1_INC = 0;
#ifdef _TIGHT_	  
          get_taylor(j);
#endif /* _TIGHT_ */
	}
	break;
	
/*--------------------------------------------------------------------------*/
      default:                                                   /* default */
	/*             Die here, we screwed up     */ 

        fprintf(DIAG_OUT,"ADOL-C fatal error in " GENERATED_FILENAME " (" 
                __FILE__
                ") : no such operation %d\n", operation);
	exit(-1);
	break;
      } /* endswitch */

      /* Get the next operation */
      operation=get_op_r();
    } /* endwhile */

  end_sweep();
  return ret_c;
}


/****************************************************************************/

END_C_DECLS
