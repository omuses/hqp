/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     fo_rev.c
 Revision: $Id: fo_rev.c,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: Contains the routines :
           fos_reverse (first-order-scalar reverse mode)  : define _FOS_
           fov_reverse (first-order-vector reverse mode)  : define _FOV_

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
          20030305 andrea: change taylor_back
          19991122 olvo:   new op_codes eq_plus_prod eq_min_prod
                           for  y += x1 * x2
                           and  y -= x1 * x2  
          19981130 olvo:   last check (includes ...)
          19980929 olvo:   allow reflexive operations for
                           - vector operations: eq_mult_av_a,
                             mult_a_av, div_av_a
          19980924 olvo:   (1) Lots of small changes (Atemp, At --> Aqo)
                           (2) new macros AQO*
                           (3) deleted all int_* opcodes
                           (4) allow reflexive operations for
                               - cond_assign, cond_assign_s
                               (changed code completely)
          19980922 olvo:   (1) allow reflexive operations for
                               - div_d_a, div_a_a
          19980921 olvo:   (1) changed save-order in sin/cos
                           (2) allow reflexive operations for
                               - pow
          19980820 olvo:   new comparison strategy 
          19980721 olvo:   write of taylors in subscript and
                           subscript_m   
          19980714 olvo:   some error elimination
                           op-code mult_av_a removed  
          19980713 olvo:   debugging and optimizing
          19980710 olvo:   sin/cos writes 2 taylors
          19980709 olvo:   new operation code: neg_sign_a
                                               pos_sign_a
          19980706 olvo:   new operation code: int_adb_d_one
                                               int_adb_d_zero
          19980703 olvo:   new operation code: assign_d_one
                                               assign_d_zero
          19980626 olvo:   vector operations & subscripts
          19980623 mitev/olvo: revision + stock stuff 
          19980616 olvo:   griewank's idea II
          19980615 olvo:   griewank's idea
          19980612 mitev
          
----------------------------------------------------------------------------*/

/*****************************************************************************
 
  There are four basic versions of the procedure `reverse', which
  are optimized for the cases of scalar or vector reverse sweeps
  with first or higher derivatives, respectively. In the calling
  sequence this distinction is apparent from the type of the
  parameters `lagrange' and `results'. The former may be left out
  and the integer parameters `depen', `indep', `degre', and `nrows'
  must be set or default according to the following matrix of
  calling cases. 

           no lagrange         double* lagrange     double** lagrange

double*   gradient of scalar   weight vector times    infeasible 
results   valued function      Jacobian product       combination

          ( depen = 1 ,         ( depen > 0 ,         
	    degre = 0 ,           degre = 0 ,              ------
	    nrows = 1 )           nrows = 1 )

double**  Jacobian of vector   weight vector times     weight matrix
results   valued function      Taylor-Jacobians        times Jacobian
           
	  ( 0 < depen           ( depen > 0 ,          ( depen > 0 ,
	      = nrows ,           degre > 0 ,            degre = 0 ,
	    degre = 0 )           nrows = 1 )            nrows > 0 )

double*** full family of         ------------          weigth matrix x
results   Taylor-Jacobians       ------------          Taylor Jacobians

*****************************************************************************/ 

/****************************************************************************/
/*                                                                   MACROS */
#undef _ADOLC_VECTOR_

/*--------------------------------------------------------------------------*/
#ifdef _FOS_
#define GENERATED_FILENAME "fos_reverse"              

#define RESULTS(l,indexi)  results[indexi]
#define LAGRANGE(l,indexd) lagrange[indexd] 

/*--------------------------------------------------------------------------*/
#elif _FOV_
#define GENERATED_FILENAME "fov_reverse"             

#define _ADOLC_VECTOR_

#define RESULTS(l,indexi)  results[l][indexi]
#define LAGRANGE(l,indexd) lagrange[l][indexd] 

#else
#error Error ! Define [_FOS_ | _FOV_] 
#endif

/*--------------------------------------------------------------------------*/
/*                                                     access to variables  */

#ifdef _FOS_
#define ARES       *Ares
#define AARG       *Aarg
#define AARG1      *Aarg1
#define AARG2      *Aarg2
#define AQO        *Aqo

#define ARES_INC   *Ares
#define AARG_INC   *Aarg
#define AARG1_INC  *Aarg1
#define AARG2_INC  *Aarg2
#define AQO_INC    *Aqo

#define ARES_INC_O  Ares
#define AARG_INC_O  /adAarg
#define AARG1_INC_O Aarg1
#define AARG2_INC_O Aarg2
#define AQO_INC_O   Aqo

#define ASSIGN_A(a,b)  a = &b;

#else  /* _FOV_, _HOS_, _HOV_ */
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
#endif

#define TRES       T[res]
#define TARG       T[arg]
#define TARG1      T[arg1]
#define TARG2      T[arg2]

#define ASSIGN_T(a,b)

/*--------------------------------------------------------------------------*/
/*                                                              loop stuff  */
#ifdef _ADOLC_VECTOR_
#define FOR_0_LE_l_LT_p for (l=0; l<p; l++)  
#define FOR_p_GT_l_GE_0 for (l=p-1; l>=0; l--)  
#else
#define FOR_0_LE_l_LT_p 
#define FOR_p_GT_l_GE_0  
#endif
 
#define FOR_0_LE_i_LT_k  
#define FOR_k_GT_i_GE_0  

#ifdef _HOV_
#define FOR_0_LE_l_LT_pk1 for (l=0; l<pk1; l++)  
#elif _FOV_
#define FOR_0_LE_l_LT_pk1 for (l=0; l<p; l++)  
#elif _HOS_
#define FOR_0_LE_l_LT_pk1 for (l=0; l<k1; l++)  
#else
#define FOR_0_LE_l_LT_pk1
#endif

/*--------------------------------------------------------------------------*/
/*                                                             other macros */
#define MAXINC(a,b) if ((a) < (b)) (a) = (b)
#define MINDEC(a,b) if ((a) > (b)) (a) = (b)

/* END Macros */


/****************************************************************************/
/*                                                       NECESSARY INCLUDES */
#include "interfaces.h"
#include "adalloc.h"
#include "oplate.h"
#include "taputil.h"
#include "taputil_p.h"
#include "tayutil.h"
#include "tayutil_p.h"

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

#ifdef _FOS_
/****************************************************************************/
/* First-Order Scalar Reverse Pass.                                         */
/****************************************************************************/
int fos_reverse(short   tnum,       /* tape id */
                int     depen,      /* consistency chk on # of deps */
                int     indep,      /* consistency chk on # of indeps */
                double  *lagrange,
                double  *results)   /*  coefficient vectors */

#elif _FOV_
/****************************************************************************/
/* First-Order Vector Reverse Pass.                                         */
/****************************************************************************/

int fov_reverse(short   tnum,        /* tape id */
                int     depen,       /* consistency chk on # of deps */
                int     indep,       /* consistency chk on # of indeps */
                int     nrows,       /* # of Jacobian rows being calculated */
                double  **lagrange,  /* domain weight vector */
                double  **results)   /* matrix of coefficient vectors */  

#endif 

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

  double coval = 0, *d = 0;

  int indexi = 0,  indexd = 0;

  /* loop indices */
#if defined(_FOS_)
  int i;
#else
  int l;
#endif
  int j, ls;

  /* other necessary variables */
  double r0, r_0;
  int buffer;
  static int rax;
#if defined(_FOV_)
  static int pax;
#endif
  int taycheck;
  int numdep,numind;
  double aTmp;

/*--------------------------------------------------------------------------*/
  /* Taylor stuff */
  static revreal *T;
 
/*--------------------------------------------------------------------------*/
  /* Adjoint stuff */
#ifdef _FOS_
  static double* A;
  double Atemp;
#else /* _FOV_, _HOS_, _HOV_ */  
  static double** A;
  static double *Atemp;
#endif
  double  *Ares, *Aarg, *Aarg1, *Aarg2, *Aqo;

/*--------------------------------------------------------------------------*/

#ifdef _ADOLC_VECTOR_
  int p = nrows;
#endif 

#ifdef _HOV_
  int pk1 = p*k1;
#endif


#ifdef DEBUG
/****************************************************************************/
/*                                                           DEBUG MESSAGES */
  fprintf(DIAG_OUT,"Call of %s(..) with tag: %d, n: %d, m %d,\n",
                   GENERATED_FILENAME, tnum, indep, depen);      
#ifdef _ADOLC_VECTOR_
  fprintf(DIAG_OUT,"                    p: %d\n\n",nrows);
#endif
  
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

/*--------------------------------------------------------------------------*/
#ifdef _FOS_                                                         /* FOS */
  if (rev_location_cnt compsize rax) 
  { if (rax)
    { free((char*) T);
      free((char*) A);
    }
    T = (revreal*) malloc(rev_location_cnt*sizeof(revreal));
    A = myalloc1(rev_location_cnt);
    if (T == NULL)
    { fprintf(DIAG_OUT,"ADOL-C error: cannot allocate %i bytes !\n",
              rev_location_cnt*sizeof(revreal));
      exit (-1);
    }
    rax = rev_location_cnt;
  }
  /* olvo 980924 is following initialization necessary ??? */
  Aqo = A;
  for (i=0; i<rev_location_cnt; i++) 
    *Aqo++ = 0.0;

/*--------------------------------------------------------------------------*/
#elif _FOV_                                                          /* FOV */
  if (rev_location_cnt compsize  rax || p compsize pax)
  { if (rax || pax)
    { free((char *) Atemp);
      free((char *) T);
      free((char *) *A); free((char*) A);
    }
    Atemp = myalloc1(p);
    T     = (revreal *)malloc(sizeof(revreal)*rev_location_cnt);
    if (T == NULL)
    { fprintf(DIAG_OUT,"ADOL-C error: cannot allocate %i bytes!\n",
              sizeof(revreal)*rev_location_cnt);
      exit (-1);
    }
    A = myalloc2(rev_location_cnt,p);
    rax = rev_location_cnt;
    pax = p;
  }
#endif


/****************************************************************************/
/*                                                    TAYLOR INITIALIZATION */
  taylor_back(tnum,T,&numdep,&numind,&taycheck);

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
        get_val_block_r(); /* Get the next val block */
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

        ret_c = 0;
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

        ASSIGN_T( Targ, T[arg])

        if (TARG == 0)
          ret_c = 0;
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
        { AARG_INC += ARES;
          ARES_INC = 0.0;
        }

	get_taylor(res);
	break;

/*--------------------------------------------------------------------------*/
      case assign_d:            /* assign an adouble variable a    assign_d */
	                        /* double value. (=) */
        res   = get_locint_r();
        coval = get_val_r();

	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
          ARES_INC = 0.0;
        
        get_taylor(res);
	break;

/*--------------------------------------------------------------------------*/
      case assign_d_zero:  /* assign an adouble variable a    assign_d_zero */
      case assign_d_one:   /* double value (0 or 1). (=)       assign_d_one */
        res   = get_locint_r();

	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
          ARES_INC = 0.0;
        
        get_taylor(res);
	break;

/*--------------------------------------------------------------------------*/
      case assign_ind:       /* assign an adouble variable an    assign_ind */
	                     /* independent double value (<<=) */
        res = get_locint_r();

	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
          RESULTS(l,indexi) = ARES_INC;
        
        get_taylor(res);
	indexi--;
	break;

/*--------------------------------------------------------------------------*/
      case assign_dep:           /* assign a float variable a    assign_dep */
	                         /* dependent adouble value. (>>=) */
        res = get_locint_r();

	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
          ARES_INC = LAGRANGE(l,indexd);

	indexd--;
	break;


/****************************************************************************/
/*                                                   OPERATION + ASSIGNMENT */

/*--------------------------------------------------------------------------*/
      case eq_plus_d:            /* Add a floating point to an    eq_plus_d */
	                         /* adouble. (+=) */
        res   = get_locint_r();
        coval = get_val_r();

	get_taylor(res);
	break;

/*--------------------------------------------------------------------------*/
      case eq_plus_a:             /* Add an adouble to another    eq_plus_a */
	                          /* adouble. (+=) */
        res = get_locint_r();
        arg = get_locint_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg]);

        FOR_0_LE_l_LT_p
          AARG_INC += ARES_INC;

	get_taylor(res);
	break;

/*--------------------------------------------------------------------------*/
      case eq_min_d:       /* Subtract a floating point from an    eq_min_d */
                           /* adouble. (-=) */
        res   = get_locint_r();
        coval = get_val_r();

	get_taylor(res);
	break;
	
/*--------------------------------------------------------------------------*/
      case eq_min_a:        /* Subtract an adouble from another    eq_min_a */
	                    /* adouble. (-=) */
        res = get_locint_r();
        arg = get_locint_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
          AARG_INC -= ARES_INC;

	get_taylor(res);
	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_d:              /* Multiply an adouble by a    eq_mult_d */
	                           /* flaoting point. (*=) */
        res   = get_locint_r();
        coval = get_val_r();

	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
          ARES_INC *= coval;

	get_taylor(res);
	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_a:       /* Multiply one adouble by another    eq_mult_a */
	                    /* (*=) */
        res = get_locint_r();
        arg = get_locint_r();

 	get_taylor(res);

	ASSIGN_A( Ares, A[res])
        ASSIGN_A( Aarg, A[arg])
        ASSIGN_T( Tres, T[res])
	ASSIGN_T( Targ, T[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
	  /* olvo 980713 nn: ARES = 0.0; */ 
	  ARES_INC =  aTmp * TARG;
	  AARG_INC += aTmp * TRES; 
        }
        break;
        
/*--------------------------------------------------------------------------*/
      case incr_a:                        /* Increment an adouble    incr_a */
      case decr_a:                        /* Increment an adouble    decr_a */
        res   = get_locint_r();

	get_taylor(res);
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
          ARES_INC = 0.0; 
          AARG1_INC += aTmp;
          AARG2_INC += aTmp;          
        }

      	get_taylor(res);
	break;

/*--------------------------------------------------------------------------*/
      case plus_d_a:             /* Add an adouble and a double    plus_d_a */
	                         /* (+) */
        res   = get_locint_r();
        arg   = get_locint_r();
        coval = get_val_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG_INC += aTmp;
        }

      	get_taylor(res);
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
          ARES_INC = 0.0; 
          AARG1_INC += aTmp;
          AARG2_INC -= aTmp;
        }

        get_taylor(res);
	break;

/*--------------------------------------------------------------------------*/
      case min_d_a:                /* Subtract an adouble from a    min_d_a */
	                           /* double (-) */
        res   = get_locint_r();
        arg   = get_locint_r();
        coval = get_val_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG_INC -= aTmp;
        }

        get_taylor(res);
	break;
	
/*--------------------------------------------------------------------------*/
      case mult_a_a:               /* Multiply two adoubles (*)    mult_a_a */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

        get_taylor(res);

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg2, A[arg2])
	ASSIGN_A( Aarg1, A[arg1])
        ASSIGN_T( Targ1, T[arg1])
	ASSIGN_T( Targ2, T[arg2])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG2_INC += aTmp * TARG1;
	  AARG1_INC += aTmp * TARG2;
        }
        break;

/*--------------------------------------------------------------------------*/
      /* olvo 991122: new op_code with recomputation */
      case eq_plus_prod:   /* increment a product of           eq_plus_prod */
                           /* two adoubles (*) */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg2, A[arg2])
	ASSIGN_A( Aarg1, A[arg1])
        ASSIGN_T( Targ1, T[arg1])
	ASSIGN_T( Targ2, T[arg2])

        /* RECOMPUTATION */
        ASSIGN_T( Tres,  T[res])
        TRES -= TARG1*TARG2;

        FOR_0_LE_l_LT_p
        { AARG2_INC += ARES    * TARG1;
	  AARG1_INC += ARES_INC * TARG2;
        }
        break;

/*--------------------------------------------------------------------------*/
      /* olvo 991122: new op_code with recomputation */
      case eq_min_prod:    /* decrement a product of            eq_min_prod */
                           /* two adoubles (*) */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg2, A[arg2])
	ASSIGN_A( Aarg1, A[arg1])
        ASSIGN_T( Targ1, T[arg1])
	ASSIGN_T( Targ2, T[arg2])

        /* RECOMPUTATION */
        ASSIGN_T( Tres,  T[res])
        TRES += TARG1*TARG2;

        FOR_0_LE_l_LT_p
        { AARG2_INC -= ARES    * TARG1;
	  AARG1_INC -= ARES_INC * TARG2;
        }
        break;

/*--------------------------------------------------------------------------*/
      case mult_d_a:         /* Multiply an adouble by a double    mult_d_a */
                             /* (*) */
        res   = get_locint_r();
        arg   = get_locint_r();
        coval = get_val_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG_INC += coval * aTmp;
        }

        get_taylor(res);
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
	ASSIGN_T( Tres,  T[res])
        ASSIGN_T( Targ2, T[arg2])

	/* olvo 980922 changed order to allow x=y/x */
	r_0 = -TRES;
        get_taylor(res);
        r0  = 1.0 / TARG2;
	r_0 *= r0;

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG1_INC += aTmp * r0;
	  AARG2_INC += aTmp * r_0;
        }

        break;

/*--------------------------------------------------------------------------*/
      case div_d_a:             /* Division double - adouble (/)    div_d_a */
        res   = get_locint_r();
        arg   = get_locint_r();
        coval = get_val_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])
	ASSIGN_T( Tres, T[res])
        ASSIGN_T( Targ, T[arg])

	/* olvo 980922 changed order to allow x=d/x */
        r0 = -TRES;
        if (arg == res)
          get_taylor(arg);
        r0 /= TARG;

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG_INC += aTmp * r0;
        }

        get_taylor(res);
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
          ARES_INC = 0.0; 
          AARG_INC += aTmp;
        }

      	get_taylor(res);
	break;

/*--------------------------------------------------------------------------*/
      case neg_sign_a:                                        /* neg_sign_a */
        res   = get_locint_r();
        arg   = get_locint_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG_INC -= aTmp;
        }

      	get_taylor(res);
	break;


/****************************************************************************/
/*                                                         UNARY OPERATIONS */

/*--------------------------------------------------------------------------*/
      case exp_op:                          /* exponent operation    exp_op */
        res = get_locint_r();
        arg = get_locint_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])
        ASSIGN_T( Tres, T[res])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG_INC += aTmp * TRES;
        }

 	get_taylor(res);
        break;

/*--------------------------------------------------------------------------*/
      case sin_op:                              /* sine operation    sin_op */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg1, A[arg1])
        ASSIGN_T( Targ2, T[arg2])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG1_INC += aTmp * TARG2;
        }

       	get_taylor(res);
       	get_taylor(arg2); /* olvo 980710 covalue */
	                  /* NOTE: A[arg2] should be 0 already */
	break;

/*--------------------------------------------------------------------------*/
      case cos_op:                            /* cosine operation    cos_op */
        res  = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg1, A[arg1])
        ASSIGN_T( Targ2, T[arg2])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG1_INC -= aTmp * TARG2;          
        }

       	get_taylor(res);
       	get_taylor(arg2); /* olvo 980710 covalue */
	                  /* NOTE A[arg2] should be 0 already */
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

        get_taylor(res);

        ASSIGN_A( Ares,  A[res])
        ASSIGN_A( Aarg1, A[arg1])
        ASSIGN_T( Targ2, T[arg2])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG1_INC += aTmp * TARG2;
        }
	break;

/*--------------------------------------------------------------------------*/
      case log_op:                                                /* log_op */
        res = get_locint_r();
        arg = get_locint_r();

        get_taylor(res);

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])
	ASSIGN_T( Targ, T[arg])

        r0 = 1.0/TARG;

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG_INC += aTmp * r0;
        }
	break;

/*--------------------------------------------------------------------------*/
      case pow_op:                                                /* pow_op */
        res   = get_locint_r();
        arg   = get_locint_r();
        coval = get_val_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])
        ASSIGN_T( Tres, T[res])
        ASSIGN_T( Targ, T[arg])

	/* olvo 980921 changed order to allow x=pow(x,n) */
        r0 = TRES;
        if (arg == res)
          get_taylor(arg);
        if (TARG == 0.0) 
          r0 = 0.0;
        else
          r0 *= coval/TARG;

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG_INC += aTmp * r0;
        }
 
        get_taylor(res);
        break;

/*--------------------------------------------------------------------------*/
      case sqrt_op:                                              /* sqrt_op */
        res = get_locint_r();
        arg = get_locint_r();

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])
        ASSIGN_T( Tres, T[res])

        if (TRES == 0.0)
          r0 = 0.0;
        else 
          r0 = 0.5 / TRES;

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG_INC += aTmp * r0;
        }
  
        get_taylor(res);
        break;

/*--------------------------------------------------------------------------*/
      case gen_quad:                                            /* gen_quad */
        res   = get_locint_r();
        arg2  = get_locint_r();
        arg1  = get_locint_r();
        coval = get_val_r();
        coval = get_val_r();

	ASSIGN_A( Ares,  A[res])
	ASSIGN_A( Aarg1, A[arg1])
        ASSIGN_T( Targ2, T[arg2])

        FOR_0_LE_l_LT_p
        { aTmp = ARES;
          ARES_INC = 0.0; 
          AARG1_INC += aTmp * TARG2;
        }
  
      	get_taylor(res);
	break;

/*--------------------------------------------------------------------------*/
      case min_op:                                                /* min_op */
        res   = get_locint_r();
        arg2  = get_locint_r();
        arg1  = get_locint_r();
        coval = get_val_r();
  
        get_taylor(res);

        ASSIGN_A( Aarg1, A[arg1])
        ASSIGN_A( Aarg2, A[arg2])
        ASSIGN_A( Ares,  A[res])
        ASSIGN_T( Targ1, T[arg1])
        ASSIGN_T( Targ2, T[arg2])

        if (TARG1 > TARG2)
          FOR_0_LE_l_LT_p
          { aTmp = ARES;
            ARES_INC = 0.0; 
            if ((coval) && (aTmp))
              MINDEC(ret_c,2);
            AARG2_INC += aTmp;
          }
        else 
          if (TARG1 < TARG2)
            FOR_0_LE_l_LT_p
            { aTmp = ARES;
              ARES_INC = 0.0; 
              if ((!coval) && (aTmp))
                MINDEC(ret_c,2);
              AARG1_INC += aTmp;
            }
          else
          { /* both are equal */
            FOR_0_LE_l_LT_p
            { aTmp = ARES / 2.0;
              ARES_INC = 0.0; 
              AARG2_INC += aTmp;
              AARG1_INC += aTmp;
            }
            if (arg1 != arg2)
              MINDEC(ret_c,1);
          }
        break;

/*--------------------------------------------------------------------------*/
      case abs_val:                                              /* abs_val */
        res   = get_locint_r();
        arg   = get_locint_r();
        coval = get_val_r();
  
        get_taylor(res);

        ASSIGN_A( Ares, A[res])
        ASSIGN_A( Aarg, A[arg])
        ASSIGN_T( Targ, T[arg])

        if (TARG < 0.0)
          FOR_0_LE_l_LT_p
          { aTmp = ARES;
            ARES_INC = 0.0; 
            if ((coval) && (aTmp))
              MINDEC(ret_c,2);
            AARG_INC -= aTmp;
          } 
        else 
          if (TARG > 0.0)
            FOR_0_LE_l_LT_p
            { aTmp = ARES;
              ARES_INC = 0.0; 
              if ((!coval) && (aTmp))
                MINDEC(ret_c,2);
              AARG_INC += aTmp;
            }
          else 
            FOR_0_LE_l_LT_p 
            { aTmp = ARES;
              ARES_INC = 0.0; 
              if (aTmp)        
                MINDEC(ret_c,1);
            }
        break;

/*--------------------------------------------------------------------------*/
      case ceil_op:                                              /* ceil_op */
        res   = get_locint_r();
        arg   = get_locint_r();
        coval = get_val_r();

        get_taylor(res);

        ASSIGN_A( Ares, A[res])
        ASSIGN_T( Targ, T[arg])

        coval = (coval != ceil(TARG) );

        FOR_0_LE_l_LT_p
        { if ((coval) && (ARES))
            MINDEC(ret_c,2);
          ARES_INC = 0.0;
        }
        break;

/*--------------------------------------------------------------------------*/
      case floor_op:                                            /* floor_op */
        res   = get_locint_r();
        arg   = get_locint_r();
        coval = get_val_r();

        get_taylor(res);

        ASSIGN_A( Ares, A[res])
        ASSIGN_A( Aarg, A[arg])
        ASSIGN_T( Targ, T[arg])

        coval = ( coval != floor(TARG1) );

        FOR_0_LE_l_LT_p
        { if ( (coval) && (ARES) )
            MINDEC(ret_c,2);
          ARES_INC = 0.0;
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
        coval  = get_val_r(); 

	get_taylor(res);

        ASSIGN_A( Aarg1, A[arg1])
	ASSIGN_A( Ares,  A[res])
        ASSIGN_A( Aarg2, A[arg2])
        ASSIGN_T( Targ,  T[arg])

	/* olvo 980924 changed code a little bit */
        if (TARG > 0.0)
        { if (res != arg1)
            FOR_0_LE_l_LT_p
            { if ((coval <= 0.0) && (ARES))
                MINDEC(ret_c,2);
              AARG1_INC += ARES;
              ARES_INC = 0.0;
	    }
          else
            FOR_0_LE_l_LT_p
              if ((coval <= 0.0) && (ARES_INC))
                MINDEC(ret_c,2);
	}
        else
        { if (res != arg2)
            FOR_0_LE_l_LT_p
            { if ((coval <= 0.0) && (ARES))
                MINDEC(ret_c,2);
              AARG2_INC += ARES;
              ARES_INC = 0.0;
	    }
          else
            FOR_0_LE_l_LT_p
              if ((coval <= 0.0) && (ARES_INC))
                MINDEC(ret_c,2);
	}
        break;

/*--------------------------------------------------------------------------*/
      case cond_assign_s:                                  /* cond_assign_s */
        res   = get_locint_r();
        arg1  = get_locint_r();
        arg   = get_locint_r(); 
        coval = get_val_r();

        get_taylor(res);

        ASSIGN_A( Aarg1, A[arg1])
	ASSIGN_A( Ares,  A[res])
        ASSIGN_T( Targ,  T[arg])
        
	/* olvo 980924 changed code a little bit */
        if (TARG > 0.0)
        { if (res != arg1)
            FOR_0_LE_l_LT_p
            { if ((coval <= 0.0) && (ARES))
                MINDEC(ret_c,2);
              AARG1_INC += ARES;
              ARES_INC = 0.0;
	    }
          else
            FOR_0_LE_l_LT_p
              if ((coval <= 0.0) && (ARES_INC))
                MINDEC(ret_c,2);
	}
        else
          if (TARG == 0.0) /* we are at the tie */
            FOR_0_LE_l_LT_p
              if (ARES_INC)
                MINDEC(ret_c,0);
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
	  { AARG_INC += ARES;
	    ARES_INC = 0.0;
	  }

	  get_taylor(res);
	  }
	break;

/*--------------------------------------------------------------------------*/
      case assign_dv:                                          /* assign_dv */
        res  = get_locint_r();
        size = get_locint_r();
        d    = get_val_v_r(size);

        res += size;
        d   += size;
	for (ls=size; ls>0; ls--)
	{ res--;                /* Location of left-hand-side */
	  coval = *(--d);       /* Value of right-hand-side   */     

          /* code for assign_d */
	  ASSIGN_A( Ares, A[res])

	  FOR_0_LE_l_LT_p
	    ARES_INC = 0.0;

	  get_taylor(res);
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
	    RESULTS(l,indexi) = ARES_INC;
	  indexi--;

	  get_taylor(res);
        }
	reset_val_r();
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
            ARES_INC = LAGRANGE(l, indexd);
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
	    AARG_INC += ARES_INC;

          get_taylor(res);
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
	    AARG_INC -= ARES_INC;

          get_taylor(res);
	}
	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_av_d:                                    /* eq_mult_av_d */
        res   = get_locint_r();
        size  = get_locint_r();
        coval = get_val_r();

        res += size;
	for (ls=size; ls>0; ls--)
        { res--;            /* Location of the left-hand-side  */
	  /* coval = fixed;    value on the right-hand-side */

	  /* code for eq_mult_d*/
	  ASSIGN_A( Ares, A[res])

          FOR_0_LE_l_LT_p
            ARES_INC *= coval;

	  get_taylor(res);
        }
	break;

/*--------------------------------------------------------------------------*/
      case eq_mult_av_a:                                    /* eq_mult_av_a */
        res  = get_locint_r();
        size = get_locint_r();
        arg  = get_locint_r();

        /* olvo 980929 new strategy to check for overwrites 
           (changes computation order) */
        if ((arg >= res) && (arg < res+size)) 
	{ /* FIRST compute the case: res==arg */
          /* simplified code for eq_mult_a*/
	  get_taylor(arg);
	
	  ASSIGN_A( Aarg, A[arg]) 
          ASSIGN_T( Targ, T[arg])

          FOR_0_LE_l_LT_p
            AARG_INC  *= 2.0 * TARG; 
	}

	res += size;
        for (ls=size; ls>0; ls--)
	{ res--;                 /* Location of the left-hand-side  */
	  /* arg    = fixed;        Location on the right-hand-side */
     
          if (res == arg) /* NOW skip this case */
            continue;

          /* code for eq_mult_a*/
	  get_taylor(res);
	
	  ASSIGN_A( Aarg, A[arg]) 
	  ASSIGN_A( Ares, A[res]) 
          ASSIGN_T( Tres, T[res])
          ASSIGN_T( Targ, T[arg])

          FOR_0_LE_l_LT_p
          { r0 = ARES;
	    /* olvo 980713 nn: ARES = 0; */
            ARES_INC  = r0 * TARG; 
            AARG_INC += r0 * TRES; 
          }     
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
            ARES_INC = 0.0; 
            AARG1_INC += aTmp;
            AARG2_INC += aTmp;          
          }

      	  get_taylor(res);
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
            ARES_INC = 0.0; 
            AARG1_INC += aTmp;
            AARG2_INC -= aTmp;          
          }

          get_taylor(res);
	}
	break;

/*--------------------------------------------------------------------------*/
      case dot_av_av:                                          /* dot_av_av */
        res  = get_locint_r();
        size = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

        get_taylor(res);

        /* save Ares to Atemp */
        ASSIGN_A( Aqo,  Atemp)
	ASSIGN_A( Ares, A[res])
        FOR_0_LE_l_LT_p
	{ AQO_INC = ARES;
          ARES_INC = 0.0; 
	}

	for (ls=0; ls<size; ls++)
	{ /* code for mult_a_a  */
	  ASSIGN_A( Aarg2, A[arg2])
	  ASSIGN_A( Aarg1, A[arg1])
          ASSIGN_A( Aqo,   Atemp)
	  ASSIGN_T( Targ1, T[arg1])
	  ASSIGN_T( Targ2, T[arg2])
         

	  FOR_0_LE_l_LT_p
          { AARG2_INC += AQO     * TARG1;
	    AARG1_INC += AQO_INC * TARG2;
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

        /* olvo 980929 new strategy to check for overwrites 
           (changes computation order) */
        if ((arg2 >= res) && (arg2 < res+size)) 
	{ /* FIRST compute the case: res==arg2 */
	  /* simplified code for mult_a_a */
    	  get_taylor(arg2);

	  ASSIGN_A( Aarg1, A[arg1+res-arg2])
	  ASSIGN_A( Aarg2, A[arg2])
	  ASSIGN_T( Targ2, T[arg2])
	  ASSIGN_T( Targ1, T[arg1+res-arg2])

          FOR_0_LE_l_LT_p
          { AARG1_INC += AARG2 * TARG2;
	    AARG2_INC *=         TARG1;
          }
	}

        res  += size;
        arg1 += size;
	for (ls=size; ls>0; ls--)
	{ arg1--;    /* Location of rght hnd side vectore[l]  */
	  res--;     /* Location of the result */
	    
          if (res == arg2) /* NOW skip this case */
            continue;

	  /* code for mult_a_a */
    	  get_taylor(res);

	  ASSIGN_A( Ares,  A[res])
	  ASSIGN_A( Aarg2, A[arg2])
	  ASSIGN_A( Aarg1, A[arg1])
	  ASSIGN_T( Targ1, T[arg1])
	  ASSIGN_T( Targ2, T[arg2])

          FOR_0_LE_l_LT_p
          { aTmp = ARES;
            ARES_INC = 0.0; 
            AARG2_INC += aTmp * TARG1;
	    AARG1_INC += aTmp * TARG2;
          }
	}
	break;

/*--------------------------------------------------------------------------*/
      case mult_d_av:                                          /* mult_d_av */
        res   = get_locint_r();
        size  = get_locint_r();
        arg   = get_locint_r();
        coval = get_val_r();

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
            ARES_INC = 0.0; 
            AARG_INC += coval * aTmp;
          }

          get_taylor(res);
        }
	break;

/*--------------------------------------------------------------------------*/
      case div_av_a:                                            /* div_av_a */
        res  = get_locint_r();
        size = get_locint_r();
        arg2 = get_locint_r();
        arg1 = get_locint_r();

        /* olvo 980929 new strategy to check for overwrites 
           (changes computation order) */
        if ((arg2 >= res) && (arg2 < res+size)) 
	{ /* FIRST compute the case: res==arg2 */
	  /* code for div_a_a */
	  ASSIGN_A( Aarg2, A[arg2])
	  ASSIGN_A( Aarg1, A[arg1+res-arg2])
	  ASSIGN_T( Targ2, T[arg2])

	  /* olvo 980922 changed order to allow x=y/x */
	  r_0 = -TARG2;
          get_taylor(arg2);
          r0  = 1.0 / TARG2;
 	  r_0 *= r0;

          FOR_0_LE_l_LT_p
          { AARG1_INC += AARG2 * r0;
	    AARG2_INC *=         r_0;
          }
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
	  ASSIGN_T( Tres,  T[res])
	  ASSIGN_T( Targ2, T[arg2])

	  /* olvo 980922 changed order to allow x=y/x */
	  r_0 = -TRES;
          get_taylor(res);
          r0  = 1.0 / TARG2;
 	  r_0 *= r0;

          FOR_0_LE_l_LT_p
          { aTmp = ARES;
            ARES_INC = 0.0; 
            AARG1_INC += aTmp * r0;
	    AARG2_INC += aTmp * r_0;
          }
	}
	break;


/****************************************************************************/
/*                                                               SUBSCRIPTS */

/*--------------------------------------------------------------------------*/
      case subscript:                                          /* subscript */
        res   = get_locint_r();
        arg1  = get_locint_r();
        arg2  = get_locint_r(); 
        coval = get_val_r();

	ASSIGN_T( Targ1, T[arg1])

        arg  = arg2 + (int)(TARG1);

	/* olvo 980721 new nl */
        get_taylor(res);

	ASSIGN_A( Aarg, A[arg])
	ASSIGN_A( Ares, A[res])

        FOR_0_LE_l_LT_p
        { if (((int)(coval) != (int)(TARG1)) && (ARES))
            MINDEC(ret_c,2);

          AARG_INC += ARES;
          if (arg != res)
            ARES_INC = 0;
#if defined(_FOV_)
          else 
            ARES_INC;
#endif
        }

	break;

/*--------------------------------------------------------------------------*/
      case subscript_l:                                      /* subscript_l */
        arg   = get_locint_r();
        arg1  = get_locint_r();
        arg2  = get_locint_r(); 
        coval = get_val_r();

	ASSIGN_T( Targ1, T[arg1])

        res   = arg2 + (int)(TARG1);

	get_taylor(res);

	ASSIGN_A( Ares, A[res])
	ASSIGN_A( Aarg, A[arg])

	FOR_0_LE_l_LT_p
        { if (((int)(coval) != (int)(TARG1)) && (ARES))
            MINDEC(ret_c,2);

          AARG_INC += ARES;
          if(arg != res)
            ARES_INC = 0;
#if defined(_FOV_)
          else
            ARES_INC;
#endif
        }
	break;
      
/*--------------------------------------------------------------------------*/
      case subscript_ld:                                    /* subscript_ld */
        arg1  = get_locint_r();
        arg2  = get_locint_r(); 
        coval = get_val_r();
        coval = get_val_r();

	ASSIGN_T( Targ1, T[arg1])

        arg = arg2 + (int)(TARG1);

        get_taylor(arg);

        if((int)(coval)!=(int)(TARG1))
          MINDEC(ret_c,2);

	ASSIGN_A( Aarg, A[arg])

        FOR_0_LE_l_LT_p
	  AARG_INC = 0.0;
	break;
       
/*--------------------------------------------------------------------------*/
      case m_subscript:                                      /* m_subscript */
        res   = get_locint_r();
        size  = get_locint_r();
        arg1  = get_locint_r();
        arg2  = get_locint_r(); 
        coval = get_val_r();

	ASSIGN_T( Targ1, T[arg1])

        arg = arg2 + ((int)(TARG1) + 1)*size;
        res += size;
        for (ls=size; ls>0; ls--)
        { res--; arg--;

   	  /* olvo 980721 new nl */
          get_taylor(res);

          ASSIGN_A( Aarg, A[arg])
	  ASSIGN_A( Ares, A[res])

          FOR_0_LE_l_LT_p
          { if (((int)(coval)!=(int)(TARG1)) && (ARES))
              MINDEC(ret_c,2);
            AARG_INC += ARES;
            if (arg != res)
              ARES_INC = 0;
#if defined(_FOV_)
            else
              ARES_INC;
#endif
          }
        }
	break;

/*--------------------------------------------------------------------------*/
      case m_subscript_l:                                  /* m_subscript_l */
        arg   = get_locint_r();
        size  = get_locint_r();
        arg1  = get_locint_r();
        arg2  = get_locint_r(); 
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
          { if (((int)(coval) != (int)(TARG1)) && (ARES))
              MINDEC(ret_c,2);
            AARG_INC += ARES;
            if (arg != res)
              ARES_INC = 0;
#if defined(_FOV_)
            else
              ARES_INC;
#endif
          } 
        } 
	break;

/*--------------------------------------------------------------------------*/
      case m_subscript_ld:                                /* m_subscript_ld */
        size   = get_locint_r();
        arg    = get_locint_r();
        arg1   = get_locint_r();
        arg2   = get_locint_r(); 
        /* olvo 980702 changed n2l */
        d      = get_val_v_r(size);
        coval  = get_val_r();

	ASSIGN_T( Targ1, T[arg1])

        if ((int)(coval) != (int)(TARG1))
          MINDEC(ret_c,2);

        res = arg2 + ((int)(TARG1) + 1)*size + arg;
        for (ls=size; ls>0; ls--)
        { res--;
          
          get_taylor(res);
        
          ASSIGN_A( Ares, A[res])
      
          FOR_0_LE_l_LT_p
            ARES_INC = 0.0;
        }
	break;


/****************************************************************************/
/*                                                          REMAINING STUFF */

/*--------------------------------------------------------------------------*/
      case take_stock_op:                                  /* take_stock_op */
        res  = get_locint_r();
        size = get_locint_r();
        d    = get_val_v_r(size);

        res += size;
	for (ls=size; ls>0; ls--)
	{ res--;

          ASSIGN_A( Ares, A[res])

          FOR_0_LE_l_LT_p
            ARES_INC = 0.0;
        }
        break;

/*--------------------------------------------------------------------------*/
      case death_not:                                          /* death_not */
        arg2 = get_locint_r();
        arg1 = get_locint_r();

	for (j=arg1;j<=arg2;j++)
	{ ASSIGN_A(Aarg1, A[j])

          FOR_0_LE_l_LT_p
            AARG1_INC = 0.0;
	  
          get_taylor(j);
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
/*                                                               THAT'S ALL */

END_C_DECLS
