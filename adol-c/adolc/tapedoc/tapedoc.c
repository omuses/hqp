/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     tapedoc/tapedoc.c
 Revision: $Id: tapedoc.c,v 1.2 2004/10/14 13:29:48 e_arnold Exp $
 Contents: Routine tape_doc(..) writes the taped operations in LaTeX-code 
           to the file tape_doc.tex

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
        19991122 olvo:   new op_codes eq_plus_prod eq_min_prod
                         for  y += x1 * x2
                         and  y -= x1 * x2  
        19981130 olvo:   last check (includes ...)
        19980930 olvo:   allow reflexive vecops
        19980924 olvo:   deleted all int_* opcodes
        19980923 olvo:   allow reflexive ops 
        19980820 olvo:   new comparison strategy 
        19980714 olvo:   removed operation code: mult_av_a
        19980709 olvo:   new operation code: neg_sign_a
                                             pos_sign_a
        19980706 olvo:   new operation code: assign_d_one
                                             assign_d_zero
                                             int_adb_d_one
                                             int_adb_d_zero
                                             incr_a
                                             decr_a
        19980623 olvo:   new operation code: take_stock_op

----------------------------------------------------------------------------*/

#include "../oplate.h"
#include "../taputil.h"
#include "../taputil_p.h"
#include "../tayutil.h"
#include "../tayutil_p.h"
#include "../adalloc.h"

#include <math.h>

BEGIN_C_DECLS

/****************************************************************************/
/*                                                                   MACROS */
#define computenumbers true


/****************************************************************************/
/*                                                         STATIC VARIABLES */

/*--------------------------------------------------------------------------*/
static short tag;

static int for_location_cnt;

static int file_cnt;
static int read_cnt;
static int pagelength;
static FILE *fp;

/*--------------------------------------------------------------------------*/
/* operation names */
static char* a[] =  {  "ignore me","death not","assign ind",
           "assign dep","assign a","assign d",
           "eq plus d","eq plus a","eq min d",
           "eq min a","eq mult d","eq mult a",
           "plus a a","plus d a","min a a",
           "min d a","mult a a","mult d a",
           "div a a","div d a","exp op",
           "cos op","sin op","atan op","log op",
           "pow op","asin op","acos op","sqrt op",
           "eq div a","eq div d","tan op","gen quad",
           "not used","not used","end of tape",
           "start of tape","end of op","end of int",
           "end of val","plus av av","plus dv av",
           "sub av av","sub dv av","sub av dv",
           "dot av av","dot dv av","mult a av",
           "mult d av","mult a dv","not used",
           "not used","assign av","assign dv",
           "assign indvec","assign depvec","eq min dv",
           "eq min av","eq plus dv","eq plus av",
           "div av a","eq mult av d","eq mult av a",
           "dot av dv","nicht belegt","nicht belegt",
           "not used","not used","not used",
           "not used","cond assign $\\longrightarrow$",
           "cond assign s $\\longrightarrow$",
           "m subscript $\\longrightarrow$",
           "m subscript l $\\longrightarrow$",
           "m subscript ld $\\longrightarrow$",
           "subscript $\\longrightarrow$",
           "subscript l $\\longrightarrow$",
           "subscript ld $\\longrightarrow$",
           "not used","not used",
           "cross av av","mult cv3 av4","not used",
           "not used","not used","not used","not used",
           "not used","not used","not used","take stock op",
           "assign d one","assign d zero","not used",
           "not used","incr a","decr a","neg sign a",
           "pos sign a","not used","min op","abs val",
           "eq zero","neq zero","le zero","gt zero",
           "ge zero","lt zero","not used","not used",
           "eq plus prod","eq min prod","not used","not used",
           "erf op","ceil op","floor op"
       };

/****************************************************************************/
/*                                                     LOCAL WRITE ROUTINES */

/*--------------------------------------------------------------------------*/
void filewrite_start( int opcode )
{ if ((fp = fopen("tape.tex","w")) == NULL) 
  { fprintf(DIAG_OUT,"cannot open file !\n");
    exit(1);
  }
  fprintf(fp,"\\documentstyle{article} \n");
  fprintf(fp,"\\headheight0cm\n");
  fprintf(fp,"\\headsep-1cm\n");
  fprintf(fp,"\\textheight25cm\n");
  fprintf(fp,"\\oddsidemargin-1cm\n");
  fprintf(fp,"\\topmargin0cm\n");
  fprintf(fp,"\\textwidth18cm\n");
  fprintf(fp,"\\begin{document}\n");
  fprintf(fp,"\\tiny\n");
#ifdef computenumbers
  fprintf(fp,"\\begin{tabular}{|r|l|r|r|r|r||r|r||r|r|r|r|} \\hline \n");
  fprintf(fp," code & op & loc & loc & loc & loc & double & double & value & value & value & value \\\\ \\hline \n");
  fprintf(fp," %i & start of tape & & & & & & & & & &  \\\\ \\hline \n",opcode);
#else
  fprintf(fp,"\\begin{tabular}{|r|l|r|r|r|r||r|r|} \\hline \n");
  fprintf(fp," code & op & loc & loc & loc & loc & double & double \\\\ \\hline \n");
  fprintf(fp," %i & start of tape & & & & & & & \\\\ \\hline \n",opcode);
#endif
  pagelength = 0;
}

/*--------------------------------------------------------------------------*/
void filewrite( unsigned short opcode, int nloc, int *loc, 
                double *val,int ncst, double* cst) 
{ 
  int i;
  
  if (pagelength == 100) 
  { fprintf(fp,"\\end{tabular}\\\\\n");
    fprintf(fp,"\\newpage\n");
#ifdef computenumbers
    fprintf(fp,"\\begin{tabular}{|r|l|r|r|r|r||r|r||r|r|r|r|} \\hline \n");
    fprintf(fp," code & op & loc & loc & loc & loc & double & double & value & value & value & value \\\\ \\hline \n");
#else
    fprintf(fp,"\\begin{tabular}{|r|l|r|r|r|r||r|r|} \\hline \n");
    fprintf(fp," code & op & loc & loc & loc & loc & double & double \\\\ \\hline \n");
#endif
    pagelength=-1;
  } 
  fprintf(fp,"%i & ",opcode);

  i=0;
  while (a[opcode][i])
  { fprintf(fp,"%c",a[opcode][i]); 
    i++;
  } 

  fprintf(fp," &");
  for(i=0; i<(4-nloc); i++)  
    fprintf(fp," &");
  for(i=0; i<nloc; i++)
    fprintf(fp," %i &",loc[i]);
#ifdef computenumbers
  for(i=0; i<(2-ncst); i++)  
    fprintf(fp," &");
  for(i=0; i<ncst; i++)
    fprintf(fp,"$ %e $&",cst[i]);
  for(i=0; i<(4-nloc); i++)  
    fprintf(fp," &");
  for(i=0; i<nloc-1; i++)
    fprintf(fp,"$ %e $&",val[i]);
  if (nloc)
    fprintf(fp,"$ %e $",val[nloc-1]);
  else
    fprintf(fp," ");
#else
  for(i=0; i<(2-ncst); i++)
    fprintf(fp," &");
  for(i=0; i<ncst-1; i++)
    fprintf(fp,"$ %e $ &",cst[i]);
  if (ncst)
    fprintf(fp,"$ %e $",val[ncst-1]);
  else
    fprintf(fp," ");
#endif
  fprintf(fp,"\\\\ \\hline \n");
  pagelength++;
}
    
/*--------------------------------------------------------------------------*/
void filewrite_end( int opcode ) 
{  
#ifdef computenumbers
  fprintf(fp," %i & end of tape & & & & & & & & & &  \\\\ \\hline \n",opcode);
#else
  fprintf(fp," %i & end of tape & & & & & & & \\\\ \\hline \n",opcode);
#endif 
  fprintf(fp,"\\end{tabular}");
  fprintf(fp,"\\end{document}"); 
  fclose(fp);                                                    
}
     

/****************************************************************************/
/*                                                             NOW THE CODE */
void tape_doc(short tnum,         /* tape id */
              int depcheck,       /* consistency chk on # of dependents */
	      int indcheck,       /* consistency chk on # of independents */
              double *basepoint,  /* independent variable values */
              double *valuepoint) /* dependent variable values */
{
/****************************************************************************/
/*                                                            ALL VARIABLES */
  unsigned char operation;
  int tape_stats[11];  /* tape stats */

  locint size = 0;
  locint res  = 0;
  locint arg  = 0;
  locint arg1 = 0;
  locint arg2 = 0;

  double coval = 0, *d = 0;

  int indexi = 0, indexd = 0;

  /* loop indices */ 
  int  l;

  /* other necessary variables */
  int buffer; 
  static int fax;

  /* Taylor stuff */
  static double* T0;
  double T0res;
 
  /* interface temporaries */
  int loc_a[4];
  double val_a[4], cst_d[2];


/****************************************************************************/
/*                                                                    INITs */
  tag = tnum;         /*tag is global which specifies which tape to look at */
  read_cnt = 0;	    
  file_cnt = 80; 
   
  tapestats(tag,tape_stats);
  for_location_cnt = tape_stats[2];
  buffer           = tape_stats[4];

  set_buf_size(buffer);

  if ((depcheck != tape_stats[1]) || (indcheck != tape_stats[0]))
  { fprintf(DIAG_OUT,"ADOL-C error: Tape_doc on tape %d  aborted!\n",tag);
    fprintf(DIAG_OUT,"Number of dependent and/or independent variables passed to Tape_doc is\ninconsistant with number recorded on tape %d \n",tag);
    exit (-1);
  }

/****************************************************************************/
/*                                                        MEMORY ALLOCATION */
  if (for_location_cnt compsize fax) {
    if (fax)
      free((char *) T0);
    T0 = myalloc1(for_location_cnt);
    fax=for_location_cnt;
  }
  
/****************************************************************************/
/*                                                       FORWARD TAPE SWEEP */

  /* Initialize the Forward Sweep */
  init_for_sweep(tag);

  operation=get_op_f();
  while (operation !=end_of_tape)
  { switch (operation){

/****************************************************************************/
/*                                                                  MARKERS */

/*--------------------------------------------------------------------------*/
      case end_of_op:                                          /* end_of_op */ 
    	filewrite(operation,0,loc_a,val_a,0,cst_d);
        get_op_block_f();
        operation=get_op_f(); 
        /* Skip next operation, it's another end_of_op */
        break;

/*--------------------------------------------------------------------------*/
      case end_of_int:                                        /* end_of_int */
        filewrite(operation,0,loc_a,val_a,0,cst_d);
        get_loc_block_f();
        break;

/*--------------------------------------------------------------------------*/
      case end_of_val:                                        /* end_of_val */
        filewrite(operation,0,loc_a,val_a,0,cst_d);
        get_val_block_f();
        break;

/*--------------------------------------------------------------------------*/
      case start_of_tape:                                  /* start_of_tape */
        filewrite_start(operation);
        break;

/*--------------------------------------------------------------------------*/
      case end_of_tape:                                      /* end_of_tape */ 
      	break;


/****************************************************************************/
/*                                                               COMPARISON */

/*--------------------------------------------------------------------------*/
      case eq_zero  :                                            /* eq_zero */
      case neq_zero :                                           /* neq_zero */
      case le_zero  :                                            /* le_zero */
      case gt_zero  :                                            /* gt_zero */
      case ge_zero  :                                            /* ge_zero */
      case lt_zero  :                                            /* lt_zero */
        arg  = get_locint_f();
        loc_a[0] = arg;
#ifdef computenumbers
        val_a[0] = T0[arg];
#endif 
        filewrite(operation,1,loc_a,val_a,0,cst_d);       
        break;


/****************************************************************************/
/*                                                              ASSIGNMENTS */

/*--------------------------------------------------------------------------*/
      case assign_a:           /* assign an adouble variable an    assign_a */
	                       /* adouble value. (=) */
        arg = get_locint_f();
        res = get_locint_f();
        loc_a[0]=arg;
        loc_a[1]=res;
#ifdef computenumbers
        val_a[0]=T0[arg];
        T0[res]= T0[arg];
        val_a[1]=T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case assign_d:            /* assign an adouble variable a    assign_d */
	                        /* double value. (=) */
        res  = get_locint_f();
        cst_d[0]=get_val_f();
        loc_a[0]=res;
#ifdef computenumbers
        T0[res]= cst_d[0];
        val_a[0]=T0[res];
#endif
        filewrite(operation,1,loc_a,val_a,1,cst_d);
	break;  

/*--------------------------------------------------------------------------*/
      case assign_d_one:    /* assign an adouble variable a    assign_d_one */
	                    /* double value. (1) (=) */
        res  = get_locint_f();
        loc_a[0]=res;
#ifdef computenumbers
        T0[res]= 1.0;
        val_a[0]=T0[res];
#endif
        filewrite(operation,1,loc_a,val_a,0,cst_d);
	break;  

/*--------------------------------------------------------------------------*/
      case assign_d_zero:  /* assign an adouble variable a    assign_d_zero */
	                   /* double value. (0) (=) */
        res  = get_locint_f();
        loc_a[0]=res;
#ifdef computenumbers
        T0[res]= 0.0;
        val_a[0]=T0[res];
#endif
        filewrite(operation,1,loc_a,val_a,0,cst_d);
	break;  

/*--------------------------------------------------------------------------*/
      case assign_ind:       /* assign an adouble variable an    assign_ind */
	                     /* independent double value (<<=) */
        res  = get_locint_f();
        loc_a[0]=res;
#ifdef computenumbers
        T0[res]= basepoint[indexi];
        cst_d[0]= basepoint[indexi];
        val_a[0]=T0[res];
        filewrite(operation,1,loc_a,val_a,1,cst_d);
#else
        filewrite(operation,1,loc_a,val_a,0,cst_d);
#endif
        indexi++;
	break;   

/*--------------------------------------------------------------------------*/
      case assign_dep:           /* assign a float variable a    assign_dep */
	                         /* dependent adouble value. (>>=) */
     	res = get_locint_f();
        loc_a[0]=res;
#ifdef computenumbers
        val_a[0]=T0[res];
        valuepoint[indexd++]=T0[res];
#endif
        filewrite(operation,1,loc_a,val_a,0,cst_d);
	break;      


/****************************************************************************/
/*                                                   OPERATION + ASSIGNMENT */

/*--------------------------------------------------------------------------*/
      case eq_plus_d:            /* Add a floating point to an    eq_plus_d */
	                         /* adouble. (+=) */
    	res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        T0[res] += coval;
        val_a[0] = T0[res];
#endif
        filewrite(operation,1,loc_a,val_a,1,cst_d);
	break;  

/*--------------------------------------------------------------------------*/
      case eq_plus_a:             /* Add an adouble to another    eq_plus_a */
	                          /* adouble. (+=) */
        arg  = get_locint_f();
        res  = get_locint_f();
        loc_a[0]=arg;
        loc_a[1]=res;
#ifdef computenumbers
        val_a[0]=T0[arg];
        T0[res]+= T0[arg];
        val_a[1]=T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,0,cst_d);
       	break; 

/*--------------------------------------------------------------------------*/
      case eq_plus_prod:    /* Add an product to an            eq_plus_prod */
	                    /* adouble. (+= x1*x2) */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        val_a[1]=T0[arg2];
        T0[res] += T0[arg1]*T0[arg2];
        val_a[2]=T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,0,cst_d);
       	break; 

/*--------------------------------------------------------------------------*/
      case eq_min_d:       /* Subtract a floating point from an    eq_min_d */
                           /* adouble. (-=) */
    	res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        T0[res] -= coval;
        val_a[0] = T0[res];
#endif
        filewrite(operation,1,loc_a,val_a,1,cst_d);
	break; 

/*--------------------------------------------------------------------------*/
      case eq_min_a:        /* Subtract an adouble from another    eq_min_a */
	                    /* adouble. (-=) */
        arg  = get_locint_f();
        res  = get_locint_f();
        loc_a[0]=arg;
        loc_a[1]=res;
#ifdef computenumbers
        val_a[0]=T0[arg];
        T0[res]-= T0[arg];
        val_a[1]=T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,0,cst_d);
       	break; 

/*--------------------------------------------------------------------------*/
      case eq_min_prod:     /* Subtract an product from an      eq_min_prod */
	                    /* adouble. (+= x1*x2) */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        val_a[1]=T0[arg2];
        T0[res] -= T0[arg1]*T0[arg2];
        val_a[2]=T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,0,cst_d);
       	break; 

/*--------------------------------------------------------------------------*/
      case eq_mult_d:              /* Multiply an adouble by a    eq_mult_d */
	                           /* flaoting point. (*=) */
    	res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        T0[res] *= coval;
        val_a[0] = T0[res];
#endif
        filewrite(operation,1,loc_a,val_a,1,cst_d);
	break; 

/*--------------------------------------------------------------------------*/
      case eq_mult_a:       /* Multiply one adouble by another    eq_mult_a */
	                    /* (*=) */
        arg  = get_locint_f();
        res  = get_locint_f();
        loc_a[0]=arg;
        loc_a[1]=res;
#ifdef computenumbers
        val_a[0]=T0[arg];
        T0[res]*= T0[arg];
        val_a[1]=T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,0,cst_d);
       	break; 

/*--------------------------------------------------------------------------*/
      case incr_a:                        /* Increment an adouble    incr_a */
    	res = get_locint_f();
        loc_a[0] = res;
#ifdef computenumbers
        T0[res]++;
        val_a[0] = T0[res];
#endif
        filewrite(operation,1,loc_a,val_a,0,cst_d);
	break;  

/*--------------------------------------------------------------------------*/
      case decr_a:                        /* Increment an adouble    decr_a */
    	res = get_locint_f();
        loc_a[0] = res;
#ifdef computenumbers
        T0[res]--;
        val_a[0] = T0[res];
#endif
        filewrite(operation,1,loc_a,val_a,0,cst_d);
	break;  


/****************************************************************************/
/*                                                        BINARY OPERATIONS */

/*--------------------------------------------------------------------------*/
      case plus_a_a:                 /* : Add two adoubles. (+)    plus a_a */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        val_a[1]=T0[arg2];
        T0[res]=T0[arg1]+T0[arg2];
        val_a[2]=T0[res];
#endif      
        filewrite(operation,3,loc_a,val_a,0,cst_d);
	break;  

/*--------------------------------------------------------------------------*/
      case plus_d_a:             /* Add an adouble and a double    plus_d_a */
	                         /* (+) */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = arg;
        loc_a[1] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        val_a[0]=T0[arg];
        T0[res]= T0[arg] + coval;
        val_a[1]=T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,1,cst_d);
       	break; 

/*--------------------------------------------------------------------------*/
      case min_a_a:              /* Subtraction of two adoubles     min_a_a */
	                         /* (-) */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        val_a[1]=T0[arg2];
        T0[res]=T0[arg1]-T0[arg2];
        val_a[2]=T0[res];
#endif      
        filewrite(operation,3,loc_a,val_a,0,cst_d);
	break;  

/*--------------------------------------------------------------------------*/
      case min_d_a:                /* Subtract an adouble from a    min_d_a */
	                           /* double (-) */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = arg;
        loc_a[1] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        val_a[0] = T0[arg];
        T0[res]  = coval - T0[arg];
        val_a[1] = T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,1,cst_d);
       	break; 

/*--------------------------------------------------------------------------*/
      case mult_a_a:               /* Multiply two adoubles (*)    mult_a_a */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        val_a[1]=T0[arg2];
        T0[res]=T0[arg1]*T0[arg2];
        val_a[2]=T0[res];
#endif      
        filewrite(operation,3,loc_a,val_a,0,cst_d);
	break;  

/*--------------------------------------------------------------------------*/
      case mult_d_a:         /* Multiply an adouble by a double    mult_d_a */
                             /* (*) */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = arg;
        loc_a[1] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        val_a[0] = T0[arg];
        T0[res]  = coval * T0[arg];
        val_a[1] = T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,1,cst_d);
       	break; 

/*--------------------------------------------------------------------------*/
      case div_a_a:           /* Divide an adouble by an adouble    div_a_a */
                              /* (/) */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        val_a[1]=T0[arg2];
        T0[res]=T0[arg1]/T0[arg2];
        val_a[2]=T0[res];
#endif      
        filewrite(operation,3,loc_a,val_a,0,cst_d);
	break;  

/*--------------------------------------------------------------------------*/
      case div_d_a:             /* Division double - adouble (/)    div_d_a */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = arg;
        loc_a[1] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        val_a[0] = T0[arg];
        T0[res]  = coval / T0[arg];
        val_a[1] = T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,1,cst_d);
       	break; 


/****************************************************************************/
/*                                                         SIGN  OPERATIONS */

/*--------------------------------------------------------------------------*/
      case pos_sign_a:                                        /* pos_sign_a */
        arg  = get_locint_f();
        res  = get_locint_f();
        loc_a[0]=arg;
        loc_a[1]=res;
#ifdef computenumbers
        val_a[0]=T0[arg];
        T0[res]= T0[arg];
        val_a[1]=T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,0,cst_d);
       	break; 

/*--------------------------------------------------------------------------*/
      case neg_sign_a:                                        /* neg_sign_a */
        arg  = get_locint_f();
        res  = get_locint_f();
        loc_a[0]=arg;
        loc_a[1]=res;
#ifdef computenumbers
        val_a[0]=T0[arg];
        T0[res]= -T0[arg];
        val_a[1]=T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,0,cst_d);
       	break; 


/****************************************************************************/
/*                                                         UNARY OPERATIONS */

/*--------------------------------------------------------------------------*/
      case exp_op:                          /* exponent operation    exp_op */
        arg  = get_locint_f();
        res  = get_locint_f(); 
        loc_a[0]=arg;
        loc_a[1]=res;
#ifdef computenumbers
        val_a[0]=T0[arg]; 
        T0[res]= exp(T0[arg]);
        val_a[1]=T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case sin_op:                              /* sine operation    sin_op */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        /* olvo 980923 changed order to allow x=sin(x) */
        val_a[0]=T0[arg1];
        T0[arg2]= cos(T0[arg1]);
        T0[res] = sin(T0[arg1]);
        val_a[1]=T0[arg2];
        val_a[2]=T0[res];
#endif      
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break; 

/*--------------------------------------------------------------------------*/
      case cos_op:                            /* cosine operation    cos_op */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        /* olvo 980923 changed order to allow x=cos(x) */
        val_a[0]=T0[arg1];
        T0[arg2]= sin(T0[arg1]);
        T0[res] = cos(T0[arg1]);
        val_a[1]=T0[arg2];
        val_a[2]=T0[res];
#endif      
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break; 

/*--------------------------------------------------------------------------*/
      case atan_op:                                              /* atan_op */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        T0[res] = atan(T0[arg1]);
        val_a[1]=T0[arg2];
        val_a[2]=T0[res];
#endif      
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break; 

/*--------------------------------------------------------------------------*/
      case asin_op:                                              /* asin_op */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        T0[res] = asin(T0[arg1]);
        val_a[1]=T0[arg2];
        val_a[2]=T0[res];
#endif      
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break; 

/*--------------------------------------------------------------------------*/
      case acos_op:                                              /* acos_op */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        T0[res] = acos(T0[arg1]);
        val_a[1]=T0[arg2];
        val_a[2]=T0[res];
#endif      
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break; 

#ifdef ATRIG_ERF

/*--------------------------------------------------------------------------*/
      case asinh_op:                                            /* asinh_op */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        T0[res] = asinh(T0[arg1]);
        val_a[1]=T0[arg2];
        val_a[2]=T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
       case acosh_op:                                           /* acosh_op */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        T0[res] = acosh(T0[arg1]);
        val_a[1]=T0[arg2];
        val_a[2]=T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case atanh_op:                                            /* atanh_op */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        T0[res] = atanh(T0[arg1]);
        val_a[1]=T0[arg2];
        val_a[2]=T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case erf_op:                                                /* erf_op */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        res  = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
        T0[res] = erf(T0[arg1]);
        val_a[1]=T0[arg2];
        val_a[2]=T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break;

#endif
/*--------------------------------------------------------------------------*/
      case log_op:                                                /* log_op */
        arg  = get_locint_f();
        res  = get_locint_f(); 
        loc_a[0]=arg;
        loc_a[1]=res;
#ifdef computenumbers
        val_a[0]=T0[arg]; 
        T0[res]= log(T0[arg]);
        val_a[1]=T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case pow_op:                                                /* pow_op */
        arg  = get_locint_f();
        res  = get_locint_f(); 
        coval   = get_val_f();
        cst_d[0]=coval;
        loc_a[0]=arg;
        loc_a[1]=res;
#ifdef computenumbers
        val_a[0]=T0[arg]; 
        T0[res] = pow(T0[arg],coval);
        val_a[1]=T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case sqrt_op:                                              /* sqrt_op */
        arg  = get_locint_f();
        res  = get_locint_f(); 
        loc_a[0]=arg;
        loc_a[1]=res;
#ifdef computenumbers
        val_a[0]=T0[arg]; 
        T0[res]= sqrt(T0[arg]);
        val_a[1]=T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case gen_quad:                                            /* gen_quad */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();		
        cst_d[0] = get_val_f();
        cst_d[1] = get_val_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        loc_a[2]=res;
#ifdef computenumbers
        val_a[0]=T0[arg1];
    	T0[res] = cst_d[1];
        val_a[1]=T0[arg2];
        val_a[2]=T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,2,cst_d);
        break;     

/*--------------------------------------------------------------------------*/
      case min_op:                                                /* min_op */
        arg1  = get_locint_f();
        arg2  = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = arg1;
        loc_a[1] = arg2;
        loc_a[2] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        val_a[0] = T0[arg1];
        val_a[1] = T0[arg2];
        if (T0[arg1] > T0[arg2])
           T0[res] = T0[arg2];
        else 
           T0[res] = T0[arg1];
        val_a[2] = T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case abs_val:                                              /* abs_val */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = arg;
        loc_a[1] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        val_a[0] = T0[arg];
        T0[res]  = fabs(T0[arg]);
        val_a[1] = T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,1,cst_d);
        break; 

/*--------------------------------------------------------------------------*/
      case ceil_op:                                              /* ceil_op */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = arg;
        loc_a[1] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        val_a[0] = T0[arg];
        T0[res]  = ceil(T0[arg]);
        val_a[1] = T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case floor_op:                 /* Compute ceil of adouble    floor_op */
        arg   = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = arg;
        loc_a[1] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        val_a[0] = T0[arg];
        T0[res]  = floor(T0[arg]);
        val_a[1] = T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,1,cst_d);
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
        loc_a[0]=arg;
        loc_a[1]=arg1;
        loc_a[2]=arg2 ;
        loc_a[3]=res;
        cst_d[0]=coval;
#ifdef computenumbers
        val_a[0]=T0[arg];
        val_a[1]=T0[arg1];
        val_a[2]=T0[arg2];
        if (T0[arg]>0)
          T0[res]=T0[arg1];
        else
          T0[res]=T0[arg2];
        val_a[3]=T0[res];
#endif
	filewrite(operation,4,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case cond_assign_s:                                  /* cond_assign_s */
        arg   = get_locint_f();
        arg1  = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();
        loc_a[0]=arg;
        loc_a[1]=arg1;
        loc_a[2]=res;
        cst_d[0]=coval;
#ifdef computenumbers
        val_a[0]=T0[arg];
        val_a[1]=T0[arg1];
        if (T0[arg]>0)
          T0[res]=T0[arg1];
        val_a[2]=T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,1,cst_d);
        break;


/****************************************************************************/
/*                                                       VECTOR ASSIGNMENTS */

/*--------------------------------------------------------------------------*/
      case assign_av:                                          /* assign_av */
        arg  = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();
        loc_a[0] = arg;
        loc_a[1] = size;
        loc_a[2] = res;
#ifdef computenumbers
        val_a[0] = T0[arg];
        val_a[1] = make_nan();
        for (l=0; l<size; l++)
          T0[res + l] = T0[arg + l];
        val_a[2] = T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break;
/*--------------------------------------------------------------------------*/
      case assign_dv:                                          /* assign_dv */
        size = get_locint_f();
        res  = get_locint_f();
        d    = get_val_v_f(size);
        loc_a[0] = size;
        loc_a[1] = res;
        cst_d[0] = d[0];
#ifdef computenumbers
        for (l=0; l<size; l++)
          T0[res + l] = d[l];  
        val_a[0] = make_nan();
        val_a[1] = T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case assign_indvec:                                  /* assign_indvec */
        size = get_locint_f();
        res  = get_locint_f();
        loc_a[0] = size;
        loc_a[1] = res;
#ifdef computenumbers
        for (l=0; l<size; l++) 
        { T0[res + l] = basepoint[indexi];
          ++indexi;
        }
        val_a[0] = make_nan();
        val_a[1] = T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case assign_depvec:                                  /* assign_depvec */
        size = get_locint_f();
        res  = get_locint_f();
        loc_a[0] = size;
        loc_a[1] = res;
#ifdef computenumbers
        for (l=0; l<size; l++)
        { valuepoint[indexd] = T0[res + l];
          indexd++;
        }
        val_a[0] = make_nan();
        val_a[1] = T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,0,cst_d);
        break;


/****************************************************************************/
/*                                            VECTOR OPERATION + ASSIGNMENT */

/*--------------------------------------------------------------------------*/
      case eq_plus_av:                                        /* eq_plus_av */
        arg  = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();
        loc_a[0] = arg;
        loc_a[1] = size;
        loc_a[2] = res;
#ifdef computenumbers
        val_a[0] = T0[arg];
        val_a[1] = make_nan();
        for (l=0; l<size; l++) 
          T0[res + l] += T0[arg + l];
        val_a[2] = T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case eq_min_av:                                          /* eq_min_av */
        arg  = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();
        loc_a[0] = arg;
        loc_a[1] = size;
        loc_a[2] = res;
#ifdef computenumbers
        val_a[0] = T0[arg];
        val_a[1] = make_nan();
        for (l=0; l<size; l++) 
          T0[res + l] -= T0[arg + l];
        val_a[2] = T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case eq_mult_av_d:                                    /* eq_mult_av_d */
        size  = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = size;
        loc_a[1] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        for (l=0; l<size; l++) 
          T0[res + l] *= coval;
        val_a[0] = make_nan();
        val_a[1] = T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case eq_mult_av_a:                                    /* eq_mult_av_a */
        arg  = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();
        loc_a[0] = arg;
        loc_a[1] = size;
        loc_a[2] = res;
#ifdef computenumbers
        val_a[0] = T0[arg];
        val_a[1] = make_nan();
        for (l=0; l<size; l++) 
          T0[res + l] *= val_a[0];
        val_a[2] = T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,0,cst_d);
        break;


/****************************************************************************/
/*                                                 BINARY VECTOR OPERATIONS */

/*--------------------------------------------------------------------------*/
      case plus_av_av:                                        /* plus_av_av */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();
        loc_a[0] = arg1;
        loc_a[1] = arg2;
        loc_a[2] = size;
        loc_a[3] = res;
#ifdef computenumbers
        val_a[0] = T0[arg1];
        val_a[1] = T0[arg2];
        val_a[2] = make_nan();
        for (l=0; l<size; l++) 
          T0[res+l] = T0[arg1+l] + T0[arg2+l];
        val_a[3] = T0[res];
#endif
        filewrite(operation,4,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case sub_av_av:                                          /* sub_av_av */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();
        loc_a[0] = arg1;
        loc_a[1] = arg2;
        loc_a[2] = size;
        loc_a[3] = res;
#ifdef computenumbers
        val_a[0] = T0[arg1];
        val_a[1] = T0[arg2];
        val_a[2] = make_nan();
        for (l=0; l<size; l++)
          T0[res+l] = T0[arg1+l] - T0[arg2+l];
        val_a[3] = T0[res];
#endif
        filewrite(operation,4,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case dot_av_av:                                          /* dot_av_av */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();
        loc_a[0] = arg1;
        loc_a[1] = arg2;
        loc_a[2] = size;
        loc_a[3] = res;
#ifdef computenumbers
        val_a[0] = T0[arg1];
        val_a[1] = T0[arg2];
        val_a[2] = make_nan();
        T0res=0;
        for (l=0; l<size; l++)
          T0res += T0[arg1+l] * T0[arg2+l];
        val_a[3] = T0[res] = T0res;
#endif
        filewrite(operation,4,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case mult_a_av:                                          /* mult_a_av */
        arg2 = get_locint_f();
        arg1 = get_locint_f();
        size = get_locint_f();
        res  = get_locint_f();
        loc_a[0] = arg2;
        loc_a[1] = arg1;
        loc_a[2] = size;
        loc_a[3] = res;
#ifdef computenumbers
        val_a[0] = T0[arg2];
        val_a[1] = T0[arg1];
        val_a[2] = make_nan();
        for (l=0; l<size; l++)
          T0[res+l] = val_a[1] * T0[arg2+l];
        val_a[3] = T0[res];
#endif
        filewrite(operation,4,loc_a,val_a,0,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case mult_d_av:                                          /* mult_d_av */
        arg   = get_locint_f();
        size  = get_locint_f();
        res   = get_locint_f();
        coval = get_val_f();
        loc_a[0] = arg;
        loc_a[1] = size;
        loc_a[2] = res; 
        cst_d[0] = coval;
#ifdef computenumbers
        val_a[0] = T0[arg];
        val_a[1] = make_nan();
        for (l=0; l<size; l++)
          T0[res+l] = T0[arg+l] * coval;
        val_a[2] = T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case div_av_a:                                            /* div_av_a */
        arg1   = get_locint_f();
        arg2   = get_locint_f();
        size   = get_locint_f();
        res    = get_locint_f();
        loc_a[0] = arg1;
        loc_a[1] = arg2;
        loc_a[2] = size;
        loc_a[3] = res;
#ifdef computenumbers
        val_a[0] = T0[arg1];
        val_a[1] = T0[arg2];
        val_a[2] = make_nan();
        for (l=0; l<size; l++)
          T0[res+l] = T0[arg1+l] / val_a[1];
        val_a[3] = T0[res];
#endif
        filewrite(operation,4,loc_a,val_a,0,cst_d);
        break;


/****************************************************************************/
/*                                                               SUBSCRIPTS */

/*--------------------------------------------------------------------------*/
      case subscript:                                          /* subscript */
        arg2 = get_locint_f(); /* Base */
        arg1 = get_locint_f(); /* pointer to variable containing offset */
        res  = get_locint_f();
        coval = get_val_f();
        loc_a[0] = arg2;
        loc_a[1] = arg1;
        loc_a[2] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        arg = arg2+(int)(T0[arg1]);
        val_a[0] = T0[arg];
        val_a[1] = T0[arg1];
        T0[res]  = T0[arg];
        val_a[2] = T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case subscript_l:                                      /* subscript_l */
        arg2 = get_locint_f();  /* Base */
	arg1 = get_locint_f();
        arg  = get_locint_f();         
        coval = get_val_f();
        loc_a[0] = arg2;
        loc_a[1] = arg1;
        loc_a[2] = arg;
        cst_d[0] = coval;
#ifdef computenumbers
        res = arg2 + (int)(T0[arg1]);
        val_a[1] = T0[arg1];
        val_a[2] = T0[arg];
        T0[res]  = T0[arg];
        val_a[0] = T0[res];
#endif
        filewrite(operation,3,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case subscript_ld:                                    /* subscript_ld */
        arg2 = get_locint_f(); /* Base */
        arg1 = get_locint_f(); /* pointer to variable containing offset */
        coval = get_val_f();
        loc_a[0] = arg2;
        loc_a[1] = arg1;
        cst_d[0] = coval;
        cst_d[1] = get_val_f();
#ifdef computenumbers
        arg = arg2 + (int)(T0[arg1]);
        val_a[1] = T0[arg1];
        T0[arg]  = coval;
        val_a[0] = T0[arg];
#endif
        filewrite(operation,2,loc_a,val_a,2,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case m_subscript:                                      /* m_subscript */
        arg2 = get_locint_f(); /* Base */
        arg1 = get_locint_f(); /* pointer to variable containing offset */
        size = get_locint_f();
        res  = get_locint_f();
        coval = get_val_f();
        loc_a[0] = arg2;
        loc_a[1] = arg1;
        loc_a[2] = size;
        loc_a[3] = res;
        cst_d[0] = coval;
#ifdef computenumbers
        arg = arg2 + (int)(T0[arg1])*size;
        val_a[0] = T0[arg];
        val_a[1] = T0[arg1];
        val_a[2] = make_nan();
        for (l=0; l<size; l++) 
          T0[res+l] = T0[arg+l];
        val_a[3] = T0[res];
#endif
        filewrite(operation,4,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case m_subscript_l:                                  /* m_subscript_l */
        arg2 = get_locint_f();  /* Base LHS */
	arg1 = get_locint_f();  /* Offset LHS */
        size = get_locint_f();
        arg  = get_locint_f();  /* RHS */
        coval = get_val_f();
        loc_a[0] = arg2;
        loc_a[1] = arg1;
        loc_a[2] = size;
        loc_a[3] = arg;
        cst_d[0] = coval;
#ifdef computenumbers
        res = arg2 + (int)(T0[arg1])*size;
        val_a[0] = T0[arg2];
        val_a[1] = T0[arg1];
        val_a[2] = make_nan();
        for (l=0; l<size; l++)
          T0[res+l] = T0[arg+l];
        val_a[3] = T0[res];
#endif
        filewrite(operation,4,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case m_subscript_ld:                                /* m_subscript_ld */
        arg2 = get_locint_f(); /* Base */
        arg1 = get_locint_f(); /* pointer to variable containing offset */
        arg  = get_locint_f(); /* offset in the vector itself */
        size = get_locint_f();
        coval = get_val_f();
        d      = get_val_v_f(size);
        loc_a[0] = arg2;
        loc_a[1] = arg1;
        loc_a[2] = arg;
        loc_a[3] = size;
        cst_d[0] = coval;
        cst_d[1] = d[0];
#ifdef computenumbers
        res = arg2 + (int)(T0[arg1])*size + arg;
        val_a[1] = T0[arg1];
        val_a[2] = make_nan();
        val_a[3] = make_nan();
        for (l=0; l<size; l++)
          T0[res+l] = d[l];
        val_a[0] = T0[res];
#endif
        filewrite(operation,4,loc_a,val_a,2,cst_d);
        break;


/****************************************************************************/
/*                                                          REMAINING STUFF */

/*--------------------------------------------------------------------------*/
      case take_stock_op:                                  /* take_stock_op */
        size = get_locint_f();
        res  = get_locint_f();
        d    = get_val_v_f(size);
        loc_a[0] = size;
        loc_a[1] = res;
        cst_d[0] = d[0];
#ifdef computenumbers
        for (l=0; l<size; l++)
          T0[res+l] = d[l];  
        val_a[0] = make_nan();
        val_a[1] = T0[res];
#endif
        filewrite(operation,2,loc_a,val_a,1,cst_d);
        break;

/*--------------------------------------------------------------------------*/
      case death_not:                                          /* death_not */
        arg1 = get_locint_f();
        arg2 = get_locint_f();
        loc_a[0]=arg1;
        loc_a[1]=arg2;
        filewrite(operation,2,loc_a,val_a,0,cst_d);
  	break; 

/*--------------------------------------------------------------------------*/
      default:                                                   /* default */
	/* Die here, we screwed up */
	fprintf(DIAG_OUT,"ADOL-C error: Fatal error in tape_doc for op %d\n",
                       operation);
	break;
	
      } /* endswitch */

      /* Read the next operation */
      operation=get_op_f();
    }  /* endwhile */

    if (operation == end_of_tape)
    {  filewrite_end(operation);
    };

  end_sweep();
} /* end tape_doc */


/****************************************************************************/
/*                                                               THAT'S ALL */

END_C_DECLS
