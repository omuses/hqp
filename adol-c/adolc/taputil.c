/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     taputil.c
 Revision: $Id: taputil.c,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: Initialization, stopage, and gets&puts of the taping process;
           as well as statistic gathering functions

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing

 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20040422 kowarz: adapted to configure - make - make install
          20030304 andrea: new global variable vs_file_name
          19991122 olvo:   new op_codes eq_plus_prod eq_min_prod
                           for  y += x1 * x2
                           and  y -= x1 * x2
                           --> new: upd_resloc_inc_prod(..)  
          19981130 olvo:   newly created by unification of taputil?.c
                           and all tape stuff
 History of taputil1.c:
          19981030 olvo:   bufsize --> BUFSIZE & TBUFSIZE
          19981019 olvo:   don't check sizeof(revreal)
          19980914 olvo:   unique size of stats block
          19980825 olvo:   #defines instead of const (C-Code!)
          19980820 olvo:   modification of statistic stuff
          19980723 olvo:   taputil3.* moved here
          19980713 olvo:   (1) no write_... routines anymore!
                           (2) statistic stuff kept here only
          19980709 olvo:   void write_pos_sign_a(..)
                           void write_neg_sign_a(..)
          19980708 olvo:   void write_upd(..)
          19980707 olvo:   void write_dot_av_av(..)
          19980706 olvo:   new operation code: incr_a
                                               decr_a
                           (void write_incr_decr_a(..) )
                                               int_adb_d_one
                                               int_adb_d_zero 
          19980703 olvo:   new operation code: assign_d_one
                                               assign_d_zero 
          19980623 olvo:   new operation code: take_stock_op

 History of taputil2.c:
          19980914 olvo:   unique size of stats block
          19980517 olvo:   griewank's idea:
                           int upd_resloc(locint, locint);

----------------------------------------------------------------------------*/
#include "taputil.h"
#include "taputil_p.h"
#include "tayutil.h"
#include "tayutil_p.h"
#include "oplate.h"

#include <malloc.h>
#include <string.h>
#include <errno.h>

BEGIN_C_DECLS

/****************************************************************************/
/*                                      GLOBAL VARIABLES (external linkage) */

/*--------------------------------------------------------------------------*/
/* Buffers for the operation tape, location tape, real tape. */
unsigned char *op_codes;
locint        *loc_tape;
double        *real_tape;

/*--------------------------------------------------------------------------*/
/* Pointers into the operation tape, location tape, real tape */
unsigned char *g_op_ptr;
locint        *g_loc_ptr;
double        *g_real_ptr;

int op_ptr;
int loc_ptr;
int real_ptr;

/*--------------------------------------------------------------------------*/
/* Statistic stuff */
int ind_ptr;
int dep_ptr;
int vs_ptr;
int revalso;

/****************************************************************************/
/*                                       LOCAL VARIABLES (internal linkage) */

/*--------------------------------------------------------------------------*/
/* Max number of tapes currently in use */ 
static int maxtapes = 0;

/*--------------------------------------------------------------------------*/
/* File Names */
static char op_file_name[20];
static char int_file_name[20];
static char val_file_name[20];

/*--------------------------------------------------------------------------*/
/* Arrays of pointers to the various tape buffers */
static unsigned char **op_tape;
static locint        **int_tape;
static double        **val_tape;

/*--------------------------------------------------------------------------*/
/* Array of pointers to the stats arrays (of size statSize) */
static int **stats;

/*--------------------------------------------------------------------------*/
static int tag;

/*--------------------------------------------------------------------------*/
/* Tape identification (ADOLC & version check) */
int adolcID[]    = { ADOLC_VERSION,
                     ADOLC_SUBVERSION,
                     ADOLC_PATCHLEVEL,
                     sizeof(locint),
                     sizeof(revreal)
                   }; 

/*--------------------------------------------------------------------------*/
/* File pointers to the operation tape, location tape, real tape */ 
static FILE *op_file_out;
static FILE *int_file_out;
static FILE *val_file_out;

/*--------------------------------------------------------------------------*/
/* Stats on operation tape, location tape, real tape */
static int op_access_ptr,
           int_access_ptr,
           val_access_ptr;

static int op_len_ptr,
           int_len_ptr,
           val_len_ptr;

/*--------------------------------------------------------------------------*/
/* Strings for the tape names (actual file names) */
static char *op_file,
            *int_file,
            *val_file;

/*--------------------------------------------------------------------------*/
/* File counts */
static long op_file_cnt,
            int_file_cnt,
            val_file_cnt;

/*--------------------------------------------------------------------------*/
/* Current buffer size */
static int buff_size;

/****************************************************************************/
/*                                                         GLOBAL VARIABLES */

/* File Name */
char vs_file_name[20];

/****************************************************************************/
/*                                                       INTERNAL FUNCTIONS */

/*--------------------------------------------------------------------------*/
void fail( int error )
{ switch (error) {
    case -1:
      fprintf(DIAG_OUT,"ADOL-C error: Malloc of memory failed!");
      exit (0);
  }
}

/*--------------------------------------------------------------------------*/
/* int2asc converts the integer num to a string, places it
   in the array string, and returns the pointer to the 
   string.  (I now that this is rather redundant, but I
   didn't write the original code for this.-- DWJ ;-)    */
char* int2asc( int num, char string[] )
{ sprintf(string,"%d",num);
  return(string);
}


/*--------------------------------------------------------------------------*/
/* The subroutine get_fstr appends to the filename fname             
   the number fnum, and puts the resulting    
   string in fstr.  */
void  get_fstr( char *fstr, short fnum, char *fname )
/**** 
  The caller of this function is responsible for allocating the appropriate 
  amount of storage for fstr [strlen(FNAME)+1 <= strlen(fstr) 
                                              <= strlen(FNAME)+5] 
****/
{ char tstr[10];

  if (fnum)
  { strcpy (fstr,fname);
    int2asc (fnum,tstr);
    strcat (fstr,tstr);
  }
  else
  { strcpy (fstr,fname);
    fstr[strlen(fstr)-1] = '\0';
  }
}

/****************************************************************************/
/*                                                           HELPFUL LOCALs */

/*--------------------------------------------------------------------------*/
static void init_stat_space( short tnum )
{ unsigned char **t1;          /* t1,t2,t3 and t4 are temporaries */
  double **t2;
  locint **t3;
  int    **t4;
  int jj;

  /* Set up space for */ 
  if (maxtapes == 0) /*this is only done at first call to start_trace or
                       init_stat_space */
  { maxtapes = 10;
    if (tnum >= maxtapes)
      maxtapes = tnum + 10;
    if ((op_tape = (unsigned char **)malloc(maxtapes*sizeof(unsigned char*))) == 0) 
      fail(-1);
    if ((int_tape = (locint **)malloc(maxtapes*sizeof(locint *))) == 0)
      fail(-1);
    if ((val_tape = (double **)malloc(maxtapes*sizeof(double *))) == 0)
      fail(-1);
      
    if ((stats = (int**)malloc(maxtapes*sizeof(int*))) == 0)
      fail(-1);
    for (jj=0; jj<maxtapes; jj++)
    { op_tape[jj]  = 0;
      int_tape[jj] = 0;
      val_tape[jj] = 0;
      stats[jj]    = 0;
    }
  }
  
  if (tnum >= maxtapes)
  { int newtapes = tnum + 10;
    t1 = op_tape;
    t3 = int_tape;
    t2 = val_tape;
    t4 = stats;
    if ((op_tape =(unsigned char**)malloc(newtapes*sizeof(unsigned char*))) == 0) 
      fail(-1);
    if ((int_tape = (locint **)malloc(newtapes*sizeof(locint *))) == 0)
      fail(-1);
    if ((val_tape = (double **)malloc(newtapes*sizeof(double *))) == 0)
      fail(-1);
    if ((stats = (int**)malloc(newtapes*sizeof(int*))) == 0)
      fail(-1);
      
    for (jj=0; jj<maxtapes; jj++)
    { op_tape[jj]  = t1[jj];
      int_tape[jj] = t3[jj];
      val_tape[jj] = t2[jj];
      stats[jj]    = t4[jj];
    }
    free((char *)t1);free((char *)t2);free((char *)t3);free((char *)t4);

    for(jj=maxtapes; jj<newtapes; jj++)
    { op_tape[jj]  = 0;
      int_tape[jj] = 0;
      val_tape[jj] = 0;
      stats[jj]    = 0;
    }
    maxtapes = newtapes;
  }
}

/*--------------------------------------------------------------------------*/
static void set_up_buffers( short tag, int buffer_size )
{ /* Return old memory ... if used */
  if (op_tape[tag]) 
    free((char*)op_tape[tag]);
  if (int_tape[tag])
    free((char*)int_tape[tag]);
  if (val_tape[tag])
    free((char*)val_tape[tag]);
  if (stats[tag])
    free((char*)stats[tag]);
  
  op_tape[tag]  = (unsigned char *)malloc(buffer_size*sizeof(unsigned char));
  int_tape[tag] = (locint *)malloc(buffer_size*sizeof(locint));
  val_tape[tag] = (double *)malloc(buffer_size*sizeof(double));
  stats[tag]    = (int*)malloc(statSize*sizeof(int));
  if ((op_tape[tag] == NULL) || (int_tape[tag]==NULL) 
      || (val_tape[tag] == NULL) || (stats[tag] == NULL)) 
  { fprintf(DIAG_OUT,"ADOL-C error: cannot allocate tape buffers!\n");
    exit (-1);
  }
}

/*--------------------------------------------------------------------------*/
static void read_tape_stats( short tag, int *stats )
{ char int_file[20];
  FILE *int_file_in;
  int version[adolcIDSize];

  get_fstr(int_file,tag,FNAME1);
  
  if ((int_file_in = fopen(int_file,"rb")) == 0) 
  { fprintf(DIAG_OUT,"ADOL-C error: Error reading integer tape number %d \n",
                   tag);
    fprintf(DIAG_OUT,"Fopen returned error number %d \n",tag);
    exit(-1);
  }
  if (fread((char *)stats,statSize*sizeof(int),1,int_file_in) != 1)
  { fprintf(DIAG_OUT,"ADOL-C error: Error reading integer tape number %d \n",
                   tag);
    fprintf(DIAG_OUT,"Fread returned error number %d \n",tag);
    exit(-1);
  }
  /* olvo 980820 version check */ 
  if (fread((char *)version,adolcIDSize*sizeof(int),1,int_file_in) != 1)
  { fprintf(DIAG_OUT,"ADOL-C error: Error reading integer tape number %d \n",
                   tag);
    fprintf(DIAG_OUT,"Fread returned error number %d \n",tag);
    exit(-1);
  }
  if (   (version[0] < adolcID[0]) 
      || ((version[0] == adolcID[0]) &&  (version[1] < adolcID[1])))
  { fprintf(DIAG_OUT,"ADOL-C error: Used tape (%d) was written with ADOL-C"
                     " version older than %d.%d\n",
                   tag,ADOLC_VERSION, ADOLC_SUBVERSION);
    fprintf(DIAG_OUT,"This is ADOL-C %d.%d.%d\n", ADOLC_VERSION,
                      ADOLC_SUBVERSION, ADOLC_PATCHLEVEL);
    exit(-1); 
  }
  if (version[3] != adolcID[3]) 
  { fprintf(DIAG_OUT,"ADOL-C error: Used tape (%d) was written with locints"
                     " of size %d, size %d required.\n",
                   tag, version[3], adolcID[3]);
    exit(-1); 
  }
  /* 981019 olvo: don't check sizeof(revreal)
  if (version[4] != adolcID[4]) 
  { fprintf(DIAG_OUT,"ADOL-C error: Used tape (%d) was written with revreals"
                     " of size %d, size %d required.\n",
                   tag, version[4], adolcID[4]);
    exit(-1); 
  } */
  fclose(int_file_in);
}

/****************************************************************************/
/*                                                          STATS functions */

/*--------------------------------------------------------------------------*/
/* Tapestats:                                                               */
/* Returns statistics on the tape tag.  The array tape_stat is assumed to   */
/* contain at least 11 elements.  The elements of the array are the         */
/* following.                                                               */
/* tape_stat[0] = # of independent variables.                               */
/* tape_stat[1] = # of dependent variables.                                 */
/* tape_stat[2] = max # of live variables.                                  */
/* tape_stat[3] = value stack size.                                         */
/* tape_stat[4] = buffer size (# of chars, # of doubles, # of locints)      */
/* tape_stat[5] = # of operations.                                          */ 
/* tape_stat[6] = operation file access flag (1 = file in use, 0 otherwise) */
/* tape_stat[7] = # of saved locations.                                     */ 
/* tape_stat[8] = location file access flag (1 = file in use, 0 otherwise)  */
/* tape_stat[9] = # of saved constant values.                               */ 
/* tape_stat[10]= value file access flag (1 = file in use, 0 otherwise)     */
/*                                                                          */
/*--------------------------------------------------------------------------*/
void tapestats( short tag,int *tape_stat )
{ int i;

  /* Make sure that there is tape access */
  init_stat_space(tag);

  if (stats[tag] == 0) 
  { /* Tape number does not exist , so read in tape data */
    read_tape_stats(tag,tape_stat);
    set_up_buffers(tag,tape_stat[4]);

    /* Copy data to stats for future use */

    for (i=0; i<statSize; i++)
      stats[tag][i] = tape_stat[i];
  }
  else
    for (i=0; i<statSize; i++)
      tape_stat[i] = stats[tag][i];
}
 
/*--------------------------------------------------------------------------*/
void get_op_stats( int tag, char **ret_op_file, int *ret_op_len, 
                   int *ret_op_access, unsigned char **ret_op_tape )
{ get_fstr(op_file_name,tag,FNAME);
  *ret_op_file   = op_file_name;
  *ret_op_len    = stats[tag][5];
  *ret_op_access = stats[tag][6];
  *ret_op_tape   = op_tape[tag];
}

/*--------------------------------------------------------------------------*/
void get_int_stats( int tag, char **ret_int_file, int *ret_int_len, 
                    int *ret_int_access, locint **ret_int_tape )
{ get_fstr(int_file_name,tag,FNAME1);
  *ret_int_file   = int_file_name;
  *ret_int_len    = stats[tag][7];
  *ret_int_access = stats[tag][8];
  *ret_int_tape   = int_tape[tag];
}

/*--------------------------------------------------------------------------*/
void get_val_stats( int tag, char **ret_val_file, int *ret_val_len, 
                    int *ret_val_access,double **ret_val_tape)
{ get_fstr(val_file_name,tag,FNAME2);
  *ret_val_file   = val_file_name;
  *ret_val_len    = stats[tag][9];
  *ret_val_access = stats[tag][10];
  *ret_val_tape   = val_tape[tag];
}


/****************************************************************************/
/*                                                                  TRACING */

/****************************************************************************/
/* start_trace: (part of trace_on)                                          */
/* Initialization for the taping process.  Sets up the arrays op_tape,      */
/* int_tape, val_tape, and stats.  Op_tape, int_tape, val_tape are arrays   */
/* of pointers to individual buffers for operations, integers (locints),    */
/* and values (doubles).  Also initializes buffers for this tape, sets      */
/* files names, and calls appropriate setup routines.                       */
/****************************************************************************/
void start_trace( short tnum, int revals )
{ int kk;
  double** dum = 0;
  int degree = 0;
  
  revalso = revals;
  tag     = tnum;
  
  /* Set buffer size to be the default in usrparms.h */
  set_buf_size(BUFSIZE);

  get_fstr(op_file_name,tag,FNAME);
  get_fstr(int_file_name,tag,FNAME1);
  get_fstr(val_file_name,tag,FNAME2);
  get_fstr(vs_file_name,tag,FNAME3);
  
  init_stat_space(tag);

  /* Return old memory ... if used */
  if(op_tape[tag])  free((char*)op_tape[tag]);
  if(int_tape[tag]) free((char*)int_tape[tag]);
  if(val_tape[tag]) free((char*)val_tape[tag]);
  if(stats[tag])    free((char*)stats[tag]);
  
  op_tape[tag]  = (unsigned char *)malloc(BUFSIZE*sizeof(unsigned char));
  int_tape[tag] = (locint *)malloc(BUFSIZE*sizeof(locint));
  val_tape[tag] = (double *)malloc(BUFSIZE*sizeof(double));
  stats[tag]    = (int*)malloc(statSize*sizeof(int));
  if ((op_tape[tag]  == NULL) || 
      (int_tape[tag] == NULL) || 
      (val_tape[tag] == NULL) ||
      (stats[tag]    == NULL)) 
  { fprintf(DIAG_OUT,"ADOL-C error: cannot allocate tape buffers!\n");
    exit (-1);
  }
  
  ind_ptr   = 0;
  dep_ptr   = 0;
  vs_ptr = 0;      
 
  /* Initialize Tapes */
  set_buffers(op_file_name,op_tape[tag],
	      int_file_name,int_tape[tag],
	      val_file_name,val_tape[tag]);
  
  /* Put operation denoting the start_of_the tape */ 
  put_op(start_of_tape);
  /* Leave space for the stats */
  /* olvo 980914 unique size of stats block */
  for (kk=0; kk<sizeof(int)/sizeof(locint)*statSpace; kk++)
    put_locint(0);
   
  if (revalso)
    taylor_begin(tag,TBUFSIZE,dum,degree);
}


/*************************************************************************/
/* Stop Tracing.  Clean up, and turn off trace_flag.                    **/
/*************************************************************************/
void stop_trace(int locations, int flag)
{ int tape_stats[statSize];
  int i,sizer;

  int loc_ptr;
  /* int op_cnt,loc_cnt,val_cnt,access_ptr; */

  loc_ptr = locations;     
  put_op(end_of_tape);        /* Mark end of tape. */

  tape_stats[0] = ind_ptr;
  tape_stats[1] = dep_ptr;
  tape_stats[2] = loc_ptr;
  tape_stats[3] = vs_ptr;
  tape_stats[4] = BUFSIZE;
  close_tape(tape_stats,flag); /** closes the tape, files up stats, and
                                   writes the tape stats to the integer
                                   tape. **/
  
  sizer = sizeof(revreal)*(vs_ptr);
  if (revalso) 
    taylor_close(sizer/(1+sizer/TBUFSIZE),dep_ptr,ind_ptr);
  for (i=0; i<statSize; i++)
    stats[tag][i]=tape_stats[i];
}


/****************************************************************************/
/*                                                          DEBUG FUNCTIONS */
#ifdef HARDDEBUG

/*--------------------------------------------------------------------------*/
unsigned char get_op_f( void )
{ unsigned char temp;
  temp= *g_op_ptr++;
  fprintf(DIAG_OUT,"f_op: %i\n",temp-'\0');
  return temp;
}

/*--------------------------------------------------------------------------*/
unsigned char get_op_r( void )
{ unsigned char temp;
  temp= *(--g_op_ptr);
  fprintf(DIAG_OUT,"r_op: %i\n",temp-'\0');
  return temp;
}

/*--------------------------------------------------------------------------*/
locint get_locint_f( void )
{ locint temp;
  temp= *g_loc_ptr++;
  fprintf(DIAG_OUT,"f_loc: %i\n",temp);
  return temp;
}

/*--------------------------------------------------------------------------*/
locint get_locint_r( void )
{ unsigned char temp;
  temp= *(--g_loc_ptr);
  fprintf(DIAG_OUT,"r_loc: %i\n",temp);
  return temp;
}

/*--------------------------------------------------------------------------*/
double get_val_f( void )
{ double temp;
  temp= *g_real_ptr++;
  fprintf(DIAG_OUT,"f_val: %e\n",temp);
  return temp;
}

/*--------------------------------------------------------------------------*/
double get_val_r( void )
{ double temp;
  temp= *(--g_real_ptr);
  fprintf(DIAG_OUT,"r_val: %e\n",temp);
  return temp;
}

#endif


/****************************************************************************/
/*                                                               LOCAL PUTs */
/****************************************************************/
/** Put_Block puts a block of tape to the disk.  I assume this **/
/** is called only during a first forward pass or during the   **/
/** the taping itself. Its purpose is to record all of the     **/
/** computations.                                              **/
/****************************************************************/
void put_op_block( int buffer_size )
{ int n;
  if (op_access_ptr == 0)
  { op_file_out = fopen(op_file,"rb");
    if (op_file_out != 0)
    {
#ifdef DEBUG
      fprintf(DIAG_OUT,"ADOL-C debug: old tapefile %s exists and deleted\n",op_file);
#endif
      fclose(op_file_out);
      if (remove(op_file))  /*  Complies with ANSI C standard */
      /* if(unlink(op_file))      works on some UNIX systems */
	fprintf(DIAG_OUT,"ADOL-C error: unable to remove old tapefile\n");
      op_file_out = fopen(op_file,"wb");
    }
    else
    { op_file_out = fopen(op_file,"wb");
      errno =0; /* Clear Out the Error */
    }
    op_access_ptr = 1;
  }
  op_len_ptr += buffer_size;

  if ((n = fwrite((char *)op_codes,buffer_size,1,op_file_out)) !=1 )
  { fprintf(DIAG_OUT,"ADOL-C error: Fatal error-doing a write %d--- error %d\n",n,errno);
    switch (errno) {
      case 28: /* ENOSPC */
	fprintf(DIAG_OUT,"No space left on device-contact sys. manager\n");
	break;
      case 27: /* EFBIG */
	fprintf(DIAG_OUT,"File too big-- tape space exhausted.\n");
	break;
      default:
	fprintf(DIAG_OUT,"Unexpected unix file error-- %d.\n",errno);
	break;
    }
    exit(-3);
  }
  op_ptr = 0;
  errno  = 0;
}

void put_locint_block( int buffer_size )
{ int n;
  if (int_access_ptr == 0)
  { int_file_out = fopen(int_file,"rb");
    if (int_file_out != 0)
    { 
#ifdef DEBUG
      fprintf(DIAG_OUT,"ADOL-C debug: old tapefile %s exists and deleted\n",int_file);
#endif
      fclose(int_file_out);
      if (remove(int_file))  /*    Complies with ANSI C standard */
      /* if(unlink(int_file))        works on some UNIX systems    */
	fprintf(DIAG_OUT,"ADOL-C error: unable to remove old tapefile\n");
      int_file_out = fopen(int_file,"wb");
    }
    else
    { int_file_out = fopen(int_file,"wb");
      errno =0; /* Clear Out the Error */
    }
    int_access_ptr = 1;
  }
  int_len_ptr += buffer_size;

  if ((n = fwrite((locint *)loc_tape,buffer_size*sizeof(locint),1,int_file_out)) != 1)
  { fprintf(DIAG_OUT,"ADOL-C error: Fatal error-doing a write %d--- error %d\n",n,errno);
    switch (errno) {
      case 28: /* ENOSPC */
	fprintf(DIAG_OUT,"No space left on device-contact sys. manager\n");
	break;
      case 27: /* EFBIG */
	fprintf(DIAG_OUT,"File too big-- tape space exhausted.\n");
	break;
      default:
	fprintf(DIAG_OUT,"Unexpected unix file error-- %d.\n",errno);
	break;
    }
    exit(-3);
  }
  loc_ptr = 0;
  errno   = 0;
}

void put_val_block( int buffer_size )
{ int n;
  if (val_access_ptr == 0)
  { val_file_out = fopen(val_file,"rb");
    if (val_file_out != 0)
    {
#ifdef DEBUG
      fprintf(DIAG_OUT,"ADOL-C debug: old tapefile %s exists and deleted\n",val_file);
#endif
      fclose(val_file_out);
      if (remove(val_file))   /* Complies with ANSI C standard */
      /* if(unlink(val_file))      works on some UNIX systems    */
	fprintf(DIAG_OUT,"ADOL-C error: unable to remove old tapefile\n");
      val_file_out = fopen(val_file,"wb");
    }
    else
    { val_file_out = fopen(val_file,"wb");
      errno =0; /* Clear Out the Error */
    }
    val_access_ptr = 1;
  }
  val_len_ptr += buffer_size;

  if ((n = fwrite((double *)real_tape,buffer_size*sizeof(double),1,val_file_out)) != 1)
  { fprintf(DIAG_OUT,"ADOL-C error: Fatal error-doing a write %d--- error %d\n",n,errno);
    switch (errno) {
      case 28: /* ENOSPC */
	fprintf(DIAG_OUT,"No space left on device-contact sys. manager\n");
	break;
      case 27: /* EFBIG */
	fprintf(DIAG_OUT,"File too big-- tape space exhausted.\n");
	break;
      default:
	fprintf(DIAG_OUT,"Unexpected unix file error-- %d.\n",errno);
	break;
    }
    exit(-3);
  }
  real_ptr = 0;
  errno    = 0;
}


/****************************************************************************/
/*                                        CONTROL STUFF (inits, ends, etc.) */

void init_for_sweep( int tag )
{ get_op_stats(tag,&op_file,&op_len_ptr,&op_access_ptr,&op_codes);
  if (op_access_ptr)
  { op_file_out = fopen(op_file,"rb");
    op_ptr = min(buff_size,op_len_ptr); 
    fread((char *)op_codes,op_ptr,1,op_file_out);
    op_len_ptr -= op_ptr;
  } 
  g_op_ptr = op_codes;

  get_int_stats(tag,&int_file,&int_len_ptr,&int_access_ptr,&loc_tape);
  if (int_access_ptr)
  { int_file_out = fopen(int_file,"rb");
    loc_ptr = min(buff_size,int_len_ptr);
    fread((locint *)loc_tape,sizeof(locint),loc_ptr,int_file_out);
    int_len_ptr -= loc_ptr;
  } 

  /* olvo 980914 unique size of stats block */
  g_loc_ptr = loc_tape+ sizeof(int)/sizeof(locint)*statSpace; 

  /* loc_tape = (loc_tape+loc_ptr); */
  get_val_stats(tag,&val_file,&val_len_ptr,&val_access_ptr,&real_tape);
  if (val_access_ptr)
  { val_file_out = fopen(val_file,"rb");
    real_ptr = min(val_len_ptr,buff_size);
    fread((char *)real_tape,real_ptr*sizeof(double),1,val_file_out);
    val_len_ptr -= real_ptr;
  } 
  g_real_ptr = real_tape;
}

void init_rev_sweep(int tag)
{ get_op_stats(tag,&op_file,&op_len_ptr,&op_access_ptr,&op_codes);
  if (op_access_ptr)
  { op_file_out = fopen(op_file,"rb");
    op_ptr = op_len_ptr % buff_size;
    fseek(op_file_out,0,2);
    op_file_cnt =  ftell(op_file_out);
    op_file_cnt -= op_ptr*sizeof(unsigned char);
    fseek(op_file_out,op_file_cnt,0);
    fread((char *)op_codes,op_ptr,1,op_file_out);
    op_file_cnt -= buff_size*sizeof(unsigned char);
    g_op_ptr = op_codes + op_ptr;
  } 
  else 
    g_op_ptr = op_codes + op_len_ptr;

  get_int_stats(tag,&int_file,&int_len_ptr,&int_access_ptr,&loc_tape);
  if (int_access_ptr)
  { int_file_out = fopen(int_file,"rb");
    loc_ptr = int_len_ptr % buff_size;
    fseek(int_file_out,0,2);
    int_file_cnt =  ftell(int_file_out);
    int_file_cnt -= loc_ptr*sizeof(locint);
    fseek(int_file_out,int_file_cnt,0);
    fread((char *)loc_tape,loc_ptr*sizeof(locint),1,int_file_out);
    int_file_cnt -= buff_size*sizeof(locint);
    g_loc_ptr = loc_tape + loc_ptr;
  } 
  else 
    g_loc_ptr = loc_tape + int_len_ptr;

  get_val_stats(tag,&val_file,&val_len_ptr,&val_access_ptr,&real_tape);
  if (val_access_ptr)
  { val_file_out = fopen(val_file,"rb");
    real_ptr = val_len_ptr % buff_size;
    fseek(val_file_out,0,2);
    val_file_cnt =  ftell(val_file_out);
    val_file_cnt -= real_ptr*sizeof(double);
    fseek(val_file_out,val_file_cnt,0);
    fread((char *)real_tape,real_ptr*sizeof(double),1,val_file_out);
    val_file_cnt -= buff_size*sizeof(double);
    g_real_ptr = real_tape + real_ptr;
  } 
  else 
    g_real_ptr = real_tape + val_len_ptr;
}

void set_buf_size( int size )
{ buff_size = size;
}

void set_buffers( char *file1, unsigned char *op_addr,
		  char *file2, locint *int_addr,
		  char *file3, double *real_addr )
{ op_codes  = op_addr;
  loc_tape  = int_addr;
  real_tape = real_addr;
  op_file  = file1;
  int_file = file2;
  val_file = file3;
  op_ptr        = loc_ptr        = real_ptr       = 0;
  op_access_ptr = int_access_ptr = val_access_ptr = 0;
  op_len_ptr    = int_len_ptr    = val_len_ptr    = 0;
}

void close_tape( int *stats, int flag )
{ int i;
  int access = (flag || op_access_ptr || int_access_ptr ||val_access_ptr);
  if (access)
  { if (op_ptr != 0)
      put_op_block(op_ptr);
    fclose(op_file_out);
  }
  else 
    op_len_ptr = op_ptr;
  stats[5] = op_len_ptr;
  stats[6] = op_access_ptr;

  if (access)
  { if (real_ptr != 0)
      put_val_block(real_ptr);
      if (val_file_out!=NULL) fclose(val_file_out);
  }
  else
    val_len_ptr = real_ptr;
  stats[9]  = val_len_ptr;
  stats[10] = val_access_ptr;

  if (access)
  { if (loc_ptr != 0)
      put_locint_block(loc_ptr);
    stats[7] = int_len_ptr;
    stats[8] = int_access_ptr;
    fseek(int_file_out,0,0);
    fwrite(stats,statSize*sizeof(int),1,int_file_out);
    /* olvo 980820 new: write ADOL-C version */
    fwrite(adolcID,adolcIDSize*sizeof(int),1,int_file_out); 
    fclose(int_file_out);
  }
  else
  { int_len_ptr = loc_ptr;
    stats[7] = int_len_ptr;
    stats[8] = int_access_ptr;
    for(i=0; i<statSize; i++)
      loc_tape[i] = stats[i];
  }
}

void end_sweep(void)
{ if (op_access_ptr)
    fclose(op_file_out); 
  if (int_access_ptr)
    fclose(int_file_out);
  if (val_access_ptr)
    fclose(val_file_out);
}


/****************************************************************************/
/*                                                                     PUTs */

/*--------------------------------------------------------------------------*/
/* Locations */
void put_locint( locint loc )
{ /*if (loc_ptr == buff_size) put_locint_block(buff_size); */
  loc_tape[loc_ptr++] = loc;
}

/*--------------------------------------------------------------------------*/
/* Operations */
void put_op( unsigned char op )
{ if (loc_ptr > buff_size-5) /* every operation writes <5 locations */
  { loc_tape[buff_size-1]=buff_size-loc_ptr;
    put_locint_block(buff_size);
    /* olvo 980720 old: put_to_op(end_of_int); */
    if (op_ptr == buff_size-1) /* every operation writes 1 opcode */
    { op_codes[op_ptr] = end_of_op;
      put_op_block(buff_size);
      op_codes[op_ptr++] = end_of_op;
    }
    op_codes[op_ptr++] = end_of_int;
  }
  if (real_ptr > buff_size-5) /* every operation writes <5 constants */
  {                           /*  3 should be sufficient */
    put_locint(buff_size-real_ptr);
    put_val_block(buff_size);
    /* olvo 980720 old: put_to_op(end_of_val); */
    if (op_ptr == buff_size-1) /* every operation writes 1 opcode */
    { op_codes[op_ptr] = end_of_op;
      put_op_block(buff_size);
      op_codes[op_ptr++] = end_of_op;
    }
    op_codes[op_ptr++] = end_of_val;
  }
  /* olvo 980720 old: put_to_op(op); */
  if (op_ptr == buff_size-1) /* every operation writes 1 opcode */
  { op_codes[op_ptr] = end_of_op;
    put_op_block(buff_size);
    op_codes[op_ptr++] = end_of_op;
  }
  op_codes[op_ptr++] = op;
}

/*--------------------------------------------------------------------------*/
/* Values */
void put_val( double r_val )
{ /* if (real_ptr == buff_size) put_val_block(buff_size); */
  real_tape[real_ptr++] = r_val;
}

void put_vals_p( double *r_val, int size )
{ int j;
  for (j=0; j<size; j++)
    real_tape[real_ptr++] = r_val[j];
  put_locint(buff_size-real_ptr);
  put_val_block(buff_size);
  /* olvo 980720 old: put_to_op(end_of_val); */
  if (op_ptr == buff_size-1) /* every operation writes 1 opcode */
  { op_codes[op_ptr] = end_of_op;
    put_op_block(buff_size);
    op_codes[op_ptr++] = end_of_op;
  }
  op_codes[op_ptr++] = end_of_val;
}

void put_vals_r( double *r_val, int size )
{ int j;
  for (j=0; j<size; j++)
    real_tape[real_ptr++] = r_val[j];
}

/*--------------------------------------------------------------------------*/
/* Update/correction of values or locations */
void reset_val_r( void )
{ if (g_real_ptr == real_tape)
    get_val_block_r();
}

int upd_resloc( locint temp, locint lhs )
{ int ret = 0;
  if (temp == loc_tape[loc_ptr-1])
  { loc_tape[loc_ptr-1]= lhs;
    ret = 1;
  }
  return ret;
}

/* olvo 991122: new routine */
int upd_resloc_inc_prod( locint temp, locint newlhs, unsigned char newop )
{ int ret = 0;
  if (   (temp == loc_tape[loc_ptr-1]) 
      && (mult_a_a == op_codes[op_ptr-1]) 
      && (newlhs != loc_tape[loc_ptr-2])   /* skipping recursive case */
      && (newlhs != loc_tape[loc_ptr-3]))
  { loc_tape[loc_ptr-1]= newlhs;
    op_codes[op_ptr-1] = newop;
    ret = 1;
  }
  return ret;
}


/****************************************************************************/
/*                                                                    GETs  */

/*--------------------------------------------------------------------------*/
/* Operations */
void get_op_block_f( void )
{ op_ptr = min(buff_size,op_len_ptr); 
  fread((char *)op_codes,op_ptr,1,op_file_out);
  op_len_ptr-= op_ptr;
  g_op_ptr = op_codes;
}

void get_op_block_r( void )
{ fseek(op_file_out,op_file_cnt,0);
  fread((char *)op_codes,buff_size,1,op_file_out);
  op_file_cnt -= buff_size*sizeof(unsigned char);
  g_op_ptr = op_codes + buff_size;
}

/*--------------------------------------------------------------------------*/
/* Locations */
void get_loc_block_f( void )
{ loc_ptr = min(buff_size,int_len_ptr);
  fread((char *)loc_tape,loc_ptr*sizeof(locint),1,int_file_out);
  int_len_ptr -= loc_ptr;
  g_loc_ptr = loc_tape;
}

void get_loc_block_r( void )
{ fseek(int_file_out,int_file_cnt,0);
  fread((char *)loc_tape,buff_size*sizeof(locint),1,int_file_out);
  int_file_cnt -= buff_size*sizeof(locint);
  g_loc_ptr = loc_tape + buff_size-loc_tape[buff_size-1];
}

/*--------------------------------------------------------------------------*/
/* Values */
int get_val_space( void )
{ if ((buff_size-real_ptr-5) < 0) 
  { put_locint(buff_size-real_ptr);
    put_val_block(buff_size);
    /* olvo980720 old: put_to_op(end_of_val); */
    if (op_ptr == buff_size-1) /* every operation writes 1 opcode */
    { op_codes[op_ptr] = end_of_op;
      put_op_block(buff_size);
      op_codes[op_ptr++] = end_of_op;
    }
    op_codes[op_ptr++] = end_of_val;
  } 
  return (buff_size - real_ptr - 5);
}

void get_val_block_f( void )
{ real_ptr = min(val_len_ptr,buff_size);
  fread((char *)real_tape,real_ptr*sizeof(double),1,val_file_out);
  val_len_ptr -= real_ptr;
  g_real_ptr = real_tape;
  g_loc_ptr++; /* get_locint_f(); value used in reverse only */
}
    
void get_val_block_r( void )
{ locint temp;
  fseek(val_file_out,val_file_cnt,0);
  fread((char *)real_tape,buff_size*sizeof(double),1,val_file_out);
  val_file_cnt -= buff_size*sizeof(double);
  temp = *(--g_loc_ptr);   /*get_locint_r();*/ 
  g_real_ptr = real_tape+buff_size-temp;
}

double *get_val_v_f( locint size )
{ double *temp = g_real_ptr;
  g_real_ptr += size;
  return temp;
}

double *get_val_v_r( locint size )
{ g_real_ptr -= size;
  return g_real_ptr;
}

/****************************************************************************/
/*                                                                    UTILs */
double make_nan( void )
{ double a,b;
#ifdef inf_num
  a = non_num;
  b = non_den;
#endif 
  return a/b;
} 

double make_inf( void )
{ double a,b;
#ifdef inf_num
  a = inf_num;
  b = inf_den;
#endif 
  return a/b;
} 

END_C_DECLS
