/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     drivers/drivers.h
 Revision: $Id: tayutil.c,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
 Contents: Taylor series utilities - primarily called from the module
           hos_forward.c (--- a forward pass generates taylor 
           coefficients which need to be saved when a variable dies, or 
           is overwritten) and the reverse modules (-- to retrieve these 
           taylor coefficients to calculate the adjoints.

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
   20040717 kowarz: "to many opened files" bug fixed (taylor_begin)
   20040607 kowarz: bug fixed in begin_taylor => ChangeLog
   20040417 kowarz: adapted to configure - make - make install
   20030317 andrea: pointer to current vs_data
   20030305 andrea: clean up for vs_data
   20030304 andrea: identify value stack by tag
   20010720 andrea: corrections is get_taylors_p(..) and write_taylors(..)
   20010715 andrea: add get_taylors_p(..)
   20010706 andrea: add write_taylors(..)
   19991122 olvo:   new op_codes eq_plus_prod eq_min_prod
                    for  y += x1 * x2
                    and  y -= x1 * x2
                    --> new: delete_scaylor(..)  
   19990816 olvo:   ec in get_taylors
   19990714 olvo:   performance tuning (get_taylors)
   19981201 olvo:   changed file name to tayutilc.c
   19980921 olvo:   new interface of void overwrite_scaylor(..) to
                    allow correction of old overwrite in store
   19980708 olvo:   (1) changed order in write_scaylor(..)
                        and              write_taylor(..)
                    (2) new:  void overwrite_scaylor(..)

----------------------------------------------------------------------------*/

#include "taputil.h"
#include "tayutil.h"
#include "taputil_p.h"
#include "tayutil_p.h"

#include <malloc.h>
#include <errno.h>

BEGIN_C_DECLS

typedef struct{
 int tag;
 int numdep;
 int numind;
 int T_file_access;
 FILE* temp2_file;
 char vs_file_name[20];
 int taylor_cnt;
 revreal * save_taylor;
 int T_write_cnt;
 int T_blocks, 
           T_tail, 
           T_buf_size, 
           T_length; 
 double  **T;
 revreal **Tr;
 revreal  *Trs;
 int degsave;
} data_array;
   

/****************************************************************************/
/*                                       LOCAL VARIABLES (internal linkage) */

data_array vs_data[TBUFNUM];
data_array *cur_vs_data;
static int num_vs_data=0;
static int index;

/****************************************************************************/
/*                                                          ACCESS ROUTINES */

/*--------------------------------------------------------------------------*/
/* Has a taylor file been written? */
int taylor_access()
{ 
   if (cur_vs_data!=NULL) return cur_vs_data->T_file_access;
   else return 0;
}

/*--------------------------------------------------------------------------*/
/* Close any open taylor file. */
void close_taylor()
{ fclose(cur_vs_data->temp2_file);
}

/****************************************************************************/
/*                                                           LOCAL ROUTINES */

/*--------------------------------------------------------------------------*/
/* T_Put_Block puts a block of tape to the disk.  I assume this 
   is called only during a successive forward pass, computation. */
static void T_put_block( int nitems )
{ int n;
  if (cur_vs_data->T_file_access == 0)
    cur_vs_data->temp2_file = fopen(cur_vs_data->vs_file_name,"w+b");
  if (cur_vs_data->T_write_cnt == 0)
    fseek(cur_vs_data->temp2_file,0,0); 
  cur_vs_data->T_file_access = 1;
  cur_vs_data->taylor_cnt = 0;
  if ((n=fwrite((char *)cur_vs_data->save_taylor,sizeof(revreal)*nitems,
                                         1,cur_vs_data->temp2_file)) != 1)
  { fprintf(DIAG_OUT,"ADOL-C error: fatal error-doing a write %d--- error %d\n",
                   n,errno);
    switch (errno) {
      case 28: /* ENOSPC */
	fprintf(DIAG_OUT,"No space left on device-contact sys. manager\n");
	break;
      case 27: /* EFBIG */
	fprintf(DIAG_OUT,"File to big-- Taylor-tape space exhausted.\n");
	break;
      default:
	fprintf(DIAG_OUT,"Unexpected error %d .\n",errno);
	break;
    }
    exit(-1);
  }
  cur_vs_data->T_write_cnt++;
}

/*--------------------------------------------------------------------------*/
/* Static function T_prev_block                                          
   called by taylor_back, taylor_back2, get_taylor, get_taylors    
   Gets the next (previous block) of size nitems  */
static int T_prev_block( int nitems )
{ int n;
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug: prev %d =nitems %d T_write_cnt \n", nitems, T_write_cnt);
#endif
  if (cur_vs_data->T_file_access)
  { if (cur_vs_data->T_write_cnt == 0)
      return 0;
    cur_vs_data->T_write_cnt--;
    fseek(cur_vs_data->temp2_file,cur_vs_data->T_buf_size*cur_vs_data->T_write_cnt*sizeof(revreal),0);
    n=fread((char *)cur_vs_data->save_taylor,sizeof(revreal),nitems,cur_vs_data->temp2_file);
    if (n != nitems)
    { fprintf(DIAG_OUT,"ADOL-C error: Read error on taylor file n= %d\n",n);
      return 0;
    }
    cur_vs_data->taylor_cnt = nitems;
    return 1;
  }
  return 0;
}


/****************************************************************************/
/*                                                            CONTROL STUFF */

/*--------------------------------------------------------------------------*/
/* Close the taylor file, reset data. */
void taylor_close( int buffer, int dep, int ind)
{ int n;
  if (buffer == -1)
    cur_vs_data->degsave = -1; /* enforces failure of reverse */
  cur_vs_data->numdep = dep;
  cur_vs_data->numind = ind;
  /* olvo 980708 changed to: ++.. */
  cur_vs_data->T_tail   = ++cur_vs_data->taylor_cnt;
  cur_vs_data->T_length = cur_vs_data->T_buf_size*cur_vs_data->T_write_cnt+cur_vs_data->taylor_cnt;
  if (cur_vs_data->T_write_cnt)
  { if (cur_vs_data->T_tail>0 ) 
      T_put_block(cur_vs_data->T_tail);
    free((char *)cur_vs_data->save_taylor);
    cur_vs_data->save_taylor = 0;
  }
  cur_vs_data->T_blocks = cur_vs_data->T_write_cnt;
  if ((cur_vs_data->T_blocks) && (cur_vs_data->T_length*sizeof(revreal) <= buffer))
  { cur_vs_data->save_taylor = (revreal *) malloc(cur_vs_data->T_length*sizeof(revreal));
    if (cur_vs_data->save_taylor == NULL)
    { fprintf(DIAG_OUT,"ADOL-C error: cannot allocate taylor buffer!\n");
      exit(-1);
    }
    fseek(cur_vs_data->temp2_file,0,0);
    n = fread((char *)cur_vs_data->save_taylor,sizeof(revreal),cur_vs_data->T_length,cur_vs_data->temp2_file);
    if ( n != cur_vs_data->T_length)
    { fprintf(DIAG_OUT,"ADOL-C error: read error in taylor_close n= %d\n",n);
      exit(-2);
    }
    cur_vs_data->T_tail = cur_vs_data->T_length;
    cur_vs_data->T_blocks = 0;
  }

#ifdef DEBUG
  if (cur_vs_data->T_blocks)
    fprintf(DIAG_OUT,"\n ADOL-C debug: taylor file of length %d bytes completed\n", cur_vs_data->T_length*sizeof(revreal));
  else
    fprintf(DIAG_OUT,"\n ADOL-C debug: taylor array of length %d bytes completed\n", cur_vs_data->T_length*sizeof(revreal));
#endif
}

/*--------------------------------------------------------------------------*/
/* Set up statics for writing taylor data */ 
void taylor_begin( short tag, int buffer, double** Tg, int degree )
{ 
  index = 0;
  if (num_vs_data != 0)
   {
     while((vs_data[index].tag != tag) && (index<num_vs_data))
       index++;
   }
  if (index >= num_vs_data)
    num_vs_data++;
  if (index >= TBUFNUM)
  { fprintf(DIAG_OUT,"ADOL-C error: to many taylor buffers!\n");
    fprintf(DIAG_OUT,"              Increase TBUFNUM\n");
    exit (-1);
  }
  cur_vs_data = &vs_data[index];
  cur_vs_data->tag = tag;
  if (cur_vs_data->save_taylor)
    free((char *)cur_vs_data->save_taylor);
  cur_vs_data->T_file_access = 0;
  /* 20040717 kowarz: added to fix problem "to many opened files" */
  if (cur_vs_data->temp2_file!=NULL) fclose(cur_vs_data->temp2_file);
  cur_vs_data->temp2_file = NULL;
  cur_vs_data->save_taylor = 0;
  cur_vs_data->T = Tg;
  get_fstr(cur_vs_data->vs_file_name,tag,FNAME3);
  cur_vs_data->T_buf_size = 1+buffer/sizeof(revreal);
  cur_vs_data->save_taylor = (revreal *)malloc(sizeof(revreal)*cur_vs_data->T_buf_size);
  if (cur_vs_data->save_taylor == NULL)
  { fprintf(DIAG_OUT,"ADOL-C error: cannot allocate taylor buffer!\n");
    exit (-1);
  }
  cur_vs_data->T_write_cnt = 0;
  cur_vs_data->T_length    = 0;
  cur_vs_data->taylor_cnt  = 0; 
  cur_vs_data->degsave     = degree;
}


/*--------------------------------------------------------------------------*/
void taylor_back2(int tag, revreal** Trg, int* dep, int* ind, int* degree)
{ 
  index = 0;
  if (num_vs_data != 0)
   {
     while((vs_data[index].tag != tag) && (index<num_vs_data))
       index++;
   }
  if (index >= num_vs_data)
  { fprintf(DIAG_OUT,"ADOL-C error: no taylor buffer for this tag \n");
    exit(-2);
  }
  cur_vs_data = &vs_data[index];
  *dep    = cur_vs_data->numdep;
  *ind    = cur_vs_data->numind;
  *degree = cur_vs_data->degsave;
  cur_vs_data->Tr = Trg;
  cur_vs_data->T_write_cnt = cur_vs_data->T_blocks;
  cur_vs_data->taylor_cnt  = cur_vs_data->T_tail;
  if (cur_vs_data->T_blocks == 0 && cur_vs_data->save_taylor == 0 )
  { fprintf(DIAG_OUT,"ADOL-C error: no temp file or array for reverse sweep \n");
    exit(-2);
  }
  if (cur_vs_data->T_blocks)
  { if (cur_vs_data->save_taylor)
      free((char*) cur_vs_data->save_taylor);
    cur_vs_data->save_taylor = (revreal*) malloc(cur_vs_data->T_buf_size*sizeof(revreal));
    if (cur_vs_data->save_taylor == NULL)
    { fprintf(DIAG_OUT,"ADOL-C error: cannot allocate taylor buffer!\n");
      exit (-1);
    } 
    if (T_prev_block(cur_vs_data->T_tail) == 0) 
      fprintf(DIAG_OUT,"ADOL-C error: problems in taylorback2 \n");
  }
}

/*--------------------------------------------------------------------------*/
void taylor_back( int tag, revreal* Trg, int* dep, int* ind, int* degree)
{ 
  index = 0;
  if (num_vs_data != 0)
   {
     while((vs_data[index].tag != tag) && (index<num_vs_data))
       index++;
   }
  if (index >= num_vs_data)
  { fprintf(DIAG_OUT,"ADOL-C error: no taylor buffer for this tag \n");
    exit(-2);
  }
  cur_vs_data = &vs_data[index];
  *dep    = cur_vs_data->numdep;
  *ind    = cur_vs_data->numind;
  *degree = cur_vs_data->degsave;
  cur_vs_data->Trs = Trg;
  cur_vs_data->T_write_cnt = cur_vs_data->T_blocks;
  cur_vs_data->taylor_cnt  = cur_vs_data->T_tail;
  if (cur_vs_data->T_blocks == 0 && cur_vs_data->save_taylor == 0 )
  { fprintf(DIAG_OUT,"ADOL-C error: no temp file or array for reverse sweep \n");
    exit(-2);
  }
  if (cur_vs_data->T_blocks)
  { if (cur_vs_data->save_taylor)
      free((char*) cur_vs_data->save_taylor);
    cur_vs_data->save_taylor = (revreal*) malloc(cur_vs_data->T_buf_size*sizeof(revreal));
    if (cur_vs_data->save_taylor == NULL)
    { fprintf(DIAG_OUT,"ADOL-C error: cannot allocate taylor buffer!\n");
      exit (-1);
    } 
    if (T_prev_block(cur_vs_data->T_tail) == 0)
      fprintf(DIAG_OUT,"ADOL-C error: problems in taylor_back \n");
  }
}


/****************************************************************************/
/*                                                                   WRITEs */

/*--------------------------------------------------------------------------*/
/* Write_taylor writes the block of size depth of taylor coefficients  
   from point loc to the taylor buffer.  If the buffer is filled, then 
   it is written to the taylor tape (T_put_block). */
void write_taylor( locint loc, int depth )
{ int i;
  double* Tloc = cur_vs_data->T[loc];
  for (i=0;i<depth;i++)
  { /* olvo 980708 changed order */
    if ((++cur_vs_data->taylor_cnt) == cur_vs_data->T_buf_size) 
      T_put_block(cur_vs_data->T_buf_size);
    cur_vs_data->save_taylor[cur_vs_data->taylor_cnt]=*Tloc++;          
                                 /* In this assignment the */
                                 /* precision will be sacrificed if the */
                                 /* type revreal is defined as float.   */
  }
}

/*--------------------------------------------------------------------------*/
/* Write_taylor writes the block of size depth of taylor coefficients  
   from point loc to the taylor buffer.  If the buffer is filled, then 
   it is written to the taylor tape (T_put_block). */
void write_taylors( locint loc, int depth, int k, int nrows )
{ int i,j;
  double* Tloc = cur_vs_data->T[loc];
  for (j=0;j<nrows;j++)
  {
   for (i=0;i<depth;i++)
   { /* olvo 980708 changed order */
     if ((++cur_vs_data->taylor_cnt) == cur_vs_data->T_buf_size) 
       T_put_block(cur_vs_data->T_buf_size);
     cur_vs_data->save_taylor[cur_vs_data->taylor_cnt]=*Tloc++;          
                                  /* In this assignment the */
                                  /* precision will be sacrificed if the */
                                  /* type revreal is defined as float.   */
   }
   for(i=depth;i<k;i++)
   {  
       *Tloc++;
   }
  }
}

/*--------------------------------------------------------------------------*/
/* Overwrite_scaylor overwrites the last (single) element (x) of the       
   taylor buffer.  New by olvo 980708;
   changed interface since 980921 to allow correction of
   old overwrite in store */
void overwrite_scaylor( revreal newVal, revreal* oldVal )
{ *oldVal = cur_vs_data->save_taylor[cur_vs_data->taylor_cnt];
  cur_vs_data->save_taylor[cur_vs_data->taylor_cnt] = newVal;
}

/*--------------------------------------------------------------------------*/
/* Delete_scaylor deletes the last (single) element (x) of the
   taylor buffer.  New by olvo 981122 */
void delete_scaylor( revreal* oldVal )
{ *oldVal = cur_vs_data->save_taylor[cur_vs_data->taylor_cnt--];
}

/*--------------------------------------------------------------------------*/
/* Write_scaylor writes a single element (x) to the taylor buffer.  If full
   the buffer is written out. */
void write_scaylor( revreal x )
{ /* olvo 980708 changed order */
  if ((++cur_vs_data->taylor_cnt) == cur_vs_data->T_buf_size) 
    T_put_block(cur_vs_data->T_buf_size);
  cur_vs_data->save_taylor[cur_vs_data->taylor_cnt]= x;
}


/*--------------------------------------------------------------------------*/
/* Write_scaylors writes # size elements from x to the taylor buffer.
   If full, the buffer is written out. */
void write_scaylors(double *x, int size)
{ int i;
  for(i=0; i<size; i++)
  { /* olvo 980708 changed order */
    if ((++cur_vs_data->taylor_cnt) == cur_vs_data->T_buf_size)
      T_put_block(cur_vs_data->T_buf_size);
    cur_vs_data->save_taylor[cur_vs_data->taylor_cnt]= x[i];
  }
}


/****************************************************************************/
/*                                                                     GETs */

/*--------------------------------------------------------------------------*/
void get_taylors( locint loc, int depth)
{ int i;
  revreal* Trloc = cur_vs_data->Tr[loc];
  for (i=depth-1;i >= 0;i--)
  { if (cur_vs_data->taylor_cnt == 0)
      if (!T_prev_block(cur_vs_data->T_buf_size))
      { fprintf(stderr,"ADOL-C error: Fatal Error in get_taylors ");
	exit(-1);
      }
    Trloc[i] = cur_vs_data->save_taylor[--cur_vs_data->taylor_cnt];
  }
}

/*--------------------------------------------------------------------------*/
void get_taylors_p( locint loc, int depth, int p)
{ int i,j,cnt,base;
  revreal* Trloc = cur_vs_data->Tr[loc];
  
  cnt=p*depth-1;
  base = cur_vs_data->taylor_cnt-cnt-1;
  for (j=p;j > 0;j--)
  {
   for (i=depth-1;i > 0;i--)
   { if (cur_vs_data->taylor_cnt == 0)
        if (!T_prev_block(cur_vs_data->T_buf_size))
        { fprintf(stderr,"ADOL-C error: Fatal Error in get_taylors ");
          exit(-1);
        }
    Trloc[cnt--] = cur_vs_data->save_taylor[--cur_vs_data->taylor_cnt];  /* directions */
   }
   cnt--;
  }
  if (cur_vs_data->taylor_cnt == 0)  
    if (!T_prev_block(cur_vs_data->T_buf_size))
     { fprintf(stderr,"ADOL-C error: Fatal Error in get_taylors ");
       exit(-1);
     }
   Trloc[0] = cur_vs_data->save_taylor[--cur_vs_data->taylor_cnt]; /* base point */
   for(j=1;j<p ;j++)
    Trloc[j*depth] = cur_vs_data->save_taylor[cur_vs_data->taylor_cnt];
   
}

/*--------------------------------------------------------------------------*/
void get_taylor( locint loc )
{ if (cur_vs_data->taylor_cnt == 0)
  { if (!T_prev_block(cur_vs_data->T_buf_size))
    { fprintf(stderr,"ADOL-C error: Fatal Error in get_taylor ");
      exit(-1);
    }
  }
  cur_vs_data->Trs[loc] = cur_vs_data->save_taylor[--cur_vs_data->taylor_cnt];
}

/****************************************************************************/
/*                                                       DEALLOCATE VS_DATA */

void clean_vs_data( int tag )
{
  int i;
  index = 0;
  if (num_vs_data != 0)
   {
     while((vs_data[index].tag != tag) && (index<num_vs_data))
       index++;
   }
  if (index >= num_vs_data)
  { fprintf(DIAG_OUT,"ADOL-C error: no taylor buffer for this tag \n");
    exit(-2);
  }
  printf(" tag = %d index = %d \n",tag,index);
 free((char*) vs_data[index].save_taylor);
 for(i=index;i<num_vs_data-1;i++)
   vs_data[index] = vs_data[index+1];
 free((char*) vs_data[num_vs_data-1].save_taylor);
 num_vs_data--;
}

/****************************************************************************/
/*                                                               THAT'S ALL */

END_C_DECLS
