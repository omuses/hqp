/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     avector.cpp
 Revision: $Id: avector.cpp,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: Avector.C contains the necessary routines for vector operations       
           that are defined in avector.h.  Note: avector.h is included 
           automatically by adouble.h, and hence does not need to be 
           included here again.
           
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
          19981130 olvo:   last check (includes ...)
	                        NOTICE: I think everything concerning vectors 
                                   has to be checked again in detail!
          19980930 olvo:   allow overwrites in av*=a, av*a, a*av
          19980924 olvo:   changed all int_* opcodes
          19980721 olvo:   write of taylors in m_subscript
          19980714 olvo:   debugging vector - matrix stuff
          19980713 olvo:   elimination of "writes" from taputil1.c completed
          19980707 olvo:   (1) used void write_dot_av_av(..)
                           (2) taping with keep

----------------------------------------------------------------------------*/

#include "adouble.h"
#include "oplate.h"
#include "taputil.h"
#include "taputil_p.h"
#include "tayutil.h"
#include "tayutil_p.h"

#include <math.h>

/****************************************************************************/
/*                                                   GLOBAL VARS & ROUTINES */

/*--------------------------------------------------------------------------*/
extern double* store;
extern int trace_flag;

/*--------------------------------------------------------------------------*/
extern locint next_loc(int size);
extern locint next_loc();
extern locint free_loc(locint,int);


/****************************************************************************/
/*                                                      VECTOR CONSTRUCTORS */
  
/*--------------------------------------------------------------------------*/
adoublev::adoublev( int n )
{
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug:Declaring active vector\n");
#endif
  size = n;
  start_loc = next_loc(size); 
}

/*--------------------------------------------------------------------------*/
adoublev::adoublev( const adoublev &arg )
{
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug:Declaring active vector and"
                   " initializing from adoublev\n");
#endif
  size = arg.size;
  start_loc = next_loc(size);

  if (trace_flag) 
  { // old: write_intvec_assign_av(size,start_loc,arg.start_loc);
    put_op(assign_av);
    put_locint(arg.start_loc); // = arg
    put_locint(size);
    put_locint(start_loc);     // = res

    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }

  for (int i=0; i<size; i++)
    store[start_loc+i] = store[arg.start_loc+i];    
}

/*--------------------------------------------------------------------------*/
adoublev::adoublev( const adubv& a)
{ /* olvo 980713 what about size? */
  size = a.sz();
  start_loc = next_loc(size);

  if (trace_flag)
  { // old: write_intvec_assign_av(size,start_loc,a.loc());
    put_op(assign_av);
    put_locint(a.loc());   // = arg
    put_locint(size);
    put_locint(start_loc); // = res

    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }

  for (int i=0; i<size; i++)
    store[start_loc+i] = store[a.loc()+i];
}


/****************************************************************************/
/*                                                       VECTOR ASSIGNMENTS */

/*--------------------------------------------------------------------------*/
adoublev& adoublev::operator = ( double* coval ) 
{ if (trace_flag) 
  { // old: write_assign_vec_dv(size,start_loc,coval);
    locint space_left = get_val_space(), 
           vals_left  = size, 
           loc        = start_loc;
    double *d = coval;

    while (space_left < vals_left)
    { put_op(assign_dv);
      put_locint(space_left);
      put_locint(loc);
      put_vals_p(d,space_left);
      d         += space_left;
      vals_left -= space_left;
      loc       += space_left;
      space_left = get_val_space();
    } 

    if (vals_left > 0)
    { put_op(assign_dv);
      put_locint(vals_left);
      put_locint(loc);
      put_vals_r(d,vals_left);
    }

    vs_ptr += size;
    if (revalso)
      write_scaylors((store+start_loc),size);
  }

  for (int i=0; i<size; i++)
    store[start_loc+i] = coval[i];
  return *this;
}

/*--------------------------------------------------------------------------*/
adoublev& adoublev::operator = ( double coval )
{
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug:In adoublev=double\n");
#endif
  /* olvo 980713 very tricky */
  if (trace_flag)
  { vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }

  for (int i=0; i<size; i++)
     store[start_loc+i] = coval;

  if (trace_flag)
  { // old: write_assign_vec_dv(size,start_loc,store+start_loc);
    locint space_left = get_val_space(), 
           vals_left  = size,
           loc        = start_loc;
    double *d         = store + start_loc;

    while (space_left < vals_left)
    { put_op(assign_dv);
      put_locint(space_left);
      put_locint(loc);
      put_vals_p(d,space_left);
      d         += space_left;
      vals_left -= space_left;
      loc       += space_left;
      space_left = get_val_space();
    } 

    if (vals_left > 0)
    { put_op(assign_dv);
      put_locint(vals_left);
      put_locint(loc);
      put_vals_r(d,vals_left);
    } 
  }

  return *this;
}

/*--------------------------------------------------------------------------*/
adoublev& adoublev::operator = ( const badoublev& x ) 
{ if(start_loc != x.loc())
  /* test this to avoid  adjoint(x)=0 for x=x in reverse */
  { if (trace_flag) 
    { // old: write_assign_av(size,start_loc,x.loc());
      put_op(assign_av);
      put_locint(x.loc());   // = arg
      put_locint(size);
      put_locint(start_loc); // = res

      vs_ptr += size;
      if (revalso) 
        write_scaylors((store+start_loc),size);
    }
    
    for (int i=0; i<size; i++)
      store[start_loc+i] = store[x.loc()+i];
  } 
  return *this;
}

/*--------------------------------------------------------------------------*/
adoublev& adoublev::operator = ( const adoublev& x ) 
{ if(start_loc != x.start_loc)
  /* test this to avoid  adjoint(x)=0 for x=x in reverse */
  { if (trace_flag) 
    { // old: write_assign_av(size,start_loc,x.start_loc);
      put_op(assign_av);
      put_locint(x.start_loc); // = arg
      put_locint(size);
      put_locint(start_loc);   // = res

      vs_ptr += size;
      if (revalso) 
        write_scaylors((store+start_loc),size);
    }

    for (int i=0; i<size; i++)
      store[start_loc+i] = store[x.start_loc+i];
  } 
  return *this;
}

/*--------------------------------------------------------------------------*/
adoublev& adoublev::operator = ( const adubv& a )
{ if (start_loc != a.loc())
  /* test this to avoid  adjoint(x)=0 for x=x in reverse */
  { if (trace_flag)
    { // old: write_assign_av(size,start_loc,a.loc());
      put_op(assign_av);
      put_locint(a.loc());   // = arg
      put_locint(size);
      put_locint(start_loc); // = res

      vs_ptr += size;
      if (revalso) 
        write_scaylors((store+start_loc),size);
    }

    for (int i=0; i<size; i++)
      store[start_loc+i] = store[a.loc()+i];
  } 
  return *this;
}


/****************************************************************************/
/*                                                              DESTRUCTORS */
#ifdef overwrite

/*--------------------------------------------------------------------------*/
adoublev::~adoublev()
{
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug:Destructing active vector\n");
#endif
  free_loc(start_loc,size);
}

/*--------------------------------------------------------------------------*/
adubv::~adubv()
{ free_loc(start_loc,size);
}

/*--------------------------------------------------------------------------*/
adoublem::~adoublem()
{
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug:Destructing active matrix\n");
#endif
  delete[] index;
}

#endif


/****************************************************************************/
/*                                                           INPUT / OUTPUT */

/*--------------------------------------------------------------------------*/
std::ostream& operator << ( std::ostream& out, const badoublev &arg )
{ out << "(";
  for (int i=0; i<arg.size-1; i++)
    out << store[arg.start_loc+i] << ", ";
  out << store[(arg.start_loc+arg.size)-1] << ")(a)";
  return out;
}  


/****************************************************************************/
/*                                                                    INDEX */

/*--------------------------------------------------------------------------*/
badouble badoublev::operator[](int i) const  
{ /* Used so can access the vector like an array with the [] */
  /* Check if out of range */
  if (i<0 || i>=size)
  { fprintf (DIAG_OUT,"ADOL-C error: adoublev index out of range.\n");
    exit(-3);
  }

  return start_loc+i;
}


/****************************************************************************/
/*                                                              ASSIGNMENTS */

/*--------------------------------------------------------------------------*/
badoublev& badoublev::operator = ( const badoublev &arg )
{
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug:In badoublev = badoublev\n");
#endif
  if (start_loc != arg.start_loc)
  /* test this to avoid  adjoint(x)=0 for x=x in reverse */
  { if (trace_flag)
    { // old: write_assign_av(size,start_loc,arg.start_loc);
      put_op(assign_av);
      put_locint(arg.start_loc); // = arg
      put_locint(size);
      put_locint(start_loc);     // = res

      vs_ptr += size;
      if (revalso) 
        write_scaylors((store+start_loc),size);
    }
  
    for (int i=0; i<arg.size; i++)
      store[start_loc+i] = store[arg.start_loc+i];
  }
  return *this;
}

/*--------------------------------------------------------------------------*/
badoublev& badoublev::operator = ( const adubv &arg )
{ locint arg_start_loc = arg.loc();
  locint arg_size      = arg.sz();

  if (start_loc != arg_start_loc)
  /* test this to avoid  adjoint(x)=0 for x=x in reverse */
  { if (trace_flag)
    { // old: write_assign_av(size,start_loc,arg_start_loc);
      put_op(assign_av);
      put_locint(arg_start_loc); // = arg
      put_locint(size);
      put_locint(start_loc);     // = res

      vs_ptr += size;
      if (revalso) 
        write_scaylors((store+start_loc),size);
    }
  
    for (locint i=0; i<arg_size; i++)
      store[start_loc+i] = store[arg_start_loc+i];
  }
  return *this;
}

/*--------------------------------------------------------------------------*/
badoublev& badoublev::operator = ( const adoublev& x ) 
{ locint x_start_loc=x.loc();

  if (start_loc != x_start_loc)
  /* test this to avoid  adjoint(x)=0 for x=x in reverse */
  { if (trace_flag) 
    { // old: write_assign_av(size,start_loc,x_start_loc);
      put_op(assign_av);
      put_locint(x_start_loc); // = arg
      put_locint(size);
      put_locint(start_loc);   // = res

      vs_ptr += size;
      if (revalso) 
        write_scaylors((store+start_loc),size);
    }

    for (int i=0; i<size; i++)
      store[start_loc+i] = store[x_start_loc+i];
  }
  return *this;
}   

/*--------------------------------------------------------------------------*/
/* Assign an adouble vector an independent float vector */
adoublev& adoublev::operator <<= ( double* coval ) 
{
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug:IND EQ double*\n");
#endif

  if (trace_flag)
  { // old: write_assign_indvec(size,start_loc,coval);
    ind_ptr += size;
    put_op(assign_indvec);
    put_locint(size);
    put_locint(start_loc); // = res

    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }

  for (int i=0; i<size; i++)
    store[(start_loc)+i] = coval[i];
  return *this;
}   

/*--------------------------------------------------------------------------*/
/* Assign a float vector a dependent adouble vector */
adoublev& adoublev::operator >>= ( double* coval ) 
{  
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug:DEP EQ double* operator\n");
#endif
  if (trace_flag) 
  { // old: write_assign_depvec(size,start_loc);
    dep_ptr += size;
    put_op(assign_depvec);
    put_locint(size);
    put_locint(start_loc); // = res
  }

  for (int i=0; i<size; i++)
    coval[i] = double (store[(start_loc)+i]);
  return *this;
}   


/****************************************************************************/
/*                                            VECTOR OPERATION + ASSIGNMENT */

/*--------------------------------------------------------------------------*/
badoublev& badoublev::operator -= ( const badoublev& y ) 
{ if (trace_flag)
  { // old: write_av_same_arg(eq_min_av,size,start_loc,y.start_loc);
    put_op(eq_min_av);
    put_locint(y.start_loc); // = arg
    put_locint(size);
    put_locint(start_loc);   // = res

    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }

  for (int i=0; i<size; i++)
    store[start_loc+i] -= store[y.start_loc+i];
  return *this;
}

/*--------------------------------------------------------------------------*/
badoublev& badoublev::operator += ( const badoublev& y ) 
{ if (trace_flag)
  { // old: write_av_same_arg(eq_plus_av,size,start_loc,y.start_loc);
    put_op(eq_plus_av);
    put_locint(y.start_loc); // = arg
    put_locint(size);
    put_locint(start_loc);   // = res

    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }
  
  for (int i=0; i<size; i++)
    store[start_loc+i] += store[y.start_loc+i];
  return *this;
}

/*--------------------------------------------------------------------------*/
badoublev& badoublev::operator *= ( double coval ) 
{ if (trace_flag)
  { // old: write_samearg_av_d(eq_mult_av_d,size,start_loc,coval);
    put_op(eq_mult_av_d);
    put_locint(size);
    put_locint(start_loc); // = res
    put_val(coval);

    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }

  for (int i=0; i<size; i++)
    store[start_loc+i] *= coval;
  return *this;
}

/*--------------------------------------------------------------------------*/
badoublev& badoublev::operator *= ( const badouble& y ) 
{ int loc = y.loc();

  if (trace_flag)
  { // old: write_av_same_arg(eq_mult_av_a,size,start_loc,loc);
    put_op(eq_mult_av_a);
    put_locint(loc);       // = arg
    put_locint(size);
    put_locint(start_loc); // = res

    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }

  /* olvo 980930 use tempory to allow overwrites */
  double tmpVal = store[loc];
  for (int i=0; i<size; i++)
    store[start_loc+i] *= tmpVal;
  return *this;
}

/*--------------------------------------------------------------------------*/
badoublev& badoublev::operator /= ( double coval ) 
{ *this = *this / coval;
  return *this;
}

/*--------------------------------------------------------------------------*/
badoublev& badoublev::operator /= ( const badouble& y ) 
{*this = *this * (1.0/y);
  return *this;
}


/****************************************************************************/
/*                                                 BINARY VECTOR OPERATIONS */

/*--------------------------------------------------------------------------*/
adubv operator + ( const badoublev &arg1, const badoublev &arg2 )
{ locint size      = arg1.size;
  locint start_loc = next_loc(size);

#ifdef DEBUG
  if (arg1.size != arg2.size)
  { fprintf(DIAG_OUT,"ADOL-C error: Can not add vectors as not same size\n");
    exit(-3);
  }
#endif

  if (trace_flag) 
  { // old: write_two_av_rec(plus_av_av,size,start_loc,
    //			     arg1.start_loc,arg2.start_loc);
    put_op(plus_av_av);
    put_locint(arg1.start_loc); // = arg1
    put_locint(arg2.start_loc); // = arg2
    put_locint(size);
    put_locint(start_loc);      // = res

    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }

  for (locint i=0; i<size; i++)
    store[start_loc+i] =   store[arg1.start_loc+i]
                         + store[arg2.start_loc+i];
  return adubv(start_loc,size);
}

/*--------------------------------------------------------------------------*/
adubv operator * ( const badoublev &arg, double coval )
{ locint size = arg.size;
  locint start_loc = next_loc(size);

  if (trace_flag)
  { // old:  write_args_d_av(mult_d_av,size,start_loc,coval,arg.start_loc);
    put_op(mult_d_av);
    put_locint(arg.start_loc); // = arg
    put_locint(size);
    put_locint(start_loc);     // = res
    put_val(coval);            // = coval

    vs_ptr += size;
    if (revalso)
       write_scaylors((store+start_loc),size);
  }

  for (locint i=0; i<size; i++)
    store[start_loc+i] = store[arg.start_loc+i]*coval;
  return adubv(start_loc,size);
}

/*--------------------------------------------------------------------------*/
adubv operator * ( double coval, const badoublev &arg )
{ locint size = arg.size;
  locint start_loc = next_loc(size);

  if (trace_flag)
  { // old: write_args_d_av(mult_d_av,size,start_loc,coval,arg.start_loc);
    put_op(mult_d_av);
    put_locint(arg.start_loc); // = arg
    put_locint(size);
    put_locint(start_loc);     // = res
    put_val(coval);            // = coval

    vs_ptr += size;
    if (revalso)
       write_scaylors((store+start_loc),size);
  }    

  for (locint i=0; i<size; i++)
    store[start_loc+i] = store[arg.start_loc+i]*coval;  
  return adubv(start_loc,size);
}

/*--------------------------------------------------------------------------*/
adub operator* ( const badoublev &arg1, const badoublev &arg2 ) 
{ double x = 0;
  locint locat = next_loc();
  
#ifdef DEBUG
  if (arg1.size!=arg2.size)
  { fprintf(DIAG_OUT,"ADOL-C error: Can not take dot product,"
                     " vectors are not same size\n");
    exit(-3);
  }
#endif

  if (trace_flag)
  { // old: write_dot_av_av(arg1.size,locat,arg1.start_loc,arg2.start_loc);
    put_op(dot_av_av);
    put_locint(arg1.start_loc); // = arg1
    put_locint(arg2.start_loc); // = arg2
    put_locint(arg1.size);
    put_locint(locat);      // = res

    vs_ptr++;
    if (revalso)
      write_scaylor(store[locat]);
  }
      
  for (int i=0; i<arg1.size; i++)
     x += store[arg1.start_loc+i] * store[arg2.start_loc+i];
  store[locat] = x;
  return locat;
}

/*--------------------------------------------------------------------------*/
adubv operator / ( const badoublev &x, const badouble &y )
{ int loc  = y.loc();
  int size = x.size;
  locint start_loc = next_loc(size);

  if (trace_flag)
  { // old: write_av_a_rec(div_av_a,size,start_loc,x.start_loc,loc);
    put_op(div_av_a);
    put_locint(x.start_loc); // = arg1
    put_locint(loc);         // = arg2
    put_locint(size);
    put_locint(start_loc);   // = res

    vs_ptr += size;
    if (revalso)
      write_scaylors((store+start_loc),size);
  }
  
  for (int i=0; i<size; i++)
    store[start_loc+i] = store[x.start_loc+i]*(1.0/store[loc]);
  return adubv(start_loc,size);
}

/*--------------------------------------------------------------------------*/
adubv operator - ( const badoublev &arg1, const badoublev &arg2 )
{ locint size = arg1.size;
  locint start_loc = next_loc(size);
#ifdef DEBUG
  if (arg1.size != arg2.size)
  { fprintf(DIAG_OUT,"ADOL-C error: Can not add vectors as not same size\n");
    exit(-3);
  }
#endif
      
  if (trace_flag) 
  { // old: write_two_av_rec(sub_av_av,size,start_loc,
    //			     arg1.start_loc,arg2.start_loc);
    put_op(sub_av_av);
    put_locint(arg1.start_loc); // = arg1
    put_locint(arg2.start_loc); // = arg2
    put_locint(size);
    put_locint(start_loc);      // = res

    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }
      
  for (locint i=0; i<size; i++)
    store[start_loc+i] =  store[arg1.start_loc+i]
                        - store[arg2.start_loc+i];
  return adubv(start_loc,size);
}

/*--------------------------------------------------------------------------*/
adubv operator * ( const badoublev &arg, const badouble &n )
{ int loc = n.loc();
  int size = arg.size;
  locint start_loc = next_loc(size);

  if (trace_flag)
  { // old: write_av_a_rec(mult_av_a,size,start_loc,arg.start_loc,loc);
    put_op(mult_a_av);
    put_locint(arg.start_loc); // = arg1
    put_locint(loc);           // = arg2
    put_locint(size);
    put_locint(start_loc);     // = res

    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }

  /* olvo 980930 use tempory to allow overwrites */
  double tmpVal = store[loc];
  for (int i=0; i<size; i++)
    store[start_loc+i] = store[arg.start_loc+i]*tmpVal;    
  return adubv(start_loc,size);
}

/*--------------------------------------------------------------------------*/
adubv operator * ( const badouble &n, const badoublev &arg )
{ int loc = n.loc();
  int size = arg.size;
  locint start_loc = next_loc(size);

  if (trace_flag)
  { // old: write_av_a_rec(mult_a_av,size,start_loc,arg.start_loc,loc);
    put_op(mult_a_av);
    put_locint(arg.start_loc); // = arg1
    put_locint(loc);           // = arg2
    put_locint(size);
    put_locint(start_loc);     // = res

    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);
  }

  /* olvo 980930 use tempory to allow overwrites */
  double tmpVal = store[loc];
  for (int i=0; i<size; i++)
    store[start_loc+i] = store[arg.start_loc+i]*tmpVal;    
  return adubv(start_loc,size);
}

/****************************************************************************/
/*                                                             MATRIX STUFF */

/*--------------------------------------------------------------------------*/
adoublem::adoublem(int row, int col)
{ m = row;
  n = col;
  index = new adoublev[m];
  for (int i=0; i<m; i++)
  { index[i].size = n;
    index[i].start_loc = next_loc(n);
  }
}

/*--------------------------------------------------------------------------*/
adoublem::adoublem(const adoublem &arg)
{
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug:Declaring active matrix and initializing"
                   " from adoublem\n");
#endif
  m = arg.m;
  n = arg.n;
  index = new adoublev[m];
  for (int i=0; i < m; i++)
  { index[i].size = n;
    index[i].start_loc = next_loc(n);

    if (trace_flag)
    { /* old: write_intvec_assign_av(n, index[i].start_loc,
         arg.index[i].start_loc); */
      put_op(assign_av);
      put_locint(arg.index[i].start_loc); // = arg
      put_locint(n);                      // = size
      put_locint(index[i].start_loc);     // = res

      vs_ptr += n;
      if (revalso) 
        write_scaylors((store+(index[i].start_loc)),n);
    }

    for (int j=0; j < n; j++)
      store[index[i].start_loc+j] = store[arg.index[i].start_loc+j];
  } 
}

/*--------------------------------------------------------------------------*/
adoublev& adoublem::operator[]( int i ) 
{ if (i<0 || i>=m)
  { fprintf (DIAG_OUT,"ADOL-C error: adoublem index out of range.\n");
    exit(-3);
  }
  return index[i];
}

/*--------------------------------------------------------------------------*/
asub badoublev::operator[]( const along &i )  const
{ int j=(int)(store[i.loc()]);
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug:In along overloaded []\n");
#endif
  /* Used so can access the vector like an array with the [] */
  /* Check if out of range */
  if ((j<0) || (j>=size))
    fprintf (DIAG_OUT,"ADOL-C warning:: adoublev index out of range.\n");
  return asub(start_loc,i.loc());
}

/****************************************************************************/
/*                                                         ASUBV OPERATIONS */

#ifdef overwrite
/*--------------------------------------------------------------------------*/
asubv::~asubv()
{
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug:Destructing active subscript vector\n");
#endif
  free_loc(start_loc,size);
}

#endif

/*--------------------------------------------------------------------------*/
asubv::asubv( adoublev* start, locint index )
{
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug: Constructing an asubv with 3 arguments\n");
#endif
  begin  = (start[0]).loc(); /* start of matrix */
  base   = (start[(int)store[index]]).loc(); /* start of the i-th row */
  offset = index;
  size   = (start[(int)store[index]]).sz(); /* size of the row-vector */
  start_loc = next_loc(size);

  if (trace_flag)
  { // old: write_associating_vector(m_subscript,start_loc,begin,offset,size);
    put_op(m_subscript);
    put_locint(begin);
    put_locint(offset);
    put_locint(size);
    put_locint(start_loc);
    put_val(store[offset]);

    /* olvo 980721 new n3l */
    vs_ptr += size;
    if (revalso)
      write_scaylors(store+start_loc,size); 
  }

  for(int i=0;i<size;i++)
    store[start_loc+i] = store[base+i];
}


/*--------------------------------------------------------------------------*/
asubv adoublem::operator[]( const along &i )
{ int j = (int)(store[i.loc()]);
#ifdef DEBUG
  fprintf(DIAG_OUT,"ADOL-C debug: In along overloaded []\n");
#endif
  /* Used so can access the vector like an array with the [] */
  /* Check if out of range */
  if (j<0 || j>=n)
    fprintf (DIAG_OUT,"ADOL-C warning:: adoublem index out of range.\n");
  return asubv(index,i.loc());
}

/*--------------------------------------------------------------------------*/
asubv& asubv::operator = ( const adubv& a )
{ if (trace_flag)
  { // old: write_associating_vector(m_subscript_l,a.loc(),begin,offset,size);
    put_op(m_subscript_l);
    put_locint(begin);
    put_locint(offset);
    put_locint(size);
    put_locint(a.loc());
    put_val(store[offset]);

    vs_ptr+=size;
    if (revalso)
      write_scaylors((store+(begin+(int)store[offset])),size); 
  }

  for(int i=0;i<size;i++)
    store[base+i] = store[a.loc()+i] ;
  return *this;
}

/*--------------------------------------------------------------------------*/
asubv& asubv::operator = ( const badoublev& x )
{ if (trace_flag)
  { //old: write_associating_vector(m_subscript_l,x.loc(),begin,offset,size);
    put_op(m_subscript_l);
    put_locint(begin);
    put_locint(offset);
    put_locint(size);
    put_locint(x.loc());
    put_val(store[offset]);

    vs_ptr+=size;
    if (revalso)
      write_scaylors((store+(begin+(int)store[offset])),size); 
  }

  for(int i=0;i<size;i++)
    store[base+i] = store[x.loc()+i];
  return *this;
}  

/*--------------------------------------------------------------------------*/
asubv& asubv::operator <<= (double* y)
{ if (trace_flag)
  { // old: write_assign_indvec(size,start_loc,y);
    ind_ptr += size;
    put_op(assign_indvec);
    put_locint(size);
    put_locint(start_loc);
  
    vs_ptr += size;
    if (revalso) 
      write_scaylors((store+start_loc),size);

    /* old: write_associating_vector(m_subscript_l,start_loc,
       begin,offset,size); */
    put_op(m_subscript_l);
    put_locint(begin);
    put_locint(offset);
    put_locint(size);
    put_locint(start_loc);
    put_val(store[offset]);

    vs_ptr+=size;
    if (revalso)
      write_scaylors((store+(begin+(int)store[offset])),size); 
  }
  for(int i=0;i<size;i++)
    store[base+i] = y[i];
  return *this;
}

/*--------------------------------------------------------------------------*/
asubv& asubv::operator = (double* x)
{ if (trace_flag)
  { // old: write_associating_vector_ld(x,begin,offset,size);
    locint space_left = get_val_space(), 
           vals_left  = size, 
           loc        = 0;
    double *d = x;    
  
    while (space_left < vals_left)
    { put_op(m_subscript_ld);
      put_locint(begin);
      put_locint(offset);
      put_val(store[offset]);
      put_locint(loc);
      put_locint(space_left);
      put_vals_p(d,space_left);
      d         += space_left;
      vals_left -= space_left;
      loc       += space_left;
      space_left=get_val_space();
    }

    if (vals_left > 0)
    { put_op(m_subscript_ld);
      put_locint(begin);
      put_locint(offset);
      put_val(store[offset]);
      put_locint(loc);
      put_locint(vals_left);
      put_vals_r(d,vals_left);
    }

    vs_ptr += size;
    if (revalso)
      write_scaylors((store+(begin+(int)store[offset])),size);
  }

  for(int i=0;i<size;i++)
    store[base+i] = x[i];
  return *this;
}

/*--------------------------------------------------------------------------*/
asubv& asubv::operator = ( const asubv& x )
{ if (trace_flag)
  { // old: write_associating_vector(m_subscript_l,x.loc(),begin,offset,size);
    put_op(m_subscript_l);
    put_locint(begin);
    put_locint(offset);
    put_locint(size);
    put_locint(x.loc());
    put_val(store[offset]);

    vs_ptr+=size;
    if (revalso)
      write_scaylors((store+(begin+(int)store[offset])),size); 
  }

  for(int i=0;i<size;i++)
    store[base+i] = store[x.loc()+i];
  return *this;
}

/*--------------------------------------------------------------------------*/
asubv& asubv::operator += ( const badoublev& x )
{ *this = *this + x;
  return *this;
}

/*--------------------------------------------------------------------------*/
asubv& asubv::operator -= ( const badoublev& x )
{ *this = *this - x;
  return *this;
}

/*--------------------------------------------------------------------------*/
asubv& asubv::operator *= ( const badouble& x )
{ *this = *this * x;
  return *this;
}

/*--------------------------------------------------------------------------*/
asubv& asubv::operator *= ( double coval )
{ *this = *this * coval;
  return *this;
}

/*--------------------------------------------------------------------------*/
asubv& asubv::operator /= ( const badouble& x )
{ *this = *this * (1.0/x);
  return *this;
}

/*--------------------------------------------------------------------------*/
asubv& asubv::operator /= ( double coval )
{ *this = *this / coval;
  return *this;
}

/*--------------------------------------------------------------------------*/
adubv operator+ ( const badoublev& x )
{ return x * (1.0);
}
  
/*--------------------------------------------------------------------------*/
adubv operator- ( const badoublev& x )
{ return x * (-1.0);
}

/****************************************************************************/
/*                                                                THAT'S ALL*/
