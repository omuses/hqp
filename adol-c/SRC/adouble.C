#define _ADOUBLE_CPP_
#define _ADOLC_SRC_
/*
   --------------------------------------------------------------
   File adouble.C of ADOL-C version 1.8.5         as of Dec/10/99
   --------------------------------------------------------------
   adouble.C contains that definitions of procedures used to 
   define various badouble, adub, asub and adouble operations. 
   These operations actually have two purposes.
   The first purpose is to actual compute the function, just as 
   the same code written for double precision (single precision -
   complex - interval) arithmetic would.  The second purpose is 
   to write a transcript of the computation for the reverse pass 
   of automatic differentiation.

   Last changes:
      991210 olvo  checking the changes
      991122 olvo  new op_codes eq_plus_prod eq_min_prod
                   for  y += x1 * x2
                   and  y -= x1 * x2  
      981130 olvo  last check (includes ...)
      981119 olvo  changed tanh as J.M. Aparicio suggested
      981020 olvo  skip upd_resloc(..) if no tracing performed
      980924 olvo  changed all int_* opcodes
      980921 olvo  (1) changed save-order in sin/cos
                   (2) new interface in call to
                       void overwrite_scaylor(..) which
                       allows correction of old overwrite in store
      980820 olvo  new comparison strategy
      980721 olvo: write of taylors in subscript
      980713 olvo: elimination of "writes" from taputil1.c completed
      980710 olvo  sin/cos writes 2 taylors
      980709 olvo: elimination of "writes" from taputil1.c
      980708 olvo: new: write_upd(..)
      980707 olvo: taping with keep
      980706 olvo: new operation code: incr_a
                                       decr_a
      980623 olvo: griewank's idea -- stock manipulation
      980518 olvo: griewank's idea -- write_death(..) killed
      980517 olvo: griewank's idea -- operator:
                   adouble& adouble::operator = (const adub& a)

   --------------------------------------------------------------
*/


/****************************************************************************/
/*                                                                 INCLUDES */

#include "dvlparms.h"
#include "usrparms.h"
#include "adouble.h"
#include "oplate.h"
#include "taputil.h"
#include "tayutil.h"

#include <stdio.h>
#include <malloc.h>
#include <math.h>


/****************************************************************************/
/*                                                        HELPFUL FUNCTIONS */

/*--------------------------------------------------------------------------*/
void condassign( double &res, const double &cond, 
                 const double &arg1, const double &arg2 )
{ res = cond ? arg1 : arg2;
}

/*--------------------------------------------------------------------------*/
void condassign( double &res, const double &cond, 
                 const double &arg)
{ res = cond ? arg : res;
}

/*--------------------------------------------------------------------------*/
double fmax( const double &x, const double &y )
{ if (y > x) 
    return y;  
  else 
    return x; 
}

/*--------------------------------------------------------------------------*/
double fmin( const double &x, const double &y )
{ if (y < x) 
    return y;
  else 
    return x;
}

/****************************************************************************/
/*                                                              GLOBAL VARS */

/*--------------------------------------------------------------------------*/
double* store;                      // = double stack
int trace_flag             = 0;

/*--------------------------------------------------------------------------*/
static locint maxloc       = sizeof(locint) ==2 ? 65535 : 2147483647;
static locint current_top  = 0;     // = largest live location + 1
static locint location_cnt = 0;     // = maximal # of lives so far
static locint maxtop       = 0;     // = current size of store
static locint maxtop2;
static locint dealloc      = 0;     // = # of locations to be freed
static locint deminloc     = 0;     // = lowest loc to be freed 



/****************************************************************************/
/*                                                         MEMORY MANAGMENT */

/*--------------------------------------------------------------------------*/
/* Return the next free location in "adouble" memory */
locint next_loc()
{ /* First deallocate dead adoubles if they form a contiguous tail: */
#ifdef overwrite
  if (dealloc && dealloc+deminloc == current_top)    
  { /* olvo 980518 deleted write_death (see version griewank) */
    // if (trace_flag) 
    //   write_death(deminloc, current_top - 1);
    current_top = deminloc ;
    dealloc     = 0; 
    deminloc    = maxloc;
  }
#endif
  if (current_top == location_cnt)
    ++location_cnt;
  if (location_cnt > maxtop) 
  { maxtop2 = ++maxtop*2 > maxloc ? maxloc : 2*maxtop;
    if (maxtop2 == maxloc)
    { fprintf(DIAG_OUT,"\nADOL-C error:\n");
      fprintf(DIAG_OUT,"maximal number (%d) of live active variables exceeded\n\n", maxloc);
      fprintf(DIAG_OUT,"Possible remedies :\n\n");
      fprintf(DIAG_OUT," 1. Use more automatic local variables and \n");
      fprintf(DIAG_OUT,"    allocate/deallocate adoubles on free store\n");
      fprintf(DIAG_OUT,"     in a strictly last in first out fashion\n\n");
      fprintf(DIAG_OUT," 2. Extend the range by redefining the type of \n");
      fprintf(DIAG_OUT,"    locint (currently %d byte) from unsigned short (%d byte)  or int \n", sizeof(locint),sizeof(unsigned short));
      fprintf(DIAG_OUT,"    to int (%d byte) or long (%d byte). \n",sizeof(int),sizeof(long));
      exit(-3);
    }
    else
    { maxtop = maxtop2;
      if (maxtop == 2)
      { store = (double *)malloc(maxtop*sizeof(double));
	deminloc = maxloc;
      }
      else
        store = (double *)realloc((char *)store,maxtop*sizeof(double));
      if (store == 0) 
      { fprintf(DIAG_OUT,"\nADOL-C error:\n");
  	fprintf(DIAG_OUT,"Failure to reallocate storage for adouble values\n");
	fprintf(DIAG_OUT,"Possible remedies :\n\n");
       	fprintf(DIAG_OUT," 1. Use more automatic local variables and \n");
	fprintf(DIAG_OUT,"    allocate/deallocate adoubles on free store\n");
	fprintf(DIAG_OUT,"    in a strictly last in first out fashion\n");
	fprintf(DIAG_OUT," 2. Enlarge your system stacksize limit\n");
	exit(-3);
      }
    }
  }
  return current_top++;
}


/*--------------------------------------------------------------------------*/
/* Return the next #size free locations in "adouble" memory */
locint next_loc( int size )
{ /* First deallocate dead adoubles if they form a contiguous tail: */
#ifdef overwrite
  if (dealloc && dealloc+deminloc == current_top)    
  { /* olvo 980518 deleted write_death (see version griewank) */
    // if (trace_flag)
    //   write_death(deminloc, current_top - 1);
    current_top = deminloc ;
    dealloc     = 0; 
    deminloc    = maxloc;
  }
#endif
  if ((current_top+size) >= location_cnt) 
    location_cnt = current_top+size+1;
  while (location_cnt > maxtop) 
  { maxtop2 = ++maxtop*2 > maxloc ? maxloc : 2*maxtop;
    if (maxtop2 == maxloc)
    { fprintf(DIAG_OUT,"\nADOL-C error:  \n");
      fprintf(DIAG_OUT,"maximal number (%d) of live active variables exceeded\n\n", maxloc);
      fprintf(DIAG_OUT,"Possible remedies :\n\n");
      fprintf(DIAG_OUT," 1. Use more automatic local variables and \n");
      fprintf(DIAG_OUT,"    allocate/deallocate adoubles on free store\n");
      fprintf(DIAG_OUT,"    in a strictly last in first out fashion\n\n");
      fprintf(DIAG_OUT," 2. Extend the range by redefining the type of \n");
      fprintf(DIAG_OUT,"    locint (currently %d byte) from unsigned short (%d byte)  or int \n", sizeof(locint),sizeof(unsigned short));
      fprintf(DIAG_OUT,"    to int (%d byte) or long (%d byte). \n",sizeof(int),sizeof(long));
      exit(-3);
    }
    else
    { maxtop = maxtop2;
      if (maxtop == 2)
      { store = (double *)malloc(maxtop*sizeof(double));
	deminloc = maxloc;
      }
      else
      { /* Allocate the storage */
	double *temp;
	temp = (double *)malloc(maxtop*sizeof(double));
	if(temp == NULL)
        { fprintf(DIAG_OUT,"\nADOL-C error: cannot allocate %i bytes\n",maxtop*sizeof(double));
          exit (-1);
        }
        /* Copy over storage */
	for (int i=0; i<current_top; i++)
	  temp[i]=store[i];
        free((char*) store);
        store = temp;
      }
      if (store == 0) 
      { fprintf(DIAG_OUT,"\nADOL-C error:\n");
  	fprintf(DIAG_OUT,"Failure to reallocate storage for adouble values\n");
	fprintf(DIAG_OUT,"Possible remedies :\n\n");
       	fprintf(DIAG_OUT," 1. Use more automatic local variables and\n");
	fprintf(DIAG_OUT,"    allocate/deallocate adoubles on free store\n");
	fprintf(DIAG_OUT,"    in a strictly last in first out fashion\n");
	fprintf(DIAG_OUT," 2. Enlarge your system stacksize limit\n");
	exit(-3);
      }
    }
  }
#ifdef DEBUG
  fprintf (DIAG_OUT,"ADOL-C debug: Top is: %d\n ",current_top+size);
#endif
  locint return_val = current_top;
  current_top += size;
  return return_val;
}

/*--------------------------------------------------------------------------*/
/* Free a location in "adouble" memory */
inline void free_loc( locint old_loc )
{ ++dealloc;
  if (old_loc < deminloc)
    deminloc = old_loc;
}

/*--------------------------------------------------------------------------*/
/* Free #size locations in "adouble" memory */
void free_loc( int old_loc, int size )
{ dealloc+=size;
  if (old_loc < deminloc)
    deminloc = old_loc ;
}
  
/****************************************************************************/
/*                                                       STOCK MANIPULATION */

/*--------------------------------------------------------------------------*/
/* olvo 980623 version griewank */
void take_stock()
{
#ifdef TAPE_DOC /* olvo 980709 ??? this case might be useless */
  for (int res =0; res< current_top; res++) 
  { double coval = store[res];   // Avoid I/O of NaN's !
    if (coval == coval) 
    { // old: write_int_assign_d(res,coval);
      if (coval == 0)
      { put_op(assign_d_zero);
        put_locint(res);
      }
      else
        if (coval == 1.0)
        { put_op(assign_d_one);
          put_locint(res);
        }
        else
        { put_op(assign_d);
          put_locint(res);
          put_val(coval);
        }
 
      ++vs_ptr; 
      if (revalso)
        write_scaylor(store[res]);
    }
  }
#else /* usual case */
  // old: write_take_stock(current_top,store);
  locint space_left,
         vals_left = current_top,
         loc       = 0;
  double *vals     = store;
  space_left       = get_val_space();
  while (space_left < vals_left)
  { put_op(take_stock_op);
    put_locint(space_left);
    put_locint(loc);
    put_vals_p(vals,space_left);
    vals      += space_left;
    vals_left -= space_left;
    loc       += space_left;
    space_left = get_val_space();
  }
  if (vals_left > 0)
  { put_op(take_stock_op);
    put_locint(vals_left);
    put_locint(loc);
    put_vals_r(vals,vals_left);
  }
#endif
  trace_flag = 1;
}

/*--------------------------------------------------------------------------*/
/* olvo 980623 version griewank */
locint keep_stock()
{ if (location_cnt > 0) 
  { // old: write_death(0,location_cnt - 1);
    locint loc2 = location_cnt - 1;

    put_op(death_not);
    put_locint(0);
    put_locint(loc2);

    vs_ptr += location_cnt;
    if (revalso) 
      do 
        write_scaylor(store[loc2]);
      while(loc2-- > 0);
  }
  trace_flag = 0;
  return location_cnt;
}

/*----------------------------------------------------------------*/
/* The remaining routines define the badouble,adub,and adouble    */
/* routines.                                                      */
/*----------------------------------------------------------------*/

/****************************************************************************/
/*                                                             CONSTRUCTORS */

/*--------------------------------------------------------------------------*/
/* just a comment:
adub::adub( double coval )
{ location = next_loc();

  if (trace_flag) 
  { // old: write_int_assign_d(location,coval);
    if (coval == 0)
    { put_op(assign_d_zero);
      put_locint(location);   // = res
    }
    else
      if (coval == 1.0)
      { put_op(assign_d_one);
        put_locint(location); // = res
      }
      else
      { put_op(assign_d);
        put_locint(location); // = res
        put_val(coval);       // = coval
      }

    ++vs_ptr; 
    if (revalso)  
      write_scaylor(store[location]);
  }

  store[res] = coval;
}
*/

/*--------------------------------------------------------------------------*/
adouble::adouble()
{ location = next_loc();
}

/*--------------------------------------------------------------------------*/
adouble::adouble( double coval )
{ location = next_loc();

  if (trace_flag) 
  { // old:  write_int_assign_d(location,coval);
    if (coval == 0)
    { put_op(assign_d_zero);
      put_locint(location);   // = res
    }
    else
      if (coval == 1.0)
      { put_op(assign_d_one);
        put_locint(location); // = res
      }
      else
      { put_op(assign_d);
        put_locint(location); // = res
        put_val(coval);       // = coval
      }

    ++vs_ptr; 
    if (revalso)  
      write_scaylor(store[location]);
  }

  store[location] = coval;
}

/*--------------------------------------------------------------------------*/
adouble::adouble( const adouble& a )
{ location = next_loc();

  if (trace_flag) 
  { // old: write_int_assign_a(location,a.location);
    put_op(assign_a);
    put_locint(a.location);   // = arg
    put_locint(location);     // = res
  
    ++vs_ptr; 
    if (revalso)
      write_scaylor(store[location]);
  }

  store[location] = store[a.location];
}

/*--------------------------------------------------------------------------*/
adouble::adouble( const adub& a )
{ location = next_loc();

  if (trace_flag) 
  { // old:  write_int_assign_a(location,a.loc());
    put_op(assign_a);
    put_locint(a.loc());  // = arg
    put_locint(location); // = res
  
    ++vs_ptr; 
    if (revalso)
      write_scaylor(store[location]);
  }

  store[location] = store[a.loc()];
}

/*--------------------------------------------------------------------------*/
adouble::adouble( const along& a )
{ location = next_loc();

  if (trace_flag) 
  { // old: write_int_assign_a(location,a.loc());
    put_op(assign_a);
    put_locint(a.loc());  // = arg
    put_locint(location); // = res
  
    ++vs_ptr; 
    if (revalso)
      write_scaylor(store[location]);
  }

  store[location] = store[a.loc()];
}

/* rf, 03/17/01: add new constructor for use with avector::operator[]       */
/*--------------------------------------------------------------------------*/
adouble::adouble( const badouble &ba )
{ location = next_loc();
  operator= (ba);
}

/****************************************************************************/
/*                                                              DESTRUCTORS */

#ifdef overwrite
/*--------------------------------------------------------------------------*/
adouble::~adouble()
{ ++dealloc; 
  if (location < deminloc) 
    deminloc = location;
}

/*--------------------------------------------------------------------------*/
adub::~adub()
{ ++dealloc; 
  if (location < deminloc)
    deminloc = location;
}

/*--------------------------------------------------------------------------*/
asub::~asub()
{ ++dealloc;
  if (location < deminloc)
    deminloc = location;
}

/*--------------------------------------------------------------------------*/
along::~along()
{ ++dealloc; 
  if (location < deminloc)
    deminloc = location;
}
#endif


/****************************************************************************/
/*                                                              ASSIGNMENTS */

/*--------------------------------------------------------------------------*/
/* Assign an adouble variable a constant value. */
badouble& badouble::operator = ( double coval ) 
{ if (trace_flag) 
  { // old:  write_assign_d(location,coval);
    if (coval == 0)
    { put_op(assign_d_zero);
      put_locint(location);   // = res
    }
    else
      if (coval == 1.0)
      { put_op(assign_d_one);
        put_locint(location); // = res
      }
      else
      { put_op(assign_d);
        put_locint(location); // = res
        put_val(coval);       // = coval
      }

    ++vs_ptr; 
    if (revalso)
      write_scaylor(store[location]);
  }

  store[location] = coval;
  return *this;
}   

/*--------------------------------------------------------------------------*/
/* Assign an adouble variable a constant value. */
adouble& adouble::operator = ( double coval ) 
{ (*this).badouble::operator=(coval);
  return (*this);
}

/*--------------------------------------------------------------------------*/
/* Assign an adouble variable to an independent value. */
badouble& badouble::operator <<= ( double coval ) 
{ if (trace_flag) 
  { // old:  write_assign_ind(location);
    ind_ptr++;

    put_op(assign_ind);
    put_locint(location); // = res

    ++vs_ptr;
    if (revalso)
      write_scaylor(store[location]);
  }
 
  store[location] = coval;
  return *this;
}   

/*--------------------------------------------------------------------------*/
/* Assign a float variable from a dependent adouble value. */
badouble& badouble::operator >>= ( double& coval ) 
{ if (trace_flag) 
  { // old:  write_assign_dep(location);
    dep_ptr++;

    put_op(assign_dep);
    put_locint(location); // = res
  }

  coval = double (store[location]);
  return *this;
}   

/*--------------------------------------------------------------------------*/
/* Assign an Badouble variable an Badouble value. */
badouble& badouble::operator = ( const badouble& x ) 
{ locint x_loc = x.loc();
  if (location!=x_loc)  
  /* test this to avoid for x=x statements adjoint(x)=0 in reverse mode */
  { if (trace_flag) 
    { // old:  write_assign_a(location,x.location);
      put_op(assign_a);
      put_locint(x_loc);    // = arg
      put_locint(location);   // = res

      ++vs_ptr;
      if (revalso)
        write_scaylor(store[location]);
    }

    store[location]=store[x_loc];
  } 
  return *this;
}  

/*--------------------------------------------------------------------------*/
/* Assign an Badouble variable an Badouble value. */
adouble& adouble::operator = ( const badouble& x ) 
{ (*this).badouble::operator=(x);
  return (*this);
}

/*--------------------------------------------------------------------------*/
/* Assign an adouble an adub */
/* olvo 980517 new version griewank */
badouble& badouble::operator = ( const adub& a )
{ locint a_loc = a.loc();
  int upd = 0;
  /* 981020 olvo  skip upd_resloc(..) if no tracing performed */
  if (trace_flag)
    upd = upd_resloc(a_loc,location);
  if (upd)
  { /* olvo 980708 new n2l & 980921 changed interface */
    revreal tempVal = store[a_loc];
    if (revalso) 
      overwrite_scaylor(store[location],&store[a_loc]); 
    if (a_loc == current_top-1)
    { current_top--;     // The temporary will die in a minute and
      dealloc--;         // by reducing dealloc and current_top 
    }                    // we neutralize that effect
    store[location] = tempVal;
  }
  else
  { if (trace_flag)
    { // old: write_assign_a(location,a_loc);
      put_op(assign_a);
      put_locint(a_loc);    // = arg
      put_locint(location); // = res

      ++vs_ptr;
      if (revalso) 
        write_scaylor(store[location]);
    }
    store[location] = store[a_loc];
  }

  return *this;
}

/*--------------------------------------------------------------------------*/
/* Assign an adouble an adub */
/* olvo 980517 new version griewank */
adouble& adouble::operator = ( const adub& a )
{ (*this).badouble::operator=(a);
  return (*this);
}


/****************************************************************************/
/*                                                           INPUT / OUTPUT */

/*--------------------------------------------------------------------------*/
/* Output an adouble value !!! No tracing of this action */
ostream& operator << ( ostream& out, const badouble& y )
{ return out << store[y.location] << "(a)" ;
}

/*--------------------------------------------------------------------------*/
/* Input adouble value */
istream& operator >> ( istream& in, const badouble& y )
{ double coval;
  in >> coval;
  if (trace_flag) 
  { // old: write_assign_d(y.location,coval);
    if (coval == 0)
    { put_op(assign_d_zero);
      put_locint(y.location);   // = res
    }
    else
      if (coval == 1.0)
      { put_op(assign_d_one);
	put_locint(y.location); // = res
      }
      else
      { put_op(assign_d);
      put_locint(y.location);   // = res
        put_val(coval);         // = coval
      }

    ++vs_ptr; 
    if (revalso)  
      write_scaylor(store[y.location]);
  }

  store[y.location] = coval;
  return in;
}
  
/****************************************************************************/
/*                                                    INCREMENT / DECREMENT */

/*--------------------------------------------------------------------------*/
/* Postfix increment */
adub adouble::operator++( int ) 
{ locint locat = next_loc();

  if (trace_flag) 
  { // old: write_assign_a(locat,location);
    put_op(assign_a);
    put_locint(location); // = arg
    put_locint(locat);    // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat]=store[location];

  if (trace_flag) 
  { // old: write_incr_decr_a(incr_a,location);
    put_op(incr_a);
    put_locint(location); // = res

    ++vs_ptr;
    if (revalso)
      write_scaylor(store[location]);
  }

  store[location]++;
  return locat;
}

/*--------------------------------------------------------------------------*/
 /* Postfix decrement */
adub adouble::operator--( int )
{ locint locat = next_loc();

  if (trace_flag) 
  { // old: write_assign_a(locat,location);
    put_op(assign_a);
    put_locint(location); // = arg
    put_locint(locat);    // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat]=store[location];
  if (trace_flag) 
  { // old: write_incr_decr_a(decr_a,location);
    put_op(decr_a);
    put_locint(location); // = res

    ++vs_ptr;
    if (revalso)
      write_scaylor(store[location]);
  }

  store[location]--;
  return locat;
}

/*--------------------------------------------------------------------------*/
 /* Prefix increment */
badouble& adouble::operator++()
{ if (trace_flag) 
  { // old: write_incr_decr_a(incr_a,location);
    put_op(incr_a);
    put_locint(location); // = res

    ++vs_ptr;
    if (revalso)
      write_scaylor(store[location]);
  }

  store[location]++;
  return *this;
}

/*--------------------------------------------------------------------------*/
/* Prefix decrement */
badouble& adouble::operator--()
{ if (trace_flag) 
  { // old: write_incr_decr_a(decr_a,location);
    put_op(decr_a);
    put_locint(location); // = res

    ++vs_ptr;
    if (revalso)
      write_scaylor(store[location]);
  }

  store[location]--;
  return *this;
}

/****************************************************************************/
/*                                                   OPERATION + ASSIGNMENT */

/*--------------------------------------------------------------------------*/
/* Adding a floating point to an adouble */
badouble& badouble::operator += ( double coval ) 
{ if (trace_flag) 
  { // old: write_d_same_arg(eq_plus_d,location,coval);
    put_op(eq_plus_d);
    put_locint(location); // = res
    put_val(coval);       // = coval

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[location]);
  }

  store[location] += coval;
  return *this; 
} 


/*--------------------------------------------------------------------------*/
/* Subtracting a floating point from an adouble */
badouble& badouble::operator -= ( double coval ) 
{ if (trace_flag) 
  { // old: write_d_same_arg(eq_min_d,location,coval);
    put_op(eq_min_d);
    put_locint(location); // = res
    put_val(coval);       // = coval

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[location]);
  }

  store[location] -= coval;
  return *this;
}

/*--------------------------------------------------------------------------*/
/* Add an adouble to another adouble */
badouble& badouble::operator += ( const badouble& y ) 
{ locint y_loc = y.loc();
  if (trace_flag) 
  { // old: write_a_same_arg(eq_plus_a,location,y.location);
    put_op(eq_plus_a);
    put_locint(y_loc); // = arg
    put_locint(location);   // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[location]);
  }

  store[location] += store[y_loc];
  return *this;
}

/*--------------------------------------------------------------------------*/
/* olvo 991122 new version for y += x1 * x2; */
badouble& badouble::operator += ( const adub& a )
{ locint a_loc = a.loc();
  int upd = 0;
  if (trace_flag)
    upd = upd_resloc_inc_prod(a_loc,location,eq_plus_prod);
  if (upd)
  { store[location] += store[a_loc];
    if (revalso) 
      delete_scaylor(&store[a_loc]); 
    if (a_loc == current_top-1)
    { current_top--;     // The temporary will die in a minute and
      dealloc--;         // by reducing dealloc and current_top 
    }                    // we neutralize that effect
    --vs_ptr;
  }
  else
  { if (trace_flag)
    { // old: write_assign_a(location,a_loc);
      put_op(eq_plus_a);
      put_locint(a_loc);    // = arg
      put_locint(location); // = res

      ++vs_ptr;
      if (revalso) 
        write_scaylor(store[location]);
    }
    store[location] += store[a_loc];
  }

  return *this;
}

/*--------------------------------------------------------------------------*/
/* Subtract an adouble from another adouble */
badouble& badouble::operator -= ( const badouble& y ) 
{ locint y_loc = y.loc();
  if (trace_flag) 
  { // old: write_a_same_arg(eq_min_a,location,y.location);
    put_op(eq_min_a);
    put_locint(y_loc); // = arg
    put_locint(location);   // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[location]);
  }

  store[location] -= store[y_loc];
  return *this;
}

/*--------------------------------------------------------------------------*/
/* olvo 991122 new version for y -= x1 * x2; */
badouble& badouble::operator -= ( const adub& a )
{ locint a_loc = a.loc();
  int upd = 0;
  if (trace_flag)
    upd = upd_resloc_inc_prod(a_loc,location,eq_min_prod);
  if (upd)
  { store[location] -= store[a_loc];
    if (revalso) 
      delete_scaylor(&store[a_loc]); 
    if (a_loc == current_top-1)
    { current_top--;     // The temporary will die in a minute and
      dealloc--;         // by reducing dealloc and current_top 
    }                    // we neutralize that effect
    --vs_ptr;
  }
  else
  { if (trace_flag)
    { // old: write_assign_a(location,a_loc);
      put_op(eq_min_a);
      put_locint(a_loc);    // = arg
      put_locint(location); // = res

      ++vs_ptr;
      if (revalso) 
        write_scaylor(store[location]);
    }
    store[location] -= store[a_loc];
  }

  return *this;
}

/*--------------------------------------------------------------------------*/
/* Multiply an adouble by a floating point */
badouble& badouble::operator *= ( double coval ) 
{ if (trace_flag) 
  { // old: write_d_same_arg(eq_mult_d,location,coval);
    put_op(eq_mult_d);
    put_locint(location); // = res
    put_val(coval);       // = coval

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[location]);
  }

  store[location] *= coval;
  return *this;
}

/*--------------------------------------------------------------------------*/
/* Multiply one adouble by another adouble*/
badouble& badouble::operator *= ( const badouble& y ) 
{ locint y_loc = y.loc();
  if (trace_flag) 
  { // old: write_a_same_arg(eq_mult_a,location,y.location);
    put_op(eq_mult_a);
    put_locint(y_loc); // = arg
    put_locint(location);   // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[location]);
  }

  store[location] *= store[y_loc];
  return *this;
}

/*--------------------------------------------------------------------------*/
badouble& badouble::operator /= (double y) 
{ *this = *this/y;
  return *this;
}

/*--------------------------------------------------------------------------*/
badouble& badouble::operator /= (const badouble& y) 
{ *this = *this * (1.0/y);
  return *this;
}

/****************************************************************************/
/*                                                               COMPARISON */
/* olvo 980819 NOTE: new comparison strategy !!! */

/*--------------------------------------------------------------------------*/
/*   The Not Equal Operator (!=) */
int operator != ( const badouble& v, double coval )
{ if (coval)
  /* rf, 07/18/01: introduce additional local variable to circumvent
     bug of gcc 2.9?? with code optimization. Original:
     return (-coval+v != 0); */
  {
    adouble a = -coval+v;
    return (a != 0);
  }
  else
  { if (trace_flag)
    { put_op(store[v.location] ? neq_zero : eq_zero);
      put_locint(v.location);
    }
    return (store[v.location] != 0);
  }
}

/*--------------------------------------------------------------------------*/
/*   The Equal Operator (==) */
int operator == ( const badouble& v, double coval)
{ if (coval)
  /* rf, 07/18/01: introduce additional local variable to circumvent
     bug of gcc 2.9?? with code optimization. Original:
     return (-coval+v == 0); */
  {
    adouble a = -coval+v;
    return (a == 0);
  }
  else
  { if (trace_flag)
    { put_op(store[v.location] ? neq_zero : eq_zero);
      put_locint(v.location);
    }
    return (store[v.location] == 0);
  }
}

/*--------------------------------------------------------------------------*/
/*   The Less than or Equal Operator (<=)      */
int operator <= ( const badouble& v, double coval )
{ if (coval)
  /* rf, 07/18/01: introduce additional local variable to circumvent
     bug of gcc 2.9?? with code optimization. Original:
     return (-coval+v <= 0); */
  {
    adouble a = -coval+v;
    return (a <= 0);
  }
  else
  { int b = (store[v.location] <= 0);
    if (trace_flag)
    { put_op(b ? le_zero : gt_zero);
      put_locint(v.location);
    }
    return b;
  }
}

/*--------------------------------------------------------------------------*/
/*   The Greater than or Equal Operator (>=)      */
int operator >= ( const badouble& v, double coval )
{ if (coval)
  /* rf, 07/18/01: introduce additional local variable to circumvent
     bug of gcc 2.9?? with code optimization. Original:
     return (-coval+v >= 0); */
  {
    adouble a = -coval+v;
    return (a >= 0);
  }
  else
  { int b = (store[v.location] >= 0);
    if (trace_flag)
    { put_op(b ? ge_zero : lt_zero);
      put_locint(v.location);
    }
    return b;
  }
}

/*--------------------------------------------------------------------------*/
/*   The Greater than Operator (>)      */
int operator > ( const badouble& v, double coval )
{ if (coval)
  /* rf, 07/18/01: introduce additional local variable to circumvent
     bug of gcc 2.9?? with code optimization. Original:
     return (-coval+v > 0); */
  {
    adouble a = -coval+v;
    return (a > 0);
  }
  else
  { int b = (store[v.location] > 0);
    if (trace_flag)
    { put_op(b ? gt_zero : le_zero);
      put_locint(v.location);
    }
    return b;
  }
}

/*--------------------------------------------------------------------------*/
/*   The Less than Operator (<)      */
int operator < ( const badouble& v, double coval )
{ if (coval)
  /* rf, 07/18/01: introduce additional local variable to circumvent
     bug of gcc 2.9?? with code optimization. Original:
     return (-coval+v < 0); */
  {
    adouble a = -coval+v;
    return (a < 0);
  }
  else
  { int b = (store[v.location] < 0);
    if (trace_flag)
    { put_op(b ? lt_zero : ge_zero);
      put_locint(v.location);
    }
    return b;
  }
}


/****************************************************************************/
/*                                                          SIGN  OPERATORS */

/*--------------------------------------------------------------------------*/
/* olvo 980709 modified positive sign operator 
   ??? possibly there is a better way */
adub operator + ( const badouble& x )
{ locint locat = next_loc();

  if (trace_flag)
  { // old: write_pos_sign_a(locat,x.location);
    put_op(pos_sign_a);
    put_locint(x.location); // = arg
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = store[x.location];  
  return locat;
} 

/*--------------------------------------------------------------------------*/
/* olvo 980709 modified negative sign operator */
adub operator - ( const badouble& x )
{ locint locat = next_loc();

  if (trace_flag)
  { // old: write_neg_sign_a(locat,x.location);
    put_op(neg_sign_a);
    put_locint(x.location); // = arg
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = -store[x.location];  
  return locat;
} 


/****************************************************************************/
/*                                                         BINARY OPERATORS */

/* NOTE: each operator calculates address of temporary  and returns
         an adub */

/*--------------------------------------------------------------------------*/
/* Adding two adoubles */
adub operator + ( const badouble& x, const badouble& y ) 
{ locint locat = next_loc();

  if (trace_flag) 
  { // old: write_two_a_rec(plus_a_a,locat,x.location,y.location);
    put_op(plus_a_a);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = store[x.location] + store[y.location];
  return locat;
}

/*--------------------------------------------------------------------------*/
/* Adding a adouble and a floating point */
adub operator + ( double coval, const badouble& y ) 
{ locint locat = next_loc(); 

  /* olvo 980708 test coval to be zero */
  if (coval)
  { if (trace_flag)
    { // old: write_args_d_a(plus_d_a,locat,coval,y.location);
      put_op(plus_d_a);
      put_locint(y.location); // = arg
      put_locint(locat);      // = res
      put_val(coval);         // = coval

      ++vs_ptr;
      if (revalso)  
        write_scaylor(store[locat]);
    }

    store[locat] = coval + store[y.location];
  } 
  else 
  { if (trace_flag)
    { // old: write_pos_sign_a(locat,y.location);
      put_op(pos_sign_a);
      put_locint(y.location); // = arg
      put_locint(locat);      // = res

      ++vs_ptr;
      if (revalso) 
        write_scaylor(store[locat]);
    }

    store[locat] = store[y.location];
  }

  return locat;
}

/*--------------------------------------------------------------------------*/
adub operator + ( const badouble& y, double coval) 
{ locint locat = next_loc(); 

  /* olvo 980708 test coval to be zero */
  if (coval)
  { if (trace_flag)
    { // old: write_args_d_a(plus_d_a,locat,coval,y.location);
      put_op(plus_d_a);
      put_locint(y.location); // = arg
      put_locint(locat);      // = res
      put_val(coval);         // = coval

      ++vs_ptr;
      if (revalso)  
        write_scaylor(store[locat]);
    }

    store[locat] = coval + store[y.location];
  } 
  else 
  { if (trace_flag)
    { // old: write_pos_sign_a(locat,y.location);
      put_op(pos_sign_a);
      put_locint(y.location); // = arg
      put_locint(locat);      // = res

      ++vs_ptr;
      if (revalso) 
        write_scaylor(store[locat]);
    }

    store[locat] = store[y.location];
  }

  return locat;
}

/*--------------------------------------------------------------------------*/
/* Subtraction of two adoubles */
adub operator - ( const badouble& x, const badouble& y )
{ locint locat = next_loc();

  if (trace_flag) 
  { // old: write_two_a_rec(min_a_a,locat,x.location,y.location);
    put_op(min_a_a);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = store[x.location] - store[y.location];
  return locat;
}


/*--------------------------------------------------------------------------*/
/* Subtract an adouble from a floating point */
adub operator - ( double coval, const badouble& y )
{ locint locat = next_loc();

  /* olvo 980708 test coval to be zero */
  if (coval)
  { if (trace_flag) 
    { // old: write_args_d_a(min_d_a,locat,coval,y.location);
      put_op(min_d_a);
      put_locint(y.location); // = arg
      put_locint(locat);      // = res
      put_val(coval);         // = coval

      ++vs_ptr;
      if (revalso)  
        write_scaylor(store[locat]);
    }

    store[locat] = coval - store[y.location];
  }
  else
  { if (trace_flag)
    { // old: write_neg_sign_a(locat,y.location);
      put_op(neg_sign_a);
      put_locint(y.location); // = arg
      put_locint(locat);      // = res

      ++vs_ptr;
      if (revalso) 
        write_scaylor(store[locat]);
    }

    store[locat] = -store[y.location];  
  }

  return locat; 
}

/*--------------------------------------------------------------------------*/
/* Multiply two adoubles */
adub operator * ( const badouble& x, const badouble& y )
{ locint locat = next_loc();

  if (trace_flag) 
  { // old: write_two_a_rec(mult_a_a,locat,x.location,y.location);
    put_op(mult_a_a);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = store[x.location] * store[y.location];
  return locat;
}

/*--------------------------------------------------------------------------*/
/* Multiply an adouble by a floating point */ 
/* olvo 980709 modified */
adub operator * ( double coval, const badouble& y )
{ locint locat = next_loc();

  if ( coval == 1.0 )
  { if (trace_flag)
    { // old: write_pos_sign_a(locat,y.location);
      put_op(pos_sign_a);
      put_locint(y.location); // = arg
      put_locint(locat);      // = res

      ++vs_ptr;
      if (revalso) 
        write_scaylor(store[locat]);
    }

    store[locat] = store[y.location];  
  }
  else
    if ( coval == -1.0 )
    { if (trace_flag)
      { // old: write_neg_sign_a(locat,y.location);
        put_op(neg_sign_a);
        put_locint(y.location); // = arg
        put_locint(locat);      // = res

        ++vs_ptr;
        if (revalso) 
          write_scaylor(store[locat]);
      }

      store[locat] = -store[y.location];  
    }
    else
    { if (trace_flag) 
      { // old: write_args_d_a(mult_d_a,locat,coval,y.location);
        put_op(mult_d_a);
        put_locint(y.location); // = arg
        put_locint(locat);      // = res
        put_val(coval);         // = coval

        ++vs_ptr;
        if (revalso)  
          write_scaylor(store[locat]);
      }
 
      store[locat] = coval * store[y.location];
    }
  return locat;
}

/*--------------------------------------------------------------------------*/
/* Divide an adouble by another adouble */
adub operator / ( const badouble& x, const badouble& y )
{ locint locat = next_loc();

  if (trace_flag)
  { // old: write_two_a_rec(div_a_a,locat,x.location,y.location);
    put_op(div_a_a);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = store[x.location] / store[y.location];
  return locat;
}

/*--------------------------------------------------------------------------*/
/* Division floating point - adouble */
adub operator / ( double coval, const badouble& y )
{ locint locat = next_loc();

  if (trace_flag)
  { // old: write_args_d_a(div_d_a,locat,coval,y.location);
    put_op(div_d_a);
    put_locint(y.location); // = arg
    put_locint(locat);      // = res
    put_val(coval);         // = coval

    ++vs_ptr;
    if (revalso)  
      write_scaylor(store[locat]);
  }

  store[locat] = coval  / store[y.location];
  return locat;
}


/****************************************************************************/
/*                                                        SINGLE OPERATIONS */

/*--------------------------------------------------------------------------*/
/* Compute exponential of adouble */
adub exp ( const badouble& x )
{ locint locat = next_loc();

  if (trace_flag)
  { // old: write_single_op(exp_op,locat,x.location);
    put_op(exp_op);
    put_locint(x.location); // = arg
    put_locint(locat);      // = res   

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = exp(store[x.location]);
  return locat; 
}

/*--------------------------------------------------------------------------*/
/* Compute logarithm of adouble */
adub log ( const badouble& x )
{ locint locat = next_loc();

  if (trace_flag)
  { // old: write_single_op(log_op,locat,x.location);
    put_op(log_op);
    put_locint(x.location); // = arg
    put_locint(locat);      // = res   

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = log(store[x.location]);
  return locat; 
}

/*--------------------------------------------------------------------------*/
/* Compute sqrt of adouble */
adub sqrt ( const badouble& x )
{ locint locat = next_loc();

  if (trace_flag)
  { // old: write_single_op(sqrt_op,locat,x.location);
    put_op(sqrt_op);
    put_locint(x.location); // = arg
    put_locint(locat);      // = res   

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = sqrt(store[x.location]);
  return locat; 
}

/****************************************************************************/
/*                                                          QUAD OPERATIONS */

/*--------------------------------------------------------------------------*/
/* Compute sin of adouble
   !!! Sin and Cos are always evaluated together
*/
adub sin ( const badouble& x ) 
{ locint locat = next_loc(); 

  adouble y;

  if (trace_flag) 
  { // old: write_quad(sin_op,locat,x.location,y.location);
    put_op(sin_op);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    vs_ptr += 2;
    if (revalso) 
    { /* olvo 980921 changed order */
      write_scaylor(store[y.location]);
      write_scaylor(store[locat]);
    }
  }

  store[locat]      = sin(store[x.location]);
  store[y.location] = cos(store[x.location]);
  return locat;
}

/*--------------------------------------------------------------------------*/
/* Compute cos of adouble */
adub cos ( const badouble& x )
{ locint locat = next_loc();

  adouble y;

  if (trace_flag)
  { // old: write_quad(cos_op, locat,x.location,y.location);
    put_op(cos_op);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    vs_ptr += 2;
    if (revalso) 
    { /* olvo 980921 changed order */
      write_scaylor(store[y.location]);
      write_scaylor(store[locat]);
    }
  }

  store[locat]      = cos(store[x.location]);
  store[y.location] = sin(store[x.location]);
  return locat; 
}

/*--------------------------------------------------------------------------*/
/* Compute tan of adouble */
adub tan ( const badouble& x ) 
{ return sin(x) / cos(x);
}

/*--------------------------------------------------------------------------*/
/* Asin value -- really a quadrature */
adub asin ( const badouble& x )
{ locint locat = next_loc();

  adouble y = 1.0 / sqrt(1.0 - x*x);

  if (trace_flag)
  { // old:  write_quad(asin_op,locat,x.location,y.location);
    put_op(asin_op);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = asin(store[x.location]);
  return locat; 
}

/*--------------------------------------------------------------------------*/
/* Acos value -- really a quadrature */
adub acos ( const badouble& x )
{ locint locat = next_loc();

  adouble y = -1.0 / sqrt(1.0 - x*x);

  if (trace_flag)
  { // old: write_quad(acos_op,locat,x.location,y.location);
    put_op(acos_op);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = acos(store[x.location]);
  return locat; 
}

/*--------------------------------------------------------------------------*/
/* Atan value -- really a quadrature */
adub atan ( const badouble& x )
{ locint locat = next_loc();

  adouble y = 1.0 / (1.0 + x*x);

  if (trace_flag)
  { // old: write_quad(atan_op,locat,x.location,y.location);
    put_op(atan_op);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = atan(store[x.location]);
  return locat; 
}

/*--------------------------------------------------------------------------*/
adouble atan2( const badouble& y, const badouble& x)
{ adouble a1, a2, ret, sy;
  const double pihalf = asin(1.0);
  /* y+0.0 is a hack since condassign is currently not defined for 
     badoubles */
  condassign( sy,  y+0.0,  1.0 , -1.0 ); 
  condassign( a1,  x+0.0, (adouble) atan(y/x), 
                           (adouble)( atan(y/x)+sy*2*pihalf));
  condassign( a2,  (adouble) fabs(y), (adouble) (sy*pihalf-atan(x/y)),
                                      (adouble) 0.0 );
  condassign( ret, (adouble) (fabs(x) - fabs(y)), a1, a2 );
  return ret;
}

/*--------------------------------------------------------------------------*/
/* power value -- adouble ^ floating point */
adub pow ( const badouble& x, double coval )
{ locint locat = next_loc();

  if (trace_flag)
  { // old: write_args_d_a(pow_op,locat,cocval,x.location);
    put_op(pow_op);
    put_locint(x.location); // = arg
    put_locint(locat);      // = res
    put_val(coval);         // = coval

    ++vs_ptr;
    if (revalso)  
       write_scaylor(store[locat]);
  }

  store[locat] = pow(store[x.location],coval);
  return locat;
}

/*--------------------------------------------------------------------------*/
/* power value --- floating point ^ adouble */
adouble pow ( double coval, const badouble& y )
{ adouble ret;

  if (coval <= 0)
  { fprintf(DIAG_OUT,"\nADOL-C message:  exponent at zero/negative constant basis deactivated\n");
  }

  condassign (ret, coval, exp(y*log(coval)), pow(coval,value(y)) );

  return ret;
}

/*--------------------------------------------------------------------------*/
/* power value --- adouble ^ adouble */
adouble pow ( const badouble& x, const badouble& y)  
{ adouble a1, a2, ret;
  double vx = value(x);
  double vy = value(y);

  if (!(vx > 0))
    if (vx < 0 || vy >= 0)
      fprintf(DIAG_OUT,"\nADOL-C message: exponent of zero/negative basis deactivated\n");
    else 
      fprintf(DIAG_OUT,"\nADOL-C message: negative exponent and zero basis deactivated\n"); 

  condassign(a1,-y,pow(vx,vy),pow(x,vy));
  condassign(a2,fabs(x),pow(x, vy),a1);
  condassign(ret,x+0.0,exp(y*log(x)),a2);

  return ret;  
}

/*--------------------------------------------------------------------------*/
/* log base 10 of an adouble */
adub log10 ( const badouble& x ) 
{ return log(x) / log(10.0);
}

/*--------------------------------------------------------------------------*/
/* Hyperbolic Sine of an adouble */
/* 981119 olvo changed as J.M. Aparicio suggested */
adub sinh ( const badouble& x ) 
{ if (value(x) < 0.0)
  { adouble temp = exp(x);
    return  0.5*(temp - 1.0/temp);
  }
  else
  { adouble temp = exp(-x);
    return 0.5*(1.0/temp - temp);
  }
}

/*--------------------------------------------------------------------------*/
/* Hyperbolic Cosine of an adouble */
/* 981119 olvo changed as J.M. Aparicio suggested */
adub cosh ( const badouble& x ) 
{ adouble temp = (value(x) < 0.0) ? exp(x) : exp(-x);
  return 0.5*(temp + 1.0/temp);
}

/*--------------------------------------------------------------------------*/
/*
  Hyperbolic Tangent of an adouble value.
*/
/* 981119 olvo changed as J.M. Aparicio suggested */
adub tanh ( const badouble& x ) 
{ if (value(x) < 0.0)
  { adouble temp = exp(2.0*x);
    return (temp - 1.0)/(temp + 1.0);
  }
  else
  { adouble temp = exp((-2.0)*x);
    return (1.0 - temp)/(temp + 1.0);
  }
}

/*--------------------------------------------------------------------------*/
/* Ceiling function (NOTE: This function is nondifferentiable) */
adub ceil ( const badouble& x ) 
{ locint locat=next_loc();
  
  double coval = ceil(store[x.location]);

  if (trace_flag)
  { // old: write_args_d_a(ceil_op,locat,coval,x.location);
    put_op(ceil_op);
    put_locint(x.location); // = arg
    put_locint(locat);      // = res
    put_val(coval);         // = coval

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = coval;
  return locat;
}

/*--------------------------------------------------------------------------*/
/* Floor function (NOTE: This function is nondifferentiable) */
adub floor ( const badouble& x ) 
{ locint locat=next_loc();

  double coval = floor(store[x.location]);

  if (trace_flag)
  { // old: write_args_d_a(floor_op,locat,coval,x.location);
    put_op(floor_op);
    put_locint(x.location); // = arg
    put_locint(locat);      // = res
    put_val(coval);         // = coval

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = coval;
  return locat;
}

#ifdef ATRIG_ERF
/* NOTE: enable if your compiler knows asinh, acosh, atanh, erf */

/*--------------------------------------------------------------------------*/
/* Asinh value -- really a quadrature */
adub asinh ( const badouble& x )
{ locint locat = next_loc();

  adouble y = 1.0 / sqrt(1.0 + x*x);

  if (trace_flag) 
  { // old: write_quad(asinh_op,locat,x.location,y.location);
    put_op(asinh_op);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = asinh(store[x.location]);
  return locat; 
} 

/*--------------------------------------------------------------------------*/
/* Acosh value -- really a quadrature */
adub acosh ( const badouble& x )
{ locint locat = next_loc();

  adouble y = 1.0 / sqrt(1.0 - x*x);

  if (trace_flag)
  { // old: write_quad(acosh_op,locat,x.location,y.location);
    put_op(acosh_op);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = acosh(store[x.location]);
  return locat; 
}

/*--------------------------------------------------------------------------*/
/* Atanh value -- really a quadrature */
adub atanh ( const badouble& x )
{ locint locat = next_loc();

  adouble y = 1.0 / (1.0 - x*x);

  if (trace_flag)
  { // old: write_quad(atanh_op,locat,x.location,y.location);
    put_op(atanh_op);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat] = atanh(store[x.location]);
  return locat; 
}

/*--------------------------------------------------------------------------*/
/*  The error function erf */
adub erf( const badouble& x ) 
{ locint locat = next_loc();

  adouble y = exp(-x*x);

  if (trace_flag)
  { // old: write_quad(erf_op,locat,x.location,y.location);
    put_op(erf_op);
    put_locint(x.location); // = arg1
    put_locint(y.location); // = arg2
    put_locint(locat);      // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat]=erf(store[x.location]);
  return locat;
} 

#endif

/*--------------------------------------------------------------------------*/
/* Fabs Function (NOTE: This function is also nondifferentiable at x=0) */
adub fabs ( const badouble& x )
{ locint locat = next_loc();
  
  double coval = 1.0;
  double temp  = fabs(store[x.location]);
  if (temp != store[x.location])
    coval = 0.0;

  if (trace_flag)
  { /*  write_args_d_a(abs_val,locat,coval,x.location); */
    put_op(abs_val);
    put_locint(x.location);   /* arg */
    put_locint(locat);        /* res */ 
    put_val(coval);           /* coval */
 
    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }
  store[locat] = temp;
  return locat;
}

/*--------------------------------------------------------------------------*/
/* max and min functions  (changed : 11/15/95) */
adub fmin ( const badouble& x, const badouble& y )
{ /* olvo 980702 tested: return 0.5*fabs(x+y-fabs(x-y)); */
  locint locat = next_loc();

  if (store[y.location] < store[x.location])
  { if (trace_flag) 
    { // old: write_min_op(x.location,y.location,locat,0.0);
      put_op(min_op);
      put_locint(x.location); // = arg1
      put_locint(y.location); // = arg2
      put_locint(locat);      // = res
      put_val(0.0);           // = coval

      ++vs_ptr;
      if (revalso) 
        write_scaylor(store[locat]);
    }

    store[locat]=store[y.location];
  } 
  else 
  { if (trace_flag) 
    { // old: write_min_op(x.location,y.location,locat,1.0);
      put_op(min_op);
      put_locint(x.location); // = arg1
      put_locint(y.location); // = arg2
      put_locint(locat);      // = res
      put_val(1.0);           // = coval

      ++vs_ptr;
      if (revalso) 
        write_scaylor(store[locat]);
    }

    store[locat]=store[x.location];
  }
  return locat;
}

/*--------------------------------------------------------------------------*/
/*21.8.96*/
adub fmin ( double d, const badouble& y )
{ adouble x = d;
  return (fmin (x,y));
}

/*--------------------------------------------------------------------------*/
adub fmin ( const badouble& x, double d )
{ adouble y = d;
  return (fmin (x,y));
}

/*--------------------------------------------------------------------------*/
adub fmax ( const badouble& x, const badouble& y )
{ return (-fmin(-x,-y));
}

/*--------------------------------------------------------------------------*/
/*21.8.96*/
adub fmax ( double d, const badouble& y )
{ adouble x = d;
  return (-fmin(-x,-y));
}

/*--------------------------------------------------------------------------*/
adub fmax ( const badouble& x, double d )
{ adouble y = d;
  return (-fmin(-x,-y));
}

/*--------------------------------------------------------------------------*/
/* Ldexp Function */
adub ldexp ( const badouble& x, int exp ) 
{ return x*ldexp(1.0,exp);
}

/*--------------------------------------------------------------------------*/
/* Macro for user defined quadratures, example myquad is below.*/
/* the forward sweep tests if the tape is executed exactly at  */
/* the same argument point otherwise it stops with a returnval */
#define extend_quad(func,integrand)\
adouble func ( const badouble& arg )\
{  adouble temp; \
    adouble val; \
    integrand; \
    if (trace_flag) \
    { put_op(gen_quad); \
      put_locint(arg.location); \
      put_locint(val.location); \
      put_locint(temp.location); \
      ++vs_ptr; \
      if (revalso) \
        write_scaylor(store[temp.location]); \
    } \
    store[temp.location]=func(store[arg.location]); \
    if (trace_flag) \
    { put_val(store[arg.location]); \
      put_val(store[temp.location]); \
    } \
    return temp; }

double myquad(double& x)
{
  double res;
  res = log(x);
  return res;
}

/* This defines the natural logarithm as a quadrature */

extend_quad(myquad,val = 1/arg)


/****************************************************************************/
/*                                                             CONDITIONALS */

/* For the time being condassign is defined using adoubles in two 
   versions with adouble and along as left hand side.  This implies 
   some problems when badoubles are used as arguments, e.g. inside 
   the pow definition. For later versions we will replace this with
   complete definition for all parameter type constellations */

/*--------------------------------------------------------------------------*/
void condassign( adouble &res,        const adouble &cond, 
                 const adouble &arg1, const adouble &arg2 ) 
{ if (trace_flag)
  { // old: write_condassign(res.location,cond.location,arg1.location,
    //		     arg2.location);
    put_op(cond_assign);
    put_locint(cond.location); // = arg
    put_val(store[cond.location]);
    put_locint(arg1.location); // = arg1
    put_locint(arg2.location); // = arg2
    put_locint(res.location);  // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[res.location]);
  }

  if (store[cond.location] > 0)
    store[res.location] = store[arg1.location];
  else
    store[res.location] = store[arg2.location];
}

/*--------------------------------------------------------------------------*/
void condassign( adouble &res, const adouble &cond, const adouble &arg ) 
{ if (trace_flag)		
  { // old: write_condassign2(res.location,cond.location,arg.location);
    put_op(cond_assign_s);
    put_locint(cond.location); // = arg
    put_val(store[cond.location]);
    put_locint(arg.location);  // = arg1
    put_locint(res.location);  // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[res.location]);
  }

  if (store[cond.location] > 0)
    store[res.location] = store[arg.location];
}

/*--------------------------------------------------------------------------*/
void condassign( along &res, const adouble &cond, 
                 const adouble &arg1, const adouble &arg2 ) 
{ if (trace_flag)
  { // old: write_condassign(res.location,cond.location,arg1.location,
    //		     arg2.location);
    put_op(cond_assign);
    put_locint(cond.location); // = arg
    put_val(store[cond.location]);
    put_locint(arg1.location); // = arg1
    put_locint(arg2.location); // = arg2
    put_locint(res.location);  // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[res.location]);
  }

  if (store[cond.location] > 0)
    store[res.location] = store[arg1.location];
  else
    store[res.location] = store[arg2.location];
}

/*--------------------------------------------------------------------------*/
void condassign( along &res, const adouble &cond, const adouble &arg) 
{ if (trace_flag)		
  { // old: write_condassign2(res.location,cond.location,arg.location);
    put_op(cond_assign_s);
    put_locint(cond.location); // = arg
    put_val(store[cond.location]);
    put_locint(arg.location);  // = arg1
    put_locint(res.location);  // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[res.location]);
  }

  if (store[cond.location] > 0)
    store[res.location] = store[arg.location];
}

/****************************************************************************/
/*                                   SUBSCRIPTS (CONSTRUCTOR / ASSIGNMENTS) */

/*--------------------------------------------------------------------------*/
asub::asub(locint start, locint index)
{
#ifdef DEBUG
  fprintf(DIAG_OUT,"\nADOL-C debug: Constructing an asub with 2 arguments\n");
#endif  
  base   = start;
  offset = index;

  location = next_loc();

  if (trace_flag)
  { // old:write_associating_value(subscript,location,base,offset);
    put_op(subscript);
    put_locint(base);
    put_locint(offset);
    put_locint(location);
    put_val(store[offset]);

    /* olvo 980721 new n3l */
    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[location]);
  }

  store[location] = store[base+(int)store[offset]];  
}

/*--------------------------------------------------------------------------*/
asub& asub::operator <<= ( double coval ) 
{ locint res = base+(int)store[offset];
  
  if (trace_flag)
  { // old: write_assign_ind(location);
    ind_ptr++;
    put_op(assign_ind);
    put_locint(location); // = res

    ++vs_ptr;
    if (revalso)
      write_scaylor(store[location]);

  /* olvo 980711 necessary ??? 
  }
  store[location] = coval;  
  if (trace_flag)
  { */

    // old: write_associating_value(subscript_l,location,base,offset);
    put_op(subscript_l);
    put_locint(base);
    put_locint(offset);
    put_locint(location);
    put_val(store[offset]);

    ++vs_ptr;
    if (revalso) /* olvo 980711 ??? next line */
      write_scaylor(store[location]); 
      /* this is correct since location is the location of the copy which 
         already contains the value to be stored */
  }

  store[res] = coval;
  return *this;
}   

/*--------------------------------------------------------------------------*/
asub& asub::operator = ( double coval )
{ locint res = base+(int)store[offset];

  if (trace_flag)
  { // old: write_associating_value_ld(subscript_ld,coval,base,offset);
    put_op(subscript_ld);
    put_val(coval);
    put_locint(base);
    put_locint(offset);
    put_val(store[offset]);

   ++vs_ptr;
   if (revalso) 
     write_scaylor(store[res]);
  }

  store[res] = coval;
  return *this;
}

/*--------------------------------------------------------------------------*/
asub& asub::operator = ( const badouble& x ) 
{ locint res = base+(int)store[offset];

  if (trace_flag)
  { // old: write_associating_value(subscript_l,x.loc(),base,offset);
    put_op(subscript_l);
    put_locint(base);
    put_locint(offset);
    put_locint(x.loc());
    put_val(store[offset]);

    ++vs_ptr;
    if (revalso) /* olvo 980711 ??? next line */
      write_scaylor(store[x.loc()]);
      /* this is correct since location is the location of the copy which 
         already contains the value to be stored */
  }

  store[res]=store[x.loc()];
  return *this;
}   

/*--------------------------------------------------------------------------*/
/* 20.08.96 */
asub& asub::operator = ( const asub& x ) 
{ locint res = base+(int)store[offset];

  if (trace_flag)
  { // : write_associating_value(subscript_l,x.loc(),base,offset);
    put_op(subscript_l);
    put_locint(base);
    put_locint(offset);
    put_locint(x.loc());
    put_val(store[offset]);

    ++vs_ptr;
    if (revalso) /* olvo 980711 ??? next line */
      write_scaylor(store[x.loc()]);
      /* this is correct since location is the location of the copy which 
         already contains the value to be stored */
  }

  store[res]=store[x.loc()];
  return *this;
} 


/****************************************************************************/
/*                                      SUBSCRIPTS (OPERATION + ASSIGNMENT) */

/* olvo 980713 !!! seems to be a temporary version */
 
/*--------------------------------------------------------------------------*/
/* Sep/01/96 */
asub& asub::operator += ( double x )
{ *this = *this + x;
  return *this;
}

/*--------------------------------------------------------------------------*/
asub& asub::operator += ( const badouble& x )
{ *this = *this + x;
  return *this;
}

/*--------------------------------------------------------------------------*/
asub& asub::operator -= ( double x )
{ *this = *this - x;
  return *this;
}

/*--------------------------------------------------------------------------*/
asub& asub::operator -= ( const badouble& x )
{ *this = *this - x;
  return *this;
}

/*--------------------------------------------------------------------------*/
asub& asub::operator *= ( double x )
{ *this = *this * x;
  return *this;
}

/*--------------------------------------------------------------------------*/
asub& asub::operator *= ( const badouble& x )
{ *this = *this * x;
  return *this;
}

/*--------------------------------------------------------------------------*/
asub& asub::operator /= ( double x )
{ *this = *this / x;
  return *this;
}

/*--------------------------------------------------------------------------*/
asub& asub::operator /= ( const badouble& x )
{ *this = *this / x;
  return *this;
}

/****************************************************************************/
/*                                       SUBSCRIPTS (INCREMENT / DECREMENT) */

/*--------------------------------------------------------------------------*/
/* postfix increment */
adub asub::operator++( int )
{ locint locat = next_loc();

  if (trace_flag)
  { // old: write_assign_a(locat,location);
    put_op(assign_a);
    put_locint(location); // = arg
    put_locint(locat);    // = res

    ++vs_ptr;
    if (revalso)
      write_scaylor(store[locat]);
  }

  store[locat]=store[location]; /* location is the local copy of the
  asub for which this definition is invoked */
  *this = *this + 1;
  return locat ;
}

/*--------------------------------------------------------------------------*/
/* postfix decrement */
adub asub::operator--( int )
{ locint locat = next_loc();

  if (trace_flag)
  { // old: write_assign_a(locat,location);
    put_op(assign_a);
    put_locint(location); // = arg
    put_locint(locat);    // = res

    ++vs_ptr;
    if (revalso)
      write_scaylor(store[locat]);
  }

  store[locat]=store[location]; /* location is the local copy of the 
  asub for which this definition is invoked */
  *this = *this - 1;
  return locat ;
}

/*--------------------------------------------------------------------------*/
/* prefix increment */
asub& asub::operator++() 
{ *this = *this + 1;
  return *this;
}

/*--------------------------------------------------------------------------*/
/* prefix decrement */
asub& asub::operator--() 
{ *this = *this - 1;
  return *this;
}

/****************************************************************************/
/*                                                              ALONG STUFF */

/*--------------------------------------------------------------------------*/
along::along()
{ location = next_loc();
}

/*--------------------------------------------------------------------------*/
along& along::operator = ( int coval ) 
{ if (trace_flag)
  { // old: write_assign_d(location,coval);
    if (coval == 0)
    { put_op(assign_d_zero);
      put_locint(location);   // = res
    }
    else
      if (coval == 1.0)
      { put_op(assign_d_one);
        put_locint(location); // = res
      }
      else
      { put_op(assign_d);
        put_locint(location); // = res
        put_val(coval);
      }

    ++vs_ptr;
    if (revalso)  
      write_scaylor(store[location]);
  }

  store[location]=coval;
  return *this;
}   

/*--------------------------------------------------------------------------*/
along& along::operator = ( const badouble& x ) 
{ if (location != x.loc())  
  /* test this to avoid for x=x statements adjoint(x)=0 in reverse mode */
  { if (trace_flag)
    { // old:   write_assign_a(location,x.loc());
      put_op(assign_a);
      put_locint(x.loc());  // = arg
      put_locint(location); // = res

      ++vs_ptr;
      if (revalso)
        write_scaylor(store[location]);
    }

    store[location]=store[x.loc()];
  } 
  return *this;
}

/*--------------------------------------------------------------------------*/
along& along::operator = ( const along& x ) 
{ if (location != x.location) 
  /* test this to avoid for x=x statements adjoint(x)=0 in reverse mode */
  { if (trace_flag)
    { // old: write_assign_a(location,x.location);
      put_op(assign_a);
      put_locint(x.location); // = arg
      put_locint(location);   // = res

      ++vs_ptr;
      if (revalso)
        write_scaylor(store[location]);
    }

    store[location]=store[x.location];
  } 
  return *this;
}

/*--------------------------------------------------------------------------*/
along& along::operator = ( const adub& a )
{ if (location != a.loc())  
  /* test this to avoid for x=x statements adjoint(x)=0 in reverse mode */
  { if (trace_flag)
    { // old: write_assign_a(location,a.loc());
      put_op(assign_a);
      put_locint(a.loc());  // = arg
      put_locint(location); // = res

      ++vs_ptr;
      if (revalso)
        write_scaylor(store[location]);
    }

    store[location]=store[a.loc()] ;
  } 
  return *this;
}

/*--------------------------------------------------------------------------*/
along::along( int coval )
{ location = next_loc();

  if (trace_flag)
  { // old: write_int_assign_d(location,coval);
    if (coval == 0)
    { put_op(assign_d_zero);
      put_locint(location);   // = res
    }
    else
      if (coval == 1.0)
      { put_op(assign_d_one);
        put_locint(location); // = res
      }
      else
      { put_op(assign_d);
        put_locint(location); // = res
        put_val(coval);
      }

    ++vs_ptr; 
    if (revalso)  
      write_scaylor(store[location]);
  }

  store[location] = coval;
}

/*--------------------------------------------------------------------------*/
along::along( const along& a )
{ location = next_loc();

  if (trace_flag)
  { // old: write_int_assign_a(location,a.location);
    put_op(assign_a);
    put_locint(a.location); // = arg
    put_locint(location);   // = res
 
    ++vs_ptr; 
    if (revalso)
      write_scaylor(store[location]);
  }

  store[location]=store[a.location];
}

/*--------------------------------------------------------------------------*/
along::along( const adub& a )
{ location = next_loc();

  if (trace_flag)
  { // old: write_int_assign_a(location,a.loc());
    put_op(assign_a);
    put_locint(a.loc());  // = arg
    put_locint(location); // = res
 
    ++vs_ptr; 
    if (revalso)
      write_scaylor(store[location]);
  }

  store[location]=store[a.loc()];
}

/*--------------------------------------------------------------------------*/
/* postfix increment */
adub along::operator++( int ) 
{ locint locat = next_loc();

  if (trace_flag) 
  { // old: write_assign_a(locat,location);
    put_op(assign_a);
    put_locint(location); // = arg
    put_locint(locat);    // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat]=store[location];

  if (trace_flag) 
  { // old: write_incr_decr_a(incr_a,location); 
    put_op(incr_a);
    put_locint(location);

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[location]);
  }

  store[location]++;
  return locat ;
}

/*--------------------------------------------------------------------------*/
/* postfix decrement */
adub along::operator--( int ) 
{ locint locat = next_loc();

  if (trace_flag)
  { // old: write_assign_a(locat,location);
    put_op(assign_a);
    put_locint(location); // = arg
    put_locint(locat);    // = res

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[locat]);
  }

  store[locat]=store[location];

  if (trace_flag)
  { // old: write_incr_decr_a(decr_a,location);
    put_op(decr_a);
    put_locint(location);

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[location]);
  }

  store[location]--;
  return locat ;
}

/*--------------------------------------------------------------------------*/
/* prefix increment */
along& along::operator++()
{ if (trace_flag)
  { // old: write_incr_decr_a(incr_a,location);
    put_op(incr_a);
    put_locint(location);

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[location]);
  }

  store[location]++;
  return *this;
}

along& along::operator--() /* prefix decrement */
{ if (trace_flag) 
  { // old: write_incr_decr_a(decr_a,location);
    put_op(decr_a);
    put_locint(location);

    ++vs_ptr;
    if (revalso) 
      write_scaylor(store[location]);
  }

  store[location]--;
  return *this;
}


/****************************************************************************/
/*                                                                THAT'S ALL*/
#undef _ADOLC_SRC_
#undef _ADOUBLE_CPP_




