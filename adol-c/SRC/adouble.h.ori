#ifndef _ADOUBLE_H_
#define _ADOUBLE_H_
/*
   --------------------------------------------------------------
   File adouble.h of ADOL-C version 1.8.6         as of Jan/07/00
   --------------------------------------------------------------
   adouble.h contains the basis for the class of adouble
   included here are all the possible functions defined on
   the adouble class.  Notice that, as opposed to ealier versions,
   both the class adub and the class adouble are derived from a base
   class (badouble).  See below for further explanation.

   Last changes:
    20000107 olvo  iostream.h instaed of stream.h  
      991210 olvo  checking the changes
      991122 olvo  new op_codes eq_plus_prod eq_min_prod
                   for  y += x1 * x2
                   and  y -= x1 * x2  
      981201 olvo   last check: 
                    - taputil things changed, includes 
      980820 olvo   new comparison strategy & some inlines
      980709 olvo   modified sign operators

   --------------------------------------------------------------
*/


/****************************************************************************/
/*
  NOTICE that the purpose of the class adub is merely to avoid the 
  generation and recording of an extra return adouble for each elementary 
  operation and function call. The same result can be achieved much
  more elegantly with GNUs named return variables, which would also 
  achieve the desired last in first out pattern for adouble construction 
  and destruction.
*/


/****************************************************************************/
/*      D I S C L A I M E R 

The ADOL-C Software is provided under the following disclaimer:

NO WARRANTY.  The software was created in the course of a research
endeavor. It is not a commercial package.  The present version is
still in development, and is distributed "AS IS, WITH ALL DEFECTS."
By using the software, each user agrees to assume all responsibility
for any and all such use.  The authors and Argonne National Laboratory
are not aware that the software or the use thereof infringe any
proprietary right belonging to a third party.  However, NO WARRANTY,
CONDITION, OR REPRESENTATION OF ANY KIND, EXPRESS OR IMPLIED, is made
about the software, including without limitation any warranty of title,
noninfringement, merchantability, or fitness for a particular purpose,
by the authors or their affiliated institutions.

NO CONSEQUENTIAL DAMAGES.  Independent of the foregoing disclaimer
of warranties, each person that uses the software thereby agrees, that
NEITHER ARGONNE NATIONAL LABORATORY NOR THE AUTHORS OR THEIR AFFILIATED
INSTITUTIONS SHALL BE LIABLE FOR ANY INCIDENTAL OR CONSEQUENTIAL DAMAGES
IN CONNECTION WITH THE USE OF THE SOFTWARE, INCLUDING WITHOUT LIMITATION
LOST PROFITS OR INJURY TO BUSINESS, WHETHER OR NOT ARGONNE NATIONAL
LABORATORY, AND THE AUTHORS AND THEIR AFFILIATED INSTITUTIONS KNOW OR
HAVE REASON TO KNOW OF THE POSSIBILITY OF SUCH DAMAGES.

INDEMNITY.  Each person that uses the software thereby agrees, to
indemnify and defend Argonne National Laboratory and the authors and
their affiliated institutions, or any of them, against any loss, expense,
claim, damage, or liability of any kind arising from or connected with
their respective uses of the software, and to hold them or any of them
harmless from any of the same, WHETHER OR NOT ARISING IN WHOLE OR IN PART
FROM THE NEGLIGENCE OR GROSS NEGLIGENCE OF ARGONNE NATIONAL LABORATORY OR
ANY OF THE AUTHORS OR THEIR AFFILIATED INSTITUTIONS.

SUPPORT. Each person that uses this software understands that the software
is not supported by the authors or by their affiliated institutions.
*/

/****************************************************************************/
/*                                                         THIS FILE IS C++ */
#ifdef __cplusplus


/****************************************************************************/
/*                                                                 INCLUDES */
#include "dvlparms.h" /* Developers Parameters */
#include "usrparms.h" /* Users Parameters */

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

/* NOTICE: There are automatic includes at the end of this file! */


/****************************************************************************/
/*                                                     FORWARD DECLARATIONS */

/*--------------------------------------------------------------------------*/
class adouble;
class adub;
class badouble;
class badoublev;
class adoublev;
class adubv;
class along;
/* class doublev;  that's history */

/*--------------------------------------------------------------------------*/
void condassign( double &res, const double &cond, 
                 const double &arg1, const double &arg2 );
void condassign( double &res, const double &cond,
                 const double &arg );

double fmin( const double &x, const double &y );
double fmax( const double &x, const double &y );


/****************************************************************************/
/*                                                           CLASS BADOUBLE */

/* 
   The class badouble contains the basic definitions for 
   the arithmetic operations, comparisons, etc. 
   This is a basic class from which the adub and adouble are 
   derived.  Notice that the constructors/destructors for 
   the class badouble are of the trivial variety.  This is the
   main difference among badoubles, adubs, and adoubles.
*/
class badouble{
  friend class badoublev;
 protected:
  locint location;
  badouble( void ){};
  badouble( const badouble& a ) {location = a.location;};
  badouble( locint lo ) {location = lo;};

 public:
/*--------------------------------------------------------------------------*/
  inline locint loc( void ) const;                         /* Helpful stuff */
  inline friend double value( const badouble& );

  /*------------------------------------------------------------------------*/
  badouble& operator >>= ( double& );                        /* Assignments */
  badouble& operator <<= ( double );
  badouble& operator = ( double );
  badouble& operator = ( const badouble& );
  badouble& operator = ( const adub& );
  /* badouble& operator = ( const adouble& ); 
     !!! olvo 991210: was the same as badouble-assignment */

/*--------------------------------------------------------------------------*/
  friend ostream& operator << ( ostream&, const badouble& );  /* IO friends */
  friend istream& operator >> ( istream&, const badouble& );

  /*------------------------------------------------------------------------*/
  badouble& operator += ( double );               /* Operation + Assignment */
  badouble& operator += ( const badouble& );
  badouble& operator -= ( double y );
  badouble& operator -= ( const badouble& );
  badouble& operator *= ( double );
  badouble& operator *= ( const badouble& );
  badouble& operator /= ( double );
  badouble& operator /= ( const badouble& );
  /* olvo 991122 n2l: new special op_codes */
  badouble& operator += ( const adub& );
  badouble& operator -= ( const adub& );

/*--------------------------------------------------------------------------*/
                                                    /* Comparison (friends) */
  inline friend int operator != ( const badouble&, const badouble& );
  inline friend int operator != ( double, const badouble& );
         friend int operator != ( const badouble&, double );
  inline friend int operator == ( const badouble&, const badouble& );
  inline friend int operator == ( double, const badouble& );
         friend int operator == ( const badouble&, double );
  inline friend int operator <= ( const badouble&, const badouble& );
  inline friend int operator <= ( double, const badouble& );
         friend int operator <= ( const badouble&, double );
  inline friend int operator >= ( const badouble&, const badouble& );
  inline friend int operator >= ( double, const badouble& );
         friend int operator >= ( const badouble&, double );
  inline friend int operator >  ( const badouble&, const badouble& );
  inline friend int operator >  ( double, const badouble& );
         friend int operator >  ( const badouble&, double );
  inline friend int operator <  ( const badouble&, const badouble& );
  inline friend int operator <  ( double, const badouble& );
         friend int operator <  ( const badouble&, double );


/*--------------------------------------------------------------------------*/
                                                /* sign operators (friends) */
         friend adub operator + ( const badouble& x ); 
         friend adub operator - ( const badouble& x ); 

/*--------------------------------------------------------------------------*/
                                              /* binary operators (friends) */
         friend adub operator + ( const badouble&, const badouble& ); 
         friend adub operator + ( double, const badouble& ); 
         friend adub operator + ( const badouble&, double ); 
         friend adub operator - ( const badouble&, const badouble& ); 
  inline friend adub operator - ( const badouble&, double ); 
         friend adub operator - ( double, const badouble& ); 
         friend adub operator * ( const badouble&, const badouble& ); 
         friend adub operator * ( double, const badouble& ); 
  inline friend adub operator * ( const badouble&, double );
  inline friend adub operator / ( const badouble&, double );
         friend adub operator / ( const badouble&, const badouble& ); 
         friend adub operator / ( double, const badouble& ); 

/*--------------------------------------------------------------------------*/
                                               /* unary operators (friends) */
         friend adub exp  ( const badouble& ); 
         friend adub log  ( const badouble& ); 
         friend adub sqrt ( const badouble& );
         friend adub sin  ( const badouble& ); 
         friend adub cos  ( const badouble& );
         friend adub tan  ( const badouble& );
         friend adub asin ( const badouble& );
         friend adub acos ( const badouble& );
         friend adub atan ( const badouble& ); 

/*--------------------------------------------------------------------------*/
                                             /* special operators (friends) */
         friend adouble atan2 ( const badouble&, const badouble& ); 
         /* no internal use of condassign: */
         friend adub    pow   ( const badouble&, double );
         /* uses condassign internally */
         friend adouble pow   ( const badouble&, const badouble& );
         friend adouble pow   ( double, const badouble& );
         friend adub    log10 ( const badouble& );
         /* User defined version of logarithm to test extend_quad macro */
         friend adouble myquad( const badouble& );

/*--------------------------------------------------------------------------*/
        /* Additional ANSI C standard Math functions Added by DWJ on 8/6/90 */
         friend adub sinh  ( const badouble& );
         friend adub cosh  ( const badouble& );
         friend adub tanh  ( const badouble& );
         friend adub asinh ( const badouble& );
         friend adub acosh ( const badouble& );
         friend adub atanh ( const badouble& );

         friend adub fabs  ( const badouble& );
         friend adub ceil  ( const badouble& );
         friend adub floor ( const badouble& );

         friend adub fmax ( const badouble&, const badouble& );
         friend adub fmax ( double, const badouble& );
         friend adub fmax ( const badouble&, double );
         friend adub fmin ( const badouble&, const badouble& );
         friend adub fmin ( double, const badouble& );
         friend adub fmin ( const badouble&, double );

         friend adub ldexp ( const badouble&, int );
         friend adub frexp ( const badouble&, int* );
         friend adub erf   ( const badouble& );

/*--------------------------------------------------------------------------*/
                                                            /* Conditionals */
         friend void condassign( adouble &res, const adouble &cond,
                          const adouble &arg1, const adouble &arg2 );
         friend void condassign( adouble &res, const adouble &cond,
                                               const adouble &arg );
         friend void condassign( along &res, const adouble &cond,
                        const adouble &arg1, const adouble &arg2 );
         friend void condassign( along &res, const adouble &cond,
                                             const adouble &arg );
};



/****************************************************************************/
/*                                                               CLASS ADUB */

/* 
   The class Adub
   ---- Basically used as a temporary result.  The address for an
        adub is usually generated within an operation.  That address
        is "freed" when the adub goes out of scope (at destruction time).
   ---- operates just like a badouble, but it has a destructor defined for it.
*/

class adub:public badouble{
  friend class adouble;
/* added Sep/01/96 */
  friend class asub;
  friend class along;
 protected:
  adub( locint lo ):badouble(lo){};
  adub( void ):badouble(0)
  { fprintf(DIAG_OUT,"ADOL-C error: illegal default construction of adub"
                     " variable\n");
    exit(-2);
  };
  adub( double ):badouble(0)
  { fprintf(DIAG_OUT,"ADOL-C error: illegal  construction of adub variable"
                     " from double\n");
    exit(-2);
  };

 public:

/*--------------------------------------------------------------------------*/
                                                /* sign operators (friends) */
         friend adub operator + ( const badouble& x ); 
         friend adub operator - ( const badouble& x ); 

/*--------------------------------------------------------------------------*/
                                              /* binary operators (friends) */
         friend adub operator + ( const badouble&, const badouble& ); 
         friend adub operator + ( double, const badouble& ); 
         friend adub operator + ( const badouble&, double ); 
         friend adub operator - ( const badouble&, const badouble& ); 
  inline friend adub operator - ( const badouble&, double ); 
         friend adub operator - ( double, const badouble& ); 
         friend adub operator * ( const badouble&, const badouble& ); 
         friend adub operator * ( double, const badouble& ); 
  inline friend adub operator * ( const badouble&, double );
  inline friend adub operator / ( const badouble&, double );
         friend adub operator / ( const badouble&, const badouble& ); 
         friend adub operator / ( double, const badouble& ); 

/*--------------------------------------------------------------------------*/
                                               /* unary operators (friends) */
         friend adub exp  ( const badouble& ); 
         friend adub log  ( const badouble& ); 
         friend adub sqrt ( const badouble& );
         friend adub sin  ( const badouble& ); 
         friend adub cos  ( const badouble& );
         friend adub tan  ( const badouble& );
         friend adub asin ( const badouble& );
         friend adub acos ( const badouble& );
         friend adub atan ( const badouble& ); 

/*--------------------------------------------------------------------------*/
                                             /* special operators (friends) */
         /* no internal use of condassign: */
         friend adub    pow   ( const badouble&, double );
         friend adub    log10 ( const badouble& );

/*--------------------------------------------------------------------------*/
        /* Additional ANSI C standard Math functions Added by DWJ on 8/6/90 */
         friend adub sinh  ( const badouble& );
         friend adub cosh  ( const badouble& );
         friend adub tanh  ( const badouble& );
         friend adub asinh ( const badouble& );
         friend adub acosh ( const badouble& );
         friend adub atanh ( const badouble& );

         friend adub fabs  ( const badouble& );
         friend adub ceil  ( const badouble& );
         friend adub floor ( const badouble& );

         friend adub fmax ( const badouble&, const badouble& );
         friend adub fmax ( double, const badouble& );
         friend adub fmax ( const badouble&, double );
         friend adub fmin ( const badouble&, const badouble& );
         friend adub fmin ( double, const badouble& );
         friend adub fmin ( const badouble&, double );

         friend adub ldexp ( const badouble&, int );
         friend adub frexp ( const badouble&, int* );
         friend adub erf   ( const badouble& );

/*--------------------------------------------------------------------------*/
	                                     /* vector operations (friends) */
         friend adub operator*( const badoublev&, double* );
         friend adub operator*( double*, const badoublev& );
         friend adub operator*( const badoublev&, const badoublev& );

#ifdef overwrite
         ~adub();
#endif
};


/****************************************************************************/
/*                                                            CLASS ADOUBLE */
/*
  The class adouble.
  ---Derived from badouble.  Contains the standard constructors/destructors.
  ---At construction, it is given a new address, and at destruction, that
     address is freed.
*/
class adouble:public badouble{
  friend class along;
 public:
  adouble( const adub& );
  adouble( const along& );
  adouble( const adouble& );
  adouble( void );
  adouble( double );
  /* adub prevents postfix operators to occur on the left 
     side of an assignment which would not work  */
  adub operator++( int );
  adub operator--( int );
  badouble& operator++( void );
  badouble& operator--( void );

#ifdef overwrite
  ~adouble();
#endif

  adouble& operator = ( double );
  adouble& operator = ( const badouble& );
  /* adouble& operator = ( const adouble& );
     !!! olvo 991210 was the same as badouble-assignment */ 
  adouble& operator = ( const adub& );
};


/****************************************************************************/
/*                                                               CLASS ASUB */
class asub:public badouble
{ locint base,
         offset;
 public:
  asub( locint start, locint index );
#ifdef overwrite
  ~asub();
#endif
  asub& operator <<= ( double );
  asub& operator =   ( double );
  /* asub& operator =   ( const adub& );
     !!! olvo 991210 is the same as normal assignment */
  asub& operator =   ( const badouble& );
  /* added Sep/01/96 */
  /* olvo 991210 seems to be a provisional version */
  asub& operator =   ( const asub& );
  asub& operator +=  ( double );
  asub& operator +=  ( const badouble& );
  asub& operator -=  ( double x );
  asub& operator -=  ( const badouble& );
  asub& operator *=  ( double x );
  asub& operator *=  ( const badouble& );
  asub& operator /=  ( double x );
  asub& operator /=  ( const badouble& );
  /* adub prevents postfix operators to occur on the left 
     side of an assignment which would not work  */
  adub operator++( int );
  adub operator--( int );
  asub& operator++( void );
  asub& operator--( void );
};


/****************************************************************************/
/*                                                              CLASS ALONG */
/* The class along was originally designed for the sole purpose of 
   allowing subscripting operations with computable ON THE TAPE.
   The current definition refers to badouble, i.e. it is a 
   floating point class. The current constuction allows to 
   use alongs in the same way as adoubles. Especially adubs can 
   can be used as subscripts, i.e. the results of "active computations".
   This useful because people want to compute the indices on the tape.
   Notice, that along does NOT perform integer arithmetic. 
   This is the major disadvantage of the current along. */
class along:public badouble
{ friend class adouble;
 public:
  along( const adub& );
  along( const along& );
  along( void );
  along( int );
#ifdef overwrite
  ~along();
#endif
  along& operator = ( int );
  /* along& operator = ( const adouble& );
     !!! olvo 991210 is the same as badouble-assignment */ 
  along& operator = ( const badouble& );
  along& operator = ( const along& );
  along& operator = ( const adub& );
  /* adub prevents postfix operators to occur on the left
     side of an assignment which would not work  */
  adub operator++( int );
  adub operator--( int );
  along&  operator++( void );
  along&  operator--( void );
};


/****************************************************************************/
/*                                                       INLINE DEFINITIONS */

/*--------------------------------------------------------------------------*/
inline locint badouble::loc( void ) const
{ return location;
} 

/*--------------------------------------------------------------------------*/
inline double value(const badouble& x) 
{ extern double* store;
  return store[x.location];
}

/*--------------------------------------------------------------------------*/
                                                              /* Comparison */
inline int operator != ( const badouble& u, const badouble& v )
{ return (u-v != 0);
}

inline int operator != ( double coval, const badouble& v)
{ if (coval)
    return (-coval+v != 0);
  else
    return (v != 0);
}

inline int operator == ( const badouble& u, const badouble& v )
{ return (u-v == 0);
}

inline int operator == ( double coval, const badouble& v)
{ if (coval)
    return (-coval+v == 0);
  else
    return (v == 0);
}

inline int operator <= ( const badouble& u, const badouble& v )
{ return (u-v <= 0);
}

inline int operator <= ( double coval, const badouble& v )
{ if (coval)
    return (-coval+v >= 0);
  else
    return (v >= 0);
}

inline int operator >= ( const badouble& u, const badouble& v )
{ return (u-v >= 0);
}

inline int operator >= ( double coval, const badouble& v )
{ if (coval)
    return (-coval+v <= 0);
  else
    return (v <= 0);
}

inline int operator > ( const badouble& u, const badouble& v )
{ return (u-v > 0);
}

inline int operator > ( double coval, const badouble& v )
{ if (coval)
    return (-coval+v < 0);
  else
    return (v < 0);
}

inline int operator < ( const badouble& u, const badouble& v )
{ return (u-v < 0);
}

inline int operator < ( double coval, const badouble& v )
{ if (coval)
    return (-coval+v > 0);
  else
    return (v > 0);
}

/*--------------------------------------------------------------------------*/
/* Subtract a floating point from an adouble  */
inline adub operator - ( const badouble& x , double coval )
{ return (-coval) + x;
}

/*--------------------------------------------------------------------------*/
/* Multiply an adouble by a floating point */ 
inline adub operator * (const badouble& x, double coval)
{ return coval * x;
}

/*--------------------------------------------------------------------------*/
/* Divide an adouble by a floating point */
inline adub operator / (const badouble& x, double coval)
{ return (1.0/coval) * x;
}


/****************************************************************************/
/*                                                         STOCK OPERATIONS */
void   take_stock (void);
locint keep_stock (void);


/****************************************************************************/
/*                                                        AUTOMTIC INCLUDES */
#include "avector.h" /* active vector classes */
#include "taputil.h" /* trace on/off */


/****************************************************************************/
/*                                                                THAT'S ALL*/
#endif
#endif






