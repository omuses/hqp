/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     adouble.h
 Revision: $Id: adouble.h,v 1.1 2004/10/13 14:18:11 e_arnold Exp $
 Contents: adouble.h contains the basis for the class of adouble
           included here are all the possible functions defined on
           the adouble class.  Notice that, as opposed to ealier versions,
           both the class adub and the class adouble are derived from a base
           class (badouble).  See below for further explanation.

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
          20000107 olvo:   iostream.h instaed of stream.h  
          19991210 olvo:   checking the changes
          19991122 olvo:   new op_codes eq_plus_prod eq_min_prod
                           for  y += x1 * x2
                           and  y -= x1 * x2  
          19981201 olvo:   last check: 
                           - taputil things changed, includes 
          19980820 olvo:   new comparison strategy & some inlines
          19980709 olvo:   modified sign operators

----------------------------------------------------------------------------*/

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

#if !defined(ADOLC_ADOUBLE_H)
#define ADOLC_ADOUBLE_H 1

/****************************************************************************/
/*                                                         THIS FILE IS C++ */
#ifdef __cplusplus

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "common.h"
#include "taputil.h"

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
  friend std::ostream& operator << ( std::ostream&, const badouble& );  /* IO friends */
  friend std::istream& operator >> ( std::istream&, const badouble& );

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

/****************************************************************************/
/*                                                                THAT'S ALL*/

#endif

#endif
