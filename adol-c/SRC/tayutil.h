#ifndef _TAYUTIL_H_
#define _TAYUTIL_H_
/*
   --------------------------------------------------------------
   File tayutil.h of ADOL-C version 1.8.5         as of Nov/22/99
   --------------------------------------------------------------
   tayutil.h defines the prototypes for the functions from 
   tayutil.c.  See tayutil.c for an explanation of the 
   functionality of these routines.

   Last changes: 
     991122 olvo  new op_codes eq_plus_prod eq_min_prod
                  for  y += x1 * x2
                  and  y -= x1 * x2
                  --> new: delete_scaylor(..)  
     981130 olvo: automatic cleanup from utils.C moved here
     980921 olvo: new interface of void overwrite_scaylor(..) to
                  allow correction of old overwrite in store
     980708 olvo  new:  void overwrite_scaylor(..)

   --------------------------------------------------------------
*/

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                        NO PUBLIC EXPORTS */


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                  ADOL-C INTERNAL EXPORTS */
#ifdef _ADOLC_SRC_


#ifdef __cplusplus
/****************************************************************************/
/****************************************************************************/
/*                                                       Now the C++ THINGS */



/****************************************************************************/
/****************************************************************************/
/*                                                         Now the C THINGS */
extern "C" {
#endif


/****************************************************************************/
#ifndef __STDC__
   int unlink(char *); 
#endif


/****************************************************************************/
/*                                                            CONTROL STUFF */
int taylor_access();
void close_taylor();
void taylor_begin( int, double**,int );
void taylor_close( int, int, int );
void taylor_back ( revreal*, int*, int*, int* );
void taylor_back2( revreal**, int*, int*, int* );


/****************************************************************************/
/*                                                                   WRITEs */
void write_taylor( locint, int );
void write_scaylor( revreal );
/* olvo 980708 new nl */
void overwrite_scaylor( revreal, revreal* );
/* olvo 991122 new nl */
void delete_scaylor( revreal* );
void write_scaylors( double*, int );


/****************************************************************************/
/*                                                                     GETs */
void get_taylor( locint );
void get_taylors( locint, int );

/****************************************************************************/
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

#endif /* _ADOLC_SRC_ */

#endif

