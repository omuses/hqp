/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     tayutil_p.h
 Revision: $Id: tayutil_p.h,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: Management for the value stack tape (Taylors)
           (ADOL-C internal routines)

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
         20030304 andrea: identify value stack by tag
         20010719 andrea: add write_taylors(..)
                          add get_taylors_p(..)
         19991122 olvo:   new op_codes eq_plus_prod eq_min_prod
                          for  y += x1 * x2
                          and  y -= x1 * x2
                          --> new: delete_scaylor(..)  
         19981130 olvo:   automatic cleanup from utils.C moved here
         19980921 olvo:   new interface of void overwrite_scaylor(..) to
                          allow correction of old overwrite in store
         19980708 olvo:   new:  void overwrite_scaylor(..)
         
----------------------------------------------------------------------------*/

#if !defined(ADOLC_TAYUTIL_P_H)
#define ADOLC_TAYUTIL_P_H 1

#include "common.h"

BEGIN_C_DECLS

/****************************************************************************/
/*                                                         GLOBAL VARIABLES */

/*--------------------------------------------------------------------------*/
/* File Name */
extern char vs_file_name[20];

/****************************************************************************/
#if !defined(__STDC__)
   int unlink(char *); 
#endif

/****************************************************************************/
/*                                                            CONTROL STUFF */
int taylor_access();
void close_taylor();
void taylor_begin( short, int, double**,int );
void taylor_close( int, int, int );
void taylor_back ( int, revreal*, int*, int*, int* );
void taylor_back2( int, revreal**, int*, int*, int* );

/****************************************************************************/
/*                                                                   WRITEs */
void write_taylor( locint, int );
void write_taylors( locint, int, int, int );
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
void get_taylors_p( locint, int, int );

END_C_DECLS

/****************************************************************************/
#endif
