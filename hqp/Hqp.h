/**
 * @file Hqp.h
 *   declaration of external API for HQP (Huge Quadratic Programming)
 *
 * rf, 5/28/94
 */

/*
    Copyright (C) 1994--2014  Ruediger Franke

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; 
    version 2 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library (file COPYING.LIB);
    if not, write to the Free Software Foundation, Inc.,
    59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef Hqp_H
#define Hqp_H

/* include command interface */
#include <If.h>

/* include type declarations */
#include <Meschach.h>

/** define HQP_API when compiling a Dynamic Link Library (DLL) */
#ifndef HQP_API
#define HQP_API 
#endif

/** initialize Tcl module Hqp */
extern "C" HQP_API int Hqp_Init(Tcl_Interp *interp);

/** call Tcl procedure hqp_exit for signals SIGINT, SIGFPE and SIGXCPU */
extern "C" HQP_API int Hqp_InitSignalHandler();

#endif
