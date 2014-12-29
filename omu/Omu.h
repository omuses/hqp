/**
 * @file Omu.h
 *   declaration of external API for Omuses
 *
 * rf, 12/29/14
 */

/*
    Copyright (C) 1997--2014  Ruediger Franke

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

#ifndef Omu_H
#define Omu_H

/** define OMU_API when compiling a Dynamic Link Library (DLL) */
#ifndef OMU_API
#define OMU_API
#endif

/* include underlying HQP */
#include <Hqp.h>

/** initialize Tcl module Omuses */
extern "C" OMU_API int Omu_Init(Tcl_Interp *interp);

#endif
