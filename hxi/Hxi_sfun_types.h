/** 
 * @file Hxi_sfun_types.h
 *   Alternative native type definitions required to compile
 *   a Simulink(R) S-function with HQP.
 *
 * (Simulink is a registered trademark of The MathWorks, Inc.)
 *
 * rf, 05/06/2001
 */

/*
    Copyright (C) 1994--2001  Ruediger Franke

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

#if !defined(Hxi_sfun_types_H)
#define Hxi_sfun_types_H

/** Use this macro to avoid compiler warning about unused function args. */
#define UNUSED_ARG(arg) 		(arg)=(arg)

/** HXI_REAL_T can be defined before including this file (default: adouble). */
#if !defined(HXI_REAL_T)
#include <adouble.h>
#define HXI_REAL_T adouble
#endif

/** Real type used in S-function. */
typedef HXI_REAL_T real_T;
typedef real_T **InputRealPtrsType;

/** Integer type used in S-function. */
typedef int int_T;

/** Unsigned integer type used in S-function. */
typedef unsigned uint_T;

/** Character type used in S-function. */
typedef char char_T;

#endif

