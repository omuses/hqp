/**
 * @file Hxi_MEX_SFunction.h
 *   Interface to a Simulink(R) S-function given as binary MEX object.
 *   Several S-function methods are supported that redirect the call
 *   to the MEX object. Data is communicated through a SimStruct.
 *   The two functions Hxi_SimStruct_create and Hxi_SimStruct_destroy
 *   are provided for allocating and releasing a SimStruct.
 *
 * (Simulink is a registered trademark of The MathWorks, Inc.)
 *
 * rf, 07/14/2001
 */

/*
    Copyright (C) 1994--2005  Ruediger Franke

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

/** Avoid multiple inclusion */
#if !defined(Hxi_MEX_SFunction_H)
#define Hxi_MEX_SFunction_H

/** define HXI_MEX_S_FUNCTION prior to inclusion of simstruc.h,
    in order to select the MEX version of SimStruct. */
#define HXI_MEX_S_FUNCTION
#include "simstruc.h"

/** Create a MEX SimStruct */
SimStruct *Hxi_MEX_SimStruct_create();

/** Free memory of MEX SimStruct */
void Hxi_MEX_SimStruct_destroy(SimStruct *S);

/** Initialize MEX S-function, including call to mdlInitializeSizes */
void Hxi_MEX_SFunction_init(SimStruct *S);

#endif
