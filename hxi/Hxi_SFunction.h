/**
 * @file Hxi_SFunction.h
 *   Interface to a Simulink(R) S-function given as binary object.
 *   Several S-function methods are supported that redirect the call
 *   to the binary object. Data is communicated through a SimStruct.
 *   The two functions Hxi_SimStruct_create and Hxi_SimStruct_destroy
 *   are provided for allocating and releasing the SimStruct.
 *
 * (Simulink is a registered trademark of The MathWorks, Inc.)
 *
 * rf, 01/07/2005
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
#if !defined(Hxi_SFunction_H)
#define Hxi_SFunction_H

#if defined(MATLAB_MEX_FILE)
#error "Hxi_SFunction can't be used for MATLAB_MEX_FILE"
#endif

#include <simstruc.h>

/** Create SimStruct for S-function. */
SimStruct *Hxi_SimStruct_create();

/** Realease S-function and SimStruct. */
void Hxi_SimStruct_destroy(SimStruct *S);

/** Load binary object and initialize sizes of data vectors in S. 
    This function loads the shared object found under ssGetPath(S),
    calls the S-function method mdlInitializeSizes
    and allocates memory required for a level 2 S-function.
    The model parameters must have been initialized prior to calling this
    function. */
void mdlInitializeSizes(SimStruct *S);

/// @name Macros for calling S-function methods
//@{
#define mdlCheckParameters(S) 	    (*((S)->getmdlCheckParameters()))(S)
#define mdlInitializeSampleTimes(S) (*((S)->getmdlInitializeSampleTimes()))(S)
#define mdlStart(S) 		    (*((S)->getmdlStart()))(S)
#define mdlInitializeConditions(S)  (*((S)->getmdlInitializeConditions()))(S)
#define mdlOutputs(S, tid) 	    (*((S)->getmdlOutputs()))(S, tid)
#define mdlUpdate(S, tid) 	    (*((S)->getmdlUpdate()))(S, tid)
#define mdlDerivatives(S) 	    (*((S)->getmdlDerivatives()))(S)
#define mdlJacobian(S) 		    (*((S)->getmdlJacobian()))(S)
#define mdlTerminate(S) 	    (*((S)->getmdlTerminate()))(S)
//@}


#endif // !defined(Hxi_SFunction_H)
