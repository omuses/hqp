/**
 * @file Hxi_SFunction.h
 *   Interface to a Simulink(R) S-function given as binary object.
 *   Several S-function methods are supported that redirect the call
 *   to the binary object. Data is communicated through a SimStruct.
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

#include "Hxi_SimStruct.h"

/** @name Supported S-function methods. */
//@{
/** Initialize sizes of data vectors in %SimStruct. 
    This function works with the shared object found under ssGetPath(S),
    calls the S-function method mdlInitializeSizes
    and allocates memory required for a level 2 S-function.
    The model parameters must have been initialized prior to calling this
    function. */
inline void mdlInitializeSizes(SimStruct *S) {
  Hxi_mdlInitializeSizes(S);
}

/** Optional: Allocate local ressources for simulation. */
inline void mdlStart(SimStruct *S) {
  Hxi_mdlStart(S);
}

/** Initialize sample times. */
inline void mdlInitializeSampleTimes(SimStruct *S) {
  Hxi_mdlInitializeSampleTimes(S);
}

/** Optional: Compute initial conditions. */
inline void mdlInitializeConditions(SimStruct *S) {
  Hxi_mdlInitializeConditions(S);
}

/** Compute model outputs. */
inline void mdlOutputs(SimStruct *S, int_T tid) {
  Hxi_mdlOutputs(S, tid);
}

/** Optional: Update discrete-time states. */
inline void mdlUpdate(SimStruct *S, int_T tid) {
  Hxi_mdlUpdate(S, tid);
}

/** Optional: Compute derivatives for continuous-time states. */
inline void mdlDerivatives(SimStruct *S) {
  Hxi_mdlDerivatives(S);
}

/** Optional: Compute Jacobian J = d(dxc,xd,y)/d(xc,xd,u). */
inline void mdlJacobian(SimStruct *S) {
  Hxi_mdlJacobian(S);
}

/** Release resources allocated for simulation. */
inline void mdlTerminate(SimStruct *S) {
  Hxi_mdlTerminate(S);
}
//@}


#endif // !defined(Hxi_SFunction_H)
