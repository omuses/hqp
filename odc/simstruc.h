/*
 * simstruc.h:
 *   Provide approprirate SimStruct for compiling an S-function
 *   either for Simulink(R) or for HQP.
 *
 * Intended usage:
 *  This file should be included at the beginning of an S-function with
 *  #include "simstruc.h"
 *  As it is placed in the same directory as the S-function C file, it 
 *  will be found by a C preprocessor instead of the orignal simstruc.h
 *  coming with Simulink(R).
 *  The S-function can either be compiled for Simulink(R) or for HQP.
 *  * Define the macro MATLAB_MEX_FILE (e.g. when invoking the Matlab(R) 
 *  command mex) if the S-function should be compiled for Simulink(R).
 *  Then the original simstruc.h found in the include path will be 
 *  included from here via
 *  #include <simstruc.h>
 *  Note: Take care to not having -I. in your include path!
 *  * The alternative file Hxi_SimStruct.h will be included to compile
 *  the S-function for HQP.
 *
 * (Matlab and Simulink are registered trademarks of The MathWorks, Inc.)
 *
 * rf, 05/05/2001
 */

#if !defined(HXI_SIMSTRUC_H)
#define HXI_SIMSTRUC_H

#if defined(MATLAB_MEX_FILE)
// The S-function is being compiled as Matlab Mex file;
// we include the simstruc.h found in the include path.

#include <simstruc.h>

#else
// We are compiling the S-function for HQP
// using ADOL-C for automatic differentiation.

#include <adouble.h>
#define HXI_REAL_T adouble
#include <Hxi_SimStruct.h>

#endif // defined(MATLAB_MEX_FILE)

#endif // !defined(HXI_SIMSTRUC_H)
