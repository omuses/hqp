/*
 * simstruc.h: provide approprirate SimStruct for an S-function
 *
 * rf, 05/05/2001
 */

#if !defined(HXI_SIMSTRUC_H)
#define HXI_SIMSTRUC_H

#if defined(MATLAB_MEX_FILE)
// this S-function is being compiled as Matlab Mex file;
// we include the simstruc.h provided by Matlab

#include <simstruc.h>

#else
// we are compiling this S-function for HQP using ADOL-C

#include <adouble.h>

#define HXI_REAL_T adouble
#include <Hxi_SimStruct.h>

#endif // defined(MATLAB_MEX_FILE)

#endif // !defined(HXI_SIMSTRUC_H)
