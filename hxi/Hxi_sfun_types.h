/** 
 * @file Hxi_sfun_types.h
 *   Native type definitions required to compile
 *   a Simulink(R) S-function with Hqp.
 *
 * (Simulink is a registered trademark of The MathWorks, Inc.)
 *
 * rf, 05/06/2001
 */

/** Avoid multiple inclusion */
#if !defined(Hxi_sfun_types_H)
#define Hxi_sfun_types_H

/** Use this macro to avoid compiler warning about unused function args. */
#if !defined(UNUSED_ARG)
#define UNUSED_ARG(arg) 		(arg)=(arg)
#endif

/** HXI_REAL_T can be defined before including this file (default: double). */
#if !defined(HXI_REAL_T)
#define HXI_REAL_T double
#endif

/** Real type used in S-function. */
typedef HXI_REAL_T real_T;

/** Pointer type used for Real inputs of S-function. */
typedef real_T **InputRealPtrsType;

/** Integer type used in S-function. */
typedef int int_T;

/** Unsigned integer type used in S-function. */
typedef unsigned uint_T;

/** Character type used in S-function. */
typedef char char_T;

/** mxArray element types */
typedef enum mxComplexity {mxREAL=0, mxComplex} mxComplexity;

#if defined(__cplusplus)
/** ADOL-C's value function for double. It may be needed to compile
    code written for real_T adouble with real_T double. */
inline double value(double a) {return a;}
#endif


#endif
