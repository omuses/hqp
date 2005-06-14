/**
 * @file Hxi_SimStruct.h
 *    Declaration of Simulink(R) S-function methods supported by Hqp.
 *
 * (Simulink is a registered trademark of The MathWorks, Inc.)
 *
 * rf, 05/05/2001
 */

/** Avoid multiple inclusion */
#if !defined(Hxi_SimStruct_H)
#define Hxi_SimStruct_H

#include <stdio.h>
#include <stdlib.h>

#include "Hxi_sfun_types.h"

#if defined(S_FUNCTION_NAME)
/** @name Definitions expected in S-functions and not used otherwise
    (only a subset is supported) */
/*@{*/
#define SS_OPTION_EXCEPTION_FREE_CODE 		0x0001
#define SS_OPTION_RUNTIME_EXCEPTION_FREE_CODE 	0x0002
#define SS_OPTION_PORT_SAMPLE_TIMES_ASSIGNED 	0x0004
#define SS_REUSABLE_AND_LOCAL 			0x0100
#define SS_NOT_REUSABLE_AND_GLOBAL 		0x0200

#define SS_STDIO_AVAILABLE 			true

/** Print anything */
#define ssPrintf 	printf

/** Print a warning */
#define ssWarning(S, msg) \
  printf("Warning: S-function \"%s\": %s\n", ssGetPath(S), msg)
/*@}*/
#endif /* defined(S_FUNCTION_NAME) */

/** @name Definitions used outside S-functions as well */
/*@{*/
#define CONTINUOUS_SAMPLE_TIME 			0.0

/** Free memory allocated though mx methods, e.g. mxArrayToString. */
#define mxFree(p)	free(p)
/*@}*/

/* forward declare SimStruct and mxArray */
typedef struct SimStruct SimStruct;
typedef struct mxArray mxArray;

#if defined(_MSC_VER)
  /** define inline for Microsoft C compiler */
# define inline __forceinline
#endif

#if defined(HXI_INLINE_S_FUNCTION)
#  define HXI_EXTERN static
#else
#  if defined(__cplusplus)
#    define HXI_EXTERN extern "C"
#  else
#    define HXI_EXTERN extern
#  endif
#endif

/** S-function method type. */
typedef void (SFunctionMethod1_type)(SimStruct *S);
/** S-function method type with additional task id argument. */
typedef void (SFunctionMethod2_type)(SimStruct *S, int_T tid);

/** Construct a SimStruct for given path. Store path in SimStruct. */
HXI_EXTERN SimStruct *Hxi_SimStruct_create(const char *path);

/** Delete a SimStruct. */
HXI_EXTERN void Hxi_SimStruct_destroy(SimStruct *S);

#if !defined(HXI_INLINE_S_FUNCTION)
/** @name Supported S-function methods. */
/*@{*/
/** Initialize sizes of data vectors in %SimStruct. 
    This function works with the shared object found under ssGetPath(S),
    calls the S-function method mdlInitializeSizes
    and allocates memory required for a level 2 S-function.
    The model parameters must have been initialized prior to calling this
    function. */
HXI_EXTERN void Hxi_mdlInitializeSizes(SimStruct *S);

/** Optional: Allocate local ressources for simulation. */
HXI_EXTERN void Hxi_mdlStart(SimStruct *S);

/** Initialize sample times. */
HXI_EXTERN void Hxi_mdlInitializeSampleTimes(SimStruct *S);

/** Optional: Compute initial conditions. */
HXI_EXTERN void Hxi_mdlInitializeConditions(SimStruct *S);

/** Compute model outputs. */
HXI_EXTERN void Hxi_mdlOutputs(SimStruct *S, int_T tid);

/** Optional: Update discrete-time states. */
HXI_EXTERN void Hxi_mdlUpdate(SimStruct *S, int_T tid);

/** Optional: Compute derivatives for continuous-time states. */
HXI_EXTERN void Hxi_mdlDerivatives(SimStruct *S);

/** Optional: Compute Jacobian J = d(dxc,xd,y)/d(xc,xd,u). */
HXI_EXTERN void Hxi_mdlJacobian(SimStruct *S);

/** Release resources allocated for simulation. */
HXI_EXTERN void Hxi_mdlTerminate(SimStruct *S);
/*@}*/
#endif //!defined(HXI_INLINE_S_FUNCTION)

/** @name Declaration of S-function methods. */ 
/*@{*/
#define HXI_SS_NOSET0(ITEM) \
  inline void ssSet##ITEM(SimStruct *S) { \
  }
#define HXI_SS_SETGET1(ITEM, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE ARG); \
  inline TYPE ssSet##ITEM(SimStruct *S, TYPE ARG) { \
    return hssSet##ITEM(S, ARG); \
  } \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S); \
  inline TYPE ssGet##ITEM(SimStruct *S) { \
    return hssGet##ITEM(S); \
  }
#define HXI_SS_SET1(ITEM, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE ARG); \
  inline TYPE ssSet##ITEM(SimStruct *S, TYPE ARG) { \
    return hssSet##ITEM(S, ARG); \
  }
#define HXI_SS_NOSET1(ITEM, TYPE, ARG) \
  inline TYPE ssSet##ITEM(SimStruct *S, TYPE ARG) { \
    return (TYPE)0; \
  }
#define HXI_SS_GET1(ITEM, TYPE) \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S); \
  inline TYPE ssGet##ITEM(SimStruct *S) { \
    return hssGet##ITEM(S); \
  }
#define HXI_SS_IS1(ITEM) \
  HXI_EXTERN int_T hssIs##ITEM(SimStruct *S); \
  inline int_T ssIs##ITEM(SimStruct *S) { \
    return hssIs##ITEM(S); \
  }

#define HXI_SS_SETGET2(ITEM, TYPE1, ARG1, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE1 ARG1, TYPE ARG); \
  inline TYPE ssSet##ITEM(SimStruct *S, TYPE1 ARG1, TYPE ARG) { \
    return hssSet##ITEM(S, ARG1, ARG); \
  } \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S, TYPE1 ARG1); \
  inline TYPE ssGet##ITEM(SimStruct *S, TYPE1 ARG1) { \
    return hssGet##ITEM(S, ARG1); \
  }
#define HXI_SS_NOSETGET2(ITEM, TYPE1, ARG1, TYPE, ARG) \
  inline TYPE ssSet##ITEM(SimStruct *S, TYPE1 ARG1, TYPE ARG) { \
    return (TYPE)0; \
  } \
  inline TYPE ssGet##ITEM(SimStruct *S, TYPE1 ARG1) { \
    return (TYPE)0; \
  }
#define HXI_SS_SET2(ITEM, TYPE1, ARG1, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE1 ARG1, TYPE ARG); \
  inline TYPE ssSet##ITEM(SimStruct *S, TYPE1 ARG1, TYPE ARG) { \
    return hssSet##ITEM(S, ARG1, ARG); \
  }
#define HXI_SS_NOSET2(ITEM, TYPE1, ARG1, TYPE, ARG) \
  inline TYPE ssSet##ITEM(SimStruct *S, TYPE1 ARG1, TYPE ARG) { \
    return (TYPE)0; \
  }
#define HXI_SS_GET2(ITEM, TYPE1, ARG1, TYPE) \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S, TYPE1 ARG1); \
  inline TYPE ssGet##ITEM(SimStruct *S, TYPE1 ARG1) { \
    return hssGet##ITEM(S, ARG1); \
  }
#define HXI_SS_IS2(ITEM, TYPE1, ARG1) \
  HXI_EXTERN int_T hssIs##ITEM(SimStruct *S, TYPE1 ARG1); \
  inline int_T ssIs##ITEM(SimStruct *S, TYPE1 ARG1) { \
    return hssIs##ITEM(S, ARG1); \
  }
/*@}*/

/** @name Declarations of mxArray methods. */ 
/*@{*/
#define HXI_MX_CREATE1(WHAT, TYPE1, ARG1) \
  HXI_EXTERN mxArray* hmxCreate##WHAT(TYPE1 ARG1); \
  inline mxArray* mxCreate##WHAT(TYPE1 ARG1) { \
    return hmxCreate##WHAT(ARG1); \
  }
#define HXI_MX_CREATE3(WHAT, TYPE1, ARG1, TYPE2, ARG2, TYPE3, ARG3) \
  HXI_EXTERN mxArray* hmxCreate##WHAT(TYPE1 ARG1, TYPE2 ARG2, TYPE3 ARG3); \
  inline mxArray* mxCreate##WHAT(TYPE1 ARG1, TYPE2 ARG2, TYPE3 ARG3) { \
    return hmxCreate##WHAT(ARG1, ARG2, ARG3); \
  }
#define HXI_MX_DESTROY(WHAT) \
  HXI_EXTERN void hmxDestroy##WHAT(mxArray *a); \
  inline void mxDestroy##WHAT(mxArray *a) { \
    hmxDestroy##WHAT(a); \
  }
#define HXI_MX_SETGET1(ITEM, TYPE, ARG) \
  HXI_EXTERN void hmxSet##ITEM(mxArray *a, TYPE ARG); \
  inline void mxSet##ITEM(mxArray *a, TYPE ARG) { \
    hmxSet##ITEM(a, ARG); \
  } \
  HXI_EXTERN TYPE hmxGet##ITEM(const mxArray *a); \
  inline TYPE mxGet##ITEM(const mxArray *a) { \
    return hmxGet##ITEM(a); \
  }
#define HXI_MX_GET1(ITEM, TYPE) \
  HXI_EXTERN TYPE hmxGet##ITEM(const mxArray *a); \
  inline TYPE mxGet##ITEM(const mxArray *a) { \
    return hmxGet##ITEM(a); \
  }
#define HXI_MX_IS1(ITEM) \
  HXI_EXTERN int_T hmxIs##ITEM(const mxArray *a); \
  inline int_T mxIs##ITEM(const mxArray *a) { \
    return hmxIs##ITEM(a); \
  }
#define HXI_MX_TO1(ITEM, TYPE) \
  HXI_EXTERN TYPE hmxArrayTo##ITEM(const mxArray *a); \
  inline TYPE mxArrayTo##ITEM(const mxArray *a) { \
    return hmxArrayTo##ITEM(a); \
  }
/*@}*/

#include "Hxi_SimStruct_methods.h"

#if defined(HXI_INLINE_S_FUNCTION)
/* include implementation file for inline S-function */

#undef HXI_SS_NOSET0
#undef HXI_SS_SETGET1
#undef HXI_SS_SET1
#undef HXI_SS_NOSET1
#undef HXI_SS_GET1
#undef HXI_SS_IS1
#undef HXI_SS_SETGET2
#undef HXI_SS_NOSETGET2
#undef HXI_SS_SET2
#undef HXI_SS_NOSET2
#undef HXI_SS_GET2
#undef HXI_SS_IS2

#undef HXI_MX_CREATE1
#undef HXI_MX_CREATE3
#undef HXI_MX_DESTROY
#undef HXI_MX_SETGET1
#undef HXI_MX_GET1
#undef HXI_MX_IS1
#undef HXI_MX_TO1

#include "Hxi_SimStruct.C"

#endif /* defined(HXI_INLINE_S_FUNCTION) */

#endif
