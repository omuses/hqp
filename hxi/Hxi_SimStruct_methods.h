/**
 * @file Hxi_SimStruct_methods.h
 *    Simulink(R) S-function methods supported by Hqp.
 *
 * (Simulink is a registered trademark of The MathWorks, Inc.)
 *
 * rf, 01/15/2005
 */

/** @name Macros to access SimStruct (only a subset is supported) */
/*@{*/
HXI_SS_SETGET1(NumSFcnParams, int_T, np);
HXI_SS_SETGET1(SFcnParamsCount, int_T, num);
HXI_SS_SETGET2(SFcnParam, int_T, idx, /* const */ mxArray*, ptr);
HXI_SS_SETGET1(NumContStates, int_T, nc);
HXI_SS_GET1(ContStates, real_T*);
HXI_SS_GET1(dX, real_T*);
HXI_SS_SETGET1(NumDiscStates, int_T, nd);
HXI_SS_GET1(DiscStates, real_T*);
HXI_SS_GET1(RealDiscStates, real_T*);

HXI_SS_SETGET1(NumInputPorts, int_T, nports);
HXI_SS_SETGET2(InputPortWidth, int_T, port, int_T, nu);
HXI_SS_SET2(InputPortVectorDimension, int_T, port, int_T, nu);
HXI_SS_GET2(InputPortRealSignal, int_T, port, const real_T*);
HXI_SS_GET2(InputPortRealSignalPtrs, int_T, port, InputRealPtrsType);
HXI_SS_SETGET2(InputPortDirectFeedThrough, int_T, port, int_T, dft);
HXI_SS_SETGET2(InputPortRequiredContiguous, int_T, port, int_T, flag);
HXI_SS_NOSET2(InputPortOverWritable, int_T, port, int_T, val);
HXI_SS_NOSET2(InputPortSampleTime, int_T, port, int_T, val);
HXI_SS_NOSET2(InputPortOffsetTime, int_T, port, int_T, val);
HXI_SS_NOSET2(InputPortOptimOpts, int_T, port, int_T, val);

HXI_SS_SETGET1(NumOutputPorts, int_T, nports);
HXI_SS_SETGET2(OutputPortWidth, int_T, port, int_T, nu);
HXI_SS_SET2(OutputPortVectorDimension, int_T, port, int_T, nu);
HXI_SS_GET2(OutputPortRealSignal, int_T, port, real_T*);
HXI_SS_GET2(OutputPortSignal, int_T, port, void*);
HXI_SS_NOSET2(OutputPortSampleTime, int_T, port, int_T, val);
HXI_SS_NOSET2(OutputPortOffsetTime, int_T, port, int_T, val);
HXI_SS_NOSET2(OutputPortOptimOpts, int_T, port, int_T, val);

HXI_SS_SETGET1(NumSampleTimes, int_T, nst);
HXI_SS_SETGET2(SampleTime, int_T, idx, real_T, val);
HXI_SS_SETGET2(OffsetTime, int_T, idx, real_T, val);
HXI_SS_GET1(SampleHitPtr, int_T*);
HXI_SS_GET2(SampleTimeTaskID, int_T, sti, int_T);
HXI_SS_IS1(ContinuousTask, int_T, tid);
HXI_SS_IS2(SampleHit, int_T, st_index, int_T, tid);
HXI_SS_SETGET1(NumNonsampledZCs, int_T, nzcs);
HXI_SS_GET1(NonsampledZCs, real_T*);

HXI_SS_SETGET1(JacobianNzMax, int_T, nnz);
HXI_SS_GET1(JacobianPr, real_T*);
HXI_SS_GET1(JacobianIr, int_T*);
HXI_SS_GET1(JacobianJc, int_T*);

HXI_SS_SETGET1(NumRWork, int_T, nrw);
HXI_SS_GET1(RWork, real_T*);
HXI_SS_SETGET1(NumIWork, int_T, niw);
HXI_SS_GET1(IWork, int_T*);
HXI_SS_SETGET1(NumPWork, int_T, npw);
HXI_SS_GET1(PWork, void**);
HXI_SS_SETGET1(NumDWork, int_T, ndw);
HXI_SS_SETGET2(DWorkWidth, int_T, idx, int_T, width);
HXI_SS_NOSETGET2(DWorkName, int_T, idx, const char*, name);
HXI_SS_SETGET2(DWorkUsedAsDState, int_T, idx, int_T, usage);
HXI_SS_GET2(DWork, int_T, idx, void*);
HXI_SS_SETGET1(NumModes, int_T, nm);
HXI_SS_GET1(ModeVector, int_T*);
HXI_SS_SETGET1(UserData, void*, ptr);
HXI_SS_SETGET1(Options, uint_T, opts);
HXI_SS_SETGET1(T, real_T, t);
HXI_SS_SETGET1(ModelName, const char_T*, name);
HXI_SS_SETGET1(Path, const char_T*, path);
HXI_SS_SETGET1(Version, int_T, ver);
HXI_SS_SETGET1(ErrorStatus, const char_T*, msg);
HXI_SS_NOSET1(RTWGeneratedSFcn, int_T, val);
HXI_SS_NOSET1(Checksum0, uint_T, val);
HXI_SS_NOSET1(Checksum1, uint_T, val);
HXI_SS_NOSET1(Checksum2, uint_T, val);
HXI_SS_NOSET1(Checksum3, uint_T, val);
/*@}*/

/** @name Macros for specifying solver information */
/*@{*/
HXI_SS_SET1(MinorTimeStep, int_T, step);
HXI_SS_IS0(MinorTimeStep);
HXI_SS_IS0(MajorTimeStep);
HXI_SS_SET1(VariableStepSolver, int_T, val);
HXI_SS_IS0(VariableStepSolver);
HXI_SS_SETGET1(SolverMaxOrder, int_T, order);
HXI_SS_SETGET1(SolverName, const char_T*, name);
HXI_SS_NOSET0(SolverNeedsReset);
/*@}*/

/** @name Macros for accessing S-function methods */
/*@{*/
HXI_SS_SET1(mdlInitializeSizes, SFunctionMethod1_type*, m);
HXI_SS_SETGET1(mdlCheckParameters, SFunctionMethod1_type*, m);
HXI_SS_SET1(mdlInitializeSampleTimes, SFunctionMethod1_type*, m);
HXI_SS_SETGET1(mdlStart, SFunctionMethod1_type*, m);
HXI_SS_SETGET1(mdlInitializeConditions, SFunctionMethod1_type*, m);
HXI_SS_SET1(mdlOutputs, SFunctionMethod2_type*, m);
HXI_SS_SETGET1(mdlUpdate, SFunctionMethod2_type*, m);
HXI_SS_SETGET1(mdlDerivatives, SFunctionMethod1_type*, m);
HXI_SS_SETGET1(mdlJacobian, SFunctionMethod1_type*, m);
HXI_SS_SET1(mdlTerminate, SFunctionMethod1_type*, m);
/*@}*/

/** @name mxArray macros */
/*@{*/
HXI_MX_CREATE3(DoubleMatrix, int_T, m, int_T, n, mxComplexity, ty);
HXI_MX_CREATE1(String, const char_T*, s);
HXI_MX_DESTROY(Array);
HXI_MX_TO1(String, char*);
HXI_MX_GET1(NumberOfElements, int_T);
HXI_MX_SETGET1(M, int_T, m);
HXI_MX_SETGET1(N, int_T, n);
HXI_MX_GET1(Pr, real_T*);
HXI_MX_IS0(Empty);
HXI_MX_IS0(Char);
HXI_MX_IS0(Double);
HXI_MX_IS0(Sparse);
HXI_MX_IS0(Complex);
HXI_MX_IS0(Numeric);
/*@}*/
