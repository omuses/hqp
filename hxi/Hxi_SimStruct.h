/**
 * @file Hxi_SimStruct.h
 *    Native %SimStruct for compiling a Simulink(R) S-function with Hqp.
 *
 * (Simulink is a registered trademark of The MathWorks, Inc.)
 *
 * rf, 05/05/2001
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
#if !defined(Hxi_SimStruct_H)
#define Hxi_SimStruct_H

#include <assert.h>
#include <vector>

#include "Hxi_sfun_types.h"
#include "Hxi_mxArray.h"

using namespace std;

/// @name Definitions expected in S-functions (only a subset is supported)
//@{
#define SS_OPTION_EXCEPTION_FREE_CODE 		0x0001
#define SS_OPTION_RUNTIME_EXCEPTION_FREE_CODE 	0x0002
#define SS_OPTION_PORT_SAMPLE_TIMES_ASSIGNED 	0x0004
#define SS_REUSABLE_AND_LOCAL 			0x0100
#define SS_NOT_REUSABLE_AND_GLOBAL 		0x0200

#define SS_STDIO_AVAILABLE 			true
#define CONTINUOUS_SAMPLE_TIME 			0.0
#define MINOR_TIME_STEP				0
#define MAJOR_TIME_STEP				1
//@}

/// Empty implementation for some macros.
#define HXI_NOT_IMPLEMENTED

/// @name Macros to access SimStruct (only a subset is supported)
//@{
#define ssSetNumSFcnParams(S, np) 	(S)->setNumSFcnParams(np)
#define ssGetNumSFcnParams(S) 		(S)->getNumSFcnParams()
#define ssSetSFcnParamsCount(S, num) 	(S)->setSFcnParamsCount(num)
#define ssGetSFcnParamsCount(S) 	(S)->getSFcnParamsCount()
#define ssSetSFcnParam(S, idx, ptr) 	(S)->setSFcnParam(idx, ptr)
#define ssGetSFcnParam(S, idx) 		(S)->getSFcnParam(idx)
#define ssSetNumContStates(S, nc)	(S)->setNumContStates(nc)
#define ssGetNumContStates(S)		(S)->getNumContStates()
#define ssGetContStates(S) 		(S)->getContStates()
#define ssGetdX(S) 			(S)->getdX()
#define ssSetNumDiscStates(S, nd)	(S)->setNumDiscStates(nd)
#define ssGetNumDiscStates(S)		(S)->getNumDiscStates()
#define ssGetDiscStates(S) 		(S)->getDiscStates()
#define ssGetRealDiscStates(S) 		ssGetDiscStates(S)
#define ssSetNumInputPorts(S, nports) 	(S)->setNumInputPorts(nports)
#define ssGetNumInputPorts(S) 		(S)->getNumInputPorts()
#define ssSetInputPortWidth(S, port, nu) (S)->setInputPortWidth(port, nu)
#define ssSetInputPortVectorDimension(S, port, ny) \
  ssSetInputPortWidth(S, port, ny)
#define ssGetInputPortWidth(S, port) 	(S)->getInputPortWidth(port)
#define ssGetInputPortRealSignal(S, port) (S)->getInputPortRealSignal(port)
#define ssGetInputPortRealSignalPtrs(S, port) \
  (S)->getInputPortRealSignalPtrs(port)
#define ssGetInputPortSignalPtrs(S, port) \
  ssGetInputPortRealSignalPtrs(S, port) 
#define ssSetInputPortDirectFeedThrough(S, port, dft) \
  (S)->setInputPortDirectFeedThrough(port, dft)
#define ssGetInputPortDirectFeedThrough(S, port) \
  (S)->getInputPortDirectFeedThrough(port)
#define ssSetInputPortRequiredContiguous(S, port, flag) \
  (S)->setInputPortRequiredContiguous(port, flag)
#define ssGetInputPortRequiredContiguous(S, port) \
  (S)->getInputPortRequiredContiguous(port)
#define ssSetInputPortOverWritable(S, port, val) HXI_NOT_IMPLEMENTED
#define ssSetInputPortSampleTime(S, port, val)   HXI_NOT_IMPLEMENTED
#define ssSetInputPortOffsetTime(S, port, val)   HXI_NOT_IMPLEMENTED
#define ssSetInputPortOptimOpts(S, port, val)    HXI_NOT_IMPLEMENTED
#define ssSetNumOutputPorts(S, nports) 	(S)->setNumOutputPorts(nports)
#define ssGetNumOutputPorts(S) 		(S)->getNumOutputPorts()
#define ssSetOutputPortWidth(S, port, ny) (S)->setOutputPortWidth(port, ny)
#define ssSetOutputPortVectorDimension(S, port, ny) \
  ssSetOutputPortWidth(S, port, ny)
#define ssGetOutputPortWidth(S, port) 	(S)->getOutputPortWidth(port)
#define ssGetOutputPortRealSignal(S, port) (S)->getOutputPortRealSignal(port)
#define ssGetOutputPortSignal(S, port)  ssGetOutputPortRealSignal(S, port)
#define ssSetOutputPortSampleTime(S, port, val) HXI_NOT_IMPLEMENTED
#define ssSetOutputPortOffsetTime(S, port, val) HXI_NOT_IMPLEMENTED
#define ssSetOutputPortOptimOpts(S, port, val)  HXI_NOT_IMPLEMENTED
#define ssSetNumSampleTimes(S, nst) 	(S)->setNumSampleTimes(nst)
#define ssGetNumSampleTimes(S) 		(S)->getNumSampleTimes()
#define ssSetSampleTime(S, idx, val) 	(S)->setSampleTime(idx, val)
#define ssGetSampleTime(S, idx) 	(S)->getSampleTime(idx)
#define ssSetOffsetTime(S, idx, val) 	(S)->setOffsetTime(idx, val)
#define ssGetOffsetTime(S, idx) 	(S)->getOffsetTime(idx)
#define ssSetNumNonsampledZCs(S, nzcs) 	(S)->setNumNonsampledZCs(nzcs)
#define ssGetNumNonsampledZCs(S) 	(S)->getNumNonsampledZCs()
#define ssGetNonsampledZCs(S) 		(S)->getNumNonsampledZCs()
#define ssSetJacobianNzMax(S, nnz) 	(S)->setJacobianNzMax(nnz)
#define ssGetJacobianNzMax(S) 		(S)->getJacobianNzMax()
#define ssGetJacobianPr(S) 		(S)->getJacobianPr()
#define ssGetJacobianIr(S) 		(S)->getJacobianIr()
#define ssGetJacobianJc(S) 		(S)->getJacobianJc()
#define ssIsContinuousTask(S, tid) \
  (ssGetSampleTime(S, tid) == 0.0 && ssGetOffsetTime(S, tid) == 0.0)
#define ssSetNumRWork(S, nrw) 		(S)->setNumRWork(nrw)
#define ssGetNumRWork(S) 		(S)->getNumRWork()
#define ssGetRWork(S) 			(S)->getRWork()
#define ssSetNumIWork(S, niw) 		(S)->setNumIWork(niw)
#define ssGetNumIWork(S) 		(S)->getNumIWork()
#define ssGetIWork(S) 			(S)->getIWork()
#define ssSetNumPWork(S, npw) 		(S)->setNumPWork(npw)
#define ssGetNumPWork(S) 		(S)->getNumPWork()
#define ssGetPWork(S) 			(S)->getPWork()
#define ssSetNumDWork(S, ndw) 		(S)->setNumDWork(ndw)
#define ssGetNumDWork(S) 		(S)->getNumDWork()
#define ssSetDWorkWidth(S, idx, width) 	(S)->setDWorkWidth(idx, width)
#define ssGetDWorkWidth(S, idx) 	(S)->getDWorkWidth(idx)
#define ssSetDWorkName(S, idx, name) 	HXI_NOT_IMPLEMENTED
#define ssGetDWorkName(S, idx) 		HXI_NOT_IMPLEMENTED
#define ssSetDWorkUsedAsDState(S, idx, usage) \
  (S)->setDWorkUsedAsDState(idx, usage)
#define ssGetDWorkUsedAsDState(S, idx)  (S)->getDWorkUsedAsDState(idx)
#define ssGetDWork(S, idx) 		(S)->getDWork(idx)
#define ssSetNumModes(S, nm)		(S)->setNumModes(nm)
#define ssGetNumModes(S)		(S)->getNumModes()
#define ssGetModeVector(S)		(S)->getModeVector()
#define ssSetUserData(S, ptr) 		(S)->setUserData(ptr)
#define ssGetUserData(S) 		(S)->getUserData()
#define ssSetOptions(S, opts) 		(S)->setOptions(opts)
#define ssGetOptions(S) 		(S)->getOptions()
#define ssSetT(S, t)			(S)->setT(t)
#define ssGetT(S)			(S)->getT()
#define ssSetModelName(S, name) 	(S)->setModelName(name)
#define ssGetModelName(S) 		(S)->getModelName()
#define ssSetPath(S, name) 		(S)->setPath(name)
#define ssGetPath(S) 			(S)->getPath()
#define ssSetVersion(S, ver) 		(S)->setVersion(ver)
#define ssGetVersion(S) 		(S)->getVersion()
#define ssSetErrorStatus(S, msg)	(S)->setErrorStatus(msg)
#define ssGetErrorStatus(S)		(S)->getErrorStatus()
#define ssWarning(S, msg) \
  printf("Warning: S-function \"%s\": %s\n", ssGetPath(S), msg)
#define ssPrintf 			printf
#define ssSetRTWGeneratedSFcn(S, val)  	HXI_NOT_IMPLEMENTED
#define ssSetChecksum0(S, val)  	HXI_NOT_IMPLEMENTED
#define ssSetChecksum1(S, val)  	HXI_NOT_IMPLEMENTED
#define ssSetChecksum2(S, val)  	HXI_NOT_IMPLEMENTED
#define ssSetChecksum3(S, val)  	HXI_NOT_IMPLEMENTED
//@}

/// @name Macros for specifying solver information
//@{
#define ssSetSimTimeStep(S, step) 	(S)->setSimTimeStep(step)
#define ssIsMajorTimeStep(S) \
  ((S)->getSimTimeStep() == MAJOR_TIME_STEP)
#define ssSetVariableStepSolver(S, val) (S)->setVariableStepSolver(val)
#define ssIsVariableStepSolver(S) \
  ((S)->getVariableStepSolver() != 0)
#define ssSetSolverMaxOrder(S, order) 	(S)->setSolverMaxOrder(order)
#define ssGetSolverMaxOrder(S) 		(S)->getSolverMaxOrder()
#define ssSetSolverName(S, name) 	(S)->setSolverName(name)
#define ssGetSolverName(S) 		(S)->getSolverName()
#define ssSetSolverNeedsReset(S) 	HXI_NOT_IMPLEMENTED
//@}

/// @name Macros for accessing S-function methods
//@{
#define ssSetmdlInitializeSizes(S, m) 	(S)->setmdlInitializeSizes(m)
#define ssSetmdlCheckParameters(S, m) 	(S)->setmdlCheckParameters(m)
#define ssGetmdlCheckParameters(S) 	(S)->getmdlCheckParameters()
#define ssSetmdlInitializeSampleTimes(S, m) (S)->setmdlInitializeSampleTimes(m)
#define ssSetmdlStart(S, m) 		(S)->setmdlStart(m)
#define ssGetmdlStart(S) 		(S)->getmdlStart()
#define ssSetmdlInitializeConditions(S, m) (S)->setmdlInitializeConditions(m)
#define ssGetmdlInitializeConditions(S) (S)->getmdlInitializeConditions()
#define ssSetmdlOutputs(S, m) 		(S)->setmdlOutputs(m)
#define ssSetmdlUpdate(S, m) 		(S)->setmdlUpdate(m)
#define ssGetmdlUpdate(S) 		(S)->getmdlUpdate()
#define ssSetmdlDerivatives(S, m) 	(S)->setmdlDerivatives(m)
#define ssGetmdlDerivatives(S) 		(S)->getmdlDerivatives()
#define ssSetmdlJacobian(S, m) 		(S)->setmdlJacobian(m)
#define ssGetmdlJacobian(S) 		(S)->getmdlJacobian()
#define ssSetmdlTerminate(S, m) 	(S)->setmdlTerminate(m)
//@}

namespace Hxi {

/**
 *  Native %SimStruct for Hqp.
 */
class SimStruct {
protected:
  /** Argument passed to resize of vectors.
      This is done to avoid memory management problems with ADOL-C 1.8.7. */
  real_T 	 _dummy;
				
  real_T	 _t;		///< current simulation time

  int_T 	 _p_sfun_size;///< number of parameters expected in S-function
  vector<int_T>  _dwork_usage;	///< usage of each data work vector

  vector<const mxArray *> _p; 	///< parameters provided by calling program
  vector<real_T> _xc; 		///< continuous states
  vector<real_T> _dxc; 		///< derivatives of continuous states
  vector<real_T> _xd; 		///< discrete states
  vector< vector<real_T> > _u; 	///< inputs
  vector< vector<real_T *> > _uPtrs; ///< pointers to inputs
  vector<int_T>  _u_dft;      ///< mark if input port is accessed in mdlOutputs
  vector< vector<real_T> > _y; 	///< outputs
  vector< vector<real_T> > _dwork; ///< data work vectors

  vector<real_T> _st_period;	///< sample time period 
  vector<real_T> _st_offset;	///< sample time offset
  vector<real_T> _zc_signals; 	///< values of zero crossings

  vector<real_T> _jacobianPr; 	///< Jacobian elements
  vector<int_T>  _jacobianIr; 	///< Jacobian row indices
  vector<int_T>  _jacobianJc; 	///< Jacobian start indices per column

  vector<real_T> _rwork; 	///< real work array
  vector<int_T>  _iwork; 	///< int work array
  vector<void *> _pwork; 	///< pointer work array
  vector<int_T>  _modes;	///< modes array
  void 		*_userData; 	///< pointer to user data

  uint_T 	_options; 	///< option flags

  const char_T 	*_model_name; 	///< relative name of this model
  const char_T 	*_path; 	///< absolute name of this model
  int_T 	_version; 	///< version of this model
  const char_T 	*_error_msg;	///< used to report errors from S-function

  int_T		_simTimeStep; 	///< simulation time step (major or minor)
  int_T		_variableStepSolver; ///< indicate variable step solver
  int_T		_solverMaxOrder; ///< indicate variable step solver
  const char_T 	*_solverName; 	///< name of ODE solver

  /// @name S-function methods
  //@{
  typedef void (SFunctionMethod1)(SimStruct *S);
  typedef void (SFunctionMethod2)(SimStruct *S, int_T tid);
  SFunctionMethod1 *_mdlInitializeSizes;
  SFunctionMethod1 *_mdlCheckParameters;
  SFunctionMethod1 *_mdlInitializeSampleTimes;
  SFunctionMethod1 *_mdlStart;
  SFunctionMethod1 *_mdlInitializeConditions;
  SFunctionMethod2 *_mdlOutputs;
  SFunctionMethod2 *_mdlUpdate;
  SFunctionMethod1 *_mdlDerivatives;
  SFunctionMethod1 *_mdlJacobian;
  SFunctionMethod1 *_mdlTerminate;
  //@}

public:
  /*
   * Specific methods of this implementation.
   */

  /** Constructor sets default values. */
  SimStruct()
  {
    _t = 0.0;

    _p_sfun_size = 0;
    _userData = NULL;

    _options = 0;

    _model_name = "Hxi_SimStruct";
    _path = "Hxi_SimStruct";
    _version = 0;
    _error_msg = NULL;

    _simTimeStep = MINOR_TIME_STEP;
    _variableStepSolver = 0;
    _solverMaxOrder = 0;
    _solverName = "unspecified";
  }

  /*
   * General S-function methods.
   */

  /** Set number of externally provided parameters. */
  int_T setSFcnParamsCount(int_T num) {
    _p.resize(num);
    return _p.size();
  }
  /** Get number of externally provided parameters. */
  int_T getSFcnParamsCount() {
    return _p.size();
  }
  /** Set number of parameters used by S-function. */
  int_T setNumSFcnParams(int_T np) {
    return _p_sfun_size = np;
  }
  /** Get number of parameters used by S-function. */
  int_T getNumSFcnParams() {
    return _p_sfun_size;
  }
  /** Set mxArray pointing to a parameter. */
  const mxArray *setSFcnParam(int_T idx, const mxArray *val) {
    return (_p[idx] = val);
  }
  /** Get mxArray pointing to a parameter. */
  const mxArray *getSFcnParam(int_T idx) {
    return _p[idx];
  }

  /** Set number of continuous states. */
  int_T setNumContStates(int_T nc) {
    _xc.resize(nc, _dummy);
    _dxc.resize(nc, _dummy);
    return _xc.size();
  }
  /** Get number of continuous states. */
  int_T getNumContStates() {
    return _xc.size();
  }
  /** Get continuous states. */
  real_T *getContStates() {
    return &_xc[0];
  }
  /** Get derivatives of continuous states. */
  real_T *getdX() {
    return &_dxc[0];
  }

  /** Set number of discrete states. */
  int_T setNumDiscStates(int_T nd) {
    _xd.resize(nd, _dummy);
    return _xd.size();
  }
  /** Get number of discrete states. */
  int_T getNumDiscStates() {
    return _xd.size();
  }
  /** Get pointer to discrete states. */
  real_T *getDiscStates() {
    return &_xd[0];
  }

  /** Set number of input ports. */
  int_T setNumInputPorts(int_T nports) {
    _u.resize(nports);
    _uPtrs.resize(nports);
    _u_dft.resize(nports);
    for (int_T i = 0; i < nports; ++i)
      _u_dft[i] = 0;	// initialize direct feed through flags
    return _u.size();
  }
  /** Get number of input ports. */
  int_T getNumInputPorts() {
    return _u.size();
  }
  /** Set number of inputs of a port. */
  int_T setInputPortWidth(int_T port, int_T nu) {
    _u[port].resize(nu, _dummy);
    _uPtrs[port].resize(nu);
    for (int i = 0; i < nu; ++i)
      _uPtrs[port][i] = &_u[port][i];
    return _u[port].size();
  }
  /** Get number of inputs of a port. */
  int_T getInputPortWidth(int_T port) {
    return _u[port].size();
  }
  /** Get pointer to inputs of a port. */
  real_T *getInputPortRealSignal(int_T port) {
    return &_u[port][0];
  }
  /** Get pointer to pointers to inputs of a port. */
  InputRealPtrsType getInputPortRealSignalPtrs(int_T port) {
    return &_uPtrs[port][0];
  }
  /** Specify if input port is used in mdlOutputs or mdlGetTimeOfNextVarHit. */
  int_T setInputPortDirectFeedThrough(int_T port, int_T dft) {
    _u_dft[port] = dft;
    return _u_dft[port];
  }
  /** Get specification for direct feed through. */
  int_T getInputPortDirectFeedThrough(int_T port) {
    return _u_dft[port];
  }
  /** Specify that input port vector should be layed out contiguous in memory.
      Note: here this is always fulfilled. */
  int_T setInputPortRequiredContiguous(int_T port, int_T flag) {
    return 1;
  }
  /** Get specification about contiguous memory layout.
      Note: here this is always fulfilled. */
  int_T getInputPortRequiredContiguous(int_T port) {
    return 1;
  }

  /** Set number of output ports. */
  int_T setNumOutputPorts(int_T nports) {
    _y.resize(nports);
    return _y.size();
  }
  /** Get number of output ports. */
  int_T getNumOutputPorts() {
    return _y.size();
  }
  /** Set number of outputs of a port. */
  int_T setOutputPortWidth(int_T port, int_T ny) {
    _y[port].resize(ny, _dummy);
    return _y[port].size();
  }
  /** Get number of outputs of a port. */
  int_T getOutputPortWidth(int_T port) {
    return _y[port].size();
  }
  /** Get pointer to outputs of a port. */
  real_T *getOutputPortRealSignal(int_T port) {
    return &_y[port][0];
  }

  /** Set number of data vectors. */
  int_T setNumDWork(int_T num) {
    _dwork.resize(num);
    _dwork_usage.resize(num);
    return _dwork.size();
  }
  /** Get number of data work vectors. */
  int_T getNumDWork() {
    return _dwork.size();
  }
  /** Set number of elements of a data work vector. */
  int_T setDWorkWidth(int_T idx, int_T num) {
    _dwork[idx].resize(num, _dummy);
    return _dwork[idx].size();
  }
  /** Get number of elements of a data work vector. */
  int_T getDWorkWidth(int_T idx) {
    return _dwork[idx].size();
  }
  /** Specify that a data work vector is used for discrete states. */
  int_T setDWorkUsedAsDState(int_T idx, int_T usage) {
    return _dwork_usage[idx] = usage;
  }
  /** Determine whether a data work vector is used for discrete states. */
  int_T getDWorkUsedAsDState(int_T idx) {
    return _dwork_usage[idx];
  }
  /** Get pointer to a data work vector. */
  real_T *getDWork(int_T idx) {
    return &_dwork[idx][0];
  }

  /** Set number of sample times. */
  int_T setNumSampleTimes(int_T nst) {
    _st_period.resize(nst, _dummy);
    _st_offset.resize(nst, _dummy);
    return _st_period.size();
  }
  /** Get number of sample times. */
  int_T getNumSampleTimes() {
    return _st_period.size();
  }
  /** Set sample time period for given st_index. */
  real_T setSampleTime(int_T st_index, real_T period) {
    return _st_period[st_index] = period;
  }
  /** Get sample time period for given st_index. */
  real_T getSampleTime(int_T st_index) {
    return _st_period[st_index];
  }
  /** Set offset time for given st_index. */
  real_T setOffsetTime(int_T st_index, real_T offset) {
    return _st_offset[st_index] = offset;
  }
  /** Get offset time for given st_index. */
  real_T getOffsetTime(int_T st_index) {
    return _st_offset[st_index];
  }

  /** Set number of signals for which zero crossings may occur. */
  int_T setNumNonsampledZCs(int_T nzcs) {
    _zc_signals.resize(nzcs);
    return _zc_signals.size();
  }
  /** Get number of signals for which zero crossings may occur. */
  int_T getNumNonsampledZCs() {
    return _zc_signals.size();
  }
  /** Get vector of zero crossing signals. */
  real_T *getNonsampledZCs(int_T ) {
    return &_zc_signals[0];
  }

  /** Set number of non-zero elements in Jacobian. */
  int_T setJacobianNzMax(int_T nnz) {
    if (nnz < 0) {
      // nnz == -1 allocates space for a full Jacobian
      nnz = (_xc.size() + _xd.size() + _y.size())
        * (_xc.size() + _dxc.size() + _u.size());
    }
    _jacobianPr.resize(nnz, _dummy);
    _jacobianIr.resize(nnz);
    _jacobianJc.resize(nnz);
    return _jacobianPr.size();
  }
  /** Get number of non-zero elements in Jacobian. */
  int_T getJacobianNzMax() {
    return _jacobianPr.size();
  }
  /** Get pointer to first Jacobian element. */
  real_T *getJacobianPr() {
    return &_jacobianPr[0];
  }
  /** Get pointer to first Jacobian index. */
  int_T *getJacobianIr() {
    return &_jacobianIr[0];
  }
  /** Get pointer to first column start index of Jacobian. */
  int_T *getJacobianJc() {
    return &_jacobianJc[0];
  }

  /** Set size of real work array. */
  int_T setNumRWork(int_T nrw) {
    _rwork.resize(nrw, _dummy);
    return _rwork.size();
  }
  /** Get size of real work array. */
  int_T getNumRWork() {
    return _rwork.size();
  }
  /** Get pointer to first element of real work vector. */
  real_T *getRWork() {
    return &_rwork[0];
  }

  /** Set size of int work vector. */
  int_T setNumIWork(int_T niw) {
    _iwork.resize(niw);
    return _iwork.size();
  }
  /** Get size of int work vector. */
  int_T getNumIWork() {
    return _iwork.size();
  }
  /** Get pointer to first element of int work vector. */
  int_T *getIWork() {
    return &_iwork[0];
  }

  /** Set size of pointer work vector. */
  int_T setNumPWork(int_T npw) {
    _pwork.resize(npw);
    return _pwork.size();
  }
  /** Get size of pointer work vector. */
  int_T getNumPWork() {
    return _pwork.size();
  }
  /** Get pointer to first element of pointer work vector. */
  void **getPWork() {
    return &_pwork[0];
  }

  /** Set size of modes vector. */
  int_T setNumModes(int_T nm) {
    _modes.resize(nm);
    for (int_T i = 0; i < nm; i++)
      _modes[i] = 0; // initialize with zeros
    return _modes.size();
  }
  /** Get size of modes vector. */
  int_T getNumModes() {
    return _modes.size();
  }
  /** Get pointer to first element of modes vector. */
  int_T *getModeVector() {
    return &_modes[0];
  }

  /** Set pointer to user data. */
  void *setUserData(void *ptr) {
    return _userData = ptr;
  }
  /** Get pointer to user data. */
  void *getUserData() {
    return _userData;
  }

  /** Set option flags. */
  uint_T setOptions(uint_T opts) {
    return _options = opts;
  }
  /** Get option flags. */
  uint_T getOptions() {
    return _options;
  }

  /** Set current simulation time. */
  real_T setT(real_T t) {
    return _t = t;
  }
  /** Get current simulation time. */
  real_T getT() {
    return _t;
  }

  /** Set relative model name. */
  const char_T *setModelName(const char_T *name) {
    return _model_name = name;
  }
  /** Get relative model name. */
  const char_T *getModelName() const {
    return _model_name;
  }
  /** Set absolute model name. */
  const char_T *setPath(const char_T *path) {
    return _path = path;
  }
  /** Get absolute model name. */
  const char_T *getPath() const {
    return _path;
  }
  /** Set model version. */
  int_T setVersion(int_T version) {
    return _version = version;
  }
  /** Get model version. */
  int_T getVersion() const {
    return _version;
  }

  /** Set error message. */
  const char_T *setErrorStatus(const char_T *msg) {
    return _error_msg = msg;
  }
  /** Get error message. */
  const char_T *getErrorStatus() const {
    return _error_msg;
  }

  /** Set simulation time stepping mode */
  int_T setSimTimeStep(int_T val) {
    return _simTimeStep = val;
  }
  /** Get simulation time stepping mode */
  int_T getSimTimeStep() {
    return _simTimeStep;
  }

  /** Set variable step solver */
  int_T setVariableStepSolver(int_T val) {
    return _variableStepSolver = val;
  }
  /* Get variable step solver */
  int_T getVariableStepSolver() {
    return _variableStepSolver;
  }

  /** Set solver max order */
  int_T setSolverMaxOrder(int_T val) {
    return _solverMaxOrder = val;
  }
  /* Get solver max order */
  int_T getSolverMaxOrder() {
    return _solverMaxOrder;
  }

  /** Set name of ODE solver. */
  const char_T *setSolverName(const char_T *name) {
    return _solverName = name;
  }
  /** Get name of ODE solver. */
  const char_T *getSolverName() const {
    return _solverName;
  }

  /** @name Access methods for S-function methods */
  //@{
  SFunctionMethod1 *setmdlInitializeSizes(SFunctionMethod1 *m) {
    return _mdlInitializeSizes = m;
  }
  SFunctionMethod1 *getmdlInitializeSizes() const {
    return _mdlInitializeSizes;
  }

  SFunctionMethod1 *setmdlCheckParameters(SFunctionMethod1 *m) {
    return _mdlCheckParameters = m;
  }
  SFunctionMethod1 *getmdlCheckParameters() const {
    return _mdlCheckParameters;
  }

  SFunctionMethod1 *setmdlInitializeSampleTimes(SFunctionMethod1 *m) {
    return _mdlInitializeSampleTimes = m;
  }
  SFunctionMethod1 *getmdlInitializeSampleTimes() const {
    return _mdlInitializeSampleTimes;
  }

  SFunctionMethod1 *setmdlStart(SFunctionMethod1 *m) {
    return _mdlStart = m;
  }
  SFunctionMethod1 *getmdlStart() const {
    return _mdlStart;
  }

  SFunctionMethod1 *setmdlInitializeConditions(SFunctionMethod1 *m) {
    return _mdlInitializeConditions = m;
  }
  SFunctionMethod1 *getmdlInitializeConditions() const {
    return _mdlInitializeConditions;
  }

  SFunctionMethod2 *setmdlOutputs(SFunctionMethod2 *m) {
    return _mdlOutputs = m;
  }
  SFunctionMethod2 *getmdlOutputs() const {
    return _mdlOutputs;
  }

  SFunctionMethod2 *setmdlUpdate(SFunctionMethod2 *m) {
    return _mdlUpdate = m;
  }
  SFunctionMethod2 *getmdlUpdate() const {
    return _mdlUpdate;
  }

  SFunctionMethod1 *setmdlDerivatives(SFunctionMethod1 *m) {
    return _mdlDerivatives = m;
  }
  SFunctionMethod1 *getmdlDerivatives() const {
    return _mdlDerivatives;
  }

  SFunctionMethod1 *setmdlJacobian(SFunctionMethod1 *m) {
    return _mdlJacobian = m;
  }
  SFunctionMethod1 *getmdlJacobian() const {
    return _mdlJacobian;
  }

  SFunctionMethod1 *setmdlTerminate(SFunctionMethod1 *m) {
    return _mdlTerminate = m;
  }
  SFunctionMethod1 *getmdlTerminate() const {
    return _mdlTerminate;
  }
  //@}
};

}; // namespace Hxi

#if defined(HXI_INLINE_S_FUNCTION)

/** Construct a SimStruct. */
static Hxi::SimStruct *Hxi_SimStruct_create() {
  void Hxi_SimStruct_init(Hxi::SimStruct *S);
  Hxi::SimStruct *S = new Hxi::SimStruct();
  Hxi_SimStruct_init(S);
  return S;
}

/** Delete a SimStruct. */
static void Hxi_SimStruct_destroy(Hxi::SimStruct *S) {
  delete S;
}

#endif // defined(HXI_INLINE_S_FUNCTION)

#endif
