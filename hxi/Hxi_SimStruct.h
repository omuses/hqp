/*
 * Hxi_SimStruct.h: implement SimStruct for evaluating an S-function from HQP
 *
 * rf, 05/05/2001
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

#if !defined(Hxi_SimStruct_H)
#define Hxi_SimStruct_H

#include <assert.h>
#include <vector>

// HXI_REAL_T can be defined before including this file (default: double).
#if !defined(HXI_REAL_T)
#define HXI_REAL_T double
#endif

// definitions expected in S-functions (only a subset is supported)
#define SS_OPTION_EXCEPTION_FREE_CODE 	0x0001
#define CONTINUOUS_SAMPLE_TIME 		0.0
#define UNUSED_ARG(arg) 		(arg)=(arg)

// Macros to access SimStruct (only a subset is supported)
#define ssSetNumSFcnParams(S, np) 	S->setNumSFcnParams(np)
#define ssGetNumSFcnParams(S) 		S->getNumSFcnParams()
#define ssGetSFcnParamsCount(S) 	S->getSFcnParamsCount()
#define ssSetNumContStates(S, nc)	S->setNumContStates(nc)
#define ssGetNumContStates(S)		S->getNumContStates()
#define ssGetContStates(S) 		S->getContStates()
#define ssGetdX(S) 			S->getdX()
#define ssSetNumDiscStates(S, nd)	S->setNumDiscStates(nd)
#define ssGetNumDiscStates(S)		S->getNumDiscStates()
#define ssGetDiscStates(S) 		S->getDiscStates()
#define ssGetRealDiscStates(S) 		ssGetDiscStates(S)
#define ssSetNumInputPorts(S, nports) 	S->setNumInputPorts(nports)
#define ssGetNumInputPorts(S) 		S->getNumInputPorts()
#define ssSetInputPortWidth(S, port, nu) S->setInputPortWidth(port, nu)
#define ssGetInputPortWidth(S, port) 	S->getInputPortWidth(port)
#define ssGetInputPortRealSignal(S, port) S->getInputPortRealSignal(port)
#define ssGetInputPortRealSignalPtrs(S, port) \
  S->getInputPortRealSignalPtrs(port)
#define ssSetInputPortDirectFeedThrough(S, port, dft) \
  S->setInputPortDirectFeedThrough(port, dft)
#define ssGetInputPortDirectFeedThrough(S, port) \
  S->getInputPortDirectFeedThrough(port)
#define ssSetInputPortRequiredContiguous(S, port) \
  S->setInputPortRequiredContiguous(port)
#define ssGetInputPortRequiredContiguous(S, port) \
  S->getInputPortRequiredContiguous(port)
#define ssSetNumOutputPorts(S, nports) 	S->setNumOutputPorts(nports)
#define ssGetNumOutputPorts(S) 		S->getNumOutputPorts()
#define ssSetOutputPortWidth(S, port, ny) S->setOutputPortWidth(port, ny)
#define ssGetOutputPortWidth(S, port) 	S->getOutputPortWidth(port)
#define ssGetOutputPortRealSignal(S, port) S->getOutputPortRealSignal(port)
#define ssGetOutputPortRealSignalPtrs(S, port) \
  S->getOutputPortRealSignalPtrs(port)
#define ssSetNumSampleTimes(S, ns) 	S->setNumSampleTimes(ns)
#define ssGetNumSampleTimes(S) 		S->getNumSampleTimes()
#define ssSetSampleTime(S, idx, val) 	S->setSampleTime(idx, val)
#define ssGetSampleTime(S, idx) 	S->getSampleTime(idx)
#define ssSetOffsetTime(S, idx, val) 	S->setOffsetTime(idx, val)
#define ssGetOffsetTime(S, idx) 	S->getOffsetTime(idx)
#define ssSetNumRWork(S, nrw) 		S->setNumRWork(nrw)
#define ssGetNumRWork(S) 		S->getNumRWork()
#define ssGetRWork(S) 			S->getRWork()
#define ssSetNumIWork(S, niw) 		S->setNumIWork(niw)
#define ssGetNumIWork(S) 		S->getNumIWork()
#define ssGetIWork(S) 			S->getIWork()
#define ssSetNumPWork(S, npw) 		S->setNumPWork(npw)
#define ssGetNumPWork(S) 		S->getNumPWork()
#define ssGetPWork(S) 			S->getPWork()
#define ssSetNumModes(S, nm)		S->setNumModes(nm)
#define ssGetNumModes(S)		S->getNumModes()
#define ssGetModeVector(S)		S->getModeVector()
#define ssSetNumNonsampledZCs(S, nzcs) 	S->setNumNonsampledZCs(nzcs)
#define ssGetNumNonsampledZCs(S) 	S->getNumNonsampledZCs()
#define ssSetOptions(S, opts) 		S->setOptions(opts)
#define ssGetOptions(S) 		S->getOptions()
#define ssSetT(S, t)			S->setT(t)
#define ssGetT(S)			S->getT()

/** Real type used in S-function. */
typedef HXI_REAL_T real_T;
typedef real_T **InputRealPtrsType;

/** Integer type used in S-function. */
typedef int int_T;

/** Unsigned integer type used in S-function. */
typedef unsigned uint_T;

/** SimStruct for HQP */
class SimStruct {
protected:
  real_T	 _t;		// current simulation time
  vector<real_T> _p; 		// parameters
  vector<real_T> _xc; 		// continuous states
  vector<real_T> _dxc; 		// derivatives of continuous states
  vector<real_T> _xd; 		// discrete states
  vector< vector<real_T> > _u; 	// inputs
  vector< vector<real_T *> > _uPtrs; // pointers to inputs
  vector<int_T>  _u_dft;	// mark if input port is accessed in mdlOutputs
  vector< vector<real_T> > _y; 	// outputs
  vector< vector<real_T *> > _yPtrs; // pointers to outputs
  vector<real_T> _st_period; 	// sample time periods
  vector<real_T> _st_offset; 	// sample time offsets

  vector<real_T> _rwork; 	// real work array
  vector<int_T>  _iwork; 	// int work array
  vector<void *> _pwork; 	// pointer work array
  vector<int_T>  _modes;	// modes array

  uint_T 	_options; 	// option flags

public:
  /** Constructor sets default values. */
  SimStruct()
  {
    _t = 0.0;
    setNumInputPorts(1);
    setNumOutputPorts(1);
  }

  /** Set number of parameters. */
  int_T setNumSFcnParams(int_T np) {
    _p.resize(np);
    return _p.size();
  }
  /** Get number of parameters. */
  int_T getNumSFcnParams() {
    return _p.size();
  }
  /** Get number of externally provided parameters.
      Currently just the size of the parameter vector is returned. */
  int_T getSFcnParamsCount() {
    return _p.size();
  }

  /** Set number of continuous states. */
  int_T setNumContStates(int_T nc) {
    _xc.resize(nc);
    _dxc.resize(nc);
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
    _xd.resize(nd);
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
      _u_dft[i] = 0;	// initialize
    return _u.size();
  }
  /** Get number of input ports. */
  int_T getNumInputPorts() {
    return _u.size();
  }
  /** Set number of inputs of a port. */
  int_T setInputPortWidth(int_T port, int_T nu) {
    _u[port].resize(nu);
    _uPtrs[port].resize(nu);
    for (int_T i = 0; i < nu; i++)
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
  int_T setInputPortRequiredContiguous(int_T port) {
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
    _yPtrs.resize(nports);
    return _y.size();
  }
  /** Get number of output ports. */
  int_T getNumOutputPorts() {
    return _y.size();
  }
  /** Set number of outputs of a port. */
  int_T setOutputPortWidth(int_T port, int_T ny) {
    _y[port].resize(ny);
    _yPtrs[port].resize(ny);
    for (int_T i = 0; i < ny; i++)
      _yPtrs[port][i] = &_y[port][i];
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
  /** Get pointer to pointers to outputs of a port. */
  InputRealPtrsType getOutputPortRealSignalPtrs(int_T port) {
    return &_yPtrs[port][0];
  }

  /** Set number of sample times. */
  int_T setNumSampleTimes(int_T ns) {
    _st_period.resize(ns);
    _st_offset.resize(ns);
    return _st_period.size();
  }
  /** Get number of sample times. */
  int_T getNumSampleTimes() {
    return _st_period.size();
  }
  /** Set sample time period for given st_index. */
  real_T setSampleTime(int_T st_index, real_T period) {
    _st_period[st_index] = period;
    return _st_period[st_index];
  }
  /** Get sample time period for given st_index. */
  real_T getSampleTime(int_T st_index) {
    return _st_period[st_index];
  }
  /** Set offset time for given st_index. */
  real_T setOffsetTime(int_T st_index, real_T offset) {
    _st_offset[st_index] = offset;
    return _st_offset[st_index];
  }
  /** Get offset time for given st_index. */
  real_T getOffsetTime(int_T st_index) {
    return _st_offset[st_index];
  }

  /** Set size of real work vector. */
  int_T setNumRWork(int_T nrw) {
    _rwork.resize(nrw);
    return _rwork.size();
  }
  /** Get size of real work vector. */
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

  /** Set number of states for which zero crossings may occur.
      Note: state events are not supported by HQP. */
  int_T setNumNonsampledZCs(int_T nzcs) {
    assert(nzcs == 0);
    return 0;
  }
  /** Get number of states for which zero crossings may occur.
      Note: state events are not supported by HQP. */
  int_T getNumNonsampledZCs(int_T nzcs) {
    return 0;
  }

  /** Set option flags. */
  uint_T setOptions(uint_T opts) {
    _options = opts;
    return _options;
  }
  /** Get option flags. */
  uint_T getOptions() {
    return _options;
  }

  /** Set current simulation time. */
  real_T setT(real_T t) {
    _t = t;
    return _t;
  }
  /** Get current simulation time. */
  real_T getT() {
    return _t;
  }
};

#endif
