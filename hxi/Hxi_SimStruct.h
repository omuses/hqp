/*
 * Hxi_SimStruct.h:
 *   SimStruct for compiling a Simulink(R) S-function for HQP
 *
 * (Simulink is a registered trademark of The MathWorks, Inc.)
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

#include "Hxi_sfun_types.h"
#include "Hxi_mxArray.h"

// definitions expected in S-functions (only a subset is supported)
#define SS_OPTION_EXCEPTION_FREE_CODE 	0x0001
#define SS_STDIO_AVAILABLE 		true

#define CONTINUOUS_SAMPLE_TIME 		0.0

// Macros to access SimStruct (only a subset is supported)
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
#define ssGetInputPortWidth(S, port) 	(S)->getInputPortWidth(port)
#define ssGetInputPortRealSignal(S, port) (S)->getInputPortRealSignal(port)
#define ssGetInputPortRealSignalPtrs(S, port) \
  (S)->getInputPortRealSignalPtrs(port)
#define ssSetInputPortDirectFeedThrough(S, port, dft) \
  (S)->setInputPortDirectFeedThrough(port, dft)
#define ssGetInputPortDirectFeedThrough(S, port) \
  (S)->getInputPortDirectFeedThrough(port)
#define ssSetInputPortRequiredContiguous(S, port) \
  (S)->setInputPortRequiredContiguous(port)
#define ssGetInputPortRequiredContiguous(S, port) \
  (S)->getInputPortRequiredContiguous(port)
#define ssSetNumOutputPorts(S, nports) 	(S)->setNumOutputPorts(nports)
#define ssGetNumOutputPorts(S) 		(S)->getNumOutputPorts()
#define ssSetOutputPortWidth(S, port, ny) (S)->setOutputPortWidth(port, ny)
#define ssGetOutputPortWidth(S, port) 	(S)->getOutputPortWidth(port)
#define ssGetOutputPortRealSignal(S, port) (S)->getOutputPortRealSignal(port)
#define ssGetOutputPortRealSignalPtrs(S, port) \
  (S)->getOutputPortRealSignalPtrs(port)
#define ssSetNumSampleTimes(S, ns) 	(S)->setNumSampleTimes(ns)
#define ssGetNumSampleTimes(S) 		(S)->getNumSampleTimes()
#define ssSetSampleTime(S, idx, val) 	(S)->setSampleTime(idx, val)
#define ssGetSampleTime(S, idx) 	(S)->getSampleTime(idx)
#define ssSetOffsetTime(S, idx, val) 	(S)->setOffsetTime(idx, val)
#define ssGetOffsetTime(S, idx) 	(S)->getOffsetTime(idx)
#define ssSetNumRWork(S, nrw) 		(S)->setNumRWork(nrw)
#define ssGetNumRWork(S) 		(S)->getNumRWork()
#define ssGetRWork(S) 			(S)->getRWork()
#define ssSetNumIWork(S, niw) 		(S)->setNumIWork(niw)
#define ssGetNumIWork(S) 		(S)->getNumIWork()
#define ssGetIWork(S) 			(S)->getIWork()
#define ssSetNumPWork(S, npw) 		(S)->setNumPWork(npw)
#define ssGetNumPWork(S) 		(S)->getNumPWork()
#define ssGetPWork(S) 			(S)->getPWork()
#define ssSetNumModes(S, nm)		(S)->setNumModes(nm)
#define ssGetNumModes(S)		(S)->getNumModes()
#define ssGetModeVector(S)		(S)->getModeVector()
#define ssSetNumNonsampledZCs(S, nzcs) 	(S)->setNumNonsampledZCs(nzcs)
#define ssGetNumNonsampledZCs(S) 	(S)->getNumNonsampledZCs()
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

/**
 * SimStruct for HQP.
 * Note: the standard vector class is used for memory management.
 * Special treatment is required if real_T is defined as adouble,
 * as adoubles should be allocated locally during each evaluation.
 * That is why external memory can be provided for all vector<real_T>
 * to avoid their actual allocation.
 * (Unfortunately the predefined vector type adoublev can not be used
 *  as the S-function treats vectors C-like, using pointers to real_T.)
 */
class SimStruct {
protected:
  real_T	 _t;		// current simulation time

  int_T 	 _p_sfun_size;	// number of parameters expected in S-function
  int_T 	 _xc_size;	// number of continuous states
  int_T 	 _xd_size;	// number of discrete states
  vector<int_T>  _u_width;	// number of inputs per port
  vector<int_T>  _y_width;	// number of outputs per port

  vector<mxArray> _p; 		// parameters provided by calling program
  vector<real_T> _xc; 		// continuous states
  vector<real_T> _dxc; 		// derivatives of continuous states
  vector<real_T> _xd; 		// discrete states
  vector< vector<real_T> > _u; 	// inputs
  vector< vector<real_T *> > _uPtrs; // pointers to inputs
  vector<int_T>  _u_dft;	// mark if input port is accessed in mdlOutputs
  vector< vector<real_T> > _y; 	// outputs
  vector< vector<real_T *> > _yPtrs; // pointers to outputs

  int_T 	 _st_size;	// number of sample times
  real_T 	 _st_period;	// sample time period 
  real_T 	 _st_offset;	// sample time offset

  int_T 	 _rwork_size; 	// size of real work array
  vector<real_T> _rwork; 	// real work array
  vector<int_T>  _iwork; 	// int work array
  vector<void *> _pwork; 	// pointer work array
  vector<int_T>  _modes;	// modes array

  uint_T 	_options; 	// option flags

  real_T	*_xc_ext;	// externally provided memory for cont. states
  real_T	*_dxc_ext;	// externally provided memory for derivatives
  real_T	*_xd_ext;	// externally provided memory for disc. states
  real_T	*_u_ext;	// externally provided memory for inputs
  real_T	*_y_ext;	// externally provided memory for outputs
  real_T	*_rwork_ext;	// externally provided memory for work array

  const char_T 	*_model_name; 	// relative name of this model
  const char_T 	*_path; 	// absolute name of this model
  int_T 	_version; 	// version of this model
  const char_T 	*_error_msg;	// used to report errors from S-function

public:
  /*
   * Specific methods of this implementation.
   */

  /** Constructor sets default values. */
  SimStruct()
  {
    _t = 0.0;

    _p_sfun_size = 0;
    _xc_size = 0;
    _xd_size = 0;
    _st_size = 0;
    _rwork_size = 0;

    _options = 0;

    _xc_ext = NULL;
    _dxc_ext = NULL;
    _xd_ext = NULL;
    _u_ext = NULL;
    _y_ext = NULL;
    _rwork_ext = NULL;

    _model_name = "Hxi_SimStruct";
    _path = "Hxi_SimStruct";
    _version = 1;
    _error_msg = NULL;
  }

  /** Set external memory for continuous states */
  void set_xc_ext(real_T *rp) {
    _xc_ext = rp;
  }
  /** Set external memory for derivatives */
  void set_dxc_ext(real_T *rp) {
    _dxc_ext = rp;
  }
  /** Set external memory for discrete states */
  void set_xd_ext(real_T *rp) {
    _xd_ext = rp;
  }
  /** Set external memory for inputs */
  void set_u_ext(real_T *rp) {
    if (rp != _u_ext) {
      _u_ext = rp;
      // mark pointers uninitialized
      for (int j = 0; j < (int)_uPtrs.size(); ++j) {
	if (_uPtrs[j].size() > 0)
	  _uPtrs[j][0] = NULL;
      }
    }
  }
  /** Set external memory for outputs */
  void set_y_ext(real_T *rp) {
    if (rp != _y_ext) {
      _y_ext = rp;
      // mark pointers uninitialized
      for (int j = 0; j < (int)_yPtrs.size(); ++j) {
	if (_yPtrs[j].size() > 0)
	  _yPtrs[j][0] = NULL;
      }
    }
  }
  /** Set external memory for work array */
  void set_rwork_ext(real_T *rp) {
    _rwork_ext = rp;
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
  /** Set mxArray pointing to parameter. A copy of the passed mxArray
      is stored to avoid memory problems and side effects. */
  const mxArray *setSFcnParam(int_T idx, const mxArray *val) {
    return &(_p[idx] = *val);
  }
  /** Get mxArray pointing to parameter. */
  const mxArray *getSFcnParam(int_T idx) {
    return &_p[idx];
  }

  /** Set number of continuous states. */
  int_T setNumContStates(int_T nc) {
    return _xc_size = nc;
  }
  /** Get number of continuous states. */
  int_T getNumContStates() {
    return _xc_size;
  }
  /** Get continuous states. */
  real_T *getContStates() {
    if (!_xc_ext && (int_T)_xc.size() != _xc_size)
      _xc.resize(_xc_size);
    return _xc_ext? _xc_ext: &_xc[0];
  }
  /** Get derivatives of continuous states. */
  real_T *getdX() {
    if (!_dxc_ext && (int_T)_dxc.size() != _xc_size)
      _dxc.resize(_xc_size);
    return _dxc_ext? _dxc_ext: &_dxc[0];
  }

  /** Set number of discrete states. */
  int_T setNumDiscStates(int_T nd) {
    return _xd_size = nd;
  }
  /** Get number of discrete states. */
  int_T getNumDiscStates() {
    return _xd_size;
  }
  /** Get pointer to discrete states. */
  real_T *getDiscStates() {
    if (!_xd_ext && (int_T)_xd.size() != _xd_size)
      _xd.resize(_xd_size);
    return _xd_ext? _xd_ext: &_xd[0];
  }

  /** Set number of input ports. */
  int_T setNumInputPorts(int_T nports) {
    _u.resize(nports);
    _u_width.resize(nports);
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
    _uPtrs[port].resize(nu);
    if (nu > 0)
      _uPtrs[port][0] = NULL; // mark uninitialized
    return _u_width[port] = nu;
  }
  /** Get number of inputs of a port. */
  int_T getInputPortWidth(int_T port) {
    return _u_width[port];
  }
  /** Get pointer to inputs of a port. */
  real_T *getInputPortRealSignal(int_T port) {
    real_T *rp;
    if (_u_ext) {
      // externally provided inputs have linear memory layout
      rp = _u_ext;
      for (int_T j = 0; j < port; ++j)
	rp += _u[j].size();
    }
    else {
      if ((int_T)_u[port].size() != _u_width[port])
	_u[port].resize(_u_width[port]);
      rp = &_u[port][0];
    }
    return rp;
  }
  /** Get pointer to pointers to inputs of a port. */
  InputRealPtrsType getInputPortRealSignalPtrs(int_T port) {
    // initialize pointers if required
    if (_uPtrs[port].size() > 0 && _uPtrs[port][0] == NULL) {
      real_T *rp = getInputPortRealSignal(port);
      for (int i = 0; i < (int)_uPtrs[port].size(); ++i)
	_uPtrs[port][i] = rp++;
    }
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
    _y_width.resize(nports);
    _yPtrs.resize(nports);
    return _y.size();
  }
  /** Get number of output ports. */
  int_T getNumOutputPorts() {
    return _y.size();
  }
  /** Set number of outputs of a port. */
  int_T setOutputPortWidth(int_T port, int_T ny) {
    _yPtrs[port].resize(ny);
    if (ny > 0)
      _yPtrs[port][0] = NULL; // mark uninitialized
    return _y_width[port] = ny;
  }
  /** Get number of outputs of a port. */
  int_T getOutputPortWidth(int_T port) {
    return _y_width[port];
  }
  /** Get pointer to outputs of a port. */
  real_T *getOutputPortRealSignal(int_T port) {
    real_T *rp;
    if (_y_ext) {
      // externally provided outputs have linear memory layout
      rp = _y_ext;
      for (int_T j = 0; j < port; ++j)
	rp += _y[j].size();
    }
    else {
      if ((int_T)_y[port].size() != _y_width[port])
	_y[port].resize(_y_width[port]);
      rp = &_y[port][0];
    }
    return rp;
  }
  /** Get pointer to pointers to outputs of a port. */
  InputRealPtrsType getOutputPortRealSignalPtrs(int_T port) {
    // initialize pointers if required
    if (_yPtrs[port].size() > 0 && _yPtrs[port][0] == NULL) {
      real_T *rp = getOutputPortRealSignal(port);
      for (int i = 0; i < (int)_yPtrs[port].size(); ++i)
	_yPtrs[port][i] = rp++;
    }
    return &_yPtrs[port][0];
  }

  /** Set number of sample times. Currently at most one is supported. */
  int_T setNumSampleTimes(int_T ns) {
    assert(ns <= 1);
    return _st_size = ns;
  }
  /** Get number of sample times. */
  int_T getNumSampleTimes() {
    return _st_size;
  }
  /** Set sample time period for given st_index. */
  real_T setSampleTime(int_T st_index, real_T period) {
    assert(st_index == 0);
    assert(_st_size == 1);  // must be initialized
    return _st_period = period;
  }
  /** Get sample time period for given st_index. */
  real_T getSampleTime(int_T st_index) {
    assert(st_index == 0);
    assert(_st_size == 1);  // must be initialized
    return _st_period;
  }
  /** Set offset time for given st_index. */
  real_T setOffsetTime(int_T st_index, real_T offset) {
    assert(st_index == 0);
    assert(_st_size == 1);  // must be initialized
    return _st_offset = offset;
  }
  /** Get offset time for given st_index. */
  real_T getOffsetTime(int_T st_index) {
    assert(st_index == 0);
    assert(_st_size == 1);  // must be initialized
    return _st_offset;
  }

  /** Set size of real work array. */
  int_T setNumRWork(int_T nrw) {
    return _rwork_size = nrw;
  }
  /** Get size of real work array. */
  int_T getNumRWork() {
    return _rwork_size;
  }
  /** Get pointer to first element of real work vector. */
  real_T *getRWork() {
    if (!_rwork_ext && (int_T)_rwork.size() != _rwork_size)
      _rwork.resize(_rwork_size);
    return _rwork_ext? _rwork_ext: &_rwork[0];
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
};

#endif
