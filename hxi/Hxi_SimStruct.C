/**
 * @file Hxi_SimStruct.C
 *    Implementation of native %SimStruct for compiling a
 *    Simulink(R) S-function with Hqp.
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

#include <assert.h>

#include <vector>
#include <string>
using namespace std;

#define HQP_MAGIC_CODE 1994.0528 
// (creation date of first Hqp component -- the IP solver)

/*
 * Note: Don't include Hxi_SimStruct.h to separate calls to ss* from
 *       SimStruct implementation. Hxi_SimStruct.h declares ss* as functions;
 *       the implementation is as C macros for MEX and C++ methods for HXI.
 */

#if !defined(Hxi_SimStruct_H)

#if defined(HXI_WITH_MEX)
// support MEX SimStruct and Hxi::SimStruct

// first include MEX SimStruct, in addition to Hxi::SimStruct defined below
#include "Hxi_MEX_SFunction.h"

// additionally include Hxi types for Hxi::SimStruct below
namespace Hxi {
#include "Hxi_sfun_types.h"
}
using Hxi::value;

#else
// only support Hxi::SimStruct

#include "Hxi_sfun_types.h"

struct SimStruct;
struct mxArray;

#endif

// define S-function method types if this hasn't been done already
/** S-function method type. */
typedef void (SFunctionMethod1_type)(SimStruct *S);
/** S-function method type with additional task id argument. */
typedef void (SFunctionMethod2_type)(SimStruct *S, int_T tid);

// forward declaration of default S-function methods (implementation below)
static void defaultSFunctionMethod1(SimStruct *S);
static void defaultSFunctionMethod2(SimStruct *S, int_T);

#endif // !defined(Hxi_SimStruct_H)

// check for Hxi::SimStruct
inline bool isHxiSimStruct(const SimStruct *S) {
  return (value(*((real_T*)S)) == HQP_MAGIC_CODE);
}

namespace Hxi {

/**
 *  Native %SimStruct for Hqp.
 */
template <class REAL_T>
class SimStruct {
protected:
  /** Argument passed to resize of vectors.
      This is done to avoid memory management problems with ADOL-C 1.8.7. */
  REAL_T 	 _dummy;
				
  REAL_T	 _t;		///< current simulation time

  int_T 	 _p_sfun_size;///< number of parameters expected in S-function
  vector<int_T>  _dwork_usage;	///< usage of each data work vector

  vector<const ::mxArray *> _p; 	///< parameters provided by calling program
  vector<REAL_T> _xc; 		///< continuous states
  vector<REAL_T> _dxc; 		///< derivatives of continuous states
  vector<REAL_T> _xd; 		///< discrete states
  vector< vector<REAL_T> > _u; 	///< inputs
  vector< vector<REAL_T *> > _uPtrs; ///< pointers to inputs
  vector<int_T>  _u_dft;      ///< mark if input port is accessed in mdlOutputs
  vector< vector<REAL_T> > _y; 	///< outputs
  vector< vector<REAL_T> > _dwork; ///< data work vectors

  vector<REAL_T> _st_period;	///< sample time period 
  vector<REAL_T> _st_offset;	///< sample time offset
  vector<REAL_T> _zc_signals; 	///< values of zero crossings

  vector<REAL_T> _jacobianPr; 	///< Jacobian elements
  vector<int_T>  _jacobianIr; 	///< Jacobian row indices
  vector<int_T>  _jacobianJc; 	///< Jacobian start indices per column

  vector<REAL_T> _rwork; 	///< real work array
  vector<int_T>  _iwork; 	///< int work array
  vector<void *> _pwork; 	///< pointer work array
  vector<int_T>  _modes;	///< modes array
  void 		*_userData; 	///< pointer to user data

  uint_T 	_options; 	///< option flags

  const char_T 	*_model_name; 	///< relative name of this model
  const char_T 	*_path; 	///< absolute name of this model
  int_T 	_version; 	///< version of this model
  const char_T 	*_error_msg;	///< used to report errors from S-function

  int_T		_minorTimeStep; ///< minor simulation time step 
  int_T		_variableStepSolver; ///< indicate variable step solver
  int_T		_solverMaxOrder; ///< indicate variable step solver
  const char_T 	*_solverName; 	///< name of ODE solver

  /// @name S-function methods
  //@{
  SFunctionMethod1_type *_mdlInitializeSizes;
  SFunctionMethod1_type *_mdlCheckParameters;
  SFunctionMethod1_type *_mdlInitializeSampleTimes;
  SFunctionMethod1_type *_mdlStart;
  SFunctionMethod1_type *_mdlInitializeConditions;
  SFunctionMethod2_type *_mdlOutputs;
  SFunctionMethod2_type *_mdlUpdate;
  SFunctionMethod1_type *_mdlDerivatives;
  SFunctionMethod1_type *_mdlJacobian;
  SFunctionMethod1_type *_mdlTerminate;
  //@}

public:
  /*
   * Specific methods of this implementation.
   */

  /** Constructor sets default values. */
  SimStruct()
  {
    _dummy = HQP_MAGIC_CODE; 
    _t = 0.0;

    _p_sfun_size = 0;
    _userData = NULL;

    _options = 0;

    _model_name = "Hxi_SimStruct";
    _path = "Hxi_SimStruct";
    _version = 0;
    _error_msg = NULL;

    _minorTimeStep = 0;
    _variableStepSolver = 0;
    _solverMaxOrder = 0;
    _solverName = "unspecified";

    // assign defaults to required S-function methods
    // as they are not checked for NULL when called
    _mdlInitializeSizes = &defaultSFunctionMethod1;
    _mdlCheckParameters = NULL;
    _mdlInitializeSampleTimes = &defaultSFunctionMethod1;
    _mdlStart = NULL;
    _mdlInitializeConditions = NULL;
    _mdlOutputs = &defaultSFunctionMethod2;
    _mdlUpdate = NULL;
    _mdlDerivatives = NULL;
    _mdlJacobian = NULL;
    _mdlTerminate = &defaultSFunctionMethod1;
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
  /* const */ ::mxArray *setSFcnParam(int_T idx, const ::mxArray *val) {
    return (::mxArray *)(_p[idx] = val);
  }
  /** Get mxArray pointing to a parameter. */
  /* const */ ::mxArray *getSFcnParam(int_T idx) {
    return (::mxArray *)_p[idx];
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
  REAL_T *getContStates() {
    return &_xc[0];
  }
  /** Get derivatives of continuous states. */
  REAL_T *getdX() {
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
  REAL_T *getDiscStates() {
    return &_xd[0];
  }
  /** Alternative get pointer to discrete states. */
  REAL_T *getRealDiscStates() {
    return getDiscStates();
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
  /** Alternative set number of inputs of a port. */
  int_T setInputPortVectorDimension(int_T port, int_T nu) {
    return setInputPortWidth(port, nu);
  }
  /** Get number of inputs of a port. */
  int_T getInputPortWidth(int_T port) {
    return _u[port].size();
  }
  /** Get pointer to inputs of a port. */
  REAL_T *getInputPortRealSignal(int_T port) {
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
  /** Alternative set number of outputs of a port. */
  int_T setOutputPortVectorDimension(int_T port, int_T ny) {
    return setOutputPortWidth(port, ny);
  }
  /** Get number of outputs of a port. */
  int_T getOutputPortWidth(int_T port) {
    return _y[port].size();
  }
  /** Get pointer to outputs of a port. */
  REAL_T *getOutputPortRealSignal(int_T port) {
    return &_y[port][0];
  }
  /** Alternative get pointer to outputs of a port. */
  void *getOutputPortSignal(int_T port) {
    return getOutputPortRealSignal(port);
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
  void *getDWork(int_T idx) {
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
  REAL_T setSampleTime(int_T st_index, REAL_T period) {
    return _st_period[st_index] = period;
  }
  /** Get sample time period for given st_index. */
  REAL_T getSampleTime(int_T st_index) {
    return _st_period[st_index];
  }
  /** Set offset time for given st_index. */
  REAL_T setOffsetTime(int_T st_index, REAL_T offset) {
    return _st_offset[st_index] = offset;
  }
  /** Get offset time for given st_index. */
  REAL_T getOffsetTime(int_T st_index) {
    return _st_offset[st_index];
  }
  /** Test for continuous task. */
  int_T isContinuousTask(int_T tid) {
    return (getSampleTime(tid) == 0.0 && getOffsetTime(tid) == 0.0);
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
  REAL_T *getNonsampledZCs() {
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
  REAL_T *getJacobianPr() {
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
  REAL_T *getRWork() {
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
  REAL_T setT(REAL_T t) {
    return _t = t;
  }
  /** Get current simulation time. */
  REAL_T getT() {
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

  /** Set minor simulation time stepping mode */
  int_T setMinorTimeStep(int_T val) {
    return _minorTimeStep = val;
  }
  /** Test for minor simulation time stepping mode */
  int_T isMinorTimeStep() {
    return (_minorTimeStep != 0);
  }
  /** Test for major simulation time stepping mode */
  int_T isMajorTimeStep() {
    return (_minorTimeStep == 0);
  }

  /** Set variable step solver */
  int_T setVariableStepSolver(int_T val) {
    return _variableStepSolver = val;
  }
  /* Get variable step solver */
  int_T getVariableStepSolver() {
    return _variableStepSolver;
  }
  /* Test variable step solver */
  int_T isVariableStepSolver() {
    return (getVariableStepSolver() != 0);
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
  SFunctionMethod1_type *setmdlInitializeSizes(SFunctionMethod1_type *m) {
    return _mdlInitializeSizes = m;
  }
  SFunctionMethod1_type *getmdlInitializeSizes() const {
    return _mdlInitializeSizes;
  }

  SFunctionMethod1_type *setmdlCheckParameters(SFunctionMethod1_type *m) {
    return _mdlCheckParameters = m;
  }
  SFunctionMethod1_type *getmdlCheckParameters() const {
    return _mdlCheckParameters;
  }

  SFunctionMethod1_type *setmdlInitializeSampleTimes(SFunctionMethod1_type *m) {
    return _mdlInitializeSampleTimes = m;
  }
  SFunctionMethod1_type *getmdlInitializeSampleTimes() const {
    return _mdlInitializeSampleTimes;
  }

  SFunctionMethod1_type *setmdlStart(SFunctionMethod1_type *m) {
    return _mdlStart = m;
  }
  SFunctionMethod1_type *getmdlStart() const {
    return _mdlStart;
  }

  SFunctionMethod1_type *setmdlInitializeConditions(SFunctionMethod1_type *m) {
    return _mdlInitializeConditions = m;
  }
  SFunctionMethod1_type *getmdlInitializeConditions() const {
    return _mdlInitializeConditions;
  }

  SFunctionMethod2_type *setmdlOutputs(SFunctionMethod2_type *m) {
    return _mdlOutputs = m;
  }
  SFunctionMethod2_type *getmdlOutputs() const {
    return _mdlOutputs;
  }

  SFunctionMethod2_type *setmdlUpdate(SFunctionMethod2_type *m) {
    return _mdlUpdate = m;
  }
  SFunctionMethod2_type *getmdlUpdate() const {
    return _mdlUpdate;
  }

  SFunctionMethod1_type *setmdlDerivatives(SFunctionMethod1_type *m) {
    return _mdlDerivatives = m;
  }
  SFunctionMethod1_type *getmdlDerivatives() const {
    return _mdlDerivatives;
  }

  SFunctionMethod1_type *setmdlJacobian(SFunctionMethod1_type *m) {
    return _mdlJacobian = m;
  }
  SFunctionMethod1_type *getmdlJacobian() const {
    return _mdlJacobian;
  }

  SFunctionMethod1_type *setmdlTerminate(SFunctionMethod1_type *m) {
    return _mdlTerminate = m;
  }
  SFunctionMethod1_type *getmdlTerminate() const {
    return _mdlTerminate;
  }
  //@}
};

/** Native %mxArray for Hqp. */
template <class REAL_T>
class mxArray {
protected:
  REAL_T 	 	_dummy; ///< argument for ADOL-C memory management
  vector<REAL_T> 	_realData; 	///< data vector for array of reals
  string 		_charData; 	///< data for char array
  int_T 		_m;	///< number of rows (more generally first dim)
  int_T 		_n;	///< number of cols (more generally other dims)
  bool 			_isNumeric; ///< indicate numeric data

public:
  /** Constuctor for real array. */
  mxArray(int_T m, int_T n, int_T) : _realData(m*n, _dummy) {
    _isNumeric = true;
    _m = m;
    _n = n;
    for (int i = 0; i < (int)_realData.size(); i++)
      _realData[i] = 0.0;
  }

  /** Constuctor for char array. */
  mxArray(const char *s) : _charData(s) {
    _isNumeric = false;
  }

  /** Convert mxArray to char string and return it as new memory.
      The memory should be freed by the caller using mxFree. */
  char *toString() const {
    char *str = NULL;
    int i;
    if (isChar()) {
      str = (char *)malloc(_charData.size() + 1);
      for (i = 0; i < (int)_charData.size(); i++)
        str[i] = _charData[i];
      str[i] = '\0';
    }
    return str;
  }

  /** Get address of first element of real data. */
  REAL_T *getPr() const {
    return _realData.size() > 0? (REAL_T *)&_realData[0]: NULL;
  }

  /** Set number of data elements. */
  int_T setNumberOfElements(int_T num) {
    _realData.resize(num, _dummy);
    return _realData.size();
  }
  /** Get number of data elements. */
  int_T getNumberOfElements() const {
    return _realData.size();
  }
  /** Set number of rows. */
  int_T setM(int_T m) {
    return _m = m;
  }
  /** Get number of rows. */
  int_T getM() const {
    return _m;
  }
  /** Set number of columns. */
  int_T setN(int_T n) {
    return _n = n;
  }
  /** Get number of columns. */
  int_T getN() const {
    return _n;
  }

  /** Check if data object is empty. */
  bool isEmpty() const {
    return (_realData.size() == 0 && _charData.size() == 0);
  }
  /** Check if data object is char array. */
  bool isChar() const {
    return !_isNumeric;
  }
  /** Check if data object is real array. */
  bool isDouble() const {
    return _isNumeric;
  }
  /** Check if data object is sparse. */
  bool isSparse() const {
    return false;
  }
  /** Check if data object is complex. */
  bool isComplex() const {
    return false;
  }
  /** Check if data object is of numeric type. */
  bool isNumeric() const {
    return _isNumeric;
  }
};

}; // namespace Hxi

#if !defined(HXI_INLINE_S_FUNCTION)
#if !defined(_MSC_VER)
#  if defined(__cplusplus)
#    define HXI_EXTERN extern "C"
#  else
#    define HXI_EXTERN extern
#  endif
#else
#  if defined(__cplusplus)
#    define HXI_EXTERN extern "C" __declspec(dllexport)
#  else
#    define HXI_EXTERN extern __declspec(dllexport)
#  endif
#endif
#endif

#if defined(HXI_WITH_MEX)
/// @name Implementation of S-function methods accessing the SimStruct.
//@{
#define HXI_SS_NOSET0(ITEM)
#define HXI_SS_SETGET1(ITEM, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE ARG) { \
    if (isHxiSimStruct(S)) \
      return ((Hxi::SimStruct<real_T>*)S)->set##ITEM(ARG); \
    else \
      return ssSet##ITEM(S, ARG); \
  } \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S) { \
    if (isHxiSimStruct(S)) \
      return ((Hxi::SimStruct<real_T>*)S)->get##ITEM(); \
    else \
      return ssGet##ITEM(S); \
  }
#define HXI_SS_SET1(ITEM, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE ARG) { \
    if (isHxiSimStruct(S)) \
      return ((Hxi::SimStruct<real_T>*)S)->set##ITEM(ARG); \
    else \
      return ssSet##ITEM(S, ARG); \
  }
#define HXI_SS_NOSET1(ITEM, TYPE, ARG)
#define HXI_SS_GET1(ITEM, TYPE) \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S) { \
    if (isHxiSimStruct(S)) \
      return ((Hxi::SimStruct<real_T>*)S)->get##ITEM(); \
    else \
      return ssGet##ITEM(S); \
  }
#define HXI_SS_IS1(ITEM) \
  HXI_EXTERN int_T hssIs##ITEM(SimStruct *S) { \
    if (isHxiSimStruct(S)) \
      return ((Hxi::SimStruct<real_T>*)S)->is##ITEM(); \
    else \
      return ssIs##ITEM(S); \
  }

#define HXI_SS_SETGET2(ITEM, TYPE1, ARG1, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE1 ARG1, TYPE ARG) { \
    if (isHxiSimStruct(S)) \
      return ((Hxi::SimStruct<real_T>*)S)->set##ITEM(ARG1, ARG);  \
    else \
      return ssSet##ITEM(S, ARG1, ARG); \
  } \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S, TYPE1 ARG1) { \
    if (isHxiSimStruct(S)) \
      return ((Hxi::SimStruct<real_T>*)S)->get##ITEM(ARG1); \
    else \
      return (TYPE)ssGet##ITEM(S, ARG1); \
  }
#define HXI_SS_NOSETGET2(ITEM, TYPE1, ARG1, TYPE, ARG)
#define HXI_SS_SET2(ITEM, TYPE1, ARG1, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE1 ARG1, TYPE ARG) { \
    if (isHxiSimStruct(S)) \
      return ((Hxi::SimStruct<real_T>*)S)->set##ITEM(ARG1, ARG); \
    else \
      return ssSet##ITEM(S, ARG1, ARG); \
  }
#define HXI_SS_NOSET2(ITEM, TYPE1, ARG1, TYPE, ARG) 
#define HXI_SS_GET2(ITEM, TYPE1, ARG1, TYPE) \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S, TYPE1 ARG1) { \
    if (isHxiSimStruct(S)) \
      return ((Hxi::SimStruct<real_T>*)S)->get##ITEM(ARG1); \
    else \
      return ssGet##ITEM(S, ARG1); \
  }
#define HXI_SS_IS2(ITEM, TYPE1, ARG1) \
  HXI_EXTERN int_T hssIs##ITEM(SimStruct *S, TYPE1 ARG1) { \
    if (isHxiSimStruct(S)) \
      return ((Hxi::SimStruct<real_T>*)S)->is##ITEM(ARG1); \
    else \
      return ssIs##ITEM(S, ARG1); \
  }
//@}

/// @name Implementation of mxArray methods.
///       Always use MEX mxArrays if MEX is enabled.
//@{
#define HXI_MX_CREATE1(WHAT, TYPE1, ARG1) \
  HXI_EXTERN mxArray* hmxCreate##WHAT(TYPE1 ARG1) { \
    return mxCreate##WHAT(ARG1); \
  }
#define HXI_MX_CREATE3(WHAT, TYPE1, ARG1, TYPE2, ARG2, TYPE3, ARG3) \
  HXI_EXTERN mxArray* hmxCreate##WHAT(TYPE1 ARG1, TYPE2 ARG2, TYPE3 ARG3) { \
    return mxCreate##WHAT(ARG1, ARG2, ARG3); \
  }
#define HXI_MX_DESTROY(WHAT) \
  HXI_EXTERN void hmxDestroy##WHAT(mxArray *a) { \
    mxDestroy##WHAT(a); \
  }
#define HXI_MX_SETGET1(ITEM, TYPE, ARG) \
  HXI_EXTERN void hmxSet##ITEM(mxArray *a, TYPE ARG) { \
    mxSet##ITEM(a, ARG); \
  } \
  HXI_EXTERN TYPE hmxGet##ITEM(const mxArray *a) { \
    return mxGet##ITEM(a); \
  }
#define HXI_MX_NOSET1(ITEM, TYPE, ARG)
#define HXI_MX_GET1(ITEM, TYPE) \
  HXI_EXTERN TYPE hmxGet##ITEM(const mxArray *a) { \
    return mxGet##ITEM(a); \
  }
#define HXI_MX_IS1(ITEM) \
  HXI_EXTERN int_T hmxIs##ITEM(const mxArray *a) { \
    return mxIs##ITEM(a); \
  }
#define HXI_MX_TO1(ITEM, TYPE) \
  HXI_EXTERN TYPE hmxArrayTo##ITEM(const mxArray *a) { \
    return mxArrayTo##ITEM(a); \
  }
//@}

#else // defined(HXI_WITH_MEX)
/// @name Implementation of S-function methods accessing a Hxi::SimStruct.
//@{
#define HXI_SS_NOSET0(ITEM)
#define HXI_SS_SETGET1(ITEM, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE ARG) { \
    return ((Hxi::SimStruct<real_T>*)S)->set##ITEM(ARG); \
  } \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S) { \
    return ((Hxi::SimStruct<real_T>*)S)->get##ITEM(); \
  }
#define HXI_SS_SET1(ITEM, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE ARG) { \
    return ((Hxi::SimStruct<real_T>*)S)->set##ITEM(ARG); \
  }
#define HXI_SS_NOSET1(ITEM, TYPE, ARG)
#define HXI_SS_GET1(ITEM, TYPE) \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S) { \
    return ((Hxi::SimStruct<real_T>*)S)->get##ITEM(); \
  }
#define HXI_SS_IS1(ITEM) \
  HXI_EXTERN int_T hssIs##ITEM(SimStruct *S) { \
    return ((Hxi::SimStruct<real_T>*)S)->is##ITEM(); \
  }

#define HXI_SS_SETGET2(ITEM, TYPE1, ARG1, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE1 ARG1, TYPE ARG) { \
    return ((Hxi::SimStruct<real_T>*)S)->set##ITEM(ARG1, ARG); \
  } \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S, TYPE1 ARG1) { \
    return ((Hxi::SimStruct<real_T>*)S)->get##ITEM(ARG1); \
  }
#define HXI_SS_NOSETGET2(ITEM, TYPE1, ARG1, TYPE, ARG)
#define HXI_SS_SET2(ITEM, TYPE1, ARG1, TYPE, ARG) \
  HXI_EXTERN TYPE hssSet##ITEM(SimStruct *S, TYPE1 ARG1, TYPE ARG) { \
    return ((Hxi::SimStruct<real_T>*)S)->set##ITEM(ARG1, ARG); \
  }
#define HXI_SS_NOSET2(ITEM, TYPE1, ARG1, TYPE, ARG) 
#define HXI_SS_GET2(ITEM, TYPE1, ARG1, TYPE) \
  HXI_EXTERN TYPE hssGet##ITEM(SimStruct *S, TYPE1 ARG1) { \
    return ((Hxi::SimStruct<real_T>*)S)->get##ITEM(ARG1); \
  }
#define HXI_SS_IS2(ITEM, TYPE1, ARG1) \
  HXI_EXTERN int_T hssIs##ITEM(SimStruct *S, TYPE1 ARG1) { \
    return ((Hxi::SimStruct<real_T>*)S)->is##ITEM(ARG1); \
  }
//@}

/// @name Implementation of mxArray methods for Hxi::mxArray only.
//@{
#define HXI_MX_CREATE1(WHAT, TYPE1, ARG1) \
  HXI_EXTERN mxArray* hmxCreate##WHAT(TYPE1 ARG1) { \
    return (mxArray*)(new Hxi::mxArray<real_T>(ARG1)); \
  }
#define HXI_MX_CREATE3(WHAT, TYPE1, ARG1, TYPE2, ARG2, TYPE3, ARG3) \
  HXI_EXTERN mxArray* hmxCreate##WHAT(TYPE1 ARG1, TYPE2 ARG2, TYPE3 ARG3) { \
    return (mxArray*)(new Hxi::mxArray<real_T>(ARG1, ARG2, ARG3)); \
  }
#define HXI_MX_DESTROY(WHAT) \
  HXI_EXTERN void hmxDestroy##WHAT(mxArray *a) { \
    delete (Hxi::mxArray<real_T>*)a; \
  }
#define HXI_MX_SETGET1(ITEM, TYPE, ARG) \
  HXI_EXTERN void hmxSet##ITEM(mxArray *a, TYPE ARG) { \
    ((Hxi::mxArray<real_T>*)a)->set##ITEM(ARG); \
  } \
  HXI_EXTERN TYPE hmxGet##ITEM(const mxArray *a) { \
    return ((Hxi::mxArray<real_T>*)a)->get##ITEM(); \
  }
#define HXI_MX_NOSET1(ITEM, TYPE, ARG)
#define HXI_MX_GET1(ITEM, TYPE) \
  HXI_EXTERN TYPE hmxGet##ITEM(const mxArray *a) { \
    return ((Hxi::mxArray<real_T>*)a)->get##ITEM(); \
  }
#define HXI_MX_IS1(ITEM) \
  HXI_EXTERN int_T hmxIs##ITEM(const mxArray *a) { \
    return ((Hxi::mxArray<real_T>*)a)->is##ITEM(); \
  }
#define HXI_MX_TO1(ITEM, TYPE) \
  HXI_EXTERN TYPE hmxArrayTo##ITEM(const mxArray *a) { \
    return ((Hxi::mxArray<real_T>*)a)->to##ITEM(); \
  }
//@}
#endif // defined(HXI_WITH_MEX)

// actually define methods
#include "Hxi_SimStruct_methods.h"

static void defaultSFunctionMethod1(SimStruct *S)
{
  hssSetErrorStatus(S, "S-function method not initialized");
}

static void defaultSFunctionMethod2(SimStruct *S, int_T)
{
  hssSetErrorStatus(S, "S-function method not initialized");
}

HXI_EXTERN void Hxi_SimStruct_init(SimStruct *S); // see cg_sfun.h

/// @name Treatment of external binary S-function object.
/// The functions are defined with extern "C" linkage
/// to avoid MEX/Hxi SimStruct incompatibility. See Hxi_SFunction.C.
//@{
/// Open binary S-function object. 
extern "C" SimStruct *Hxi_SFunction_open(SimStruct *S);
/// Release binary S-function object.
extern "C" void Hxi_SFunction_close(SimStruct *S);
//@}

/** Construct a SimStruct. */
HXI_EXTERN SimStruct *Hxi_SimStruct_create(const char *path) {
  SimStruct *S; // SimStruct to be created
  Hxi::SimStruct<real_T> *S0 = new Hxi::SimStruct<real_T>(); // used during creation
  S0->setPath(path);
#if defined(HXI_INLINE_S_FUNCTION)
  // directly initialize SimStruct
  S = (SimStruct*)S0;
  Hxi_SimStruct_init(S);
#else
  S = Hxi_SFunction_open((SimStruct*)S0);
#if defined(HXI_WITH_MEX)
  if (!S) {
    // we are having a MEX S-function; don't need S0 anymore
    delete S0;
    S = Hxi_MEX_SimStruct_create();
    ssSetPath(S, path);
  }
#endif
#endif
  return S;
}

/** Delete a SimStruct. */
HXI_EXTERN void Hxi_SimStruct_destroy(SimStruct *S) {
#if !defined(HXI_INLINE_S_FUNCTION)
  // release S-function
  Hxi_SFunction_close(S);
#endif
#if defined(HXI_WITH_MEX)
  if (!isHxiSimStruct(S))
    Hxi_MEX_SimStruct_destroy(S);
  else
#endif
    delete (Hxi::SimStruct<real_T>*)S;
}

/** Initialize SimStruct for specific model. */
HXI_EXTERN void Hxi_mdlInitializeSizes(SimStruct *S) {
#if defined(HXI_WITH_MEX)
  if (!isHxiSimStruct(S)) {
    Hxi_MEX_SFunction_init(S);
  }
  else
#endif
    (*((Hxi::SimStruct<real_T>*)S)->getmdlInitializeSizes())(S);
}

/** Optional: Allocate local ressources for simulation. */
HXI_EXTERN void Hxi_mdlStart(SimStruct *S)
{
  assert(hssGetmdlStart(S) != NULL);
#if defined(HXI_WITH_MEX)
  if (!isHxiSimStruct(S)) {
    sfcnStart(S);
  }
  else
#endif
    (*((Hxi::SimStruct<real_T>*)S)->getmdlStart())(S);
}

/** Initialize sample times. */
HXI_EXTERN void Hxi_mdlInitializeSampleTimes(SimStruct *S)
{
#if defined(HXI_WITH_MEX)
  if (!isHxiSimStruct(S)) {
    sfcnInitializeSampleTimes(S);
  }
  else
#endif
    (*((Hxi::SimStruct<real_T>*)S)->getmdlInitializeSampleTimes())(S);
}

/** Optional: Compute initial conditions. */
HXI_EXTERN void Hxi_mdlInitializeConditions(SimStruct *S)
{
  assert(hssGetmdlInitializeConditions(S) != NULL);
#if defined(HXI_WITH_MEX)
  if (!isHxiSimStruct(S)) {
    if (ssGetVersion(S) == SIMSTRUCT_VERSION_LEVEL1)
      sfcnInitializeConditionsLevel1(ssGetX(S), S);
    else
      sfcnInitializeConditions(S);
  }
  else
#endif
    (*((Hxi::SimStruct<real_T>*)S)->getmdlInitializeConditions())(S);
}

/** Compute model outputs. */
HXI_EXTERN void Hxi_mdlOutputs(SimStruct *S, int_T tid)
{
#if defined(HXI_WITH_MEX)
  if (!isHxiSimStruct(S)) {
    if (ssGetVersion(S) == SIMSTRUCT_VERSION_LEVEL1)
      sfcnOutputsLevel1(ssGetY(S), ssGetX(S), ssGetU(S), S, tid);
    else
      sfcnOutputs(S, tid);
  }
  else
#endif
    (*((Hxi::SimStruct<real_T>*)S)->getmdlOutputs())(S, tid);
}

/** Optional: Update discrete-time states. */
HXI_EXTERN void Hxi_mdlUpdate(SimStruct *S, int_T tid)
{
  assert(hssGetmdlUpdate(S) != NULL);
#if defined(HXI_WITH_MEX)
  if (!isHxiSimStruct(S)) {
    if (ssGetVersion(S) == SIMSTRUCT_VERSION_LEVEL1)
      sfcnUpdateLevel1(ssGetX(S), ssGetU(S), S, tid);
    else
      sfcnUpdate(S, tid);
  }
  else
#endif
    (*((Hxi::SimStruct<real_T>*)S)->getmdlUpdate())(S, tid);
}

/** Optional: Compute derivatives for continuous-time states. */
HXI_EXTERN void Hxi_mdlDerivatives(SimStruct *S)
{
  assert(hssGetmdlDerivatives(S) != NULL);
#if defined(HXI_WITH_MEX)
  if (!isHxiSimStruct(S)) {
    if (ssGetVersion(S) == SIMSTRUCT_VERSION_LEVEL1)
      sfcnDerivativesLevel1(ssGetdX(S), ssGetX(S), ssGetU(S), S, 0);
    else
      sfcnDerivatives(S);
  }
  else
#endif
    (*((Hxi::SimStruct<real_T>*)S)->getmdlDerivatives())(S);
}

/** Optional: Compute Jacobian J = d(dxc,xd,y)/d(xc,xd,u). */
HXI_EXTERN void Hxi_mdlJacobian(SimStruct *S)
{
  assert(hssGetmdlJacobian(S) != NULL);
#if defined(HXI_WITH_MEX)
  if (!isHxiSimStruct(S)) {
    sfcnJacobian(S);
  }
  else
#endif
    (*((Hxi::SimStruct<real_T>*)S)->getmdlJacobian())(S);
}

/** Release resources allocated for simulation. */
HXI_EXTERN void Hxi_mdlTerminate(SimStruct *S)
{
#if defined(HXI_WITH_MEX)
  if (!isHxiSimStruct(S)) {
    sfcnTerminate(S);
  }
  else
#endif
    (*((Hxi::SimStruct<real_T>*)S)->getmdlTerminate())(S);
}


//===================================================================
