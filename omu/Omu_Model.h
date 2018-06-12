/**
 * @file Omu_Model.h
 *    Basic functionality for formulating an optimization problem for
 *    a model given as S-function or Functional Model Unit (FMU).
 *
 * rf, 7/25/00
 */

/*
    Copyright (C) 1997--2017  Ruediger Franke

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

#ifndef Omu_Model_H
#define Omu_Model_H

#include <If_List.h>
#include <Meschach.h>
#include <Hxi_SFunction.h>

#include "Omu.h"

/** Error kind used for failed S-function methods */
#define SMETHOD_ERROR E_FORMAT
/** Log the call to an S-function method. */
#define SMETHOD_LOG_CALL(method, S) { \
  if (_mdl_logging >= If_LogInfo) \
    If_Log("Info", "%s at %.3f", #method, ssGetT(S)); \
  ssSetErrorStatus(S, NULL); \
}
/** Log and throw an error from an S-function method. */
#define SMETHOD_LOG_ERROR(method, S, errnum) { \
  if (ssGetErrorStatus(S)) { \
    fprintf(stderr, "Error from " #method ": %s\n", ssGetErrorStatus(S)); \
    if (errnum >= 0) \
      m_error(errnum, ssGetErrorStatus(S)); \
  } \
}
/** Call an S-function method, including logging and error handling. */
#define SMETHOD_CALL(method, S) { \
  SMETHOD_LOG_CALL(method, S) \
  method(S); \
  SMETHOD_LOG_ERROR(method, S, SMETHOD_ERROR) \
}
/** Call an S-function method, including logging. */
#define SMETHOD_CALL_HOLD(method, S) { \
  SMETHOD_LOG_CALL(method, S) \
  method(S); \
  SMETHOD_LOG_ERROR(method, S, -1) \
}
/** Call an S-function method, including logging and error handling. */
#define SMETHOD_CALL2(method, S, tid) { \
  SMETHOD_LOG_CALL(method, S) \
  method(S, tid); \
  SMETHOD_LOG_ERROR(method, S, SMETHOD_ERROR) \
}

/**
 * Basic functionality for formulating an optimization problem for
 * a model given as S-function or Functional Model Unit (FMU).
 */
class OMU_API Omu_Model {

 private:
  If_List	_ifList_model; 	///< container for interface elements

 protected:
  int		_mdl_ncpu; 	///< number of cpus
  char 		*_mdl_name;	///< S-function or FMU name
  char 		*_mdl_path;	///< S-function or FMU path for loading
  bool 		_mdl_is_fmu; 	///< indicate a Function Model Unit
  char 		*_mdl_args;	///< S-function parameters
  SimStruct 	**_SS;		///< pointers to %SimStruct
  mxArray	**_mx_args; 	///< S-function parameters after parsing
  int 		_mdl_nargs; 	///< number of S-function arguments
  If_LogLevel 	_mdl_logging;	///< log level for debugging
  bool 		_mdl_jac; 	///< use Jacobian provided by model

  double 	_t0_setup_model;///< time used for initialization of model
  int		_mdl_np_total;	///< number of model parameters (incl. strings)
  int		_mdl_np;	///< number of numeric model parameters
  int		_mdl_nd;	///< number of discrete-time states
  int		_mdl_nx;	///< number of model states (incl. discrete)
  int		_mdl_nu;	///< number of model inputs
  int		_mdl_ny;	///< number of model outputs

  VECP		_mdl_p;		///< parameters
  VECP		_mdl_x_start;	///< start values of states
  VECP		_mdl_u_start;	///< start values of inputs

  VECP 		_mdl_x_nominal;	///< nominal states (for scaling)
  VECP 		_mdl_u_nominal;	///< nominal inputs (for scaling)
  VECP 		_mdl_y_nominal;	///< nominal outputs (for scaling)

  /// indicate that setup_model needs to be called
  /// as _mdl_name, _mdl_path, or _mdl_args changed
  bool 		_mdl_needs_setup;
  IVECP 	_mdl_needs_init; ///< indicate that initialization is needed

  IVECP 	_mdl_jac_y_active; // mark outputs for Jacobian
  IVECP 	_mdl_jac_u_active; // mark inputs for Jacobian
 private:
  /// model Jacobian in row major
  IVECP 	_mdl_jac_jc;   // col indices
  IVECP 	_mdl_jac_ir;   // start index for each row in jc
  IVECP 	_mdl_jac_deps; // help variable for dependencies
  const IVECP 	mdl_jac_deps() const {return _mdl_jac_deps;}
  void 		set_mdl_jac_deps(const IVECP v) {
    _mdl_jac_deps = iv_copy(v, _mdl_jac_deps);
  }

 protected:
  // methods
  virtual void setup_model(double t0, int ncpu = 1); 	///< load S-function

  void setup_jac();
  void setup_jac_row(const char *category, int idx,
                     int jc_offset, int ir_offset);

  /**
   * @name Helper methods for reading and writing S-function arguments.
   */
  //@{
  void read_mx_args(VECP p); ///< read _mx_args into vector p
  void write_mx_args(VECP p); ///< write vector p to _mx_args
  //@}

  /**
   * @name Helper methods for treatment of sample times
   */
  //@{
  /// enable continuous sample time;
  /// return true if a continuous sample time exists and has been enabled
  bool setContinuousTask(SimStruct *S, bool val);
  /// enable hit for all discrete sample times;
  /// return true if a discrete sample time exists and has been enabled
  bool setSampleHit(SimStruct *S, bool val);
  /// set all discrete sample times
  void setSampleTime(SimStruct *S, double val);
  //@}

 public:

  Omu_Model(int ncpu = 1); ///< constructor
  ~Omu_Model();            ///< destructor

  /**
   * @name Member access methods (no If prefix)
   */
  //@{

  /** S-function name */
  const char *mdl_name() const {return _mdl_name;}
  void set_mdl_name(const char *str);	///< set S-function name

  /** S-function path, including name, used for dynamic loading. */
  const char *mdl_path() const {return _mdl_path;}
  void set_mdl_path(const char *str);	///< set S-function path

  /** String representation of S-function arguments */
  const char *mdl_args() const {return _mdl_args;}
  void set_mdl_args(const char *str);	///< set S-function arguments

  /** number of cpus */
  int mdl_ncpu() const {return _mdl_ncpu;}

  /** Log level for debugging */
  int mdl_logging() const {return _mdl_logging;}
  /** set log level */
  void set_mdl_logging(int logging) {
    if (logging < If_LogNone)
      _mdl_logging = If_LogNone;
    else if (logging > If_LogAll)
      _mdl_logging = If_LogAll;
    else
      _mdl_logging = (If_LogLevel)logging;
  }

  /** flag about use of model provided Jacobian */
  bool mdl_jac() const {return _mdl_jac;}
  /** set flag about use of model provided Jacobian */
  void set_mdl_jac(bool val) {_mdl_jac = val;}

  /** parameters */
  const VECP mdl_p() const {return _mdl_p;}
  /** set parameters */
  void set_mdl_p(const VECP value) {v_copy_elements(value, _mdl_p);}

  /** initial states */
  const VECP mdl_x_start() const {return _mdl_x_start;}
  /** set initial states */
  void set_mdl_x_start(const VECP value) {v_copy_elements(value, _mdl_x_start);}

  /** inputs */
  const VECP mdl_u_start() const {return _mdl_u_start;}
  /** set inputs */
  void set_mdl_u_start(const VECP value) {v_copy_elements(value, _mdl_u_start);}

  /** nominal state values (for scaling) */
  const VECP mdl_x_nominal() const {return _mdl_x_nominal;}
  /** set nominal states */
  void set_mdl_x_nominal(const VECP v) {v_copy_elements(v, _mdl_x_nominal);}

  /** nominal input values (for scaling) */
  const VECP mdl_u_nominal() const {return _mdl_u_nominal;}
  /** set nominal inputs */
  void set_mdl_u_nominal(const VECP v) {v_copy_elements(v, _mdl_u_nominal);}

  /** nominal output values (for scaling) */
  const VECP mdl_y_nominal() const {return _mdl_y_nominal;}
  /** set nominal outputs */
  void set_mdl_y_nominal(const VECP v) {v_copy_elements(v, _mdl_y_nominal);}
  //@}
};  

#endif
