/**
 * @file Prg_SFunction.h
 *    Basic functionality for formulating an optimization problem for
 *    a model given as MEX S-function.
 *
 * rf, 7/25/00
 */

/*
    Copyright (C) 1997--2010  Ruediger Franke

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

#ifndef Prg_SFunction_H
#define Prg_SFunction_H

#include <Omu_Program.h>

#include <If_Command.h>

#include "Hxi_SFunction.h"

/**
 * Basic functionality for formulating an optimization problem for
 * a model given as MEX S-function.
 */
class Prg_SFunction: public Omu_Program {

 protected:
  char 		*_mdl_name;	///< S-function name
  char 		*_mdl_path;	///< full S-function path
  char 		*_mdl_args;	///< S-function parameters
  SimStruct 	*_SS;		///< pointer to %SimStruct
  mxArray	**_mx_args; 	///< S-function parameters after parsing
  int 		_mdl_nargs; 	///< number of S-function arguments

  double 	_t0_setup_model;///< time used for initialization of model
  int		_mdl_np;	///< number of model parameters
  int		_mdl_nd;	///< number of discrete-time states
  int		_mdl_nx;	///< number of model states (incl. discrete)
  int		_mdl_nu;	///< number of model inputs
  int		_mdl_ny;	///< number of model outputs

  VECP		_mdl_p;		///< parameters
  VECP		_mdl_x0;	///< initial states

  /// indicate that setup_model needs to be called
  /// as _mdl_name, _mdl_path, or _mdl_args changed
  bool 		_mdl_needs_setup;

  // methods
  virtual void setup_model(); 	///< load S-function

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
  bool setContinuousTask(bool val);
  /// enable hit for all discrete sample times;
  /// return true if a discrete sample time exists and has been enabled
  bool setSampleHit(bool val);
  //@}

 public:

  Prg_SFunction();		///< constructor
  ~Prg_SFunction();		///< destructor

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

  /** parameters */
  const VECP mdl_p() const {return _mdl_p;}
  /** set parameters */
  void set_mdl_p(const VECP value) {v_copy_elements(value, _mdl_p);}

  /** initial states */
  const VECP mdl_x0() const {return _mdl_x0;}
  /** set initial states */
  void set_mdl_x0(const VECP value) {v_copy_elements(value, _mdl_x0);}

  //@}
};  

#endif
