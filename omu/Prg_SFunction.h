/**
 * @file Prg_SFunction.h
 *    Basic functionality for formulating an optimization problem for
 *    a model given as MEX S-function.
 *
 * rf, 7/25/00
 */

/*
    Copyright (C) 1997--2003  Ruediger Franke

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

/** Select MEX interface to S-function */
#define MATLAB_MEX_FILE 1
#include <simstruc.h>

/**
 * Basic functionality for formulating an optimization problem for
 * a model given as MEX S-function.
 */
class Prg_SFunction: public Omu_Program {

 protected:
  char 		*_mdl_name;	///< S-function name
  char 		*_mdl_path;	///< full S-function path
  char 		*_mdl_args;	///< S-function parameters
  SimStruct 	*_S;		///< pointer to %SimStruct
  mxArray	**_mx_args; 	///< S-function parameters after parsing
  int 		_mdl_nargs; 	///< number of S-function arguments

  int		_mdl_nx;	///< number of model states
  int		_mdl_nu;	///< number of model inputs
  int		_mdl_ny;	///< number of model outputs

  VECP		_mdl_x0;	///< initial states

  /// indicate that setup_model needs to be called
  /// as _mdl_name, _mdl_path, or _mdl_args changed
  bool 		_mdl_needs_setup;

  // methods
  virtual void setup_model(); 	///< load S-function

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

  /** initial states */
  const VECP mdl_x0() const {return _mdl_x0;}
  void set_mdl_x0(const VECP value);	///< set initial states

  //@}
};  

#endif
