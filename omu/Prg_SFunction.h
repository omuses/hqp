/**
 * @file Prg_SFunction.h
 *    Basic functionality for formulating an optimization problem for
 *    a model given as MEX S-function.
 *
 * rf, 7/25/00
 */

/*
    Copyright (C) 1997--2002  Ruediger Franke

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

#define MATLAB_MEX_FILE 1
#include <simstruc.h>

/**
 * Basic functionality for formulating an optimization problem for
 * a model given as MEX S-function.
 */
class Prg_SFunction: public Omu_Program {

 protected:
  char 		*_mdl_name;	///< S-function name
  char 		*_mdl_args;	///< S-function parameters
  SimStruct 	*_S;		///< pointer to %SimStruct
  mxArray	*_mx_args;///< S-function parameters after parsing (cell array)

  int		_mdl_nx;	///< number of model states
  int		_mdl_nu;	///< number of model inputs
  int		_mdl_ny;	///< number of model outputs

  VECP		_mdl_x0;	///< initial states

  // methods
  virtual void setup_sfun(); 	///< load S-function

 public:

  Prg_SFunction();		///< constructor
  ~Prg_SFunction();		///< destructor

  int	mdl_name(IF_DEF_ARGS);	///< access S-function name
  int	mdl_args(IF_DEF_ARGS);	///< access S-function parameters as string
};  

#endif
