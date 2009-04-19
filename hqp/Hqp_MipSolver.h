/**
 * @file Hqp_MipSolver.h 
 *   Interface for a solver for mixed integer problems.
 *
 * rf, 4/11/09
 *
 */

/*
    Copyright (C) 1994--2009  Ruediger Franke

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

#ifndef Hqp_MipSolver_H
#define Hqp_MipSolver_H

#include <If_List.h>
#include <If_Command.h>
#include <If_Class.h>

#include "Hqp_impl.h"

class Hqp_SqpProgram;

IF_BASE_DECLARE(Hqp_MipSolver);

/**
 *  Base class for a solver for mixed integer problems.
 */
class Hqp_MipSolver {

 protected:
  If_List	_ifList; 	///< interface elements

  Hqp_SqpProgram *_prg;		///< optimization problem to work on

 public:
  Hqp_MipSolver(); 		///< constructor
  virtual ~Hqp_MipSolver(); 	///< destructor

  Hqp_SqpProgram *prg() {return _prg;}
  void	set_prg(Hqp_SqpProgram *);

  /** initialize a MIP for an optimization problem given in _prg */
  virtual int init(IF_DEF_ARGS) = 0;

  /** initialize and solve a MIP for an optimization problem given in _prg */
  virtual int solve(IF_DEF_ARGS) = 0;

  virtual char *name() = 0; 	///< name of mixed integer solver
};  


#endif
