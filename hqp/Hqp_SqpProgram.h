/*
 * Hqp_SqpProgram.h -- 
 *   - base class for formulating nonlinear programs for solving with 
 *     sequential quadratic programming
 *   - constructor has to allocate _x and sparsity structure in _qp
 *   - additional methods provide value updates for SQP iterations
 *
 * rf, 6/6/94
 *
 * rf, 1/25/97
 *   - new interface command prg_update_fbd that calls update_fbd
 */

/*
    Copyright (C) 1994--1998  Ruediger Franke

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

#ifndef Hqp_SqpProgram_H
#define Hqp_SqpProgram_H

#include <If_List.h>
#include <If_Command.h>
#include <If_Class.h>

#include "Hqp_impl.h"

class Hqp_Program;

IF_BASE_DECLARE(Hqp_SqpProgram);

//--------------------
class Hqp_SqpProgram {

 protected:

  If_List	_ifList;

  Hqp_Program	*_qp;		// quadratic approximation
  VECP		_x;		// vector of variables
  Real		_f;		// objective value

 public:

  Hqp_SqpProgram();
  virtual ~Hqp_SqpProgram();

  int	setup_cmd(IF_DEF_ARGS);
  int	init_x_cmd(IF_DEF_ARGS);
  int	test_cmd(IF_DEF_ARGS);
  int	qp_dump_cmd(IF_CMD_ARGS);
  int	update_fbd_cmd(IF_CMD_ARGS);

  // methods to define a program
  //----------------------------
  virtual void	setup()=0;	// realloc _qp, set up sparsity pattern
  virtual void	init_x()=0;	// set start variables

  virtual void update_fbd()=0;	// update values for _qp according to _x
  virtual void update(const VECP y, const VECP z)=0;	// update everything

  virtual void reinit_bd();	// re-initialize constants in constraints
  
  // member access
  // --> x(const VECP) may be overloaded for recomputing cached data
  //----------------------------------------------------------------
  virtual Hqp_Program *qp() {return _qp;}

  virtual const VECP x() {return _x;}
  virtual void x(const VECP);
  
  virtual Real f() {return _f;}
  virtual void f(Real f) {_f = f;}	// needed by Hqp_HL::init()

  virtual char *name()=0;
};  

#endif


