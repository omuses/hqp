/**
 * @file Hqp_SqpProgram.h
 *   Interface for a nonlinear program
 *
 * rf, 6/6/94
 */

/*
    Copyright (C) 1994--2014  Ruediger Franke

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
#include <If_Class.h>

#include "Hqp_impl.h"
#include "Hqp_Program.h"

IF_BASE_DECLARE(Hqp_SqpProgram);

/**
 *  Base class for formulating nonlinear programs.
 */
class Hqp_SqpProgram {

 protected:
  If_List	_ifList; 	///< interface elements

  Hqp_Program	*_qp;		///< quadratic approximation
  VECP		_x;		///< vector of optimization variables
  Real		_f;		///< objective function value

 public:
  Hqp_SqpProgram(); 		///< constructor
  virtual ~Hqp_SqpProgram(); 	///< destructor

  /** Substract finite difference approximation for gradients from current
      linear approximation in qp and return maximum deviation. The finite
      difference approximation is done at the current solution x and the
      obtained deviations are stored in _qp->c, _qp->A, and _qp->C. */
  virtual Real	test();

  /** Write current linear quadratic appoximation in _qp to the file 
      prg_qp_dump.out. */
  virtual void	qp_dump();

  /** @name Methods to define a program */
  //@{

  /** (Re)-allocate _x and _qp and set up sparsity patterns in _qp. */
  virtual void	setup() = 0;

  /** Set start variables for optimization in _x. */
  virtual void	init_x() = 0;

  /** Update values of objective and constraints in _f, _qp->b, _qp->d
      for current _x. */
  virtual void 	update_fbd() = 0;

  /** Update everything, including values and derivatives, in _f and _qp
      for current _x. */
  virtual void 	update(const VECP y, const VECP z) = 0;

  /** Re-initialize constants in constraints (_qp->b, _qp->d)
      for current _x. */
  virtual void 	reinit_bd();
  //@}
  
  /** @name Member access methods */
  //@{
  /** current linear quadratic approximation */
  virtual Hqp_Program *qp() {return _qp;}

  /** current optimization step (solution of linear-quadratic sub-problem)*/
  virtual const VECP s() const {return _qp->x;}
  /** set optimization step */
  virtual void set_s(const VECP v) {v_copy_elements(v, _qp->x);}

  /** current vector of optimization variables */
  virtual const VECP x() const {return _x;}
  virtual void set_x(const VECP); ///< set vector of optimization variables

  /** current value of objective function */
  virtual Real f() const {return _f;}
  virtual void set_f(Real f) {_f = f;} ///< set objective function value

  /** current violation of constraints */
  virtual Real norm_inf() const;

  //@}

  virtual const char *name()=0; ///< program name
};  

#endif


