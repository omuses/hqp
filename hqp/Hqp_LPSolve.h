/**
 * @file Hqp_LPSolve.h -- 
 *    Interface to lp_solve for the solution of linear mixed integer programs
 *     (see: http://lpsolve.sourceforge.net/5.5).
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

#ifndef Hqp_LPSolve_H
#define Hqp_LPSolve_H

#include "Hqp_MipSolver.h"

struct _lprec;
typedef struct _lprec lprec;

/**
 *  Interface to lp_solve (see: http://lpsolve.sourceforge.net/5.5).
 */
class Hqp_LPSolve : public Hqp_MipSolver {

 protected:

  lprec *_lp; 			///< the linear program for lp_solve
  double _gap;			///< allowed tolerance for branch&bound
  int 	 _timeout; 		///< timeout for solver
  int 	 _logging; 		///< verbose level
  char *_dump_format; 		///< format for dump files

 public:
  Hqp_LPSolve(); 		///< constructor
  virtual ~Hqp_LPSolve(); 	///< destructor

  virtual int init(IF_DEF_ARGS);
  virtual int solve(IF_DEF_ARGS);

  /** Write current linear program to the file mip_dump.lp */
  virtual void dump();

  /**
   * @name Member access methods
   */
  //@{

  /// absolute tolerance for branch&bound
  double gap() const {return _gap;}
  /// set absolute tolerance for branch&bound
  void set_gap(double value) {_gap = value;}

  /// timeout for solver in seconds
  int timeout() const {return _timeout;}
  /// set timeout
  void set_timeout(int value) {_timeout = value;}

  /// verbose level for log messages
  /// (0: neutral, 1: critical, 2: severe, 3: important, 4: normal, 5: detailed, 6: full)
  int logging() const {return _logging;}
  /// set logging
  void set_logging(int value) {_logging = value;}

  /// format for dump of lp
  const char *dump_format() const {return _dump_format;}
  /// set dump format ("lp" or "mps")
  void set_dump_format(const char *value);

  //@}

  virtual char *name() {return "LPSolve";} ///< solver name
};  


#endif
