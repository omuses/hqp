/*
 * Hqp_Client.h -- 
 *   - quadratic solver for communicating with a "compute server"
 *
 * rf, 8/12/94
 *
 * rf, 8/13/98
 *   - make Hqp_Client an exchangeable interface class
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

#ifndef Hqp_Client_H
#define Hqp_Client_H

#include "Hqp_Solver.h"

class Hqp_Client: public Hqp_Solver {

 public:

  Hqp_Client() {}
  ~Hqp_Client() {}

  // initializing a program (updating for SQP integration)
  // is empty as only complete tasks are processed that time
  int	init(IF_CMD_ARGS) {return IF_OK;}
  int	update(IF_CMD_ARGS) {return IF_OK;}

  // solving a program (hot start for SQP integration)
  int	cold_start(IF_CMD_ARGS) {return IF_OK;}
  int	hot_start(IF_CMD_ARGS) {return IF_OK;}
  int	step(IF_DEF_ARGS);
  int 	solve(IF_DEF_ARGS);

  const char *name() {return "Client";}
};  

#endif


