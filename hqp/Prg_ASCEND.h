/*
 * Prg_ASCEND.h -- 
 *   - convert an ASCEND system into a Hqp_SqpProgram
 *
 * rf, 12/17/99
 */

/*
    Copyright (C) 1999  Ruediger Franke

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

#ifndef Prg_ASCEND_H
#define Prg_ASCEND_H

#include "Hqp_SqpProgram.h"

/*
 * Includes required for ASCEND
 */

extern "C" {
#include "utilities/ascConfig.h"
#include "utilities/ascMalloc.h"
#include "utilities/set.h"
#include "general/time.h"
#include "utilities/mem.h"
#include "general/list.h"
#include "compiler/fractions.h"
#include "compiler/dimen.h"
#include "compiler/functype.h"
#include "compiler/func.h"
#include "solver/mtx.h"
#include "solver/linsol.h"
#include "solver/linsolqr.h"
#include "solver/slv_types.h"
#include "solver/var.h"
#include "solver/rel.h"
#include "solver/discrete.h"
#include "solver/conditional.h"
#include "solver/logrel.h"
#include "solver/bnd.h"
#include "solver/calc.h"
#include "solver/relman.h"
#include "solver/slv_common.h"
#include "solver/slv_client.h"
};

/*
 * class declaration
 */

class Prg_ASCEND: public Hqp_SqpProgram {

 protected:
  slv_system_t _slv_system;

  int _nvars;			// number of ASCEND variables
  int _nrels;			// number of ASCEND relations
  struct var_variable **_vars;	// ASCEND variables
  struct rel_relation **_rels;	// ASCEND relations
  struct rel_relation *_obj;	// ASCEND objective relation
  int _safe_calc;		// apply ASCEND's safe calculation

  VECP  _var_lb;	// lower bounds for ASCEND variables
  VECP  _var_ub;	// upper bounds for ASCEND variables
  IVECP _var_asc2hqp;	// map ASCEND var index to HQP

  VECP _derivatives;		// derivatives of ASCEND relation
  IVECP _var_master_idxs;	// variable indices for derivatives
  IVECP _var_solver_idxs;	// variable indices for derivatives

  double _Inf;		// infinite value to check for existence of bounds
  int _me_bounds;	// number of linear equality constraints for HQP
  int _m_bounds;	// number of linear inequality constraints for HQP

  slv_status_t _slv_status;

  void update_bounds();

 public:

  Prg_ASCEND();
  ~Prg_ASCEND();

  /*
   * Client functions called by ASCEND
   */

  static int slv_register(SlvFunctionsT *);
  static SlvClientToken slv_create(slv_system_t, int *);
  static int slv_destroy(slv_system_t, SlvClientToken);
  static void slv_get_status(slv_system_t, SlvClientToken,
			     slv_status_t *);
  static void slv_presolve(slv_system_t, SlvClientToken);
  static void slv_iterate(slv_system_t, SlvClientToken);
  static void slv_solve(slv_system_t, SlvClientToken);

  /*
   * implementation of interface defined by Hqp_SqpProgram
   */
     
  void	setup();
  void	init_x();

  void	update_fbd();
  void	update(const VECP y, const VECP z);

  char *name() {return "ASCEND";}
};  


#endif


