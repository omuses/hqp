/*
 * Hqp_ASCEND_dummies.C -- 
 *   dummy functions if ASCEND is not linked
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

#include <stdio.h>
#include <stdlib.h>

#define ASCEND_DUMMY(name, namestr) \
  extern "C" void name() \
  { \
    fprintf(stderr, \
            "Fatal: ASCEND function %s called but not linked!\n", \
            namestr); \
    exit(-1); \
  }

ASCEND_DUMMY(slv_register_client, "slv_register_client");
ASCEND_DUMMY(slv_get_num_solvers_vars, "slv_get_num_solvers_vars");
ASCEND_DUMMY(slv_get_solvers_var_list, "slv_get_solvers_var_list");
ASCEND_DUMMY(slv_get_num_solvers_rels, "slv_get_num_solvers_rels");
ASCEND_DUMMY(slv_get_solvers_rel_list, "slv_get_solvers_rel_list");
ASCEND_DUMMY(slv_get_obj_relation, "slv_get_obj_relation");
ASCEND_DUMMY(var_fixed, "var_fixed");
ASCEND_DUMMY(var_lower_bound, "var_lower_bound");
ASCEND_DUMMY(var_upper_bound, "var_upper_bound");
ASCEND_DUMMY(var_value, "var_value");
ASCEND_DUMMY(var_set_value, "var_set_value");
ASCEND_DUMMY(var_sindexF, "var_sindexF");
ASCEND_DUMMY(rel_n_incidencesF, "rel_n_incidencesF");
ASCEND_DUMMY(rel_incidence_list, "rel_incidence_list");
ASCEND_DUMMY(rel_equal, "rel_equal");
ASCEND_DUMMY(rel_less, "rel_less");
ASCEND_DUMMY(relman_eval, "relman_eval");
ASCEND_DUMMY(relman_diff_grad, "relman_diff_grad");
