/*
 * Hqp_impl.C -- definitions
 *
 * rf, 5/29/94
 *
 * rf, 2/12/97
 *   modified result strings
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

#include "Hqp_impl.h"

const char *hqp_result_strings[] = {
  "optimal",
  "feasible",
  "infeasible",
  "suboptimal",
  "degenerate"
};


//=========================================================================
