/*
 * Hqp_Solver.C -- 
 *   - class definition
 *
 * rf, 7/19/94
 *
 * rf, 8/13/98
 *   - make Hqp_Solver a visible interface class
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

#include <If_Int.h>
#include <If_Real.h>
#include <If_Method.h>

#include "Hqp_Solver.h"
#include "Hqp_Program.h"

IF_BASE_DEFINE(Hqp_Solver);

typedef If_Method<Hqp_Solver> If_Cmd;

//--------------------------------------------------------------------------
Hqp_Solver::Hqp_Solver()
{
  _qp = NULL;
  _result = Hqp_Infeasible;
  _iter = 0;
  _max_iters = 200;
  _eps = 1e-10;
  _y = VNULL;
  _z = VNULL;

  _ifList.append(new If_Int("qp_iter", &_iter));
  _ifList.append(new If_Int("qp_max_iters", &_max_iters));
  _ifList.append(new If_Real("qp_eps", &_eps));
  _ifList.append(new If_Cmd("qp_result", &Hqp_Solver::result_str, this));
}

//--------------------------------------------------------------------------
Hqp_Solver::~Hqp_Solver()
{
  v_free(_y);
  v_free(_z);
}

//--------------------------------------------------------------------------
int Hqp_Solver::result_str(int, const char *[], const char **result)
{
  *result = hqp_result_strings[_result];
  return IF_OK;
}

//--------------------------------------------------------------------------
void Hqp_Solver::qp(Hqp_Program *qp)
{
  _qp = qp;
}
