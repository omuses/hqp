/*
 * Hqp_MipSolver.C -- class definition
 *
 * rf, 4/11/09
 *
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

#include "Hqp_MipSolver.h"
#include "Hqp_SqpProgram.h"

#include <If_Method.h>
#include <If_Module.h>
#include <If_String.h>

#define GET_CB(vartype, name) \
  #name, \
  IF_GET_CB(vartype, Hqp_MipSolver, name)

typedef If_Method<Hqp_MipSolver> If_Cmd;

IF_BASE_DEFINE(Hqp_MipSolver);

extern Hqp_SqpProgram *theSqpProgram;

//--------------------------------------------------------------------------
Hqp_MipSolver::Hqp_MipSolver()
{
  _prg = theSqpProgram;
  _result = Hqp_Infeasible;

  //_ifList.append(new If_Cmd("mip_init", &Hqp_MipSolver::init, this));
  _ifList.append(new If_Cmd("mip_solve", &Hqp_MipSolver::solve, this));
  _ifList.append(new If_String(GET_CB(const char *, mip_result)));
}

//--------------------------------------------------------------------------
Hqp_MipSolver::~Hqp_MipSolver()
{
}

//--------------------------------------------------------------------------
void Hqp_MipSolver::set_prg(Hqp_SqpProgram *prg)
{
  _prg = prg;
}


//==========================================================================
