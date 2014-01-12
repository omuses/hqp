/*
 * Hqp_Client.C -- class definition
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

#include <stdlib.h>
#include <assert.h>

#include "Hqp_Client.h"
#include "Hqp_Program.h"

IF_CLASS_DEFINE("Client", Hqp_Client, Hqp_Solver);

//--------------------------------------------------------------------------
int Hqp_Client::step(int, const char *[], const char **result)
{
  *result = "Hqp_Client::step not implemented";
  return IF_ERROR;
}

//--------------------------------------------------------------------------
int Hqp_Client::solve(IF_CMD_ARGS)
{
  FILE *fp1, *fp2;

  if ((fp1 = fopen("/tmp/comm/pipe1", "w")) == NULL)
    perror("opening pipe1");
  sp_foutput(fp1, _qp->Q);
  v_foutput(fp1, _qp->c);
  sp_foutput(fp1, _qp->A);
  v_foutput(fp1, _qp->b);
  sp_foutput(fp1, _qp->C);
  v_foutput(fp1, _qp->d);
  fclose(fp1);

  if ((fp2 = fopen("/tmp/comm/pipe2", "r")) == NULL)
    perror("opening pipe2");
  _qp->x = v_finput(fp2, _qp->x);
  _y = v_finput(fp2, _y);
  _z = v_finput(fp2, _z);
  fclose(fp2);

  _result = Hqp_Optimal;

  return IF_OK;
}

