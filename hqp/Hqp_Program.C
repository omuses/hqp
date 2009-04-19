/*
 * Hqp_Program.C
 *
 * rf, 5/15/94
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Hqp_Program.h"

// flags 
const int Hqp_Program::IS_LOCAL = 	0x0001;
const int Hqp_Program::IS_SLACK = 	0x0002;

//-------------------------------------------------------------------------
Hqp_Program::Hqp_Program()
{
  x = VNULL;
  x_flags = x_int = IVNULL;
  Q = A = C = SMNULL;
  c = b = d = VNULL;
}

//-------------------------------------------------------------------------
Hqp_Program::~Hqp_Program()
{
  v_free(x);
  iv_free(x_flags);
  iv_free(x_int);
  sp_free(Q);
  v_free(c);
  sp_free(A);
  v_free(b);
  sp_free(C);
  v_free(d);
}

//-------------------------------------------------------------------------
void Hqp_Program::resize(int n, int me, int m,
			 int el_n, int el_me, int el_m)
{
  double log2;

  sp_free(Q);
  sp_free(A);
  sp_free(C);

  log2 = log(2.0);
  if (el_n < 1)
    el_n = n > 0? 1 + (int)(log((double)n) / log2): 0;
  if (el_me < 1)
    el_me = me > 0? 1 + (int)(log((double)me) / log2): 0;
  if (el_m < 1)
    el_m = m > 0? 1 + (int)(log((double)m) / log2): 0;

  x = v_resize(x, n);
  x_flags = iv_resize(x_flags, n);
  x_int = iv_resize(x_int, n);

  Q = sp_get(n, n, el_n);
  c = v_resize(c, n);

  A = sp_get(me, n, el_me);
  b = v_resize(b, me);

  C = sp_get(m, n, el_m);
  d = v_resize(d, m);
}

//-------------------------------------------------------------------------
void Hqp_Program::foutput(FILE *fp)
{
  fprintf(fp, "# Quadratic Program\n");
  fprintf(fp, "#  criterion: 1/2 x'Qx + c'x -> min\n");
  fprintf(fp, "#  equality constraints: Ax + b = 0\n");
  fprintf(fp, "#  inequality constraints: Cx + d >= 0\n");

  fprintf(fp, "#\n");
  fprintf(fp, "# Q:\n");
  sp_foutput(fp, Q);

  fprintf(fp, "#\n");
  fprintf(fp, "# c:\n");
  v_foutput(fp, c);

  fprintf(fp, "#\n");
  fprintf(fp, "# A:\n");
  sp_foutput(fp, A);

  fprintf(fp, "#\n");
  fprintf(fp, "# b:\n");
  v_foutput(fp, b);

  fprintf(fp, "#\n");
  fprintf(fp, "# C:\n");
  sp_foutput(fp, C);

  fprintf(fp, "#\n");
  fprintf(fp, "# d:\n");
  v_foutput(fp, d);
}

//=========================================================================
