/*
 * Hqp_DocpAdol.C --
 *   -- class implementation
 *
 * rf, 12/29/95
 */

/*
    Copyright (C) 1994--1998  Ruediger Franke

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

#include <malloc.h>

#include <If_Int.h>

#include "Hqp_DocpAdol.h"


//--------------------------------------------------------------------------
static short** myalloc_short(int m, int n)
{
  short* Adum = (short*)malloc(m*n*sizeof(short));
  short**  A = (short**)malloc(m*sizeof(short*));
  int i;
  for(i=0;i<m;i++)
    {
      A[i] = Adum;
      Adum += n;
    }
  return A;

  /* To deallocate an array set up by   A = myalloc2(m,n)   */
  /*   use  free((char*)*A); free((char*)A); in that order  */
}

//--------------------------------------------------------------------------
Hqp_DocpAdol::Hqp_DocpAdol()
{
  _k0 = 0;
  _kmax = 0;

  _static_struct = 1;	// analyze the sparsity pattern
  _hela = 0;		// calculate analytical Lagrangian Hessian
  _ad = 1;		// use automatic differentiation for derivatives

  adalloc(1, 1);

  _ifList.append(new If_Int("prg_static_struct", &_static_struct));
  _ifList.append(new If_Int("prg_hela", &_hela));
  _ifList.append(new If_Int("prg_ad", &_ad));
}

//--------------------------------------------------------------------------
Hqp_DocpAdol::~Hqp_DocpAdol()
{
  adfree();
}

//--------------------------------------------------------------------------
void Hqp_DocpAdol::adalloc(int m, int n)
{
  _am = m;
  _an = n;
  _X = myalloc(_an, 2);
  _Y = myalloc(_am, 2);
  _U = myalloc(1, _am);
  _Z = myalloc(_an, 2);
  _Z3 = myalloc(_am, _an, 1);
  _nz = myalloc_short(_am, _an);
}

//--------------------------------------------------------------------------
void Hqp_DocpAdol::adfree()
{
  free(*_X);
  free(_X);
  free(*_Y);
  free(_Y);
  free(*_U);
  free(_U);
  free(*_Z);
  free(_Z);
  free(**_Z3);
  free(*_Z3);
  free(_Z3);
  free(*_nz);
  free(_nz);
}

//--------------------------------------------------------------------------
void Hqp_DocpAdol::adresize(int m, int n)
{
  if (m > _am || n > _an) {
    adfree();
    adalloc(m, n);
  }
}

//--------------------------------------------------------------------------
void Hqp_DocpAdol::horizon(int k0, int kmax)
{
  _k0 = k0;
  _kmax = kmax;
  Hqp_Docp::horizon(k0, kmax);
}

//--------------------------------------------------------------------------
void Hqp_DocpAdol::setup_struct(int k,
				MATP fx, MATP fu, IVECP f_lin,
				MATP cx, MATP cu, IVECP c_lin,
				MATP Lxx, MATP Luu, MATP Lxu)
{
  if (!_static_struct || !_ad)
    return;	  // initialize matrices to be dense
  
  //    -- evaluate update_avals for recording, keep = 1
  //    -- calculate Jacobians in reverse mode
  //  additionally, if the Lagrangian Hessian is required:
  //    -- exlude some trivial elements
  //

  VECP xk = x(k);
  VECP uk = u(k);
  int xdim = xk->dim;
  int udim = uk->dim;
  int fdim = f_lin->dim;
  int cdim = c_lin->dim;
  adoublev ax(xdim);
  adoublev au(udim);
  adoublev af(fdim);
  adouble  af0;
  adoublev ac(cdim);
  int ndep, nindep;
  int i, j, l, idx;
  int maxdim;
  VECP dumvec;

  maxdim = max(xdim, cdim);
  maxdim = max(maxdim, udim);
  maxdim = max(1, maxdim);
  dumvec = v_get(maxdim);

  // obtain useful variable values
  // !! store the result of update_avals() back to x(k+1) !!
  init_solution(k, xk, uk);

  // record function evaluation

  trace_on(1, 1);	// tape 1, keep = 1 for following reverse call
  ax <<= xk->ve;
  au <<= uk->ve;

  update_avals(k, ax, au, af, af0, ac);

  af0 >>= dumvec->ve[0];
  if (k < _kmax)
    af >>= x(k+1)->ve;
  else
    af >>= dumvec->ve;
  ac >>= dumvec->ve;
  trace_off();

  // calculate first derivatives using reverse

  ndep = 1 + fdim + cdim;
  nindep = xdim + udim;
  adresize(ndep, nindep);

  reverse(1, ndep, nindep, 0, _Z3, _nz);

  /* sparsity structure for f0x, f0u not yet required
  for (j = 0; j < xdim; j++)
    f0x[j] = _nz[0][j] == 0? 0.0: 1.0;
  for (j = 0; j < udim; j++)
    f0u[j] = _nz[0][xdim + j] == 0? 0.0: 1.0;
  */

  for (i = 0; i < fdim; i++) {
    idx = 1 + i;
    for (j = 0; j < xdim; j++)
      fx[i][j] = _nz[idx][j] == 0? 0.0: 1.0;
    for (j = 0; j < udim; j++)
      fu[i][j] = _nz[idx][xdim + j] == 0? 0.0: 1.0;
  }

  for (i = 0; i < cdim; i++) {
    idx = 1 + fdim + i;
    for (j = 0; j < xdim; j++)
      cx[i][j] = _nz[idx][j] == 0? 0.0: 1.0;
    for (j = 0; j < udim; j++)
      cu[i][j] = _nz[idx][xdim + j] == 0? 0.0: 1.0;
  }

#if 0
  short nzij;

  // exclude all Lagrangian Hessian elements,
  // where functions depend not or linearly from a variable

  // for each column of _nz: store the maximum of all rows in _nz[0]
  for (j = 0; j < nindep; j++) {
    nzij = _nz[0][j];
    for (i = 1; i < ndep; i++)
      nzij = max(nzij, _nz[i][j]);
    _nz[0][j] = nzij;
  }

  for (i = 0; i < xdim; i++)
    for (j = i + 1; j < xdim; j++)
      if (_nz[0][i] < 2 || _nz[0][j] < 2)
	Lxx[i][j] = 0.0;

  for (i = 0; i < udim; i++)
    for (j = i + 1; j < udim; j++)
      if (_nz[0][xdim + i] < 2 || _nz[0][xdim + j] < 2)
	Luu[i][j] = 0.0;

  for (i = 0; i < xdim; i++)
    for (j = 0; j < udim; j++)
      if (_nz[0][i] < 2 || _nz[0][xdim + j] < 2)
	Lxu[i][j] = 0.0;

#else

  // insert all Lagrangian Hessian elements, 
  // where a function depends more than linearly from both variables
  // associated with the indices

  m_zero(Lxx);
  m_zero(Luu);
  m_zero(Lxu);

  for (l = 0; l < ndep; l++) {

    for (i = 0; i < xdim; i++)
      if (_nz[l][i] > 1)
	for (j = i; j < xdim; j++)
	  if (_nz[l][j] > 1)
	    Lxx[i][j] = 1.0;

    for (i = 0; i < udim; i++)
      if (_nz[l][xdim + i] > 1)
	for (j = i; j < udim; j++)
	  if (_nz[l][xdim + j] > 1)
	    Luu[i][j] = 1.0;

    for (i = 0; i < xdim; i++)
      if (_nz[l][i] > 1)
	for (j = 0; j < udim; j++)
	  if (_nz[l][xdim + j] > 1)
	    Lxu[i][j] = 1.0;
  }
#endif

  v_free(dumvec);
}

//--------------------------------------------------------------------------
void Hqp_DocpAdol::update_vals(int k, const VECP x, const VECP u,
			       VECP f, Real &f0, VECP c)
{
  //
  // just evaluate update_avals
  //

  adoublev ax(x->dim);
  adoublev au(u->dim);
  adoublev af(f->dim);
  adouble  af0;
  adoublev ac(c->dim);

  af0 <<= f0;

  ax <<= x->ve;
  au <<= u->ve;

  update_avals(k, ax, au, af, af0, ac);

  af0 >>= f0;
  af >>= f->ve;
  ac >>= c->ve;
}

//--------------------------------------------------------------------------
void Hqp_DocpAdol::update_stage(int k, const VECP x, const VECP u,
				VECP f, Real &f0, VECP c,
				MATP fx, MATP fu, VECP f0x, VECP f0u,
				MATP cx, MATP cu,
				const VECP rf, const VECP rc,
				MATP Lxx, MATP Luu, MATP Lxu)
{
  if (!_ad) {
    Hqp_Docp::update_stage(k, x, u,
			   f, f0, c,
			   fx, fu, f0x, f0u, cx, cu,
			   rf, rc, Lxx, Luu, Lxu);
    return;
  }

  //
  //    -- evaluate update_avals for recording, keep = 1
  //    -- calculate Jacobians in reverse mode
  //  additionally, if the Lagrangian Hessian is required:
  //   repeat for each row:
  //    -- calculate forward with keep = 2
  //    -- calculate second derivatives in reverse mode
  //

  int xdim = x->dim;
  int udim = u->dim;
  int fdim = rf->dim;
  int cdim = rc->dim;
  adoublev ax(xdim);
  adoublev au(udim);
  adoublev af(fdim);
  adouble  af0;
  adoublev ac(cdim);
  int ndep, nindep;
  int i, j, idx;
  int ninfs;
  double val;

  // record function evaluation

  af0 <<= f0;

  trace_on(1, 1);	// tape 1, keep = 1 for following reverse call
  ax <<= x->ve;
  au <<= u->ve;

  update_avals(k, ax, au, af, af0, ac);

  af0 >>= f0;
  af >>= f->ve;
  ac >>= c->ve;
  trace_off();

  // calculate first derivatives using reverse

  ndep = 1 + fdim + cdim;
  nindep = xdim + udim;
  adresize(ndep, nindep);

  reverse(1, ndep, nindep, 0, _Z3);

  for (j = 0; j < xdim; j++)
    f0x[j] = _Z3[0][j][0];
  for (j = 0; j < udim; j++)
    f0u[j] = _Z3[0][xdim + j][0];

  for (i = 0; i < fdim; i++) {
    idx = 1 + i;
    for (j = 0; j < xdim; j++)
      fx[i][j] = _Z3[idx][j][0];
    for (j = 0; j < udim; j++)
      fu[i][j] = _Z3[idx][xdim + j][0];
  }

  for (i = 0; i < cdim; i++) {
    idx = 1 + fdim + i;
    for (j = 0; j < xdim; j++)
      cx[i][j] = _Z3[idx][j][0];
    for (j = 0; j < udim; j++)
      cu[i][j] = _Z3[idx][xdim + j][0];
  }

  if (_hela) {

    // init weighting vector
    _U[0][0] = 1.0;
    j = 1;
    for (i = 0; i < fdim; i++, j++)
      _U[0][j] = -rf[i];
    for (i = 0; i < cdim; i++, j++)
      _U[0][j] = -rc[i];

    // init independent variable values
    j = 0;
    for (i = 0; i < xdim; i++, j++) {
      _X[j][0] = x[i];
      _X[j][1] = 0.0;
    }
    for (i = 0; i < udim; i++, j++) {
      _X[j][0] = u[i];
      _X[j][1] = 0.0;
    }

    // calculate Lagrangian Hessian by rows
    // (only the upper diagonal part)

    ninfs = 0;	// count the number of infinite values
    for (i = 0; i < nindep; i++) {
      _X[i][1] = 1.0;

      forward(1, ndep, nindep, 1, 2, _X, _Y);
      reverse(1, ndep, nindep, 1, _U[0], _Z);

      if (i < xdim) {
	for (j = 0; j < xdim; j++) {
	  val = _Z[j][1];
	  if (is_finite(val))
	    Lxx[i][j] = val;
	  else
	    ninfs++;
	}	    
	for (j = 0; j < udim; j++) {
	  val = _Z[xdim + j][1];
	  if (is_finite(val))
	    Lxu[i][j] = val;
	  else
	    ninfs++;
	}
      }
      else {
	idx = i - xdim;
	for (j = idx; j < udim; j++) {
	  val = _Z[xdim + j][1];
	  if (is_finite(val))
	    Luu[idx][j] = val;
	  else
	    ninfs++;
	}
      }

      _X[i][1] = 0.0;
    }

    if (ninfs > 0)
      fprintf(stderr, "%d infinite Lagrangian Hessian values!\n", ninfs);

  }
}

//========================================================================
