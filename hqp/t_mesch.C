/*
 * t_mesch.C --
 *   - 'time-dependent' Meschach data structures:
 *
 *        TVECP   - vector of VECP
 *        TIVECP  - vector of IVECP
 *        PERMP   - PERM *
 *        TPERMP  - vector of PERMP
 *        TMATP   - vector of MATP
 *
 * E. Arnold  03/07/97
 *            07/03/97: _x[k] = VNULL etc
 *            2002-04-17 free() replaced by tfree()
 *                       delete[]
 *
 */

/*
    Copyright (C) 1997--2002  Eckhard Arnold

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

#include "t_mesch.h"

//----------   TVECP   ----------

TVECP::TVECP()
{
  _kmax = -1;
  _x = NULL;
}

TVECP::~TVECP()
{
  tfree();
}

void TVECP::tfree()
{
  int k;
  for (k = 0; k <= _kmax; k++)
    V_FREE(_x[k]);
  if ( _x ) {
    delete[] _x;
    _x = NULL;
  }
  _kmax = -1;
}

void TVECP::resize(int kmax, int n)
{
  int k;
  if ( ( kmax <= 0 ) || ( n < 0 ) )
    m_error(E_BOUNDS,"TVECP::resize");
  if ( ( _kmax != -1 ) && ( _kmax != kmax ) ) 
    tfree();
  if ( _kmax == -1 ) {
    _kmax = kmax;
    _x = new VECP[_kmax+1];
    if ( _x == NULL )
      m_error(E_MEM,"TVECP::resize");
    for ( k = 0; k <= _kmax; k++ ) 
      _x[k] = VNULL;
  }
  for ( k = 0; k <= _kmax; k++ ) {
    if  ( ( n == 0) && ( (VEC *) _x[k] == VNULL ) ) 
      _x[k] = v_resize(_x[k], 1);       // it's not a bug, it's a feature!
    _x[k] = v_resize(_x[k], n);
  }
}

void TVECP::Print()
{
  int k;
  for (k = 0; k <= _kmax; k++) {
    printf("k = %d\n", k);
    v_output(_x[k]);
  }
}

//----------   TIVECP   ----------

TIVECP::TIVECP()
{
  _kmax = -1;
  _x = NULL;
}

TIVECP::~TIVECP()
{
  tfree();
}

void TIVECP::tfree()
{
  int k;
  for (k = 0; k <= _kmax; k++)
    IV_FREE(_x[k]);
  if ( _x ) {
    delete[] _x;
    _x = NULL;
  }
  _kmax = -1;
}

void TIVECP::resize(int kmax, int n)
{
  int k;
  if ( ( kmax <= 0 ) || ( n < 0 ) )
    m_error(E_BOUNDS,"TIVECP::resize");
  if ( ( _kmax != -1 ) && ( _kmax != kmax ) ) 
    tfree();
  if ( _kmax == -1 ) {
    _kmax = kmax;
    _x = new IVECP[_kmax+1];
    if ( _x == NULL )
      m_error(E_MEM,"TIVECP::resize");
    for ( k = 0; k <= _kmax; k++ ) 
      _x[k] = IVNULL;
  }
  for ( k = 0; k <= _kmax; k++ ) {
    if  ( ( n == 0) && ( (IVEC *) _x[k] == IVNULL ) ) 
      _x[k] = iv_resize(_x[k], 1);       // it's not a bug, it's a feature!
    _x[k] = iv_resize(_x[k], n);
  }
}

void TIVECP::Print()
{
  int k;
  for (k = 0; k <= _kmax; k++) {
    printf("k = %d\n", k);
    iv_output(_x[k]);
  }
}

//----------   TPERMP   ----------

TPERMP::TPERMP()
{
  _kmax = -1;
  _x = NULL;
}

TPERMP::~TPERMP()
{
  tfree();
}

void TPERMP::tfree()
{
  int k;
  for (k = 0; k <= _kmax; k++)
    PX_FREE(_x[k]);
  if ( _x ) {
    delete[] _x;
    _x = NULL;
  }
  _kmax = -1;
}

void TPERMP::resize(int kmax, int n)
{
  int k;
  if ( ( kmax <= 0 ) || ( n < 0 ) )
    m_error(E_BOUNDS,"TPERMP::resize");
  if ( ( _kmax != -1 ) && ( _kmax != kmax ) ) 
    tfree();
  if ( _kmax == -1 ) {
    _kmax = kmax;
    _x = new PERMP[_kmax+1];
    if ( _x == NULL )
      m_error(E_MEM,"TPERMP::resize");
    for ( k = 0; k <= _kmax; k++ ) 
      _x[k] = PNULL;
  }
  for ( k = 0; k <= _kmax; k++ ) {
    if  ( ( n == 0) && ( (PERM *) _x[k] == PNULL ) ) 
      _x[k] = px_resize(_x[k], 1);       // it's not a bug, it's a feature!
    _x[k] = px_resize(_x[k], n);
  }
}

void TPERMP::Print()
{
  int k;
  for (k = 0; k <= _kmax; k++) {
    printf("k = %d\n", k);
    px_output(_x[k]);
  }
}

//----------   TMATP   ----------

TMATP::TMATP()
{
  _kmax = -1;
  _x = NULL;
}

TMATP::~TMATP()
{
  tfree();
}

void TMATP::tfree()
{
  int k;
  for (k = 0; k <= _kmax; k++)
    M_FREE(_x[k]);
  if ( _x ) {
    delete[] _x;
    _x = NULL;
  }
  _kmax = -1;
}

void TMATP::resize(int kmax, int m, int n)
{
  int k;
  if ( ( kmax <= 0 ) || ( n < 0 )  || ( m < 0 ) )
    m_error(E_BOUNDS,"TMATP:resize");
  if ( ( _kmax != -1 ) && ( _kmax != kmax ) ) 
    tfree();
  if ( _kmax == -1 ) {
    _kmax = kmax;
    _x = new MATP[_kmax+1];
    if ( _x == NULL )
      m_error(E_MEM,"TMATP::resize");
    for ( k = 0; k <= _kmax; k++ ) 
      _x[k] = MNULL;
  }
  for ( k = 0; k <= _kmax; k++ ) {
    if  ( ( ( n == 0) || ( m == 0 ) ) && ( (MAT *) _x[k] == MNULL ) ) 
      _x[k] = m_resize(_x[k], 1, 1);       // it's not a bug, it's a feature!
    _x[k] = m_resize(_x[k], m, n);
  }
}

void TMATP::Print()
{
  int k;
  for (k = 0; k <= _kmax; k++) {
    printf("k = %d\n", k);
    m_output(_x[k]);
  }
}

