/*
 * Omu_Program.C --
 *   -- class definition
 *
 * rf, 16/1/97
 */

/*
    Copyright (C) 1997--2001  Ruediger Franke

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

#include <assert.h>

#include <adutils.h>

#include "Omu_Program.h"

#include <If_Int.h>
#include <If_IntVec.h>
#include <If_RealVec.h>

IF_BASE_DEFINE(Omu_Program);

//--------------------------------------------------------------------------
Omu_Program::Omu_Program()
{
  _K = 0;
  _KK = 0;
  _ks = iv_get(1);
  _ts = v_get(1);

  ad_alloc(1, 1);
  _has_low_level_continuous = true;	// assume it was overloaded

  _ifList.append(new If_Int("prg_K", &_K));
  _ifList.append(new If_Int("prg_KK", &_KK));
  _ifList.append(new If_IntVec("prg_ks", &_ks));
  _ifList.append(new If_RealVec("prg_ts", &_ts));
}

//--------------------------------------------------------------------------
Omu_Program::~Omu_Program()
{
  ad_free();
  v_free(_ts);
  iv_free(_ks);
}

//--------------------------------------------------------------------------
bool Omu_Program::has_low_level_continuous()
{
  return _has_low_level_continuous;
}

//--------------------------------------------------------------------------
void Omu_Program::ad_alloc(int m, int n)
{
  _Z = m_get(m, n);
  _U = m_get(m, m);
  m_ident(_U);

  _max_ndep = m;
  _max_nindep = n;
}

//--------------------------------------------------------------------------
void Omu_Program::ad_realloc(int ndep, int nindep)
{
  if (ndep > _max_ndep || nindep > _max_nindep) {
    ad_free();
    ad_alloc(ndep, nindep);
  }
}

//--------------------------------------------------------------------------
void Omu_Program::ad_free()
{
  m_free(_Z);
  m_free(_U);
}

//--------------------------------------------------------------------------
void Omu_Program::setup_stages()
{
  setup_stages(_ks, _ts);
}

//--------------------------------------------------------------------------
void Omu_Program::setup_stages(IVECP ks, VECP ts)
{
  // default problem without stages
  iv_zero(iv_resize(ks, 1));
  v_zero(v_resize(ts, 1));
}

//--------------------------------------------------------------------------
void Omu_Program::consistic(int kk, double t,
			 const adoublev &x, const adoublev &u,
			 adoublev &xt)
{
  xt = x;
}

//--------------------------------------------------------------------------
void Omu_Program::consistic(int kk, double t,
			    const Omu_Vector &x, const Omu_Vector &u,
			    VECP xt, MATP xtx, MATP xtu)
{
  int i, j;
  int nxt = x->dim;
  int nu = u->dim;
  int ndep = nxt;
  int nindep = nxt + nu;
  bool grds;
  double *Zi;

  ad_realloc(ndep, nindep);

  adoublev ax(nxt);
  adoublev au(nu);
  adoublev axt(nxt);

  for (i = 0; i < nxt; i++) {
    axt[i] = 0.0;
  }

  if ((MAT *)xtx != MNULL && (MAT *)xtu != MNULL) {
    grds = true;
    trace_on(4, 1);	// tape 4, keep result for subsequent reverse call
  }
  else
    grds = false;

  ax <<= x->ve;
  au <<= u->ve;

  consistic(kk, t, ax, au, axt);

  axt >>= xt->ve;

  if (grds) {
    trace_off();

    reverse(4, ndep, nindep, 0, ndep, _U->me, _Z->me);

    for (i = 0; i < nxt; i++) {
      Zi = _Z[i];
      for (j = 0; j < nxt; j++)
	xtx[i][j] = Zi[j];

      for (j = 0; j < nu; j++)
	xtu[i][j] = Zi[nxt+j];
    }
  }
}

//--------------------------------------------------------------------------
void Omu_Program::continuous(int kk, double t,
			     const adoublev &x, const adoublev &u,
			     const adoublev &xp, adoublev &F)
{
  // empty default implementation
}

//--------------------------------------------------------------------------
void Omu_Program::continuous(int kk, double t,
			     const VECP x, const VECP u,
			     const VECP xp, VECP F,
			     MATP Fx, MATP Fu, MATP Fxp)
{
  _has_low_level_continuous = false;	// i.e. not overloaded

  int i, j;
  int nxt = x->dim;
  int nu = u->dim;
  int ndep = nxt;
  int nindep = 2*nxt + nu;
  double *Zi;
  bool grds;

  ad_realloc(ndep, nindep);

  adoublev ax(nxt);
  adoublev au(nu);
  adoublev axp(nxt);
  adoublev aF(nxt);

  for (i = 0; i < nxt; i++) {
    aF[i] = 0.0;
  }
  axp <<= xp->ve;	// initialize without ADOL-C trace

  if ((MAT *)Fx != MNULL && (MAT *)Fu != MNULL) {
    grds = true;
    trace_on(4, 1);	// tape 4, keep result for subsequent reverse call
  }
  else
    grds = false;

  ax <<= x->ve;
  au <<= u->ve;
  if ((MAT *)Fxp != MNULL)
    axp <<= xp->ve;
  else
    nindep -= nxt;
  
  continuous(kk, t, ax, au, axp, aF);

  aF >>= F->ve;

  if (grds) {
    trace_off();

    reverse(4, ndep, nindep, 0, ndep, _U->me, _Z->me);

    for (i = 0; i < nxt; i++) {
      Zi = _Z[i];
      for (j = 0; j < nxt; j++)
	Fx[i][j] = Zi[j];

      for (j = 0; j < nu; j++)
	Fu[i][j] = Zi[nxt+j];
    }

    if ((MAT *)Fxp != MNULL) {
      for (i = 0; i < nxt; i++) {
	Zi = _Z[i];
	for (j = 0; j < nxt; j++)
	  Fxp[i][j] = Zi[nxt+nu+j];
      }
    }
  }
}

//--------------------------------------------------------------------------
void Omu_Program::init_simulation(int k,
				  Omu_Vector &x, Omu_Vector &u)
{
  if (k == 0)
    v_copy(x.initial, x);
  v_copy(u.initial, u);
}

//--------------------------------------------------------------------------
void Omu_Program::stages_alloc(IVECP ks, VECP ts, int K, int sps,
			       double t0, double tf)
{
  int i, KK = K * sps;

  assert(tf >= t0 && K >= 0 && sps >= 0);

  v_resize(ts, KK + 1);
  for (i = 0; i <= KK; i++)
    ts[i] = t0 + (double)i * (tf - t0) / (double)KK;

  iv_resize(ks, K + 1);
  for (i = 0; i <= K; i++)
    ks[i] = i * sps;

  _K = K;
  _KK = KK;
}    


//=========================================================================
