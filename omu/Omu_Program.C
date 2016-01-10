/*
 * Omu_Program.C --
 *   -- class definition
 *
 * rf, 16/1/97
 */

/*
    Copyright (C) 1997--2008  Ruediger Franke

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

#ifdef OMU_WITH_ADOLC
#include <adolc/interfaces.h>
#include <adolc/adalloc.h>
#include <adolc/taping.h>
#endif

#include "Omu_Program.h"

#include "Omu_Deps.h" // needed for depreciated methods

#include <If_Int.h>
#include <If_IntVec.h>
#include <If_Real.h>
#include <If_RealVec.h>

IF_BASE_DEFINE(Omu_Program);

#ifdef OMU_WITH_ADOLC
#include "adoublev.h"
#endif

#define GET_SET_CB(vartype, name) \
  "prg_"#name, \
  IF_GET_CB(vartype, Omu_Program, name), \
  IF_SET_CB(vartype, Omu_Program, set_##name)

#ifdef OMU_WITH_ADOLC
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
#endif

//--------------------------------------------------------------------------
Omu_Program::Omu_Program()
{
  _K = 0;
  _KK = 0;
  _t0 = 0.0;
  _tf = 1.0;
  _ks = iv_get(1);
  _ts = v_get(1);

  ad_alloc(1, 1);

  _has_low_level_continuous = true;	// assume it was overloaded

  _ifList.append(new If_Int(GET_SET_CB(int, K)));
  _ifList.append(new If_Int(GET_SET_CB(int, KK)));
  _ifList.append(new If_Real(GET_SET_CB(double, t0)));
  _ifList.append(new If_Real(GET_SET_CB(double, tf)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, ks)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, ts)));
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
#ifdef OMU_WITH_ADOLC
  _U2 = m_get(m, m);
  m_ident(_U2);
  _Z2 = m_get(m, n);
  _Z3 = myalloc(m, n, 1);
  _nz = myalloc_short(m, n);

  _max_ndep = m;
  _max_nindep = n;
#else
  _U2 = NULL;
  _Z2 = NULL;
  _Z3 = NULL;
  _nz = NULL;

  _max_ndep = 0;
  _max_nindep = 0;
#endif
}

//--------------------------------------------------------------------------
void Omu_Program::ad_realloc(int ndep, int nindep)
{
#ifdef OMU_WITH_ADOLC
  if (ndep > _max_ndep || nindep > _max_nindep) {
    ad_free();
    ad_alloc(ndep, nindep);
  }
#endif
}

//--------------------------------------------------------------------------
void Omu_Program::ad_free()
{
#ifdef OMU_WITH_ADOLC
  free((char *)*_nz);
  free((char *)_nz);
  myfree(_Z3);
  m_free(_Z2);
  m_free(_U2);
#endif
}

//--------------------------------------------------------------------------
void Omu_Program::setup_stages()
{
  setup_stages(_ks, _ts);
  _K = _ks->dim - 1;
  _KK = _ts->dim - 1;
}

//--------------------------------------------------------------------------
void Omu_Program::setup_stages(IVECP ks, VECP ts)
{
  // default problem without stages
  iv_zero(iv_resize(ks, 1));
  v_zero(v_resize(ts, 1));
}

//--------------------------------------------------------------------------
static void setup_Jacobian(Real ***Z3, short **nz,
			   Omu_DependentVec &y,
			   MATP J, int J_flags, int i_offs, int j_offs,
                           bool isFirstCall)
{
  int ny = y->dim;
  int ncols = J->n;
  bool is_lin;
  int i, j;

  // mark linear elements of y
  for (i = 0; i < ny; i++) {
    is_lin = y.is_linear_element(i, J_flags);
    for (j = 0; j < ncols; j++) {
      switch (nz[i_offs + i][j_offs + j]) {
      case 0: // no dependency
	break;
      case 1: // linear
	if (isFirstCall) {
	  // store actual constant value
	  J[i][j] = *Z3[i_offs + i][j_offs + j];
	}
	else if (J[i][j] != *Z3[i_offs + i][j_offs + j]) {
	  // there are different linear dependencies over sample periods
	  is_lin = false;
	  J[i][j] = 1.0;
	}
	break;
      default: // non-linear
	is_lin = false;	// mark non-linear dependency
	J[i][j] = 1.0;	// store any value for sparsity pattern
      }
    }
    y.set_linear_element(i, J_flags, is_lin);
  }

  // mark linear variables
  for (j = 0; j < ncols; j++) {
    is_lin = y.is_linear_variable(J_flags, j);
    for (i = 0; i < ny; i++) {
      switch (nz[i_offs + i][j_offs + j]) {
      case 0: // no dependency
	break;
      case 1: // linear
        break;
      default: // non-linear
	is_lin = false;	// mark non-linear dependency
      }
    }
    y.set_linear_variable(J_flags, j, is_lin);
  }
}

//--------------------------------------------------------------------------
void Omu_Program::setup_struct(int k,
			       const Omu_VariableVec &x,
			       const Omu_VariableVec &u,
			       Omu_DependentVec &xt, Omu_DependentVec &F,
			       Omu_DependentVec &f,
			       Omu_Dependent &f0, Omu_DependentVec &c)
{
#ifdef OMU_WITH_ADOLC
  int i;
  int kk, kkend = k < K()? ks(k+1): ks(k) + 1;
  int nxt = x->dim;
  int nu = u->dim;
  int nf = f->dim;
  int nc = c->dim;
  bool with_continuous = (f0.gxf->dim > 0);
  bool init;

  //  adoublev ax(nxt);
  static adoublev ax; ax.alloc(nxt);
  //  adoublev au(nu);
  static adoublev au; au.alloc(nu);
  //  adoublev axt(nxt);
  static adoublev axt; axt.alloc(nxt);
  //  adoublev adx(nxt);
  static adoublev adx; adx.alloc(nxt);
  //  adoublev aF(nxt);
  static adoublev aF; aF.alloc(nxt);
  //  adoublev af(nf);
  static adoublev af; af.alloc(nf);
  adouble af0;
  //  adoublev ac(nc);
  static adoublev ac; ac.alloc(nc);
  int ndep, nindep;

  // initialize Jacobians with zeros;
  // afterwards dependencies during all sample periods are collected
  xt.set_linear();
  m_zero(xt.Jx);
  m_zero(xt.Ju);

  if (with_continuous) {
    F.set_linear();
    m_zero(F.Jx);
    m_zero(F.Ju);
    m_zero(F.Jdx);
  }

  f.set_linear();
  m_zero(f.Jx);
  m_zero(f.Ju);
  m_zero(f.Jxf);

  c.set_linear();
  m_zero(c.Jx);
  m_zero(c.Ju);
  m_zero(c.Jxf);

  init = true;  // indicate first sample period of stage
  for (kk = ks(k); kk < kkend; kk++) {

    //
    // analyze evaluation of consistent initial values
    //

    ndep = nxt;
    nindep = nxt + nu;
    ad_realloc(ndep, nindep);

    trace_on(4, 1);	// tape 4, keep = 1 for following reverse call
    ax <<= x->ve;
    au <<= u->ve;

    consistic(kk, ts(kk), ax, au, axt);

    axt >>= xt->ve; 

    trace_off();

    reverse(4, ndep, nindep, 0, ndep, _U2->me, _Z3, _nz);

    setup_Jacobian(_Z3, _nz, xt, xt.Jx, Omu_Dependent::WRT_x, 0, 0, init);
    setup_Jacobian(_Z3, _nz, xt, xt.Ju, Omu_Dependent::WRT_u, 0, nxt, init);

    //
    // analyze residuum evaluation
    //

    if (with_continuous) {
      ndep = nxt;
      nindep = 2*nxt + nu;
      ad_realloc(ndep, nindep);

      trace_on(4, 1);	// tape 4, keep = 1 for following reverse call
      ax <<= xt->ve;
      au <<= u->ve;
      for (i = 0; i < nxt; i++)
	adx[i] <<= 0.0;

      continuous(kk, ts(kk), ax, au, adx, aF);

      aF >>= F->ve; 

      trace_off();

      reverse(4, ndep, nindep, 0, ndep, _U2->me, _Z3, _nz);

      setup_Jacobian(_Z3, _nz, F, F.Jx, Omu_Dependent::WRT_x, 0, 0, init);
      setup_Jacobian(_Z3, _nz, F, F.Ju, Omu_Dependent::WRT_u, 0, nxt, init);
      setup_Jacobian(_Z3, _nz, F, F.Jdx, Omu_Dependent::WRT_dx, 0,
		     nxt + nu, init);

      // initialize final values of integration
      for (i = 0; i < nxt; i++)
	f[i] = xt[i];
    }

    //
    // analyze update 
    //

    trace_on(4, 1);	// tape 4, keep = 1 for following reverse call
    ax <<= x->ve;
    au <<= u->ve;
    nindep = nxt + nu;
    if (with_continuous) {
      for (i = 0; i < nxt; i++)
	af[i] <<= f[i];
      nindep += nxt;
    }

    update(kk, ax, au, af, af0, ac);

    for (i = 0; i < nf; i++)
      af[i] >>= f[i]; 
    af0 >>= f0; 
    ac >>= c->ve; 
    ndep = nf + 1 + nc;

    trace_off();

    ad_realloc(ndep, nindep);
    reverse(4, ndep, nindep, 0, ndep, _U2->me, _Z3, _nz);

    setup_Jacobian(_Z3, _nz, f, f.Jx, Omu_Dependent::WRT_x, 0, 0, init);
    setup_Jacobian(_Z3, _nz, f, f.Ju, Omu_Dependent::WRT_u, 0, nxt, init);
    setup_Jacobian(_Z3, _nz, f, f.Jxf, Omu_Dependent::WRT_xf,
		   0, nxt + nu, init);

    setup_Jacobian(_Z3, _nz, c, c.Jx, Omu_Dependent::WRT_x, nf + 1, 0, init);
    setup_Jacobian(_Z3, _nz, c, c.Ju, Omu_Dependent::WRT_u, nf + 1, nxt, init);
    setup_Jacobian(_Z3, _nz, c, c.Jxf, Omu_Dependent::WRT_xf,
		   nf + 1, nxt + nu, init);

    init = false;
  }
#else
  m_error(E_NULL, "Omu_Program::setup_struct: was compiled without ADOL-C");
#endif
}

//--------------------------------------------------------------------------
void Omu_Program::update(int kk, 
			 const adoublev &x, const adoublev &u,
			 adoublev &f, adouble &f0, adoublev &c)
{
  m_error(E_NULL, "Program does not implement update!");
}

//--------------------------------------------------------------------------
void Omu_Program::update(int kk, 
			 const Omu_StateVec &x, const Omu_Vec &u,
			 const Omu_StateVec &xf,
			 Omu_DependentVec &f, Omu_Dependent &f0,
			 Omu_DependentVec  &c)
{
#ifdef OMU_WITH_ADOLC
  int nx = x->dim;
  int nu = u->dim;
  int nxf = xf->dim;
  int nf = f->dim;
  int nc = c->dim;

  //  adoublev ax(nx);
  static adoublev ax; ax.alloc(nx);
  //  adoublev au(nu);
  static adoublev au; au.alloc(nu);
  //  adoublev af(max(nxf, nf));
  static adoublev af; af.alloc(max(nxf, nf));
  adouble af0;
  //  adoublev ac(nc);
  static adoublev ac; ac.alloc(nc);

  int ndep, nindep;
  int i, j;
  bool grds;
  double *Zi;

  // perform function evaluation

  af0 <<= f0;
  ac <<= c->ve;

  if (f.is_required_J() || f0.is_required_g() || c.is_required_J()) {
    grds = true;
    trace_on(4, 1);	// tape 4, keep result for subsequent reverse call
  }
  else
    grds = false;

  nindep = nx + nu;
  ax <<= x->ve;
  au <<= u->ve;
  if (!f.Jxf.is_constant() || !f0.gxf.is_constant() || !c.Jxf.is_constant()) {
    nindep += nxf;
    for (i = 0; i < nxf; i++)
      af[i] <<= xf[i];
  }

  update(kk, ax, au, af, af0, ac);

  for (i = 0; i < nf; i++)
    af[i] >>= f[i];
  af0 >>= f0;
  ac >>= c->ve;
  ndep = nf + 1 + nc;

  if (grds) {
    trace_off();

    // calculate Jacobians using reverse
    ad_realloc(ndep, nindep);
    reverse(4, ndep, nindep, 0, ndep, _U2->me, _Z2->me);

    if (!f.Jx.is_constant()) {
      for (i = 0; i < nf; i++) {
	Zi = _Z2[i];
	for (j = 0; j < nx; j++)
	  f.Jx[i][j] = Zi[j];
      }
    }
    if (!f.Ju.is_constant()) {
      for (i = 0; i < nf; i++) {
	Zi = _Z2[i];
	for (j = 0; j < nu; j++)
	  f.Ju[i][j] = Zi[nx + j];
      }
    }
    if (!f.Jxf.is_constant()) {
      for (i = 0; i < nf; i++) {
	Zi = _Z2[i];
	for (j = 0; j < nxf; j++)
	  f.Jxf[i][j] = Zi[nx + nu + j];
      }
    }

    Zi = _Z2[nf];
    if (!f0.gx.is_constant()) {
      for (j = 0; j < nx; j++)
	f0.gx[j] = Zi[j];
    }
    if (!f0.gu.is_constant()) {
      for (j = 0; j < nu; j++)
	f0.gu[j] = Zi[nx + j];
    }
    if (!f0.gxf.is_constant()) {
      for (j = 0; j < nxf; j++)
	f0.gxf[j] = Zi[nx + nu + j];
    }

    if (!c.Jx.is_constant()) {
      for (i = 0; i < nc; i++) {
	Zi = _Z2[nf + 1 + i];
	for (j = 0; j < nx; j++)
	  c.Jx[i][j] = Zi[j];
      }
    }
    if (!c.Ju.is_constant()) {
      for (i = 0; i < nc; i++) {
	Zi = _Z2[nf + 1 + i];
	for (j = 0; j < nu; j++)
	  c.Ju[i][j] = Zi[nx + j];
      }
    }
    if (!c.Jxf.is_constant()) {
      for (i = 0; i < nc; i++) {
	Zi = _Z2[nf + 1 + i];
	for (j = 0; j < nxf; j++)
	  c.Jxf[i][j] = Zi[nx + nu + j];
      }
    }
  }
#else
  m_error(E_NULL, "Omu_Program::update: was compiled without ADOL-C");
#endif
}

//--------------------------------------------------------------------------
void Omu_Program::consistic(int kk, double t,
			    const adoublev &x, const adoublev &u,
			    adoublev &xt)
{
#ifdef OMU_WITH_ADOLC
  xt = x;
#else
  m_error(E_NULL, "Omu_Program::consistic: was compiled without ADOL-C");
#endif
}

//--------------------------------------------------------------------------
void Omu_Program::consistic(int kk, double t,
			    const Omu_StateVec &x, const Omu_Vec &u,
			    Omu_DependentVec &xt)
{
#ifdef OMU_WITH_ADOLC
  int i, j;
  int nxt = x->dim;
  int nu = u->dim;
  int ndep = nxt;
  int nindep = nxt + nu;
  bool grds;
  double *Zi;

  ad_realloc(ndep, nindep);

  //  adoublev ax(nxt);
  static adoublev ax; ax.alloc(nxt);
  //  adoublev au(nu);
  static adoublev au; au.alloc(nu);
  //  adoublev axt(nxt);
  static adoublev axt; axt.alloc(nxt);

  for (i = 0; i < nxt; i++) {
    axt[i] = 0.0;
  }

  if (xt.is_required_J()) {
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

    reverse(4, ndep, nindep, 0, ndep, _U2->me, _Z2->me);

    if (!xt.Jx.is_constant()) {
      for (i = 0; i < nxt; i++) {
	Zi = _Z2[i];
	for (j = 0; j < nxt; j++)
	  xt.Jx[i][j] = Zi[j];
      }
    }
    if (!xt.Ju.is_constant()) {
      for (i = 0; i < nxt; i++) {
	Zi = _Z2[i];
	for (j = 0; j < nu; j++)
	  xt.Ju[i][j] = Zi[nxt+j];
      }
    }
  }
#else
  m_error(E_NULL, "Omu_Program::consistic: was compiled without ADOL-C");
#endif
}

//--------------------------------------------------------------------------
void Omu_Program::consistic(int kk, double t,
			    const Omu_Vector &x, const Omu_Vector &u,
			    VECP xt, MATP xtx, MATP xtu)
{
#ifdef OMU_WITH_ADOLC
  int i, j;
  int nxt = x->dim;
  int nu = u->dim;
  int ndep = nxt;
  int nindep = nxt + nu;
  bool grds;
  double *Zi;

  ad_realloc(ndep, nindep);

  //  adoublev ax(nxt);
  static adoublev ax; ax.alloc(nxt);
  //  adoublev au(nu);
  static adoublev au; au.alloc(nu);
  //  adoublev axt(nxt);
  static adoublev axt; axt.alloc(nxt);

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

    reverse(4, ndep, nindep, 0, ndep, _U2->me, _Z2->me);

    for (i = 0; i < nxt; i++) {
      Zi = _Z2[i];
      for (j = 0; j < nxt; j++)
	xtx[i][j] = Zi[j];

      for (j = 0; j < nu; j++)
	xtu[i][j] = Zi[nxt+j];
    }
  }
#else
  m_error(E_NULL, "Omu_Program::consistic: was compiled without ADOL-C");
#endif
}

//--------------------------------------------------------------------------
void Omu_Program::continuous(int kk, double t,
			     const adoublev &x, const adoublev &u,
			     const adoublev &dx, adoublev &F)
{
  // empty default implementation
}

//--------------------------------------------------------------------------
void Omu_Program::continuous(int kk, double t,
			     const Omu_StateVec &x, const Omu_Vec &u,
			     const Omu_StateVec &dx, Omu_DependentVec &F)
{
  if (F.is_required_J()) {
    if (F.Jdx.is_constant())
      continuous(kk, t, x, u, dx, F, F.Jx, F.Ju, MNULL);
    else
      continuous(kk, t, x, u, dx, F, F.Jx, F.Ju, F.Jdx);
  }
  else
    continuous(kk, t, x, u, dx, F, MNULL, MNULL, MNULL);
}

#if 0 // new version using arbitrary order reverse
// (this version does not work with Omu_IntODE calling it, e.g. Crane!?)
//--------------------------------------------------------------------------
void Omu_Program::continuous(int kk, double t,
			     const VECP x, const VECP u,
			     const VECP dx, VECP F,
			     MATP Fx, MATP Fu, MATP Fdx)
{
#ifdef OMU_WITH_ADOLC
  _has_low_level_continuous = false;	// i.e. not overloaded

  int i, j;
  int nxt = x->dim;
  int nu = u->dim;
  int ndep = nxt;
  int nindep = 2*nxt + nu;
  double **Zi;
  bool grds;

  ad_realloc(ndep, nindep);

  //  adoublev ax(nxt);
  static adoublev ax; ax.alloc(nxt);
  //  adoublev au(nu);
  static adoublev au; au.alloc(nu);
  //  adoublev adx(nxt);
  static adoublev adx; adx.alloc(nxt);
  //  adoublev aF(nxt);
  static adoublev aF; aF.alloc(nxt);

  for (i = 0; i < nxt; i++) {
    aF[i] = 0.0;
  }
  adx <<= dx->ve;	// initialize without ADOL-C trace

  if ((MAT *)Fx != MNULL && (MAT *)Fu != MNULL) {
    grds = true;
    trace_on(4, 1);	// tape 4, keep result for subsequent reverse call
  }
  else
    grds = false;

  ax <<= x->ve;
  au <<= u->ve;
  if ((MAT *)Fdx != MNULL)
    adx <<= dx->ve;
  else
    nindep -= nxt;
  
  continuous(kk, t, ax, au, adx, aF);

  aF >>= F->ve;

  if (grds) {
    trace_off();

    reverse(4, ndep, nindep, 0, ndep, _U2->me, _Z3);

    for (i = 0; i < nxt; i++) {
      Zi = _Z3[i];
      for (j = 0; j < nxt; j++)
	Fx[i][j] = *Zi[j];

      for (j = 0; j < nu; j++)
	Fu[i][j] = *Zi[nxt+j];
    }

    if ((MAT *)Fdx != MNULL) {
      for (i = 0; i < nxt; i++) {
	Zi = _Z3[i];
	for (j = 0; j < nxt; j++)
	  Fdx[i][j] = *Zi[nxt+nu+j];
      }
    }
  }
#else
  m_error(E_NULL, "Omu_Program::continuous: was compiled without ADOL-C");
#endif
}
#else
//--------------------------------------------------------------------------
void Omu_Program::continuous(int kk, double t,
			     const VECP x, const VECP u,
			     const VECP dx, VECP F,
			     MATP Fx, MATP Fu, MATP Fdx)
{
#ifdef OMU_WITH_ADOLC
  _has_low_level_continuous = false;	// i.e. not overloaded

  int i, j;
  int nxt = x->dim;
  int nu = u->dim;
  int ndep = nxt;
  int nindep = 2*nxt + nu;
  double *Zi;
  bool grds;

  ad_realloc(ndep, nindep);

  //  adoublev ax(nxt);
  static adoublev ax; ax.alloc(nxt);
  //  adoublev au(nu);
  static adoublev au; au.alloc(nu);
  //  adoublev adx(nxt);
  static adoublev adx; adx.alloc(nxt);
  //  adoublev aF(nxt);
  static adoublev aF; aF.alloc(nxt);

  for (i = 0; i < nxt; i++) {
    aF[i] = 0.0;
  }
  adx <<= dx->ve;	// initialize without ADOL-C trace

  if ((MAT *)Fx != MNULL && (MAT *)Fu != MNULL) {
    grds = true;
    trace_on(4, 1);	// tape 4, keep result for subsequent reverse call
  }
  else
    grds = false;

  ax <<= x->ve;
  au <<= u->ve;
  if ((MAT *)Fdx != MNULL)
    adx <<= dx->ve;
  else
    nindep -= nxt;
  
  continuous(kk, t, ax, au, adx, aF);

  aF >>= F->ve;

  if (grds) {
    trace_off();

    reverse(4, ndep, nindep, 0, ndep, _U2->me, _Z2->me);

    for (i = 0; i < nxt; i++) {
      Zi = _Z2[i];
      for (j = 0; j < nxt; j++)
	Fx[i][j] = Zi[j];

      for (j = 0; j < nu; j++)
	Fu[i][j] = Zi[nxt+j];
    }

    if ((MAT *)Fdx != MNULL) {
      for (i = 0; i < nxt; i++) {
	Zi = _Z2[i];
	for (j = 0; j < nxt; j++)
	  Fdx[i][j] = Zi[nxt+nu+j];
      }
    }
  }
#else
  m_error(E_NULL, "Omu_Program::continuous: was compiled without ADOL-C");
#endif
}
#endif

//--------------------------------------------------------------------------
void Omu_Program::init_simulation(int k,
				  Omu_VariableVec &x, Omu_VariableVec &u)
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

  if (t0 <= -Inf)
    t0 = _t0;
  if (tf >= Inf)
    tf = _tf;

  assert(tf >= t0 && K >= 0 && sps >= 0);

  v_resize(ts, KK + 1);
  ts[0] = t0; // explicitly assign first value to avoid div. by zero for KK=0
  for (i = 1; i <= KK; i++)
    ts[i] = t0 + (double)i * (tf - t0) / (double)KK;

  iv_resize(ks, K + 1);
  for (i = 0; i <= K; i++)
    ks[i] = i * sps;

  _K = K;
  _KK = KK;
  _t0 = t0;
  _tf = tf;
}    

//--------------------------------------------------------------------------
void Omu_Program::consistic_grds(int kk, double t,
				 const Omu_StateVec &x, const Omu_Vec &u,
				 Omu_DependentVec &xt)
{
  bool is_required_J = xt.is_required_J();

  xt.set_required_J(false);

  if (is_required_J) {
    int i, j;
    int nx = x->dim;
    int nu = u->dim;
    int nxt = xt->dim;
    double vi_bak, dvi;
    VECP xt_bak = v_copy(xt, VNULL);

    if (!xt.Jx.is_constant()) {
      for (i = 0; i < nx; i++) {
	vi_bak = x->ve[i];
	dvi = 1e-4 * fabs(vi_bak) + 1e-6;
	x->ve[i] += dvi;

	consistic(kk, t, x, u, xt);

	v_sub(xt, xt_bak, xt);
	for (j = 0; j < nxt; j++)
	  xt.Jx[j][i] = xt[j] / dvi;

	x->ve[i] = vi_bak;
      }
    }

    if (!xt.Ju.is_constant()) {
      for (i = 0; i < nu; i++) {
	vi_bak = u->ve[i];
	dvi = 1e-4 * fabs(vi_bak) + 1e-6;
	u->ve[i] += dvi;

	consistic(kk, t, x, u, xt);

	v_sub(xt, xt_bak, xt);
	for (j = 0; j < nxt; j++)
	  xt.Ju[j][i] = xt[j] / dvi;

	u->ve[i] = vi_bak;
      }
    }

    v_copy(xt_bak, xt);
    v_free(xt_bak);
  }

  // restore requirement
  xt.set_required_J(is_required_J);
}

//--------------------------------------------------------------------------
void Omu_Program::continuous_grds(int kk, double t,
				  const Omu_StateVec &x, const Omu_Vec &u,
				  const Omu_StateVec &dx, Omu_DependentVec &F,
                                  int wrt)
{
  bool is_required_J = F.is_required_J();

  F.set_required_J(false);

  if (is_required_J) {
    int i, j;
    int nx = x->dim;
    int nu = u->dim;
    int ndx = dx->dim;
    int nF = F->dim;
    double vi_bak, dvi;
    VECP F_bak = v_copy(F, VNULL);

    if ((wrt & Omu_Dependent::WRT_x) && !F.Jx.is_constant()) {
      for (i = 0; i < nx; i++) {
	vi_bak = x->ve[i];
	dvi = 1e-4 * fabs(vi_bak) + 1e-6;
	x->ve[i] += dvi;

	continuous(kk, t, x, u, dx, F);

	v_sub(F, F_bak, F);
	for (j = 0; j < nF; j++)
	  F.Jx[j][i] = F[j] / dvi;

	x->ve[i] = vi_bak;
      }
    }

    if ((wrt & Omu_Dependent::WRT_u) && !F.Ju.is_constant()) {
      for (i = 0; i < nu; i++) {
	vi_bak = u->ve[i];
	dvi = 1e-4 * fabs(vi_bak) + 1e-6;
	u->ve[i] += dvi;

	continuous(kk, t, x, u, dx, F);

	v_sub(F, F_bak, F);
	for (j = 0; j < nF; j++)
	  F.Ju[j][i] = F[j] / dvi;

	u->ve[i] = vi_bak;
      }
    }

    if ((wrt & Omu_Dependent::WRT_dx) && !F.Jdx.is_constant()) {
      for (i = 0; i < ndx; i++) {
	vi_bak = dx->ve[i];
	dvi = 1e-4 * fabs(vi_bak) + 1e-6;
	dx->ve[i] += dvi;

	continuous(kk, t, x, u, dx, F);

	v_sub(F, F_bak, F);
	for (j = 0; j < nF; j++)
	  F.Jdx[j][i] = F[j] / dvi;

	dx->ve[i] = vi_bak;
      }
    }

    v_copy(F_bak, F);
    v_free(F_bak);
  }

  // restore requirement
  F.set_required_J(is_required_J);
}

//--------------------------------------------------------------------------
void Omu_Program::update_grds(int kk, 
			      const Omu_StateVec &x, const Omu_Vec &u,
			      const Omu_StateVec &xf,
			      Omu_DependentVec &f, Omu_Dependent &f0,
			      Omu_DependentVec  &c)
{
  bool is_required_f_J = f.is_required_J();
  bool is_required_f0_g = f0.is_required_g();
  bool is_required_c_J = c.is_required_J();

  f.set_required_J(false);
  f0.set_required_g(false);
  c.set_required_J(false);

  if (is_required_f_J || is_required_f0_g || is_required_c_J) {
    int i, j;
    int nx = x->dim;
    int nu = u->dim;
    int nxf = xf->dim;
    int nf = f->dim;
    int nc = c->dim;
    double vi_bak, dvi;
    VECP f_bak = v_copy(f, VNULL);
    Real f0_bak = f0;
    VECP c_bak = v_copy(c, VNULL);
    
    if (!f.Jx.is_constant() ||!f0.gx.is_constant() ||!c.Jx.is_constant()) {
      for (i = 0; i < nx; i++) {
	vi_bak = x->ve[i];
	dvi = 1e-4 * fabs(vi_bak) + 1e-6;
	x->ve[i] += dvi;

	update(kk, x, u, xf, f, f0, c);

	v_sub(f, f_bak, f);
	for (j = 0; j < nf; j++)
	  f.Jx[j][i] = f[j] / dvi;

	f0.gx[i] = (f0 - f0_bak) / dvi;
	
	v_sub(c, c_bak, c);
	for (j = 0; j < nc; j++)
	  c.Jx[j][i] = c[j] / dvi;

	x->ve[i] = vi_bak;
      }
    }

    if (!f.Ju.is_constant() ||!f0.gu.is_constant() ||!c.Ju.is_constant()) {
      for (i = 0; i < nu; i++) {
	vi_bak = u->ve[i];
	dvi = 1e-4 * fabs(vi_bak) + 1e-6;
	u->ve[i] += dvi;

	update(kk, x, u, xf, f, f0, c);

	v_sub(f, f_bak, f);
	for (j = 0; j < nf; j++)
	  f.Ju[j][i] = f[j] / dvi;

	f0.gu[i] = (f0 - f0_bak) / dvi;
	
	v_sub(c, c_bak, c);
	for (j = 0; j < nc; j++)
	  c.Ju[j][i] = c[j] / dvi;

	u->ve[i] = vi_bak;
      }
    }

    if (!f.Jxf.is_constant() || !f0.gxf.is_constant()
	|| !c.Jxf.is_constant()) {
      for (i = 0; i < nxf; i++) {
	vi_bak = xf->ve[i];
	dvi = 1e-4 * fabs(vi_bak) + 1e-6;
	xf->ve[i] += dvi;

	update(kk, x, u, xf, f, f0, c);

	v_sub(f, f_bak, f);
	for (j = 0; j < nf; j++)
	  f.Jxf[j][i] = f[j] / dvi;

	f0.gxf[i] = (f0 - f0_bak) / dvi;
	
	v_sub(c, c_bak, c);
	for (j = 0; j < nc; j++)
	  c.Jxf[j][i] = c[j] / dvi;

	xf->ve[i] = vi_bak;
      }
    }

    v_copy(f_bak, f);
    v_free(f_bak);
    f0 = f0_bak;
    v_copy(c_bak, c);
    v_free(c_bak);
  }

  // restore requirements
  f.set_required_J(is_required_f_J);
  f0.set_required_g(is_required_f0_g);
  c.set_required_J(is_required_c_J);
}


//=========================================================================
