/*
 * Hqp_Docp.C -- class definition
 *
 * rf, 11/12/94
 */

/*
    Copyright (C) 1994--2017  Ruediger Franke

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
#include <math.h>

#include "Hqp_Docp.h"
#include "Hqp_Program.h"
#include "Hqp_omp.h"

#include <If_Int.h>
#include <If_IntVec.h>
#include <If_Method.h>

typedef If_Method<Hqp_Docp> If_Cmd;

// store association of constraints with _qp
//------------------------------------------
class Hqp_DocpAssoc {

 public:

  int	granul;	// quantity for growing up vectors
  int	offs;	// offset into _qp vectors and matrices
  int	dim;	// dimension of vectors vals, idxs, and flags
  VECP	vals;	// values of bounds
  IVECP	idxs;	// index set to associate constraints with variables
  IVECP	start;	// start index for each k into vals, and idxs

  Hqp_DocpAssoc();
  ~Hqp_DocpAssoc();

  void reset();
  void new_k();
  void add_cns(Real val, int idx);

  void sp_insert_mat(int k, SPMAT *J, int joffs,
		     const MATP mat, int mat_offs);
  void sp_extract_mat(int k, const SPMAT *J, int joffs,
		      MATP mat, int mat_offs, const IVECP lin);
  void sp_update_mat(int k, SPMAT *J, int joffs,
		     const MATP mat, int mat_offs, const IVECP lin);

  void sv_mltadd_vec(int k, Real s, const VECP z,
		     VECP vec, int vec_offs);
};

//-------------------------------------------------------------------------
Hqp_DocpAssoc::Hqp_DocpAssoc()
{
  granul = 100;
  vals = v_get(granul);
  idxs = iv_get(granul);
  start = iv_get(granul);
  reset();
}

//-------------------------------------------------------------------------
Hqp_DocpAssoc::~Hqp_DocpAssoc()
{
  v_free(vals);
  iv_free(idxs);
  iv_free(start);
}

//-------------------------------------------------------------------------
void Hqp_DocpAssoc::reset()
{
  offs = 0;
  dim = 0;
  v_resize(vals, 0);
  iv_resize(idxs, 0);
  iv_resize(start, 1);
  start[0] = 0;
}

//-------------------------------------------------------------------------
void Hqp_DocpAssoc::new_k()
{
  // store current dim as start index for new k

  iv_expand(start, 1, granul);
  start[(int)start->dim - 1] = dim;
}

//-------------------------------------------------------------------------
void Hqp_DocpAssoc::add_cns(Real val, int idx)
{
  v_expand(vals, 1, granul);
  iv_expand(idxs, 1, granul);
  vals[dim] = val;
  idxs[dim] = idx;
  dim ++;
}

//-------------------------------------------------------------------------
void Hqp_DocpAssoc::sp_insert_mat(int k, SPMAT *J, int joffs,
				  const MATP mat, int mat_offs)
{
  int i, iend;
  int mat_idx;

  iend = start[k+1];
  for (i = start[k]; i < iend; i++) {
    mat_idx = idxs[i] - mat_offs;
    sp_insert_mrow(J, offs + i - mat_idx, joffs, mat, mat_idx);
  }
}

//-------------------------------------------------------------------------
void Hqp_DocpAssoc::sp_extract_mat(int k, const SPMAT *J, int joffs,
				   MATP mat, int mat_offs, const IVECP lin)
{
  int i, iend;
  int mat_idx;

  iend = start[k+1];
  for (i = start[k]; i < iend; i++) {
    if (!lin[idxs[i]]) {
      mat_idx = idxs[i] - mat_offs;
      sp_extract_mrow(J, offs + i - mat_idx, joffs, mat, mat_idx);
    }
  }
}

//-------------------------------------------------------------------------
void Hqp_DocpAssoc::sp_update_mat(int k, SPMAT *J, int joffs,
				  const MATP mat, int mat_offs,
				  const IVECP lin)
{
  int i, iend;
  int mat_idx;

  iend = start[k+1];
  for (i = start[k]; i < iend; i++) {
    if (!lin[idxs[i]]) {
      mat_idx = idxs[i] - mat_offs;
      sp_update_mrow(J, offs + i - mat_idx, joffs, mat, mat_idx);
    }
  }
}

//-------------------------------------------------------------------------
void Hqp_DocpAssoc::sv_mltadd_vec(int k, Real s, const VECP z,
				  VECP vec, int vec_offs)
{
  int i, iend;
  int ioffs = offs + start[k];
  Real *z_ve = z->ve;
  Real *v_ve = vec->ve;

  iend = start[k+1];
  for (i = start[k]; i < iend; i++, ioffs++) {
    v_ve[idxs[i]-vec_offs] += s * z_ve[ioffs];
  }
}

//=========================================================================

const int Hqp_Docp::Periodical = 1;

//-------------------------------------------------------------------------
Hqp_Docp::Hqp_Docp(int ncpu)
{
  _ncpu = ncpu;
  _xk = new VECP [_ncpu];
  _uk = new VECP [_ncpu];
  _fk = new VECP [_ncpu];
  _s1 = new VECP [_ncpu];
  _s2 = new VECP [_ncpu];
  _ck = new VECP [_ncpu];
  _fkx = new MATP [_ncpu];
  _fku = new MATP [_ncpu];
  _ckx = new MATP [_ncpu];
  _cku = new MATP [_ncpu];
  _Lkxx = new MATP [_ncpu];
  _Lkuu = new MATP [_ncpu];
  _Lkxu = new MATP [_ncpu];
  _vfk = new VECP [_ncpu];
  _vck = new VECP [_ncpu];
  _xk_head = new VEC [_ncpu];
  _uk_head = new VEC [_ncpu];
  _fk_head = new VEC [_ncpu];
  _s1_head = new VEC [_ncpu];
  _s2_head = new VEC [_ncpu];
  for (int tn = 0; tn < _ncpu; tn++) {
    _xk[tn] = &_xk_head[tn];
    _uk[tn] = &_uk_head[tn];
    _fk[tn] = &_fk_head[tn];
    _s1[tn] = &_s1_head[tn];
    _s2[tn] = &_s2_head[tn];
    _ck[tn] = v_get(1);
    _fkx[tn] = m_get(1, 1);
    _fku[tn] = m_get(1, 1);
    _ckx[tn] = m_get(1, 1);
    _cku[tn] = m_get(1, 1);
    _Lkxx[tn] = m_get(1, 1);
    _Lkuu[tn] = m_get(1, 1);
    _Lkxu[tn] = m_get(1, 1);
    _vfk[tn] = v_get(1);
    _vck[tn] = v_get(1);
  }
  _f0 = v_get(_ncpu);
  _k0 = 0;
  _kf = 50;
  _nxs = iv_get(1);
  _nus = iv_get(1);
  _xus_init = v_get(1);
  _xus_integer = iv_get(1);
  _granul = 100;
  _f0_lin = iv_get(1);
  _xu_eq = new Hqp_DocpAssoc;
  _xu_lb = new Hqp_DocpAssoc;
  _xu_ub = new Hqp_DocpAssoc;
  _f_lin = iv_get(1);
  _f_start = iv_get(1);
  _xu_start = iv_get(1);
  _cns_eq = new Hqp_DocpAssoc;
  _cns_lb = new Hqp_DocpAssoc;
  _cns_ub = new Hqp_DocpAssoc;
  _cns_lin = iv_get(1);
  _cns_start = iv_get(1);
  _x_type = iv_get(1);
  _fbd_evals = 0;

  _ifList.append(new If_Cmd("prg_simulate",
			    &Hqp_Docp::simulate, this));
  _ifList.append(new If_IntVec("prg_nxs",
			       IF_GET_CB(const IVECP, Hqp_Docp, nxs)));
  _ifList.append(new If_IntVec("prg_nus",
			       IF_GET_CB(const IVECP, Hqp_Docp, nus)));
  _ifList.append(new If_Int("prg_fbd_evals", 
			    IF_GET_CB(int, Hqp_Docp, fbd_evals)));
  _ifList.append(new If_Int("prg_ncpu",
                            IF_GET_CB(int, Hqp_Docp, ncpu),
                            IF_SET_CB(int, Hqp_Docp, set_ncpu)));
}  

//-------------------------------------------------------------------------
Hqp_Docp::~Hqp_Docp()
{
  iv_free(_x_type);
  delete _xu_eq;
  delete _xu_lb;
  delete _xu_ub;
  iv_free(_nus);
  iv_free(_nxs);
  v_free(_xus_init);
  iv_free(_xus_integer);
  iv_free(_xu_start);
  iv_free(_f_start);
  iv_free(_f_lin);
  delete _cns_eq;
  delete _cns_lb;
  delete _cns_ub;
  iv_free(_cns_lin);
  iv_free(_cns_start);
  iv_free(_f0_lin);
  v_free(_f0);
  for (int tn = 0; tn < _ncpu; tn++) {
    v_free(_ck[tn]);
    m_free(_fkx[tn]);
    m_free(_fku[tn]);
    m_free(_ckx[tn]);
    m_free(_cku[tn]);
    m_free(_Lkxx[tn]);
    m_free(_Lkuu[tn]);
    m_free(_Lkxu[tn]);
    v_free(_vfk[tn]);
    v_free(_vck[tn]);
  }
  delete [] _xk;
  delete [] _uk;
  delete [] _fk;
  delete [] _s1;
  delete [] _s2;
  delete [] _ck;
  delete [] _fkx;
  delete [] _fku;
  delete [] _ckx;
  delete [] _cku;
  delete [] _Lkxx;
  delete [] _Lkuu;
  delete [] _Lkxu;
  delete [] _vfk;
  delete [] _vck;
  delete [] _xk_head;
  delete [] _uk_head;
  delete [] _fk_head;
  delete [] _s1_head;
  delete [] _s2_head;
}

//-------------------------------------------------------------------------
void Hqp_Docp::horizon(int k0, int kf)
{
  _k0 = k0;
  _kf = kf;
}

//-------------------------------------------------------------------------
void Hqp_Docp::alloc_vars(VECP v, VECP v_min, VECP v_max, IVECP v_int, int n)
{
  if ((VEC *)v != VNULL)
    v_resize(v, n);
  if ((VEC *)v_min != VNULL)
    v_resize(v_min, n);
  if ((VEC *)v_max != VNULL)
    v_resize(v_max, n);
  if ((IVEC *)v_int != IVNULL)
    iv_resize(v_int, n);
  v_set(v, 0.0);
  v_set(v_min, -Inf);
  v_set(v_max, Inf);
  iv_set(v_int, 0);
}

//-------------------------------------------------------------------------
void Hqp_Docp::setup_horizon(int &, int &)
{
  // empty default implementation, i.e.
  //   -- use preset _k0 and _kf
}

//-------------------------------------------------------------------------
void Hqp_Docp::setup_struct(int, const VECP, const VECP,
			    MATP, MATP, IVECP,
			    VECP, VECP, int &,
			    MATP, MATP, IVECP,
			    MATP, MATP, MATP)
{
  // empty default implementation means:
  //   -- matrices are full
  //   -- objective and all constraints are nonlinear
}

//-------------------------------------------------------------------------
void Hqp_Docp::init_simulation(int k, VECP x, VECP u)
{
  // empty default implementation, i.e.
  //   the control parameters and initial states of the current
  //   iterate _x are used
}

//-------------------------------------------------------------------------
// parse_constr:
//   -- split constraints up into:
//      equality constraints, lower bounds, upper bounds
//
void Hqp_Docp::parse_constr(const VECP cmin, const VECP cmax,
			    int idx, Hqp_DocpAssoc *eq,
			    Hqp_DocpAssoc *lb, Hqp_DocpAssoc *ub)
{
  assert(cmin->dim == cmax->dim);

  int i, iend = cmin->dim;

  for (i = 0; i < iend; i++) {
    if (cmin[i] == cmax[i]) {
      assert(is_finite(cmin[i]));
      eq->add_cns(cmin[i], idx);
    }
    else {
      if (cmin[i] > -Inf)
	lb->add_cns(cmin[i], idx);
      if (cmax[i] < Inf)
	ub->add_cns(cmax[i], idx);
    }
    idx ++;
  }

  eq->new_k();
  lb->new_k();
  ub->new_k();
}

//-------------------------------------------------------------------------
void Hqp_Docp::setup()
{
  // reset counter for function evaluations
  _fbd_evals = 0;

  // setup optimization horizon
  setup_horizon(_k0, _kf);

  // setup and initialize variables
  setup_x();
  init_x();

  // setup the sparsity structure
  setup_qp();
}

//-------------------------------------------------------------------------
void Hqp_Docp::setup_x()
{
  int i, k, K;
  int nsp, nfsp;
  int ncnssp;
  VECP x, x_min, x_max;
  VECP u, u_min, u_max;
  IVECP x_int, u_int;
  VECP c, c_min, c_max;

  K = _kf - _k0;

  _xu_eq->reset();
  _xu_lb->reset();
  _xu_ub->reset();

  _cns_eq->reset();
  _cns_lb->reset();
  _cns_ub->reset();

  //
  // initialize variables and constraints
  //

  v_resize(_xus_init, 0);
  iv_resize(_xus_integer, 0);
  nsp = 0;
  nfsp = 0;
  iv_resize(_f_start, K+2);	// two additional
  iv_resize(_xu_start, K+1);
  iv_resize(_nxs, K+1);
  iv_resize(_nus, K+1);

  ncnssp = 0;
  iv_resize(_cns_start, K+2);	// one additional for the final values

  x = v_get(1);
  x_min = v_get(1);
  x_max = v_get(1);
  x_int = iv_get(1);
  u = v_get(1);
  u_min = v_get(1);
  u_max = v_get(1);
  u_int = iv_get(1);
  c = v_get(1);
  c_min = v_get(1);
  c_max = v_get(1);

  for (k = 0; k <= K; k++) {
    v_resize(x, 0);
    v_resize(x_min, 0);
    v_resize(x_max, 0);
    iv_resize(x_int, 0);
    v_resize(u, 0);
    v_resize(u_min, 0);
    v_resize(u_max, 0);
    iv_resize(u_int, 0);
    v_resize(c, 0);
    v_resize(c_min, 0);
    v_resize(c_max, 0);

    setup_vars(_k0+k,
               x, x_min, x_max, x_int,
               u, u_min, u_max, u_int,
               c, c_min, c_max);
    assert(x_min->dim == x->dim && x->dim == x_max->dim);
    assert(x->dim == x_int->dim);
    assert(u_min->dim == u->dim && u->dim == u_max->dim);
    assert(u->dim == u_int->dim);
    assert(c_min->dim == c->dim && c->dim == c_max->dim);

    if (k == 0) {
      //
      // treat special state constraints here
      //  -- assume equal initial/final state for xmin[i] = xmax[i] (= Inf)
      //
      iv_resize(_x_type, x_min->dim);
      iv_zero(_x_type);
      for (i = 0; i < (int)x_min->dim; i++) {
	if (x_min[i] == x_max[i] && x_min[i] == Inf) {
	  _x_type[i] = Periodical;
	  x_min[i] = -Inf;
	}
      }
    }

    // grow sizes of vectors for initial values and integer settings
    v_expand(_xus_init, x->dim + u->dim, _granul);
    iv_expand(_xus_integer, x->dim + u->dim, _granul);

    // store start index into global _x vector
    _xu_start[k] = nsp;

    // parse state variables and bounds
    for (i = 0; i < (int)x->dim; i++) {
      _xus_init[nsp + i] = x[i];
      _xus_integer[nsp + i] = x_int[i];
    }
    parse_constr(x_min, x_max, nsp, _xu_eq, _xu_lb, _xu_ub);
    nsp += x_min->dim;
    _nxs[k] = x_min->dim;
    if (k > 0) {
      _f_start[k-1] = nfsp;
      nfsp += x_min->dim;
    }

    // parse control parameters and bounds
    for (i = 0; i < (int)u->dim; i++) {
      _xus_init[nsp + i] = u[i];
      _xus_integer[nsp + i] = u_int[i];
    }
    parse_constr(u_min, u_max, nsp, _xu_eq, _xu_lb, _xu_ub);
    nsp += u_min->dim;
    _nus[k] = u_min->dim;

    // parse constraints
    _cns_start[k] = ncnssp;
    parse_constr(c_min, c_max, ncnssp, _cns_eq, _cns_lb, _cns_ub);
    ncnssp += c_min->dim;
  }

  // final state

  assert(_nus[K] == 0);
  _f_start[K] = nfsp;
  _f_start[K+1] = nfsp;
  _cns_start[K+1] = ncnssp;

  v_free(x);
  v_free(x_min);
  v_free(x_max);
  iv_free(x_int);
  v_free(u);
  v_free(u_min);
  v_free(u_max);
  iv_free(u_int);
  v_free(c);
  v_free(c_min);
  v_free(c_max);

  // adapt dimensions

  if (nfsp != (int)_f_lin->dim)
    iv_resize(_f_lin, nfsp);

  if (ncnssp != (int)_cns_lin->dim)
    iv_resize(_cns_lin, ncnssp);

  if (nsp != (int)_x->dim)
    v_resize(_x, nsp);

  if (K+1 != (int)_f0_lin->dim)
    iv_resize(_f0_lin, K+1);

  // initialize constraint offsets for large QP

  int mesp = nfsp;
  _xu_eq->offs = mesp; mesp += _xu_eq->dim;
  _cns_eq->offs = mesp;

  int misp = 0;
  _xu_lb->offs = misp; misp += _xu_lb->dim;
  _xu_ub->offs = misp; misp += _xu_ub->dim;
  _cns_lb->offs = misp; misp += _cns_lb->dim;
  _cns_ub->offs = misp;
}

//-------------------------------------------------------------------------
void Hqp_Docp::setup_qp()
{
  int nsp;		// variables
  int mesp;		// equality constraints
  int misp;		// inequality constraints
  int i, iend, offs;
  int k;
  MATP fx, fu, cx, cu;
  MATP Lxx, Luu, Lxu;
  IVEC f_lin, c_lin;
  int K = _kf - _k0;
  int nx, nu, nf, nc;
  int kf, kxu;

  //
  // assign constraints to _qp:
  //   _qp->A, _qp->b:
  //       + system equations
  //       + _xu_eq
  //       + _cns_eq
  //   _qp->C, _qp->d:
  //       + _xu_lb
  //       + _xu_ub
  //       + _cns_lb
  //       + _cns_ub
  //

  nsp = _x->dim;

  mesp = _f_lin->dim + _xu_eq->dim + _cns_eq->dim;

  // count initial and final state constraints
  iend = min(_nxs[0], _nxs[K]);
  for (i = 0; i < iend; i++) {
    if (_x_type[i] & Periodical)
      mesp++;
  }

  misp = _xu_lb->dim + _xu_ub->dim + _cns_lb->dim + _cns_ub->dim;

  //
  // allocate _qp
  //

  _qp->resize(nsp, mesp, misp, 10, 10, 2);

  //
  // init sparsity structure and set constant values
  //

  // integers
  iv_copy(_xus_integer, _qp->x_int);

  // initial and final state constraints

  offs = _f_lin->dim + _xu_eq->dim + _cns_eq->dim;
  for (i = 0; i < iend; i++) {
    if (_x_type[i] & Periodical) {
      sp_set_val(_qp->A, offs, i, 1.0);
      sp_set_val(_qp->A, offs, _x->dim - _nxs[K] + i, -1.0);
      offs++;
    }
  }

  // equality constraints for state/control variables (e.g. initial state)

  offs = _xu_eq->offs;
  iend = _xu_eq->dim;
  for (i = 0; i < iend; i++)
    sp_set_val(_qp->A, offs + i, _xu_eq->idxs[i], 1.0);

  // inequality constraints for state/control variables

  offs = _xu_lb->offs;
  iend = _xu_lb->dim;
  for (i = 0; i < iend; i++)
    sp_set_val(_qp->C, offs + i, _xu_lb->idxs[i], 1.0);

  offs = _xu_ub->offs;
  iend = _xu_ub->dim;
  for (i = 0; i < iend; i++)
    sp_set_val(_qp->C, offs + i, _xu_ub->idxs[i], -1.0);

  // discrete--time structure

  fx = m_get(1, 1);
  fu = m_get(1, 1);
  cx = m_get(1, 1);
  cu = m_get(1, 1);
  Lxx = m_get(1, 1);
  Luu = m_get(1, 1);
  Lxu = m_get(1, 1);

  kf = 0;
  kxu = 0;
  for (k = 0; k <= K; k++) {
    nx = _nxs[k];
    nu = _nus[k];
    nf = _f_start[k+1] - _f_start[k];
    nc = _cns_start[k+1] - _cns_start[k];
    v_part(_x, kxu, nx, _xk[0]);
    v_part(_x, kxu+nx, nu, _uk[0]);
    m_resize(fx, nf, nx);
    m_resize(fu, nf, nu);
    iv_part(_f_lin, _f_start[k], nf, &f_lin);
    m_resize(cx, nc, nx);
    m_resize(cu, nc, nu);
    iv_part(_cns_lin, _cns_start[k], nc, &c_lin);
    m_resize(Lxx, nx, nx);
    m_resize(Luu, nu, nu);
    m_resize(Lxu, nx, nu);

    v_part(_qp->c, kxu, nx, _s1[0]);	// f0x
    v_part(_qp->c, kxu+nx, nu, _s2[0]);	// f0u
    v_ones(_s1[0]);
    v_ones(_s2[0]);
    _f0_lin[k] = 0;

    m_ones(fx);
    m_ones(fu);
    iv_set(&f_lin, 0);

    m_ones(cx);
    m_ones(cu);
    iv_set(&c_lin, 0);

    m_ones(Lxx);
    m_ones(Luu);
    m_ones(Lxu);
  
    setup_struct(_k0+k, _xk[0], _uk[0],
		 fx, fu, &f_lin,
		 _s1[0], _s2[0], _f0_lin[k],
		 cx, cu, &c_lin,
		 Lxx, Luu, Lxu);

    sp_insert_mat(_qp->A, kf, kxu, fx);
    sp_insert_mat(_qp->A, kf, kxu+nx, fu);
    for (i = 0; i < nf; i++)
      sp_set_val(_qp->A, kf + i, kxu+nx+nu + i, -1.0);
    
    symsp_insert_symmat(_qp->Q, kxu, Lxx);
    symsp_insert_symmat(_qp->Q, kxu+nx, Luu);
    sp_insert_mat(_qp->Q, kxu, kxu+nx, Lxu);

    _cns_eq->sp_insert_mat(k, _qp->A, kxu, cx, _cns_start[k]);
    _cns_eq->sp_insert_mat(k, _qp->A, kxu+nx, cu, _cns_start[k]);

    _cns_lb->sp_insert_mat(k, _qp->C, kxu, cx, _cns_start[k]);
    _cns_lb->sp_insert_mat(k, _qp->C, kxu+nx, cu, _cns_start[k]);

    sm_mlt(-1.0, cx, cx);
    sm_mlt(-1.0, cu, cu);
    _cns_ub->sp_insert_mat(k, _qp->C, kxu, cx, _cns_start[k]);
    _cns_ub->sp_insert_mat(k, _qp->C, kxu+nx, cu, _cns_start[k]);

    kf += nf;
    kxu += nx+nu;
  }

  // no initial guess
  sp_zero(_qp->Q);

  m_free(fx);
  m_free(fu);
  m_free(cx);
  m_free(cu);
  m_free(Lxx);
  m_free(Luu);
  m_free(Lxu);
}

// for now just repeat setup of variables and update bounds
//-------------------------------------------------------------------------
void Hqp_Docp::reinit_bd()
{
  // store back current dimensions for check at end of this method
  int x_dim = _x->dim;
  int xu_eq_dim = _xu_eq->dim;
  int xu_lb_dim = _xu_lb->dim;
  int xu_ub_dim = _xu_ub->dim;
  int cns_eq_dim = _cns_eq->dim;
  int cns_lb_dim = _cns_lb->dim;
  int cns_ub_dim = _cns_ub->dim;

  setup_x();
  
  // check that dimenstions did not change
  if (x_dim != _x->dim
      || xu_eq_dim != _xu_eq->dim
      || xu_lb_dim != _xu_lb->dim
      || xu_ub_dim != _xu_ub->dim
      || cns_eq_dim != _cns_eq->dim
      || cns_lb_dim != _cns_lb->dim
      || cns_ub_dim != _cns_ub->dim) {
    m_error(E_SIZES, "Hqp_Docp::reinit_bd");
  }

  update_bounds();
}

//-------------------------------------------------------------------------
void Hqp_Docp::init_x()
{
  v_copy(_xus_init, _x);
}

//-------------------------------------------------------------------------
void Hqp_Docp::simulate()
{
  int k, kxu;
  int nx, nu, nf, nc;
  Real f0k;
  VECP c;
  int K = _kf - _k0;

  //
  // for k = 0, ... , K
  //  -- call init_solution() for xk and uk
  //  -- simulate for x{k+1}
  //  -- sum up all f0k's already here
  //

  c = v_get(1);
  kxu = 0;
  _f = 0.0;
  for (k = 0; k <= K; k++) {
    nx = _nxs[k];
    nu = _nus[k];
    nf = _f_start[k+1] - _f_start[k];
    nc = _cns_start[k+1] - _cns_start[k];
    v_part(_x, kxu, nx, _xk[0]);
    v_part(_x, kxu+nx, nu, _uk[0]);
    v_part(_x, kxu+nx+nu, nf, _fk[0]);
    init_simulation(_k0+k, _xk[0], _uk[0]);
    f0k = 0.0;
    c = v_resize(c, nc);
    update_vals(_k0+k, _xk[0], _uk[0], _fk[0], f0k, c);
    _f += f0k;
    kxu += nx+nu;
  }

  v_free(c);
}

//-------------------------------------------------------------------------
void Hqp_Docp::update_fbd()
{
  int k, K = _kf - _k0;

  v_zero(_f0);

  #pragma omp parallel for num_threads(_ncpu)
  for (k = 0; k <= K; k++) {
    int tn = omp_get_thread_num();
    VECP xk = _xk[tn];
    VECP uk = _uk[tn];
    VECP fk = _fk[tn];
    VECP ck = _ck[tn];
    double f0k = 0.0;
    int nx = _nxs[k];
    int nu = _nus[k];
    int kxu = _xu_start[k];
    int kf = _f_start[k];
    int nf = _f_start[k+1] - kf;
    int nc = _cns_start[k+1] - _cns_start[k];

    v_part(_x, kxu, nx, xk);
    v_part(_x, kxu+nx, nu, uk);
    v_part(_qp->b, kf, nf, fk);
    v_resize(ck, nc);
    v_zero(fk);

    update_vals(_k0+k, xk, uk, fk, f0k, ck);

    _f0[tn] += f0k;
    v_sub(fk, v_part(_x, kxu+nx+nu, nf, xk), fk);

    // break if f0 is not finite
    //if (!is_finite(_f0[k]))
    //  break;

    int offs = _cns_start[k];
    Real *v_ve = _qp->b->ve + _cns_eq->offs;
    int i, iend = _cns_eq->start[k+1];
    for (i = _cns_eq->start[k]; i < iend; i++)
      v_ve[i] = ck[_cns_eq->idxs[i]-offs] - _cns_eq->vals[i];

    v_ve = _qp->d->ve + _cns_lb->offs;
    iend = _cns_lb->start[k+1];
    for (i = _cns_lb->start[k]; i < iend; i++)
      v_ve[i] = ck[_cns_lb->idxs[i]-offs] - _cns_lb->vals[i];

    v_ve = _qp->d->ve + _cns_ub->offs;
    iend = _cns_ub->start[k+1];
    for (i = _cns_ub->start[k]; i < iend; i++)
      v_ve[i] = _cns_ub->vals[i] - ck[_cns_ub->idxs[i]-offs];
  }

  // take over calculated objective value
  _f = v_sum(_f0);

  // update state/control bounds and initial/final states
  update_bounds();

  _fbd_evals++;
}

//-------------------------------------------------------------------------
void Hqp_Docp::update_bounds()
{
  int i, iend, offs;
  int K = _kf - _k0;
  Real *x_ve, *v_ve;
  int *idxs_ive;
  VEC v_head;
  VECP v = &v_head;

  x_ve = _x->ve;

  v_part(_qp->b, _xu_eq->offs, _xu_eq->dim, v);
  v_ve = v->ve;
  idxs_ive = _xu_eq->idxs->ive;
  iend = _xu_eq->dim;
  for (i = 0; i < iend; i++)
    v_ve[i] = x_ve[idxs_ive[i]];
  v_sub(v, _xu_eq->vals, v);

  v_part(_qp->d, _xu_lb->offs, _xu_lb->dim, v);
  v_ve = v->ve;
  idxs_ive = _xu_lb->idxs->ive;
  iend = _xu_lb->dim;
  for (i = 0; i < iend; i++)
    v_ve[i] = x_ve[idxs_ive[i]];
  v_sub(v, _xu_lb->vals, v);

  v_part(_qp->d, _xu_ub->offs, _xu_ub->dim, v);
  v_ve = v->ve;
  idxs_ive = _xu_ub->idxs->ive;
  iend = _xu_ub->dim;
  for (i = 0; i < iend; i++)
    v_ve[i] = -x_ve[idxs_ive[i]];
  v_add(v, _xu_ub->vals, v);

  // update initial/final state constraints

  iend = min(_nxs[0], _nxs[K]);
  v_part(_x, 0, _nxs[0], _xk[0]);
  v_part(_x, _x->dim - _nxs[K], _nxs[K], _fk[0]);
  offs = _f_lin->dim + _xu_eq->dim + _cns_eq->dim;
  for (i = 0; i < iend; i++) {
    if (_x_type[i] & Periodical) {
      _qp->b->ve[offs] = _xk[0][i] - _fk[0][i];
      offs++;
    }
  }
}

//-------------------------------------------------------------------------
void Hqp_Docp::update(const VECP y, const VECP z)
{
  int k, K = _kf - _k0;
  SPMAT *A = _qp->A;

  if ((const VEC *)y == NULL || (const VEC *)z == NULL)
    m_error(E_NULL, "Hqp_Docp::update");
  if (y->dim != _qp->b->dim || z->dim != _qp->d->dim)
    m_error(E_SIZES, "Hqp_Docp::update");
  v_zero(_f0);

  #pragma omp parallel for num_threads(_ncpu)
  for (k = 0; k <= K; k++) {
    int tn = omp_get_thread_num();
    VECP xk = _xk[tn];
    VECP uk = _uk[tn];
    VECP fk = _fk[tn];
    VECP ck = _ck[tn];
    double f0k = 0.0;
    MATP fx = _fkx[tn];
    MATP fu = _fku[tn];
    MATP cx = _ckx[tn];
    MATP cu = _cku[tn];
    MATP Lxx = _Lkxx[tn];
    MATP Luu = _Lkuu[tn];
    MATP Lxu = _Lkxu[tn];
    VECP vf = _vfk[tn];
    VECP vc = _vck[tn];
    VECP s1 = _s1[tn];
    VECP s2 = _s2[tn];
    int nx = _nxs[k];
    int nu = _nus[k];
    int kxu = _xu_start[k];
    int kf = _f_start[k];
    int nf = _f_start[k+1] - kf;
    int nc = _cns_start[k+1] - _cns_start[k];
    int i;

    v_part(_x, kxu, nx, xk);
    v_part(_x, kxu+nx, nu, uk);
    v_part(_qp->b, kf, nf, fk);
    v_resize(ck, nc);

    v_part(_qp->c, kxu, nx, s1);	// f0x
    v_part(_qp->c, kxu+nx, nu, s2);	// f0u
    m_resize(fx, nf, nx);
    m_resize(fu, nf, nu);
    m_resize(cx, nc, nx);
    m_resize(cu, nc, nu);
    m_resize(Lxx, nx, nx);
    m_resize(Luu, nu, nu);
    m_resize(Lxu, nx, nu);
    v_resize(vc, nc);
    v_resize(vf, nf);

    // prepare value update
    v_zero(fk);

    // multiplier for system equations
    for (i = 0; i < nf; i++)
      vf[i] = y[kf + i];

    // multiplier for general constraints
    v_zero(vc);
    _cns_ub->sv_mltadd_vec(k, -1.0, z, vc, _cns_start[k]);
    _cns_lb->sv_mltadd_vec(k, 1.0, z, vc, _cns_start[k]);
    _cns_eq->sv_mltadd_vec(k, 1.0, y, vc, _cns_start[k]);

    // extract the current Lagrangian Hessian for optional user update
    sp_extract_mat(_qp->Q, kxu, kxu, Lxx);
    sp_extract_mat(_qp->Q, kxu, kxu+nx, Lxu);
    sp_extract_mat(_qp->Q, kxu+nx, kxu+nx, Luu);

    update_stage(_k0+k, xk, uk,
		 fk, f0k, ck,
		 fx, fu, s1, s2, cx, cu,
		 vf, vc, Lxx, Luu, Lxu);

    // value update

    _f0[tn] += f0k;
    v_sub(fk, v_part(_x, kxu+nx+nu, nf, xk), fk);

    // break if _f is not finite
    //if (!is_finite(_f0[k]))
    //  break;

    int offs = _cns_start[k];
    Real *v_ve = _qp->b->ve + _cns_eq->offs;
    int iend = _cns_eq->start[k+1];
    for (i = _cns_eq->start[k]; i < iend; i++)
      v_ve[i] = ck[_cns_eq->idxs[i]-offs] - _cns_eq->vals[i];

    v_ve = _qp->d->ve + _cns_lb->offs;
    iend = _cns_lb->start[k+1];
    for (i = _cns_lb->start[k]; i < iend; i++)
      v_ve[i] = ck[_cns_lb->idxs[i]-offs] - _cns_lb->vals[i];

    v_ve = _qp->d->ve + _cns_ub->offs;
    iend = _cns_ub->start[k+1];
    for (i = _cns_ub->start[k]; i < iend; i++)
      v_ve[i] = _cns_ub->vals[i] - ck[_cns_ub->idxs[i]-offs];

    // derivatives update

    for (i = 0; i < nf; i++)
      if (!_f_lin[kf+i]) {
	sp_update_mrow(A, kf, kxu, fx, i);
        sp_update_mrow(A, kf, kxu+nx, fu, i);
      }

    _cns_eq->sp_update_mat(k, _qp->A, kxu, cx, _cns_start[k], _cns_lin);
    _cns_lb->sp_update_mat(k, _qp->C, kxu, cx, _cns_start[k], _cns_lin);
    sm_mlt(-1.0, cx, cx);
    _cns_ub->sp_update_mat(k, _qp->C, kxu, cx, _cns_start[k], _cns_lin);

    _cns_eq->sp_update_mat(k, _qp->A, kxu+nx, cu, _cns_start[k], _cns_lin);
    _cns_lb->sp_update_mat(k, _qp->C, kxu+nx, cu, _cns_start[k], _cns_lin);
    sm_mlt(-1.0, cu, cu);
    _cns_ub->sp_update_mat(k, _qp->C, kxu+nx, cu, _cns_start[k], _cns_lin);

    sp_update_mat(_qp->Q, kxu, kxu, Lxx);
    sp_update_mat(_qp->Q, kxu+nx, kxu+nx, Luu);
    sp_update_mat(_qp->Q, kxu, kxu+nx, Lxu);
  }

  // take over calculated objective value
  _f = v_sum(_f0);

  // update state/control bounds and initial/final states
  update_bounds();
}

//-------------------------------------------------------------------------
// update_stage:
//  -- default implementation calls update_vals(),
//     update_grds(), and update_hela()
//
void Hqp_Docp::update_stage(int k, const VECP x, const VECP u,
			    VECP f, Real &f0, VECP c,
			    MATP fx, MATP fu, VECP f0x, VECP f0u,
			    MATP cx, MATP cu,
			    const VECP vf, const VECP vc,
			    MATP Lxx, MATP Luu, MATP Lxu)
{
  update_vals(k, x, u, f, f0, c);
  update_grds(k, x, u, fx, fu, f0x, f0u, cx, cu);
  update_hela(k, x, u, vf, vc, Lxx, Luu, Lxu);
}

//-------------------------------------------------------------------------
// update_grds:
//  -- default implementation approximates finite differences
//
void Hqp_Docp::update_grds(int k, const VECP x, const VECP u,
			   MATP fx, MATP fu, VECP f0x, VECP f0u,
			   MATP cx, MATP cu)
{
  int i, j;
  int nx, nu, nf, nc;
  Real vi_bak, dvi;
  Real f0, df0;
  VECP f, df;
  VECP c, dc;

  nx = x->dim;
  nu = u->dim;
  nf = fx->m;
  nc = cx->m;

  c = v_get(1);
  dc = v_get(1);
  f = v_get(1);
  df = v_get(1);

  v_resize(c, nc);
  v_resize(dc, nc);
  v_resize(f, nf);
  v_resize(df, nf);

  v_zero(f);
  f0 = 0.0;
  v_zero(c);

  update_vals(k, x, u, f, f0, c);

  for (i = 0; i < nx; i++) {
    vi_bak = x->ve[i];
    dvi = 1e-4 * fabs(vi_bak) + 1e-6;
    x->ve[i] += dvi;

    v_zero(df);
    df0 = 0.0;
    v_zero(dc);

    update_vals(k, x, u, df, df0, dc);

    v_sub(df, f, df);
    for (j = 0; j < nf; j++)
      fx->me[j][i] = df->ve[j] / dvi;

    f0x->ve[i] = (df0 - f0) / dvi;

    v_sub(dc, c, dc);
    for (j = 0; j < nc; j++)
      cx->me[j][i] = dc->ve[j] / dvi;

    x->ve[i] = vi_bak;
  }
  for (i = 0; i < nu; i++) {
    vi_bak = u->ve[i];
    dvi = 1e-4 * fabs(vi_bak) + 1e-6;
    u->ve[i] += dvi;

    v_zero(df);
    df0 = 0.0;
    v_zero(dc);

    update_vals(k, x, u, df, df0, dc);

    v_sub(df, f, df);
    for (j = 0; j < nf; j++)
      fu->me[j][i] = df->ve[j] / dvi;

    f0u->ve[i] = (df0 - f0) / dvi;

    v_sub(dc, c, dc);
    for (j = 0; j < nc; j++)
      cu->me[j][i] = dc->ve[j] / dvi;

    u->ve[i] = vi_bak;
  }

  v_free(f);
  v_free(df);
  v_free(c);
  v_free(dc);
}

//-------------------------------------------------------------------------
void Hqp_Docp::update_hela(int, const VECP, const VECP,
			   const VECP, const VECP,
			   MATP, MATP, MATP)
{
  // empty default implementation, i.e.
  //   update provided by SQP solver must be used
}


//=========================================================================
