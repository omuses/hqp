/*
 * Hqp_IpLQDOCP.C -- class definition
 *
 * E. Arnold 09/24/96
 *           06/30/97  copy _A_ori, _C_ori        
 *           07/08/97  bugs in resize(), free()
 *           08/20/99  Hqp_IpLQDOCP::Get_Dim()
 *           2002-04-17 free() replaced by myfree()
 *           2004-10-13 bug fix for assert()
 *
 */

/*
    Copyright (C) 1996--2014  Eckhard Arnold

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


#include "Hqp_IpLQDOCP.h"

#include <assert.h>
#include <math.h>

extern "C" {
#include <meschach/matrix2.h>
#include <meschach/sparse2.h>
}

#include <If_Int.h>
#include <If_Real.h>
// #include <If_Method.h>

#include "Hqp_Program.h"

IF_CLASS_DEFINE("LQDOCP", Hqp_IpLQDOCP, Hqp_IpMatrix);
// typedef If_Method<Hqp_IpLQDOCP> If_Cmd;

#undef LOG
#ifdef LOG
#define pr_(s)     printf(s)
#define pr(s,x)    printf(s,x)
#else
#define pr_(s)    
#define pr(s,x)   
#endif

//--------------------------------------------------------------------------
//   Calculation of matrix C^T D^{-1} C using sparsity structure of 
//   matrix C. Only the non-zero entries of CTC are visited.
//   CTC(i,j) = CT(i,:)*(yny(:).*CT(:,j)
//   E. Arnold  10/07/96
//              06/30/97 avoid diag access
//--------------------------------------------------------------------------
SPMAT *CTDC(const SPMAT *CT, const VEC *yny, SPMAT *CTC)
{
  int i;
  int j;
  SPROW *r;
  row_elt *elt;
  Real sum;

  if ( ( CT == SMNULL ) || ( yny == VNULL) || ( CTC == SMNULL ) )
    m_error(E_NULL, "CTDC");
  if ( ( CT->n != (int) yny->dim ) || ( CT->m != CTC->m ) || 
       ( CTC->m != CTC->n ) )
    m_error(E_SIZES, "CTDC");
  //  if ( !CTC->flag_diag )
  //    sp_diag_access(CTC);

  r = CTC->row;
  for ( i = 0; i < CTC->m; i++, r++ ) {
    if ( r->len > 0 ) {
      //      elt = r->elt+r->diag;
      elt = r->elt;
      //      for ( j = r->diag; j < r->len; j++, elt++ ) {
      for ( j = 0; j < r->len; j++, elt++ ) {
	if (elt->col < i )
	  continue;
	sum = sprow_inprod(CT->row+i, yny, CT->row+elt->col);
	sp_set_val(CTC, i, elt->col, sum); 
      }
    }
  }

  if ( !CTC->flag_diag )
    sp_diag_access(CTC);

  return CTC;
}
 
//--------------------------------------------------------------------------
Hqp_IpLQDOCP::Hqp_IpLQDOCP()
{

  //  _wz_tol = 1.0e6;
  //  _wz_tol = 1.0e8;
  _wz_tol = HUGE_VAL;
  //  _ge_tol = 1.0e-8;
  _ge_tol = 1.0e-6;
  //  _a_sparse = 1;
  _a_sparse = 0;

  _kmax = -1;
  Gxx = MNULL;
  m1 = MNULL;
  m2 = MNULL;
  Gx = VNULL;
  Gu = VNULL;
  ax = MNULL;
  au = MNULL;
  yny = VNULL;
  v1 = VNULL;
  v2 = VNULL;
  ct = VNULL;
  yb = VNULL;
  Z = MNULL;
  pivot = PNULL;
  blocks = PNULL;
  _ra = IVNULL;
  _rc = IVNULL;
  _rak = IVNULL;
  _rck = IVNULL;
  _rcka = IVNULL;
  _rca = IVNULL;
  _Q_ori = SMNULL;
  _A_ori = SMNULL;
  _C_ori = SMNULL;
  _Z_ori = VNULL;
  _W_ori = VNULL;
  _WZ_ori = VNULL;
  _CT = SMNULL;
  _CTC = SMNULL;
  _AT = SMNULL;

  db = MNULL;
  d22 = MNULL;

  //   residuum  02/21/97
  _res1 = VNULL;
  _res2 = VNULL;
  _res3 = VNULL;
  _res4 = VNULL;

  //   different number of variables per stage  02/22/97
  _nk = IVNULL;
  _mk = IVNULL;
  _nmk = IVNULL;
  _nks = IVNULL;
  _indk = IVNULL;

  //   diagonal scaling
  scx0 = VNULL;

  pr_("\nLQDOCP::Hqp_IpLQDOCP");
  //  mem_info();

  //  _sbw = -1;
  //  _ifList.append(new If_Int("mat_sbw", &_sbw));
  //  _ifList.append(new If_Real("mat_tol", &_tol));

  _logging = 0;

  _ifList.append(new If_Real("mat_wz_tol", &_wz_tol));
  _ifList.append(new If_Int("mat_a_sparse", &_a_sparse));
  _ifList.append(new If_Int("mat_logging", &_logging));
}

//--------------------------------------------------------------------------
Hqp_IpLQDOCP::~Hqp_IpLQDOCP()
{
  myfree();
  pr_("\nLQDOCP::~Hqp_IpLQDOCP");
}

//--------------------------------------------------------------------------
//   Get DOCP parameters from sparsity structure of qp->A.
//   _kmax: number of stages
//   _nk:   number of states of the stage
//   _mk:   number of controls of the stage
//   _nmk:  number of states and controls up to the stage
//   _nks:  number of states up to the stage
//   _indk: stage of qp variable
//
//   E. Arnold   10/07/96
//               02/22/97   _nk, _mk, _nmk, _indk
//--------------------------------------------------------------------------
int Hqp_IpLQDOCP::Get_Dim(const Hqp_Program *qp)
{
  int i, k, nk, il, icl, icl1, icl_;

  if ( qp->A->m == 0 )
    return 0;

  _kmax = 0;
  nk = 0;
  icl_ = -1;
  for ( i = 0; i < qp->A->m; i++ ) {
    il = qp->A->row[i].len;
    if ( il <= 1 )
      return 0;
    if ( qp->A->row[i].elt[il-1].val != -1.0 )
      return 0;
    icl = qp->A->row[i].elt[il-1].col;
    icl1 = qp->A->row[i].elt[il-2].col;
    if ( icl <= icl_ )
      return 0;
    if ( ( icl - icl_ > 1 ) || ( icl - icl1 < nk ) ) {
      _kmax++;
      nk = 1;
    } else 
      nk++;
    icl_ = icl;
    if ( icl == qp->A->n-1 )
      break;
    if ( i == qp->A->m-1 )
      return 0;
  }
  if ( _kmax == 0 ) 
    return 0;
  _kmax1 = _kmax-1;

  _nk = iv_resize(_nk, _kmax+1);
  _nk = iv_zero(_nk);
  _mk = iv_resize(_mk, _kmax);
  _mk = iv_zero(_mk);
  _nmk = iv_resize(_nmk, _kmax+1);
  _nmk = iv_zero(_nmk);
  _nks = iv_resize(_nks, _kmax+1);
  _nks = iv_zero(_nks);
  _indk = iv_resize(_indk, qp->A->n);
  _indk = iv_zero(_indk);

  k = 0;
  nk = 0;
  icl_ = -1;
  for ( i = 0; i < qp->A->m; i++ ) {
    il = qp->A->row[i].len;
    icl = qp->A->row[i].elt[il-1].col;
    icl1 = qp->A->row[i].elt[il-2].col;
    if ( ( icl - icl_ > 1 ) || ( icl - icl1 < nk ) ) {
      if ( k > 0 ) 
	iv_set_val(_nk, k, nk);
      k++;
      iv_set_val(_nmk, k, icl);
      nk = 1;
    } else 
      nk++;
    icl_ = icl;
    if ( icl == qp->A->n-1 ) {
      iv_set_val(_nk, k, nk);
      break;
    }
  }
  iv_set_val(_nk, 0, min(iv_entry(_nk, 1), iv_entry(_nmk, 1))); // 1999-08-20
  for ( k = 0; k <= _kmax; k++ ) 
    if ( ( iv_entry(_nk, k) < 0 ) || ( iv_entry(_nmk, k) < 0 ) )
      return 0;
  for ( k = 0; k < _kmax; k++ ) {
    iv_set_val(_mk, k, iv_entry(_nmk, k+1)-iv_entry(_nmk, k)-iv_entry(_nk, k));
    if ( iv_entry(_mk, k) < 0 ) 
      return 0;
  }
  for ( k = 1; k <= _kmax; k++ )
    iv_set_val(_nks, k, iv_entry(_nks, k-1)+iv_entry(_nk, k));

  for ( k = 0; k < _kmax; k++ )
    for ( i = iv_entry(_nmk, k); i < iv_entry(_nmk, k+1); i++ )
      iv_set_val(_indk, i, k);
  for ( i = iv_entry(_nmk, _kmax); i < (int) _indk->dim; i++ )
    iv_set_val(_indk, i, _kmax);

  return 1;
}

//--------------------------------------------------------------------------
//   Check structure of matrices qp->Q, qp->A, qp->C.
//   Returns ok = 1, if structure fits into previously determined
//   DOCP structure with _kmax, _n, _m, _ra, _rc.
//   Checks for  fixed initial state.
//
//   E. Arnold  10/07/96
//              02/22/97  _nk, _mk, _nmk, _indk, initial state
//--------------------------------------------------------------------------
int Hqp_IpLQDOCP::Check_Structure(const Hqp_Program *qp)
{
  int i, k, k_, ok = 1;
 
  //   check DOCP structure of matrix A   
  k_ = 0;
  for (i = 0; i < _nks[_kmax]; i++) { 
    ok = ok && (qp->A->row[i].len >= 2) &&
      (qp->A->row[i].elt[qp->A->row[i].len-1].val == -1.0);
    if ( !ok ) return ok;
    k = iv_entry(_indk, qp->A->row[i].elt[0].col);
    ok = ok && ( ( k == k_ ) || ( k == k_+1 ) ) &&
      ( iv_entry(_indk, qp->A->row[i].elt[qp->A->row[i].len-2].col) == k ) && 
      ( iv_entry(_indk, qp->A->row[i].elt[qp->A->row[i].len-1].col) == k+1 );
    k_ = k;
  }
  for (i = _nks[_kmax]; i < qp->A->m; i++) {
    ok = ok && (qp->A->row[i].len >= 1);
    if ( !ok ) return ok;
    k = iv_entry(_ra, i-_nks[_kmax]);
    ok = ok &&
      ( iv_entry(_indk, qp->A->row[i].elt[0].col) == k ) &&
      ( iv_entry(_indk, qp->A->row[i].elt[qp->A->row[i].len-1].col) == k);
  }

  //   check DOCP structure of matrix C   
  for (i = 0; i < qp->C->m; i++) {
    ok = ok && (qp->C->row[i].len >= 1);
    if ( !ok ) return ok;
    k = iv_entry(_rc, i);
    ok = ok &&
      ( iv_entry(_indk, qp->C->row[i].elt[0].col) == k ) && 
      ( iv_entry(_indk, qp->C->row[i].elt[qp->C->row[i].len-1].col) == k );
  }

  //   check DOCP structure of matrix Q   
  for (i = 0; i < qp->Q->m; i++) {
    ok = ok && (qp->Q->row[i].len >= 1);
    if ( !ok ) return ok;
    k = iv_entry(_indk, i);
    ok = ok &&
      ( iv_entry(_indk, qp->Q->row[i].elt[0].col) == k ) &&
      ( iv_entry(_indk, qp->Q->row[i].elt[qp->Q->row[i].len-1].col) == k );
  }

  //   check for fixed initial state
  _fixed_x0 = ( (int) _raki[0]->dim == iv_entry(_nk, 0) );
  if ( _fixed_x0 )
    for ( i = 0; i < (int) _raki[0]->dim; i++ ) {
      _fixed_x0 = _fixed_x0 && 
	( qp->A->row[iv_entry(_raki[0], i)].len == 1 ) && 
	( qp->A->row[iv_entry(_raki[0], i)].elt[0].col == i ) &&
	( qp->A->row[iv_entry(_raki[0], i)].elt[0].val == 1.0 );
    }
  
  return ok;
}

//--------------------------------------------------------------------------
//   Get DOCP structure of constraints. 
//   _ra:   stage of equality constraint
//   _rak:  number of equality constraints of the stage
//   _raki: indices of the equality constraints of the stage
//   _rc:   stage of inequality constraint
//   _rck:  number of inequality constraints of the stage
//   _rcki: indices of the inequality constraints of the stage
//
//   E. Arnold  10/07/96
//              02/22/97  _nk, _indk
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::Get_Constr_Dim(const Hqp_Program *qp)
{
  int i, j, k;

  //   get time step for equality constraints   
  _ra = iv_resize(_ra, qp->A->m-_nks[_kmax]);
  _rak = iv_resize(_rak, _kmax+1);
  _rak = iv_zero(_rak);
  for ( i = _nks[_kmax]; i < qp->A->m; i++ ) {
    k = iv_entry(_indk, qp->A->row[i].elt[0].col);
    iv_set_val(_ra , i-_nks[_kmax], k);
    iv_add_val(_rak, k, 1);
  }
  for ( k = 0; k <= _kmax; k++ ) {
    _raki[k] = iv_resize(_raki[k], iv_entry(_rak, k));
    for ( i = _nks[_kmax], j = 0; i < qp->A->m; i++ ) 
      if ( iv_entry(_ra, i-_nks[_kmax]) == k ) { 
	iv_set_val(_raki[k], j, i);
	j++;
      }
  }

  //   get time step for inequality constraints   
  _rc = iv_resize(_rc, qp->C->m);
  _rck = iv_resize(_rck, _kmax+1);
  _rck = iv_zero(_rck);
  for ( i = 0; i < qp->C->m; i++ ) {
    k = iv_entry(_indk, qp->C->row[i].elt[0].col);
    iv_set_val(_rc, i, k);
    iv_add_val(_rck, k, 1);
  }
  for ( k = 0; k <= _kmax; k++ ) {
    _rcki[k] = iv_resize(_rcki[k], iv_entry(_rck, k));
    for ( i = j = 0; i < qp->C->m; i++ ) 
      if ( iv_entry(_rc, i) == k ) {
	iv_set_val(_rcki[k], j, i);
	j++;
      }
  }
}

//--------------------------------------------------------------------------
//   Allocate memory.
//
//   E. Arnold  10/07/96
//              02/22/97  _nk, _mk
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::resize()
{
  assert(_kmax >= 1);

  Vxx.resize(_kmax);
  Gxu.resize(_kmax1);
  Guu.resize(_kmax1);
  CH_Guu.resize(_kmax1);
  CH_Guu_p.resize(_kmax1);
  CH_Guu_b.resize(_kmax1);
  Vx.resize(_kmax);
  cb.resize(_kmax);
  S.resize(_kmax1);
  cbx.resize(_kmax);
  ctx.resize(_kmax1);
  x.resize(_kmax);
  u.resize(_kmax1);
  p.resize(_kmax1);
  gx.resize(_kmax);
  gu.resize(_kmax1);
  f.resize(_kmax1);
  Ru.resize(_kmax1);
  if ( _a_sparse ) {
    fx.tfree();
    fu.tfree();
  } else {
    fx.resize(_kmax1);  
    fu.resize(_kmax1);
  }
  Rux.resize(_kmax1);
  ax = m_resize(ax, 1, 1);
  au = m_resize(au, 1, 1);
  a.resize(_kmax);
  y.resize(_kmax);
  Ryx.resize(_kmax1);
  Ry.resize(_kmax1);
  PQ.resize(_kmax1);

  m1 = m_resize(m1, 1, 1);
  m2 = m_resize(m2, 1, 1);
  Gxx = m_resize(Gxx, 1, 1);
  v1 = v_resize(v1, 1);
  v2 = v_resize(v2, 1);
  Gx = v_resize(Gx, 1);
  Gu = v_resize(Gu, 1);
  ct = v_resize(ct, 1);
  yb = v_resize(yb, 1);
  Z = m_resize(Z, 1, 1);

  _raki.resize(_kmax);
  _rcki.resize(_kmax);
  _rckai.resize(_kmax);

  _Q_ori = SMNULL;
  if ( _a_sparse )     //  06/30/97  copy _A_ori, _C_ori
    SP_FREE(_A_ori);
  else
    _A_ori = SMNULL;
  SP_FREE(_C_ori);
  _Z_ori = v_resize(_Z_ori, 1);
  _W_ori = v_resize(_W_ori, 1);
  _WZ_ori = v_resize(_WZ_ori, 1);
  SP_FREE(_CT);
  SP_FREE(_CTC);
  SP_FREE(_AT);

  yny = v_resize(yny, 1);
  db = m_resize(db, 1, 1);
  d11.resize(_kmax1);
  d12.resize(_kmax1);
  d22 = m_resize(d22, 1, 1);
  Ruyb.resize(_kmax1);
  Ryyb.resize(_kmax1);
  CD.resize(_kmax1);
  CD_p.resize(_kmax1);

  //   diagonal scaling
  sc.resize(_kmax1);
}

//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::myfree()
{
  Vxx.tfree();
  Gxu.tfree();
  Guu.tfree();
  CH_Guu.tfree();
  CH_Guu_p.tfree();
  CH_Guu_b.tfree();
  Vx.tfree();
  cb.tfree();
  S.tfree();
  cbx.tfree();
  ctx.tfree();
  x.tfree();
  u.tfree();
  p.tfree();
  gx.tfree();
  gu.tfree();
  f.tfree();
  Ru.tfree();
  fx.tfree();  
  fu.tfree();
  Rux.tfree();
  M_FREE(ax);
  M_FREE(au);
  a.tfree();
  y.tfree();
  Ryx.tfree();
  Ry.tfree();
  PQ.tfree();

  M_FREE(m1);
  M_FREE(m2);
  M_FREE(Gxx);
  V_FREE(v1);
  V_FREE(v2);
  V_FREE(Gx);
  V_FREE(Gu);
  V_FREE(ct);
  V_FREE(yb);
  M_FREE(Z);
  PX_FREE(pivot);
  PX_FREE(blocks);

  IV_FREE(_ra);
  IV_FREE(_rc);
  IV_FREE(_rak);
  IV_FREE(_rck);
  IV_FREE(_rcka);
  IV_FREE(_rca);
  _raki.tfree();
  _rcki.tfree();
  _rckai.tfree();

  //   _Q_ori is a reference!
  _Q_ori = SMNULL;
  if ( _a_sparse )  //   06/30/97  copy _A_ori, _C_ori
    SP_FREE(_A_ori);
  else
    _A_ori = SMNULL;
  SP_FREE(_C_ori);
  V_FREE(_Z_ori);
  V_FREE(_W_ori);
  V_FREE(_WZ_ori);
  SP_FREE(_CT);
  SP_FREE(_CTC);
  SP_FREE(_AT);

  V_FREE(yny);
  M_FREE(db);
  d11.tfree();
  d12.tfree();
  M_FREE(d22);
  Ruyb.tfree();
  Ryyb.tfree();
  CD.tfree();
  CD_p.tfree();

  //   residuum  02/21/97
  V_FREE(_res1);
  V_FREE(_res2);
  V_FREE(_res3);
  V_FREE(_res4);

  //   different number of variables per stage  02/22/97
  IV_FREE(_nk);
  IV_FREE(_mk);
  IV_FREE(_nmk);
  IV_FREE(_nks);
  IV_FREE(_indk);

  //   diagonal scaling
  V_FREE(scx0);
  sc.tfree();
}

//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::dump(const char *fname)
{
  int k;
  FILE *fp;
  const char *s1 = "==========";

  fp = fopen(fname, "w");
  fprintf(fp, "%s%s   Hqp_IpLQDOCP::dump  %s%s\n", s1, s1, s1, s1);
  fprintf(fp, "_kmax = %d\n", _kmax);
  fprintf(fp, "_nk: ");   iv_foutput(fp, _nk);
  fprintf(fp, "_mk: ");   iv_foutput(fp, _mk);
  fprintf(fp, "_nmk: ");  iv_foutput(fp, _nmk);
  fprintf(fp, "_nks: ");  iv_foutput(fp, _nks);
  fprintf(fp, "_indk: "); iv_foutput(fp, _indk);

  for ( k = 0; k <= _kmax; k++) {
    fprintf(fp, "\n%s%s   k = %d   %s%s\n\n", s1, s1, k, s1, s1);
    fprintf(fp, "\n%s   gx[%d]: ", s1, k); v_foutput(fp, gx[k]); 
    if ( k < _kmax ) {
      fprintf(fp, "\n%s   gu[%d]: ", s1, k); v_foutput(fp, gu[k]);      
      if ( ! _a_sparse ) {
	fprintf(fp, "\n%s   fx[%d]: ", s1, k); m_foutput(fp, fx[k]);   
	fprintf(fp, "\n%s   fu[%d]: ", s1, k); m_foutput(fp, fu[k]);   
      }
      fprintf(fp, "\n%s   f[%d]: ", s1, k); v_foutput(fp, f[k]);   
    }

    fprintf(fp, "\n%s   _raki[%d]: ", s1, k); iv_foutput(fp, _raki[k]);   
    fprintf(fp, "\n%s   _rcki[%d]: ", s1, k); iv_foutput(fp, _rcki[k]);
    fprintf(fp, "\n%s   _rckai[%d]: ", s1, k); iv_foutput(fp, _rckai[k]);

    fprintf(fp, "\n%s   a[%d]: ", s1, k); v_foutput(fp, a[k]);      

    if ( k < _kmax ) {
      fprintf(fp, "\n%s   Gxu[%d]: ", s1, k); m_foutput(fp, Gxu[k]);   
      fprintf(fp, "\n%s   Guu[%d]: ", s1, k); m_foutput(fp, Guu[k]);   
      fprintf(fp, "\n%s   S[%d]: ", s1, k); m_foutput(fp, S[k]);   
      fprintf(fp, "\n%s   PQ[%d]: ", s1, k); m_foutput(fp, PQ[k]);   
      fprintf(fp, "\n%s   CH_Guu[%d]: ", s1, k); m_foutput(fp, CH_Guu[k]);   
      fprintf(fp, "\n%s   CH_Guu_p[%d]: ", s1, k); 
      px_foutput(fp, CH_Guu_p[k]);   
      fprintf(fp, "\n%s   CH_Guu_b[%d]: ", s1, k); 
      px_foutput(fp, CH_Guu_b[k]);   
      fprintf(fp, "\n%s   ctx[%d]: ", s1, k); m_foutput(fp, ctx[k]);   
      fprintf(fp, "\n%s   d11[%d]: ", s1, k); m_foutput(fp, d11[k]);   
      fprintf(fp, "\n%s   d12[%d]: ", s1, k); m_foutput(fp, d12[k]);   
      fprintf(fp, "\n%s   CD[%d]: ", s1, k); m_foutput(fp, CD[k]);   
      fprintf(fp, "\n%s   CD_p[%d]: ", s1, k); px_foutput(fp, CD_p[k]);   

      fprintf(fp, "\n%s   Rux[%d]: ", s1, k); m_foutput(fp, Rux[k]);   
      fprintf(fp, "\n%s   Ruyb[%d]: ", s1, k); m_foutput(fp, Ruyb[k]);   
      fprintf(fp, "\n%s   Ru[%d]: ", s1, k); v_foutput(fp, Ru[k]);   
      fprintf(fp, "\n%s   Ryx[%d]: ", s1, k); m_foutput(fp, Ryx[k]); 
      fprintf(fp, "\n%s   Ry[%d]: ", s1, k); v_foutput(fp, Ry[k]);
      fprintf(fp, "\n%s   Ryyb[%d]: ", s1, k); m_foutput(fp, Ryyb[k]);        
    }
    
    fprintf(fp, "\n%s   Vxx[%d]: ", s1, k); m_foutput(fp, Vxx[k]);   
    fprintf(fp, "\n%s   Vx[%d]: ", s1, k); v_foutput(fp, Vx[k]);   
    
    fprintf(fp, "\n%s   cbx[%d]: ", s1, k); m_foutput(fp, cbx[k]);   
    fprintf(fp, "\n%s   cb[%d]: ", s1, k); v_foutput(fp, cb[k]);   
    
    fprintf(fp, "\n%s   x[%d]: ", s1, k); v_foutput(fp, x[k]);
    fprintf(fp, "\n%s   y[%d]: ", s1, k); v_foutput(fp, y[k]);      
    if ( k < _kmax ) {
      fprintf(fp, "\n%s   u[%d]: ", s1, k); v_foutput(fp, u[k]);   
      fprintf(fp, "\n%s   p[%d]: ", s1, k); v_foutput(fp, p[k]);   
    }
  }

  fprintf(fp, "\n%s   Sparse   %s\n", s1, s1);
  fprintf(fp, "\n%s   _Q_ori (symsp): ", s1); sp_foutput(fp, _Q_ori);   
  fprintf(fp, "\n%s   _A_ori: ", s1); sp_foutput(fp, _A_ori);   
  fprintf(fp, "\n%s   _C_ori: ", s1); sp_foutput(fp, _C_ori);   
  fprintf(fp, "\n%s   _CTC (symsp): ", s1); sp_foutput(fp, _CTC);   
  fprintf(fp, "\n%s   _Z_ori: ", s1); v_foutput(fp, _Z_ori);   
  fprintf(fp, "\n%s   _W_ori: ", s1); v_foutput(fp, _W_ori);   
  fprintf(fp, "\n%s   _WZ_ori: ", s1); v_foutput(fp, _WZ_ori);   

  fprintf(fp, "\n%s   Static   %s\n", s1, s1);
  fprintf(fp, "\n%s   _rak: ", s1); iv_foutput(fp, _rak);   
  fprintf(fp, "\n%s   _rck: ", s1); iv_foutput(fp, _rck);   
  fprintf(fp, "\n%s   _rcka: ", s1); iv_foutput(fp, _rcka);   
  fprintf(fp, "\n%s   m1: ", s1); m_foutput(fp, m1);   
  fprintf(fp, "\n%s   m2: ", s1); m_foutput(fp, m2);   
  fprintf(fp, "\n%s   v1: ", s1); v_foutput(fp, v1);   
  fprintf(fp, "\n%s   v2: ", s1); v_foutput(fp, v2);   
  fclose(fp);
}

//--------------------------------------------------------------------------
//   Initialize LQDOCP, called in the first SQP step.
//   Check for DOCP structure.
//   Allocate memory.
//   Allocate memory for  C'*D*C.
//
//   E. Arnold
//               02/22/97  _nk, _mk
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::init(const Hqp_Program *qp)
{
  int _kmax_old, ret;

  pr_("\nLQDOCP::init ");

  _kmax_old = _kmax;
  ret = Get_Dim(qp);
  assert(ret);
  if ( _kmax != _kmax_old )
    resize(); 
  Get_Constr_Dim(qp);
  ret = Check_Structure(qp);  
  assert(ret);  

  pr("\nLQDOCP::init   _kmax = %d\n", _kmax);

  update(qp);
}

//--------------------------------------------------------------------------
//   Update LQDOCP data structure, called once per SQP step.
//   Check DOCP structure.
// 
//   E. Arnold
//               02/22/97  _nk, _mk
//               06/30/97  copy _A_ori, _C_ori
//                         structure of _CTC 
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::update(const Hqp_Program *qp)
{
  int i, j, k, knm, kn, ret;
  double sum;

  pr_("LQDOCP::update\n");

  ret = Check_Structure(qp);
  assert(ret);

  if ( _fixed_x0 ) 
    pr_("fixed initial state\n");

  //   do not really copy the matrix!
  _Q_ori = qp->Q;

  if ( _a_sparse ) {
    _A_ori = sp_copy3(qp->A, _A_ori);
    sp_compact(_A_ori, 0.0);
    //   transpose _A_ori
    _AT = sp_transp(_A_ori, _AT);
  } else {
    _A_ori = qp->A;
    SP_FREE(_AT);

    //   extract  coefficient matrices fx[k], fu[k] from qp->A
    for ( k = 0; k < _kmax; k++ ) {
      knm = _nmk[k];
      kn = _nks[k];
      fx[k] = m_resize(fx[k], _nk[k+1], _nk[k]);
      sp_extract_mat(qp->A, kn, knm, fx[k]);   
      fu[k] = m_resize(fu[k], _nk[k+1], _mk[k]);
      sp_extract_mat(qp->A, kn, knm+_nk[k], fu[k]);   
    }
  }

  _C_ori = sp_copy3(qp->C, _C_ori);
  sp_compact(_C_ori, 0.0);

  //   sparsity structure of matrix _CTC (only elements (i,j), j>=i)
  if ( ! _CTC )
    _CTC = sp_get(qp->Q->m, qp->Q->n, 10);
  _CT = sp_transp(_C_ori, _CT);

  //   check DOCP structure of matrix _CT
  pr_("check _C_ori");
  check_sparse(_C_ori);
  pr_(" check _CT");
  check_sparse(_CT);

  _CT = sp_ones(_CT);
  v1 = v_resize(v1, _CT->n);
  v1 = v_ones(v1);
  for ( k = 0; k <= _kmax; k++ ) {
    knm = _nmk[k];
    for ( i = 0; i < (k<_kmax?_nk[k]+_mk[k]:_nk[k]); i++ )
      for ( j = i; j < (k<_kmax?_nk[k]+_mk[k]:_nk[k]); j++ ) {
	sum = sprow_inprod(&(_CT->row[knm+i]), v1, 
			   &(_CT->row[knm+j]));
	if ( sum != 0.0 ) 
	  sp_set_val(_CTC, knm+i, knm+j, sum);   // _CTC is symsp!
      }
  }
  _CT = sp_transp(_C_ori, _CT);
  pr_("update finished\n"); 
}

//--------------------------------------------------------------------------
//   LQDOCP factorization routine, called once per QP step.
//   Form _WZ_ori, check for 'active' inequalities, form _CTC.
//
//   E. Arnold
//               02/22/97  _nk, _mk
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::factor(const Hqp_Program *qp, 
			  const VEC *y_v, const VEC *ny_v)
{
  int i, j, k;

  assert(((int)y_v->dim == _C_ori->m) && (ny_v->dim == y_v->dim));
  pr_("LQDOCP::factor:\n"); 
  //  mem_info();

  _Z_ori = v_copy1(y_v, _Z_ori);
  _W_ori = v_copy1(ny_v, _W_ori);

  //   _WZ_ori = W^{-1}Z
  _WZ_ori = v_resize(_WZ_ori, _Z_ori->dim);
  _WZ_ori = v_slash(_W_ori, _Z_ori, _WZ_ori);

  //   'active' constraints
  _rca = iv_resize(_rca, _WZ_ori->dim);
  for ( i = 0; i < (int) _WZ_ori->dim; i++ )
    if ( v_entry(_WZ_ori, i)*sprow_norm2(&(_C_ori->row[i])) > _wz_tol )
      //    if ( v_entry(_WZ_ori, i) > _wz_tol )
      iv_set_val(_rca, i, 1);
    else
      iv_set_val(_rca, i, 0);

  //   _CTC = C^T(W^{-1}Z)C
  v1 = v_copy1(_WZ_ori, v1);
  for ( i = 0; i < (int) v1->dim; i++ )
    if ( iv_entry(_rca, i) )
      v_set_val(v1, i, 0.0);
  _CTC = CTDC(_CT, v1, _CTC);

  //   set diagonal access paths for _Q_ori and _CTC
  if ( !_Q_ori->flag_diag )
    sp_diag_access(_Q_ori);
  if ( !_CTC->flag_diag )
    sp_diag_access(_CTC);

  //   indices of 'active' inequality contraints
  _rcka = iv_resize(_rcka, _kmax+1);
  _rcka = iv_zero(_rcka);
  for ( k = 0; k <= _kmax; k++ ) 
    for ( i = 0; i < iv_entry(_rck, k); i++ ) 
      if ( iv_entry(_rca, iv_entry(_rcki[k], i)) )
	iv_add_val(_rcka, k, 1);
  for ( k = 0; k <= _kmax; k++ ) {
    _rckai[k] = iv_resize(_rckai[k], iv_entry(_rcka, k));
    for ( i = 0, j = 0; i < iv_entry(_rck, k); i++ )
      if ( iv_entry(_rca, iv_entry(_rcki[k], i)) ) {
	iv_set_val(_rckai[k], j, iv_entry(_rcki[k], i));
	j++;
      }
  }
  //  printf("vor ExRiccatiFactor\n");
  if ( _wz_tol != HUGE_VAL )
    ExRiccatiFactor();
  else
    ExRiccatiFactorSc();
  
  for ( k = 0, i = 0, j = 0; k <= _kmax; k++ ) {
    i += iv_entry(_rcka, k);
    if ( k < _kmax ) j = ((int)PQ[k]->m>j)? PQ[k]->m: j;
    //    if ( k < _kmax ) j = ((int)CH_Guu[k]->m>j)? CH_Guu[k]->m: j;
  }
  //  printf(" %d 'active' inequalities, max(dim(PQ)) = %d\n", i, j); 
  //  printf(" %d 'active' inequalities, max(dim(CH_Guu)) = %d\n", i, j); 
}

//--------------------------------------------------------------------------
//   Solve the IP-Newton system.
//   E. Arnold
//              2/21/97 iterative improvement
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::step(const Hqp_Program *qp, 
			const VEC *y_v, const VEC *ny_v,
			const VEC *r1, const VEC *r2, const VEC *r3,
			const VEC *r4, VEC *dx, VEC *dy, VEC *dz, VEC *dw)
{
  int k, i, knm, kn;
  double res;

  assert((int)r1->dim == _Q_ori->m && (int)dx->dim == _Q_ori->m);
  assert((int)r2->dim == _A_ori->m && (int)dy->dim == _A_ori->m);
  assert((int)r3->dim == _C_ori->m && (int)dz->dim == _C_ori->m);
  assert((int)r4->dim == _C_ori->m && (int)dw->dim == _C_ori->m);

  pr_("LQDOCP::solve  ");

  //   gx, gu, f from r1, r2
  for ( k = 0; k <= _kmax; k++ ) {
    knm = _nmk[k]; 
    kn = _nks[k];
    gx[k] = v_resize(gx[k], _nk[k]);
    gx[k] = v_move(r1, knm, _nk[k], gx[k], 0);
    gx[k] = sv_mlt(-1.0, gx[k], gx[k]);
    a[k] = v_move_iv(r2, _raki[k], a[k]);
    if ( k < _kmax ) {
      gu[k] = v_resize(gu[k], _mk[k]);
      gu[k] = v_move(r1, knm+_nk[k], _mk[k], gu[k], 0);
      gu[k] = sv_mlt(-1.0, gu[k], gu[k]);
      f[k] = v_resize(f[k], _nk[k+1]);
      f[k] = v_move(r2, kn, _nk[k+1], f[k], 0);
    }
  }

  //   modify gx, gu for reduced system
  v2 = v_resize(v2, r3->dim);
  v2 = v_star(_Z_ori, r3, v2);
  v2 = v_add(v2, r4, v2);
  v2 = v_slash(_W_ori, v2, v2);
  for ( i = 0; i < (int) v2->dim; i++ )
    //    if (v_entry(_WZ_ori, i) > _wz_tol )
    if ( iv_entry(_rca, i) )
      v_set_val(v2, i, 0.0);
  v1 = sp_vm_mlt(_C_ori, v2, v1);
  for ( k = 0; k <= _kmax; k++ ) {
    knm = _nmk[k]; 
    v2 = v_resize(v2, _nk[k]);
    v2 = v_move(v1, knm, _nk[k], v2, 0);
    gx[k] = v_add(gx[k], v2, gx[k]);
    if ( k < _kmax) {
      v2 = v_resize(v2, _mk[k]);
      v2 = v_move(v1, knm+_nk[k], _mk[k], v2, 0);
      gu[k] = v_add(gu[k], v2, gu[k]);
    }
  }

  //   append 'active' inequalities to equalities
  v1 = v_resize(v1, r3->dim);
  v1 = v_slash(_Z_ori, r4, v1);
  v1 = v_add(r3, v1, v1);
  for ( k = 0; k <= _kmax; k++ ) {
    v2 = v_move_iv(v1, _rckai[k], v2);
    a[k] = v_concat(a[k], v2, a[k]);
  }

  if ( _wz_tol != HUGE_VAL )
    ExRiccatiSolve();
  else
    ExRiccatiSolveSc();

  // dx, dy
  for ( k = 0; k <= _kmax; k++ ) {
    knm = _nmk[k]; //k*nm;
    kn = _nks[k];  //k*_n;
    dx = v_move(x[k], 0 , _nk[k], dx, knm);
    v1 = v_resize(v1, iv_entry(_rak, k));
    v1 = v_move(y[k], 0, v1->dim, v1, 0);
    dy = v_set_iv(v1, _raki[k], dy);
    if ( k < _kmax ) {
      dx = v_move(u[k], 0 , _mk[k], dx, knm+_nk[k]);
      dy = v_move(p[k], 0 , _nk[k+1], dy, kn);
    }
  }

  // change sign for dx
  dx = sv_mlt(-1.0, dx, dx);

  //  dw = C*dx-r3, dz = W^(-1)(-Z*dw+r4)
  v1 = v_resize(v1, r3->dim);
  v1 = sp_mv_mlt(_C_ori, dx, v1);
  dw = v_sub(v1, r3, dw);
  v1 = v_star(_Z_ori, dw, v1);
  v1 = v_sub(r4, v1, v1);
  dz = v_slash(_W_ori, v1, dz);

  // dz for 'active' constraints
  for ( k = 0; k <= _kmax; k++ ) {
    v1 = v_resize(v1, iv_entry(_rcka, k));
    v1 = v_move(y[k], iv_entry(_rak, k), v1->dim, v1, 0);
    dz = v_set_iv(v1, _rckai[k], dz);
  }

  //   check residuum of IP-Newton system
  Residuum(r1, r2, r3, r4, dx, dy, dz, dw, 0);

  res = v_norm_inf(_res1)+v_norm_inf(_res2)+v_norm_inf(_res3)+
    v_norm_inf(_res4);

  pr("Residuum: %g\n", res);
}

//--------------------------------------------------------------------------
//   Calculate residuum components _res1, _res2, _res3, _res4
//   of the IP-Newton system.
//   If dump != 0, dump to file fname.
//   E. Arnold 2/21/97
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::Residuum(const VEC *r1, const VEC *r2, const VEC *r3,
			    const VEC *r4, const VEC *dx, const VEC *dy, 
			    const VEC *dz, const VEC *dw, const int dump)
{
  FILE *fp = NULL;
  const char *fname = "dump_residuum.dat";

  if ( dump ) 
    fp = fopen(fname, "w");

  //   residuum of Q*dx+A'*dy+C'*dw = r1
  v1 = v_resize(v1, r1->dim);
  v1 = sp_mv_symmlt(_Q_ori, dx, v1); 
  v2 = sp_vm_mlt(_A_ori, dy, v2); 
  v1 = v_sub(v2, v1, v1);
  v2 = sp_vm_mlt(_C_ori, dz, v2); 
  v1 = v_add(v1, v2, v1);
  v1 = v_sub(v1, r1, v1);         
  _res1 = v_copy1(v1, _res1);

  if ( dump ) {
    fprintf(fp, "===== Residuum1: %g\n", v_norm_inf(v1)); 
    v_foutput(fp, v1);
    v1 = v_resize(v1, r1->dim);
    v1 = sp_mv_symmlt(_Q_ori, dx, v1); 
    fprintf(fp, "===== _Q_ori*dx:"); 
    v_foutput(fp, v1);
    v2 = sp_vm_mlt(_A_ori, dy, v2); 
    fprintf(fp, "===== _A_ori^T*dy:"); 
    v_foutput(fp, v2);
    v2 = sp_vm_mlt(_C_ori, dz, v2); 
    fprintf(fp, "===== _C_ori^T*dz:"); 
    v_foutput(fp, v2);
    fprintf(fp, "===== r1:");  
    v_foutput(fp, r1);
  }

  //   residuum of A*dx = r2
  v1 = v_resize(v1, r2->dim);
  v1 = sp_mv_mlt(_A_ori, dx, v1);
  v1 = v_sub(v1, r2, v1);
  _res2 = v_copy1(v1, _res2);

  if ( dump ) {
    fprintf(fp, "\n\n===== Residuum2: %g\n", v_norm2(v1)); 
    v_foutput(fp, v1);
    fprintf(fp, "===== r2: "); 
    v_foutput(fp, r2);
  }

  //   residuum of C*dx-dw = r3
  v1 = v_resize(v1, r3->dim);
  v1 = sp_mv_mlt(_C_ori, dx, v1);
  v1 = v_sub(v1, dw, v1);
  v1 = v_sub(v1, r3, v1);
  _res3 = v_copy1(v1, _res3);

  if ( dump ) {
    fprintf(fp, "\n\n===== Residuum3: %g\n", v_norm2(v1)); 
    v_foutput(fp, v1);
    fprintf(fp, "===== r3: "); 
    v_foutput(fp, r3);
  }

  //   residuum of W*dz+Z*dw = r4
  v1 = v_resize(v1, r4->dim);
  v1 = v_star(_W_ori, dz, v1); 
  v2 = v_star(_Z_ori, dw, v2); 
  v1 = v_add(v1, v2, v1);
  v1 = v_sub(v1, r4, v1); 
  _res4 = v_copy1(v1, _res4);

  if ( dump ) {
    fprintf(fp, "\n\n===== Residuum4: %g\n", v_norm2(v1)); 
    v_foutput(fp, v1);
    fprintf(fp, "===== dz: "); 
    v_foutput(fp, dz);
    fprintf(fp, "===== dw: "); 
    v_foutput(fp, dw);
    fprintf(fp, "===== r4: "); 
    v_foutput(fp, r4);
  }

  if ( dump )
    fclose(fp);
}

//--------------------------------------------------------------------------
//   Form Gxx = Hxx + fx'*Vxx*fx, Gxu = Hxu + fx'*Vxx*fu, 
//   Guu = Huu + fu'*Vxx*fu for stage k.
//   Dense version.
//   E. Arnold   12/02/96
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::FormGxx(int k)
{
  int knm = _nmk[k];

  //   Gxx = Hxx + fx'*Vxx*fx
  m1 = m_resize(m1, _nk[k], _nk[k]);
  Gxx = m_resize(Gxx, _nk[k], _nk[k]);
  symsp_extract_mat(_Q_ori, knm, Gxx);   
  symsp_extract_mat(_CTC, knm, m1);   
  Gxx = m_add(Gxx, m1, Gxx);
  m2 = mtrm_mlt(fx[k], Vxx[k+1], m2);
  m1 = m_mlt(m2, fx[k], m1);
  Gxx = m_add(Gxx, m1, Gxx);
  Gxx = m_symm(Gxx, Gxx);

  //   Gxu = Hxu + fx'*Vxx*fu
  Gxu[k] = m_resize(Gxu[k], _nk[k], _mk[k]);
  sp_extract_mat(_Q_ori, knm, knm+_nk[k], Gxu[k]);   
  m1 = m_resize(m1, _nk[k], _mk[k]);
  sp_extract_mat(_CTC, knm, knm+_nk[k], m1);   
  Gxu[k] = m_add(Gxu[k], m1, Gxu[k]);
  m1 = m_mlt(m2, fu[k], m1);
  Gxu[k] = m_add(Gxu[k], m1, Gxu[k]);

  //   Guu = Huu + fu'*Vxx*fu
  Guu[k] = m_resize(Guu[k], _mk[k], _mk[k]);
  symsp_extract_mat(_Q_ori, knm+_nk[k], Guu[k]); 
  m1 = m_resize(m1, _mk[k], _mk[k]);
  symsp_extract_mat(_CTC, knm+_nk[k], m1);   
  Guu[k] = m_add(Guu[k], m1, Guu[k]);
  m2 = mtrm_mlt(fu[k], Vxx[k+1], m2);
  m1 = m_mlt(m2, fu[k], m1);
  Guu[k] = m_add(Guu[k], m1, Guu[k]);
  Guu[k] = m_symm(Guu[k], Guu[k]);
}

//--------------------------------------------------------------------------
//   Form Gxx = Hxx + fx'*Vxx*fx, Gxu = Hxu + fx'*Vxx*fu, 
//   Guu = Huu + fu'*Vxx*fu for stage k.
//   Sparse version.   
//   E. Arnold   12/02/96
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::FormGxxSp(int k)
{
  int knm = _nmk[k]; 
  int kn = _nks[k];
  int k1 = k+1;
  int i, j, col, ii, j_idx;
  SPROW *row;
  row_elt *elt;
  double val;

  //   initialize Gxx, Gxu, Guu 
  Gxx = m_resize(Gxx, _nk[k], _nk[k]);
  Gxx = m_zero(Gxx);
  Gxu[k] = m_resize(Gxu[k], _nk[k], _mk[k]);
  Gxu[k] = m_zero(Gxu[k]);
  Guu[k] = m_resize(Guu[k], _mk[k], _mk[k]);
  Guu[k] = m_zero(Guu[k]);

  //   build m1 = [f_x^T; f_u^T]*Vxx
  m1 = m_resize(m1, _nk[k]+_mk[k], _nk[k+1]);
  m1 = m_zero(m1);
  row = _AT->row + knm;
  for ( i = 0; i < _nk[k]+_mk[k]; i++, row++) {
    for ( j_idx = (((k==0)||(i>=_nk[k]))?0:1); j_idx < row->len; j_idx++ ) {
      j = row->elt[j_idx].col-kn;
      if ( j >= _nk[k] )
	break;
      val = row->elt[j_idx].val;
      __mltadd__(m1->me[i], Vxx[k1]->me[j], row->elt[j_idx].val, _nk[k+1]);
    }
  }

  //   build upper triangel of m1*[fx fu]
  row = _AT->row + knm;
  for ( j = 0; j < _nk[k]+_mk[k]; j++, row++ ) {
    for ( j_idx = (((k==0)||(j>=_nk[k])) ? 0: 1); j_idx < row->len; j_idx++ ) {
      i = row->elt[j_idx].col-kn;
      if ( i >= _nk[k+1] )
	break;
      val = row->elt[j_idx].val;
      for ( ii = 0; ii <= j; ii++ ) {
	if ( ii < _nk[k] )
	  if ( j < _nk[k] )
	    m_add_val(Gxx, ii, j, m_entry(m1, ii, i)*val);
	  else
	    m_add_val(Gxu[k], ii, j-_nk[k], m_entry(m1, ii, i)*val);
	else
	  m_add_val(Guu[k], ii-_nk[k], j-_nk[k], m_entry(m1, ii, i)*val);
      }
    }
  }

  //   add _Q_ori
  row = _Q_ori->row+knm;
  for ( i = 0; i < _nk[k]+_mk[k]; i++, row++) {
    if ( row->len > 0 ) {
      if ( ( row->diag < 0 ) || ( row->diag >= row->len ) )
	m_error(E_INTERN,
		"LQDOCP::FormGxxSp: invalid diagonal access path in FormGxx");
      elt = row->elt+row->diag;
      for ( j = row->diag; j < row->len; j++, elt++ ) {
	col = elt->col-knm;
	if ( i < _nk[k] )
	  if ( col < _nk[k] )
	    m_add_val(Gxx, i, col, elt->val); 
	  else
	    m_add_val(Gxu[k], i, col-_nk[k], elt->val); 
	else
	  m_add_val(Guu[k], i-_nk[k], col-_nk[k], elt->val); 
      }
    }
  }

  //   add _CTC
  row = _CTC->row+knm;
  for ( i = 0; i < _nk[k]+_mk[k]; i++, row++) {
    if ( row->len > 0 ) {
      if ( ( row->diag < 0 ) || ( row->diag >= row->len ) )
	m_error(E_INTERN,
		"LQDOCP::FormGxxSp: invalid diagonal access path in FormGxx");
      elt = row->elt+row->diag;
      for ( j = row->diag; j < row->len; j++, elt++ ) {
	col = elt->col-knm;
	if ( i < _nk[k] )
	  if ( col < _nk[k] )
	    m_add_val(Gxx, i, col, elt->val); 
	  else
	    m_add_val(Gxu[k], i, col-_nk[k], elt->val); 
	else
	  m_add_val(Guu[k], i-_nk[k], col-_nk[k], elt->val); 
      }
    }
  }

  //   make Gxx and Guu symmetric
  for ( i = 0; i < _nk[k]; i++ )
    for ( j = i+1; j < _nk[k]; j++ )
      m_set_val(Gxx, j, i, m_entry(Gxx, i, j));
  for ( i = 0; i < _mk[k]; i++ )
    for ( j = i+1; j < _mk[k]; j++ )
      m_set_val(Guu[k], j, i, m_entry(Guu[k], i, j));
}

//--------------------------------------------------------------------------
//   Form Gx = gx + fx'*(Vx + Vxx*f), Gu = gu + fu'*Vx + Vxx*f)
//   for stage k. Dense version.
//   E. Arnold  10/08/96 
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::FormGx(int k)
{
  //   Gx = gx + fx'*(Vx + Vxx*f)
  v1 = mv_mlt(Vxx[k+1], f[k], v1);
  v1 = v_add(v1, Vx[k+1], v1);
  v2 = vm_mlt(fx[k], v1, v2);
  Gx = v_add(gx[k], v2, Gx);
  
  //   Gu = gu + fu'*(Vx + Vxx*f)
  v2 = vm_mlt(fu[k], v1, v2);
  Gu = v_add(gu[k], v2, Gu);
}

//--------------------------------------------------------------------------
//   Form Gx = gx + fx'*(Vx + Vxx*f), Gu = gu + fu'*Vx + Vxx*f)
//   for stage k. Sparse version.
//   E. Arnold  12/12/96 
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::FormGxSp(int k)
{
  int knm = _nmk[k]; 
  int kn = _nks[k];
  int i, j, j_idx;
  SPROW *row;

  //   initialize Gx, Gu
  Gx = v_copy1(gx[k], Gx);
  Gu = v_copy1(gu[k], Gu);

  //   v1 = (Vx + Vxx*f)
  v1 = mv_mlt(Vxx[k+1], f[k], v1);
  v1 = v_add(v1, Vx[k+1], v1);

  //   [Gx; Gu] = [fx'; fu']*v1
  row = _AT->row + knm;
  for ( i = 0; i < _nk[k]+_mk[k]; i++, row++ ) {
    for ( j_idx = (((k==0)||(i>=_nk[k]))?0:1); j_idx < row->len; j_idx++ ) {
      j = row->elt[j_idx].col-kn;
      if ( j >= _nk[k+1] )
	break;
      if ( i < _nk[k] )
	v_add_val(Gx, i, row->elt[j_idx].val*v_entry(v1, j));
      else
	v_add_val(Gu, i-_nk[k], row->elt[j_idx].val*v_entry(v1, j));
    }
  }
}

//--------------------------------------------------------------------------
//   Get ax, au, yny from _A_ori, _C_ori, _W_ori, _Z_ori
//   E. Arnold   10/24/96
//               02/22/97  _nk, _mk
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::Formax(int k)
{
  int i, j_idx;

  //   ax
  ax = m_resize(ax, iv_entry(_rak, k), _nk[k]);
  sp_extract_mat_iv(_A_ori, _raki[k], _nmk[k], ax);
  m1 = m_resize(m1, iv_entry(_rcka, k), _nk[k]);
  sp_extract_mat_iv(_C_ori, _rckai[k], _nmk[k], m1);
  ax = m_concatc(ax, m1, ax);

  //   au
  if ( k < _kmax ) {
    au = m_resize(au, iv_entry(_rak, k), _mk[k]);
    sp_extract_mat_iv(_A_ori, _raki[k], _nmk[k]+_nk[k], au);
    m1 = m_resize(m1, iv_entry(_rcka, k), _mk[k]);
    sp_extract_mat_iv(_C_ori, _rckai[k], _nmk[k]+_nk[k], m1);
    au = m_concatc(au, m1, au);
  }

  //   yny
  yny = v_resize(yny, iv_entry(_rak, k)+iv_entry(_rcka, k));
  yny = v_zero(yny);
  for ( i = 0; i < iv_entry(_rcka, k); i++ ) {
    j_idx = iv_entry(_rckai[k], i);
    v_set_val(yny, iv_entry(_rak, k)+i, 
	      (v_entry(_W_ori, j_idx))/(v_entry(_Z_ori, j_idx)));
  }

  //   remove fixed initial states
  if ( ( k == 0 ) && ( _fixed_x0 ) ) {
    m1 = m_resize(m1, ax->m-_nk[k], ax->n);
    m1 = m_move(ax, _nk[k], 0, m1->m, m1->n, m1, 0, 0);
    ax = m_copy1(m1, ax);
    m1 = m_resize(m1, au->m-_nk[k], au->n);
    m1 = m_move(au, _nk[k], 0, m1->m, m1->n, m1, 0, 0);
    au = m_copy1(m1, au);
    v1 = v_resize(v1, yny->dim-_nk[k]);
    v1 = v_move(yny, _nk[k], v1->dim, v1, 0);
    yny = v_copy1(v1, yny);
  }
}

//--------------------------------------------------------------------------
//   Extended Riccati solution: factorization routine.
//   E. Arnold   10/08/96
//               02/22/97  _nk, _mk
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::ExRiccatiFactor(void)
{
  int k, i, k1, kn, knm;

  //   Final stage: set Vxx, cbx, db for _kmax.

  Vxx[_kmax] = m_resize(Vxx[_kmax], _nk[_kmax], _nk[_kmax]);
  symsp_extract_mat(_Q_ori, _nmk[_kmax], Vxx[_kmax]);   
  m1 = m_resize(m1, _nk[_kmax], _nk[_kmax]);
  symsp_extract_mat(_CTC, _nmk[_kmax], m1);   
  Vxx[_kmax] = m_add(Vxx[_kmax], m1, Vxx[_kmax]);

  Formax(_kmax);
  cbx[_kmax] = m_copy1(ax, cbx[_kmax]);

  db = m_resize(db, cbx[_kmax]->m, cbx[_kmax]->m);
  db = m_zero(db);
  if ( v_norm1(yny) > 0.0 )
    for ( i = 0; i < (int) yny->dim; i++ )
      m_set_val(db, i, i, -v_entry(yny, i));

  //   Backward run: calculate Vxx, Rux, Ryx for k=_kmax-1, ... ,0.

  for (k = _kmax-1; k >= 0; k--) {
    k1 = k+1; 
    kn = _nks[k];  
    knm = _nmk[k]; 

    //   Gxx, Gxu, Guu
    if ( !_a_sparse )
      FormGxx(k);
    else
      FormGxxSp(k);

    //    printf("vor Formax\n"); 
    //   ax, au, yny
    Formax(k);

    //   m2 = cbu
    //    printf("vor cbu\n"); 
    if ( _a_sparse )
      m1 = mbsptr_mlt(cbx[k1], _AT, kn, knm+_nk[k], _mk[k], m1);
    else
      m1 = m_mlt(cbx[k1], fu[k], m1);
    m2 = m_concatc(au, m1, m2);

    // Generalized elimination 
    if ( (MAT *)Z == MNULL )
      Z = m_resize(Z, 1, 1);     // Z must be allocated for GE_QP()!
    GE_QP(m2, S[k], Z, PQ[k], _ge_tol);

    //   cbx
    if ( _a_sparse )
      m1 = mbsptr_mlt(cbx[k1], _AT, kn, knm, _nk[k], m1);
    else
      m1 = m_mlt(cbx[k1], fx[k], m1);
    m2 = m_concatc(ax, m1, m2);
    cbx[k] = mtrm_mlt(PQ[k], m2, cbx[k]);

    //   d11, d12, d22  for 'active' inequalities
    d11[k] = m_resize(d11[k], S[k]->n, S[k]->n);
    d12[k] = m_resize(d12[k], S[k]->n, cbx[k]->m-S[k]->n);
    d22 = m_resize(d22, cbx[k]->m-S[k]->n, cbx[k]->m-S[k]->n);
    if ( ( m_norm1(db) > 0.0 ) || ( v_norm1(yny) > 0.0 ) ) {
      m1 = m_resize(m1, cbx[k]->m, cbx[k]->m);
      m1 = m_zero(m1);
      if ( v_norm1(yny) > 0.0 )
	for ( i = 0; i < (int) yny->dim; i++ )
	  m_set_val(m1, i, i, -v_entry(yny, i));
      if ( m_norm1(db) > 0.0 )
	m1 = m_move(db, 0, 0, db->m, db->n, m1, yny->dim, yny->dim);
      m2 = m_mlt(m1, PQ[k], m2);
      m1 = mtrm_mlt(PQ[k], m2, m1);
      d11[k] = m_move(m1, 0, 0, d11[k]->m, d11[k]->n, d11[k], 0, 0);
      d12[k] = m_move(m1, 0, S[k]->n, d12[k]->m, d12[k]->n, d12[k], 0, 0);
      d22 = m_move(m1, S[k]->n, S[k]->n, d22->m, d22->n, d22, 0, 0);
    } else {

      //   d11, d12, d22 for no 'active' inequalities
      d11[k] = m_zero(d11[k]);
      d12[k] = m_zero(d12[k]);
      d22 = m_zero(d22);
    }

    if ( S[k]->n == 0 ) {

      //   CH_Guu for unconstrained u: CH_Guu is CHfactor!
      CH_Guu[k] = m_copy1(Guu[k], CH_Guu[k]);
      //      printf("CH_Guu[%d]: ", k); m_output(CH_Guu[k]);
      //      CH_Guu[k] = CHfactor(CH_Guu[k]);
      CH_Guu_p[k] = px_resize(CH_Guu_p[k], CH_Guu[k]->m);
      CH_Guu_b[k] = px_resize(CH_Guu_b[k], CH_Guu[k]->m);
      CH_Guu[k] = BKPfactor(CH_Guu[k], CH_Guu_p[k], CH_Guu_b[k]);

      //   Rux, Ryx, Ruyb, Ryyb, ctx for unconstrained u
      //      Rux[k] = CHsolveMT(CH_Guu[k], Gxu[k], Rux[k]);
      Rux[k] = BKPsolveMT(CH_Guu[k], CH_Guu_p[k], CH_Guu_b[k], Gxu[k], Rux[k]);
      Ryx[k] = m_resize(Ryx[k], 0, _nk[k]);
      Ruyb[k] = m_resize(Ruyb[k], _mk[k], d22->m);
      Ruyb[k] = m_zero(Ruyb[k]);
      Ryyb[k] = m_resize(Ryyb[k], 0, d22->m);
      ctx[k] = m_resize(ctx[k], 0, _nk[k]);
    } else {

      //   ctx
      ctx[k] = m_resize(ctx[k], S[k]->n, _nk[k]);
      ctx[k] = m_move(cbx[k], 0, 0, S[k]->n, _nk[k], ctx[k], 0, 0);

      //   Rux, m2 = cttilx
      m1 = mmtr_mlt(d11[k], S[k], m1);
      m2 = mmtr_mlt(m1, Gxu[k], m2);
      m2 = m_sub(ctx[k], m2, m2);
      Rux[k] = m_mlt(S[k], m2, Rux[k]);
      Ruyb[k] = m_mlt(S[k], d12[k], Ruyb[k]);

      //   cbx
      //      i = cbx[k]->m-S[k]->n;
      m1 = m_resize(m1, cbx[k]->m-S[k]->n, _nk[k]);
      m1 = m_move(cbx[k], S[k]->n, 0, m1->m, _nk[k], m1, 0, 0);
      cbx[k] = m_copy1(m1, cbx[k]);

      if ( (int) S[k]->n < _mk[k] ) {

	//   CH_Guu for constrained u
	m1 = m_mlt(Guu[k], Z, m1);
	CH_Guu[k] = mtrm_mlt(Z, m1, CH_Guu[k]);
	//	CH_Guu[k] = CHfactor(CH_Guu[k]);
	//	m2 = CHsolveMT(CH_Guu[k], Z, m2);
	pivot = px_resize(pivot, CH_Guu[k]->m);
	blocks = px_resize(blocks, CH_Guu[k]->m);
	CH_Guu[k] = BKPfactor(CH_Guu[k], pivot, blocks);
	m2 = BKPsolveMT(CH_Guu[k], pivot, blocks, Z, m2);
        CH_Guu[k] = m_mlt(Z, m2, CH_Guu[k]);

	//   Rux for constrained u
	m1 = m_mlt(Guu[k], Rux[k], m1);
	m2 = m_transp(Gxu[k], m2);
	m1 = m_sub(m2, m1, m1);
	m2 = m_mlt(CH_Guu[k], m1, m2);
	Rux[k] = m_add(m2, Rux[k], Rux[k]);

	//   Ruyb for constrained u
	m1 = m_mlt(Guu[k], Ruyb[k], m1);
	m2 = m_mlt(CH_Guu[k], m1, m2);
	Ruyb[k] = m_sub(Ruyb[k], m2, Ruyb[k]);
      } 

      //   correction of Rux, Ruyb
      if ( ( m_norm1(db) > 0.0 ) || ( v_norm1(yny) > 0.0 ) ) {
	m1 = mtrm_mlt(S[k], Guu[k], m1);
	m2 = m_mlt(d11[k], m1, m2);
	CD[k] = m_mlt(S[k], m2, CD[k]);
	if ( (int) S[k]->n < _mk[k] ) {
	  m1 = m_mlt(Guu[k], CD[k], m1);
	  m2 = m_mlt(CH_Guu[k], m1, m2);
	  CD[k] = m_sub(CD[k], m2, CD[k]);
	}
	m1 = m_resize(m1, _mk[k], _mk[k]);
	m1 = m_ident(m1);
	CD[k] = m_sub(m1, CD[k], CD[k]);
	//	printf("CD[%d]: ", k); m_output(CD[k]);
	CD_p[k] = px_resize(CD_p[k], CD[k]->m);
	CD[k] = LUfactor(CD[k], CD_p[k]);
	{
	  Real cond = LUcondest(CD[k], CD_p[k]);
	  if ( _logging && cond > 1.0e10 )
	    fprintf(stderr, "LQDOCP: cond(CD[%d]): %g\n", k, cond);
	}
	Rux[k] = LUsolveM(CD[k], CD_p[k], Rux[k], Rux[k]);
	Ruyb[k] = LUsolveM(CD[k], CD_p[k], Ruyb[k], Ruyb[k]);
      } else
	CD[k] = m_resize(CD[k], 0, 0);

      //   Ryx
      m1 = m_transp(Gxu[k], m1);
      m2 = m_mlt(Guu[k], Rux[k], m2);
      m1 = m_sub(m1, m2, m1);
      Ryx[k] = mtrm_mlt(S[k], m1, Ryx[k]);

      //   Ryyb
      m1 = m_mlt(Guu[k], Ruyb[k], m1);
      Ryyb[k] = mtrm_mlt(S[k], m1, Ryyb[k]);
      Ryyb[k] = sm_mlt(-1.0, Ryyb[k], Ryyb[k]);
    }

    //   Vxx
    m1 = m_mlt(Gxu[k], Rux[k], m1);
    Vxx[k] = m_sub(Gxx, m1, Vxx[k]);
    if (Ryx[k]->m > 0) {
      m1 = mtrm_mlt(ctx[k], Ryx[k], m1);
      Vxx[k] = m_sub(Vxx[k], m1, Vxx[k]);
    }

    //   check for symmetry of Vxx
    if ( _logging && rel_symm(Vxx[k])  > 1e-8 ) {
      fprintf(stderr,
	      "LQDOCP: symmetry check Vxx[%d]: %g\n", k, rel_symm(Vxx[k]));
      if ( rel_symm(Vxx[k])  > 1e-3 ) {
	m_foutput(stderr, Vxx[k]);
	//	dump();
	//	m_error(E_POSDEF,"ExRiccatiFactor");
      }
    }
    Vxx[k] = m_symm(Vxx[k], Vxx[k]);

    if ( S[k]->n > 0 ) {

      //   check for symmetry of cbx
      m1 = m_mlt(Gxu[k], Ruyb[k], m1);
      m2 = m_transp(m1, m2);
      m1 = mtrm_mlt(Ryyb[k], ctx[k], m1);
      m2 = m_add(m1, m2, m2);
      m2 = m_sub(cbx[k], m2, m2);
      
      //   cbx
      m1 = mtrm_mlt(d12[k], Ryx[k], m1);
      cbx[k] = m_sub(cbx[k], m1, cbx[k]);
      m1 = m_sub(cbx[k], m2, m1);
      if ( m_norm1(m1) > 1e-8 ) {
	fprintf(stderr, "LQDOCP: symmetry check (k=%d):\n", k); 
	fprintf(stderr, "cbx1: "); m_output(cbx[k]);
	fprintf(stderr, "cbx2: "); m_output(m2);
	dump();
	m_error(E_POSDEF,"LQDOCP::ExRiccatiFactor");
      }
      
      //  db 
      m1 = mtrm_mlt(d12[k], Ryyb[k], m1);
      db = m_sub(d22, m1, db);

      //   check for symmetry of db
      m1 = m_transp(db, m1);
      m1 = m_sub(m1, db, m1);
      if ( _logging && m_norm1(m1) > 1e-8 ) {
	fprintf(stderr, "LQDOCP: symmetry check db[%d]:\n", k); 
	m_foutput(stderr, db);
	dump();
	m_error(E_POSDEF,"LQDOCP::ExRiccatiFactor");
      }
    } else
      db = m_copy1(d22, db);

    //    printf("Vxx[%d]", k); m_output(Vxx[k]);
    //    printf("Rux[%d]", k); m_output(Rux[k]);

    //    printf("k = %d, VEC: %d, MAT: %d\n",k,mem_info_numvar(TYPE_VEC,0),
    //   	   mem_info_numvar(TYPE_MAT,0));
  }

  //   Initial state
  if ( _fixed_x0 ) {
    m1 = m_copy1(db, m1);
  } else {
    m1 = m_resize(m1, _nk[0]+cbx[0]->m, _nk[0]+cbx[0]->m);
    m1 = m_zero(m1);
    m1 = m_move(Vxx[0], 0, 0, _nk[0], _nk[0], m1, 0, 0);
    m1 = m_move(cbx[0], 0, 0, cbx[0]->m, _nk[0], m1, _nk[0], 0);
    m2 = m_transp(cbx[0], m2);
    m1 = m_move(m2, 0, 0, _nk[0], cbx[0]->m, m1, 0, _nk[0]);
    m1 = m_move(db, 0, 0, db->m, db->n, m1, _nk[0], _nk[0]);
  }
  pivot = px_resize(pivot, m1->m);
  blocks = px_resize(blocks, m1->m);
  m1 = BKPfactor(m1, pivot, blocks);
  //   do not change m1, pivot, blocks!

}

//--------------------------------------------------------------------------
//   Extended Riccati solution: solve routine.
//   E. Arnold   10/08/96
//               02/22/97  _nk, _mk
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::ExRiccatiSolve(void)
{
  int k, k1, knm, kn;

  //   Backward run for Vx, Ru

  Vx[_kmax] = v_copy1(gx[_kmax], Vx[_kmax]);
  cb[_kmax] = v_copy1(a[_kmax], cb[_kmax]);

  for (k = _kmax-1; k >= 0; k--) {
    k1 = k+1;

    //   printf("k = %d, start\n", k);
    //   Gx, Gu
    if ( !_a_sparse )
      FormGx(k);
    else
      FormGxSp(k);

    //    printf("k = %d, cb\n", k);
    //   cb = PQ*[a; cb*f] (remove fixed initial states)
    v1 = mv_mlt(cbx[k1], f[k], v1);
    v1 = v_add(cb[k1], v1, v1);
    v2 = v_concat(a[k], v1, v2);
    if ( ( k == 0 ) && ( _fixed_x0 ) ) {
      v1 = v_resize(v1, v2->dim-_nk[k]);
      v1 = v_move(v2, _nk[k], v1->dim, v1, 0);
      v2 = v_copy1(v1, v2);
    }
    cb[k] = vm_mlt(PQ[k], v2, cb[k]);

    if ( S[k]->n == 0 ) {
      //      printf("k = %d, Ru,Ry\n", k);
      //   Ru, Ry for unconstrained u
      //      Ru[k] = CHsolve(CH_Guu[k], Gu, Ru[k]);
      Ru[k] = BKPsolve(CH_Guu[k], CH_Guu_p[k], CH_Guu_b[k], Gu, Ru[k]);
      Ry[k] = v_resize(Ry[k], 0);
    } else {

      //   ct
      //      i = cb[k]->dim-S[k]->n;
      ct = v_resize(ct, S[k]->n);
      ct = v_move(cb[k], 0, S[k]->n, ct, 0);

      //      printf("k = %d, ct\n", k);
      //   v1 = cttil
      v1 = vm_mlt(S[k], Gu, v1);
      v2 = mv_mlt(d11[k], v1, v2);
      v1 = v_sub(ct, v2, v1);

      //      printf("k = %d, Ru\n", k);
      //   Ru for constrained u
      Ru[k] = mv_mlt(S[k], v1, Ru[k]);

      //      printf("k = %d, cb\n", k);
      //   cb
      v1 = v_resize(v1, cb[k]->dim-S[k]->n);
      v1 = v_move(cb[k], S[k]->n, v1->dim, v1, 0);
      cb[k] = v_copy1(v1, cb[k]);

      if ( (int) S[k]->n < _mk[k] ) {
	//	printf("k = %d, Ru2\n", k);
	//   Ru for constrained u
	v1 = mv_mlt(Guu[k], Ru[k], v1);
	v1 = v_sub(Gu, v1, v1);
	v2 = mv_mlt(CH_Guu[k], v1, v2);
	Ru[k] = v_add(v2, Ru[k], Ru[k]);
      } 

      //   correction of Ru
      if ( CD[k]->n > 0 )
	Ru[k] = LUsolve(CD[k], CD_p[k], Ru[k], Ru[k]);

      //      printf("k = %d, Ry\n", k);
      //   Ry for constrained u
      v1 = mv_mlt(Guu[k], Ru[k], v1);
      v1 = v_sub(Gu, v1, v1);
      Ry[k] = vm_mlt(S[k], v1, Ry[k]);
    }

    //    printf("k = %d, Vx\n", k);
    //   Vx
    v1 = mv_mlt(Gxu[k], Ru[k], v1);
    Vx[k] = v_sub(Gx, v1, Vx[k]);
    if ( Ry[k]->dim > 0 ) {
      v1 = vm_mlt(ctx[k], Ry[k], v1);
      Vx[k] = v_sub(Vx[k], v1, Vx[k]);
    }

    //    printf("k = %d, cb\n", k);
    //   cb
    v1 = vm_mlt(d12[k], Ry[k], v1);
    cb[k] = v_sub(cb[k], v1, cb[k]);
  }

  //   Initial state
  x[0] = v_resize(x[0], _nk[0]);
  if ( _fixed_x0 ) {
    x[0] = v_move(a[0], 0, _nk[0], x[0], 0);
    x[0] = sv_mlt(-1.0, x[0], x[0]);
    if ( cb[0]->dim > 0 ) {   
      v1 = v_resize(v1, cb[0]->dim);
      v1 = mv_mltadd(cb[0], x[0], cbx[0], 1.0, v1);
      v1 = sv_mlt(-1.0, v1, v1);
      yb = BKPsolve(m1, pivot, blocks, v1, yb);
    } else
      yb = v_resize(yb, 0);
  } else {
    v1 = v_concat(Vx[0], cb[0], v1);
    v1 = sv_mlt(-1.0, v1, v1);
    v1 = BKPsolve(m1, pivot, blocks, v1, v1);
    x[0] = v_move(v1, 0, _nk[0], x[0], 0);
    yb = v_resize(yb, v1->dim-_nk[0]);
    yb = v_move(v1, _nk[0], v1->dim-_nk[0], yb, 0);
  }

  //   Forward run for x, u, p

  for (k = 0; k < _kmax; k++) {

    x[k+1] = v_resize(x[k+1], _nk[k+1]);
    u[k] = v_resize(u[k], _mk[k]);
    p[k] = v_resize(p[k], _nk[k+1]);

    //   u = -Rux*x -Ru -Ruyb*yb
    v1 = mv_mlt(Rux[k], x[k], v1);
    v1 = v_add(v1, Ru[k], v1);
    if ( yb->dim > 0 ) {
      v2 = mv_mlt(Ruyb[k], yb, v2);
      v1 = v_add(v1, v2, v1);
    }
    u[k] = sv_mlt(-1.0, v1, u[k]);

    //   x[k+1] = fx*x+fu*u+f
    k1 = k+1;
    knm = _nmk[k]; 
    kn = _nks[k];  
    if ( _a_sparse )
      v1 = bspv_mlt(_A_ori, kn, knm, _nk[k+1], x[k], v1);
    else
      v1 = mv_mlt(fx[k], x[k], v1);
    x[k1] = v_add(f[k], v1, x[k1]);
    if ( _a_sparse )
      v1 = bspv_mlt(_A_ori, kn, knm+_nk[k], _nk[k+1], u[k], v1);
    else
      v1 = mv_mlt(fu[k], u[k], v1);
    x[k1] = v_add(x[k1], v1, x[k1]);

    //   [y; yb] = -(Ryx*x + Ry + Ryyb*yb)
    if ( ( Ry[k]->dim+yb->dim > 0 ) || ( ( k == 0 ) && ( _fixed_x0 ) ) ) {
      v1 = mv_mlt(Ryx[k], x[k], v1);
      v1 = v_add(v1, Ry[k], v1);
      if ( yb->dim > 0 ) {
	v2 = mv_mlt(Ryyb[k], yb, v2);
	v1 = v_add(v1, v2, v1);
      }
      v1 = sv_mlt(-1.0, v1, v1);
      v2 = v_concat(v1, yb, v2);
      v1 = mv_mlt(PQ[k], v2, v1);
      if ( ( k == 0 ) && ( _fixed_x0 ) ) {
	v2 = v_resize(v2, _nk[k]);
	v2 = vm_mltadd(Vx[0], yb, cbx[0], 1.0, v2);
	v2 = mv_mltadd(v2, x[0], Vxx[0], 1.0, v2);
	v2 = sv_mlt(-1.0, v2, v2);
	v2 = v_concat(v2, v1, v2);
	v1 = v_copy1(v2, v1);
      }
      y[k] = v_resize(y[k], a[k]->dim);
      y[k] = v_move(v1, 0, a[k]->dim, y[k], 0);
      yb = v_resize(yb, v1->dim-y[k]->dim);
      yb = v_move(v1, a[k]->dim, v1->dim-a[k]->dim, yb, 0);
    }

    //   p = Vxx*x + Vx + cbx*yb
    p[k] = mv_mlt(Vxx[k1], x[k1], p[k]);
    p[k] = v_add(p[k], Vx[k1], p[k]);
    if ( yb->dim > 0 ) {
      v1 = vm_mlt(cbx[k1], yb, v1);
      p[k] = v_add(p[k], v1, p[k]);
    }
  }

  //   y[_kmax] = yb
  y[_kmax] = v_resize(y[_kmax], yb->dim);
  y[_kmax] = v_copy1(yb, y[_kmax]);
}

//--------------------------------------------------------------------------
//   Extended Riccati solution: factorization routine.
//   Diagonal scaled version without inequalities (_wz_tol = HUGE_VAL).
//
//   E. Arnold   02/22/97  
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::ExRiccatiFactorSc(void)
{
  int k, i, k1, kn, knm;

  //   Final stage: set Vxx, cbx, db for _kmax.

  Vxx[_kmax] = m_resize(Vxx[_kmax], _nk[_kmax], _nk[_kmax]);
  symsp_extract_mat(_Q_ori, _nmk[_kmax], Vxx[_kmax]);   
  m1 = m_resize(m1, _nk[_kmax], _nk[_kmax]);
  symsp_extract_mat(_CTC, _nmk[_kmax], m1);   
  Vxx[_kmax] = m_add(Vxx[_kmax], m1, Vxx[_kmax]);

  Formax(_kmax);
  cbx[_kmax] = m_copy1(ax, cbx[_kmax]);

  //   Backward run: calculate Vxx, Rux, Ryx for k=_kmax-1, ... ,0.

  for (k = _kmax-1; k >= 0; k--) {
    k1 = k+1; 
    kn = _nks[k];  
    knm = _nmk[k]; 

    //   Gxx, Gxu, Guu
    if ( !_a_sparse )
      FormGxx(k);
    else
      FormGxxSp(k);

    //    printf("vor Formax\n"); 
    //   ax, au, yny
    Formax(k);

    if ( cbx[k1]->m + au->m ) {
      //   m2 = cbu
      //      printf("vor cbu\n");
      if ( _a_sparse )
	m1 = mbsptr_mlt(cbx[k1], _AT, kn, knm+_nk[k], _mk[k], m1);
      else
	m1 = m_mlt(cbx[k1], fu[k], m1);
      m2 = m_concatc(au, m1, m2);

      // Generalized elimination 
      if ( (MAT *)Z == MNULL )
	Z = m_resize(Z, 1, 1);     // Z must be allocated for GE_QP()!
      GE_QP(m2, S[k], Z, PQ[k], _ge_tol);
      
      //   cbx
      if ( _a_sparse )
	m1 = mbsptr_mlt(cbx[k1], _AT, kn, knm, _nk[k], m1);
      else
	m1 = m_mlt(cbx[k1], fx[k], m1);
      m2 = m_concatc(ax, m1, m2);
      cbx[k] = mtrm_mlt(PQ[k], m2, cbx[k]);
    } else {
      S[k] = m_resize(S[k], _mk[k], 0);
      Z = m_resize(Z, _mk[k], 0);
      PQ[k] = m_resize(PQ[k], 0, 0);
      cbx[k] = m_resize(cbx[k], 0, _nk[k]);
    }

    if ( S[k]->n == 0 ) {

      //   CH_Guu for unconstrained u: CH_Guu is CHfactor!
      CH_Guu[k] = m_copy1(Guu[k], CH_Guu[k]);

      //   diagonal scaling
      sc[k] = v_diag(CH_Guu[k], sc[k]);
      for ( i = 0; i < (int) sc[k]->dim; i++ )
	if ( v_entry(sc[k], i) > 1.0 )
	  v_set_val(sc[k], i, 1.0/sqrt(v_entry(sc[k], i)));
	else
	  v_set_val(sc[k], i, 1.0);
      CH_Guu[k] = dm_mlt(CH_Guu[k], sc[k], CH_Guu[k]);
      CH_Guu[k] = md_mlt(CH_Guu[k], sc[k], CH_Guu[k]);
      //      printf("CH_Guu[%d]: ", k); m_output(CH_Guu[k]);
      //      CH_Guu[k] = CHfactor(CH_Guu[k]);
      CH_Guu_p[k] = px_resize(CH_Guu_p[k], CH_Guu[k]->m);
      CH_Guu_b[k] = px_resize(CH_Guu_b[k], CH_Guu[k]->m);
      CH_Guu[k] = BKPfactor(CH_Guu[k], CH_Guu_p[k], CH_Guu_b[k]);

      //   Rux, Ryx, Ruyb, Ryyb, ctx for unconstrained u
      //      Rux[k] = CHsolveMT(CH_Guu[k], Gxu[k], Rux[k]);
      //      printf("vor Rux\n"); 
      m1 = md_mlt(Gxu[k], sc[k], m1);
      //      m1 = m_copy1(Gxu[k], m1);
      Rux[k] = BKPsolveMT(CH_Guu[k], CH_Guu_p[k], CH_Guu_b[k], m1, Rux[k]);
      Rux[k] = dm_mlt(Rux[k], sc[k], Rux[k]);
      Ryx[k] = m_resize(Ryx[k], 0, _nk[k]);
      ctx[k] = m_resize(ctx[k], 0, _nk[k]);
    } else {

      //   ctx
      ctx[k] = m_resize(ctx[k], S[k]->n, _nk[k]);
      ctx[k] = m_move(cbx[k], 0, 0, S[k]->n, _nk[k], ctx[k], 0, 0);

      //   Rux, m2 = cttilx
      //      m1 = mmtr_mlt(d11[k], S[k], m1);
      //      m2 = mmtr_mlt(m1, Gxu[k], m2);
      //      m2 = m_sub(ctx[k], m2, m2);
      Rux[k] = m_mlt(S[k], ctx[k], Rux[k]);
      //      Ruyb[k] = m_mlt(S[k], d12[k], Ruyb[k]);

      //   cbx
      //      i = cbx[k]->m-S[k]->n;
      m1 = m_resize(m1, cbx[k]->m-S[k]->n, _nk[k]);
      m1 = m_move(cbx[k], S[k]->n, 0, m1->m, _nk[k], m1, 0, 0);
      cbx[k] = m_copy1(m1, cbx[k]);

      if ( (int) S[k]->n < _mk[k] ) {

	//   CH_Guu for constrained u
	m1 = m_mlt(Guu[k], Z, m1);
	CH_Guu[k] = mtrm_mlt(Z, m1, CH_Guu[k]);
	//   diagonal scaling
	sc[k] = v_diag(CH_Guu[k], sc[k]);
	for ( i = 0; i < (int) sc[k]->dim; i++ )
	  if ( v_entry(sc[k], i) > 1.0 )
	    v_set_val(sc[k], i, 1.0/sqrt(v_entry(sc[k], i)));
	  else
	    v_set_val(sc[k], i, 1.0);
	CH_Guu[k] = dm_mlt(CH_Guu[k], sc[k], CH_Guu[k]);
	CH_Guu[k] = md_mlt(CH_Guu[k], sc[k], CH_Guu[k]);
	Z = md_mlt(Z, sc[k], Z);
	//	CH_Guu[k] = CHfactor(CH_Guu[k]);
	//	m2 = CHsolveMT(CH_Guu[k], Z, m2);
	pivot = px_resize(pivot, CH_Guu[k]->m);
	blocks = px_resize(blocks, CH_Guu[k]->m);
	CH_Guu[k] = BKPfactor(CH_Guu[k], pivot, blocks);
	m2 = BKPsolveMT(CH_Guu[k], pivot, blocks, Z, m2);
        CH_Guu[k] = m_mlt(Z, m2, CH_Guu[k]);

	//   Rux for constrained u
	m1 = m_mlt(Guu[k], Rux[k], m1);
	m2 = m_transp(Gxu[k], m2);
	m1 = m_sub(m2, m1, m1);
	m2 = m_mlt(CH_Guu[k], m1, m2);
	Rux[k] = m_add(m2, Rux[k], Rux[k]);
      } 

      //   Ryx
      m1 = m_transp(Gxu[k], m1);
      m2 = m_mlt(Guu[k], Rux[k], m2);
      m1 = m_sub(m1, m2, m1);
      Ryx[k] = mtrm_mlt(S[k], m1, Ryx[k]);
    }

    //   Vxx
    m1 = m_mlt(Gxu[k], Rux[k], m1);
    Vxx[k] = m_sub(Gxx, m1, Vxx[k]);
    if (Ryx[k]->m > 0) {
      m1 = mtrm_mlt(ctx[k], Ryx[k], m1);
      Vxx[k] = m_sub(Vxx[k], m1, Vxx[k]);
    }

    //   check for symmetry of Vxx
    //    printf("vor rel_symm \n"); 
    if ( rel_symm(Vxx[k])  > 1e-8 ) {
      //      printf("symmetry check Vxx[%d]: %g\n", k, rel_symm(Vxx[k]));
      if ( _logging && rel_symm(Vxx[k])  > 1e-3 ) {
	fprintf(stderr,
		"LQDOCP: symmetry check Vxx[%d]: %g\n", k, rel_symm(Vxx[k]));
	//	m_output(Vxx[k]);
	//	dump();
	//	m_error(E_POSDEF,"ExRiccatiFactor");
      }
    }
    //    printf("vor m_symm \n"); 
    Vxx[k] = m_symm(Vxx[k], Vxx[k]);

    //    printf("Vxx[%d]", k); m_output(Vxx[k]);
    //    printf("Rux[%d]", k); m_output(Rux[k]);

    //    printf("k = %d, VEC: %d, MAT: %d\n",k,mem_info_numvar(TYPE_VEC,0),
    //       	   mem_info_numvar(TYPE_MAT,0)); 

  }

  //   Initial state
  if ( _fixed_x0 ) {
    //    m1 = m_copy1(db, m1);
    m1 = m_resize(m1, 0, 0);
  } else {
    m1 = m_resize(m1, _nk[0]+cbx[0]->m, _nk[0]+cbx[0]->m);
    m1 = m_zero(m1);
    m1 = m_move(Vxx[0], 0, 0, _nk[0], _nk[0], m1, 0, 0);
    m1 = m_move(cbx[0], 0, 0, cbx[0]->m, _nk[0], m1, _nk[0], 0);
    m2 = m_transp(cbx[0], m2);
    m1 = m_move(m2, 0, 0, _nk[0], cbx[0]->m, m1, 0, _nk[0]);
    //    m1 = m_move(db, 0, 0, db->m, db->n, m1, _nk[0], _nk[0]);
  }
  //   diagonal scaling
  scx0 = v_diag(m1, scx0); 
  for ( i = 0; i < (int) scx0->dim; i++ )
    if ( v_entry(scx0, i) > 1.0 )
      v_set_val(scx0, i, 1.0/sqrt(v_entry(scx0, i)));
    else 
      v_set_val(scx0, i, 1.0);
  m1 = dm_mlt(m1, scx0, m1);
  m1 = md_mlt(m1, scx0, m1);

  pivot = px_resize(pivot, m1->m);
  blocks = px_resize(blocks, m1->m);
  m1 = BKPfactor(m1, pivot, blocks);
  //   do not change m1, pivot, blocks!

}

//--------------------------------------------------------------------------
//   Extended Riccati solution: solve routine.
//   Diagonal scaled version without inequalities (_wz_tol = HUGE_VAL).
//
//   E. Arnold   02/22/97  
//--------------------------------------------------------------------------
void Hqp_IpLQDOCP::ExRiccatiSolveSc(void)
{
  int k, k1, knm, kn;

  //   Backward run for Vx, Ru

  Vx[_kmax] = v_copy1(gx[_kmax], Vx[_kmax]);
  cb[_kmax] = v_copy1(a[_kmax], cb[_kmax]);

  for (k = _kmax-1; k >= 0; k--) {
    k1 = k+1;

    //   printf("k = %d, start\n", k);
    //   Gx, Gu
    if ( !_a_sparse )
      FormGx(k);
    else
      FormGxSp(k);

    //    printf("k = %d, cb\n", k);
    //   cb = PQ*[a; cb*f] (remove fixed initial states)
    v1 = mv_mlt(cbx[k1], f[k], v1);
    v1 = v_add(cb[k1], v1, v1);
    v2 = v_concat(a[k], v1, v2);
    if ( ( k == 0 ) && ( _fixed_x0 ) ) {
      v1 = v_resize(v1, v2->dim-_nk[k]);
      v1 = v_move(v2, _nk[k], v1->dim, v1, 0);
      v2 = v_copy1(v1, v2);
    }
    cb[k] = vm_mlt(PQ[k], v2, cb[k]);

    if ( S[k]->n == 0 ) {
      //      printf("k = %d, Ru,Ry\n", k);
      //   Ru, Ry for unconstrained u
      //      Ru[k] = CHsolve(CH_Guu[k], Gu, Ru[k]);

      v1 = v_star(sc[k], Gu, v1);
      Ru[k] = BKPsolve(CH_Guu[k], CH_Guu_p[k], CH_Guu_b[k], v1, Ru[k]);
      Ru[k] = v_star(Ru[k], sc[k], Ru[k]);
      Ry[k] = v_resize(Ry[k], 0);
    } else {

      //   ct
      //      i = cb[k]->dim-S[k]->n;
      ct = v_resize(ct, S[k]->n);
      ct = v_move(cb[k], 0, S[k]->n, ct, 0);

      //      printf("k = %d, ct\n", k);
      //   v1 = cttil
      //      v1 = vm_mlt(S[k], Gu, v1);
      //      v2 = mv_mlt(d11[k], v1, v2);
      //      v1 = v_sub(ct, v2, v1);

      //      printf("k = %d, Ru\n", k);
      //   Ru for constrained u
      Ru[k] = mv_mlt(S[k], ct, Ru[k]);

      //      printf("k = %d, cb\n", k);
      //   cb
      v1 = v_resize(v1, cb[k]->dim-S[k]->n);
      v1 = v_move(cb[k], S[k]->n, v1->dim, v1, 0);
      cb[k] = v_copy1(v1, cb[k]);

      if ( (int) S[k]->n < _mk[k] ) {
	//	printf("k = %d, Ru2\n", k);
	//   Ru for constrained u
	v1 = mv_mlt(Guu[k], Ru[k], v1);
	v1 = v_sub(Gu, v1, v1);
	v2 = mv_mlt(CH_Guu[k], v1, v2);
	Ru[k] = v_add(v2, Ru[k], Ru[k]);
      } 

      //      printf("k = %d, Ry\n", k);
      //   Ry for constrained u
      v1 = mv_mlt(Guu[k], Ru[k], v1);
      v1 = v_sub(Gu, v1, v1);
      Ry[k] = vm_mlt(S[k], v1, Ry[k]);
    }

    //    printf("k = %d, Vx\n", k);
    //   Vx
    v1 = mv_mlt(Gxu[k], Ru[k], v1);
    Vx[k] = v_sub(Gx, v1, Vx[k]);
    if ( Ry[k]->dim > 0 ) {
      v1 = vm_mlt(ctx[k], Ry[k], v1);
      Vx[k] = v_sub(Vx[k], v1, Vx[k]);
    }
  }

  //   Initial state
  x[0] = v_resize(x[0], _nk[0]);
  if ( _fixed_x0 ) {
    x[0] = v_move(a[0], 0, _nk[0], x[0], 0);
    x[0] = sv_mlt(-1.0, x[0], x[0]);
    if ( cb[0]->dim > 0 ) {   
      v1 = v_resize(v1, cb[0]->dim);
      v1 = mv_mltadd(cb[0], x[0], cbx[0], 1.0, v1);
      v1 = sv_mlt(-1.0, v1, v1);
      v1 = v_star(v1, scx0, v1);
      yb = BKPsolve(m1, pivot, blocks, v1, yb);
      yb = v_star(yb, scx0, yb);
    } else
      yb = v_resize(yb, 0);
  } else {
    v1 = v_concat(Vx[0], cb[0], v1);
    v1 = sv_mlt(-1.0, v1, v1);
    v1 = v_star(v1, scx0, v1);
    v1 = BKPsolve(m1, pivot, blocks, v1, v1);
    v1 = v_star(v1, scx0, v1);
    x[0] = v_move(v1, 0, _nk[0], x[0], 0);
    yb = v_resize(yb, v1->dim-_nk[0]);
    yb = v_move(v1, _nk[0], v1->dim-_nk[0], yb, 0);
  }

  //   Forward run for x, u, p

  for (k = 0; k < _kmax; k++) {

    x[k+1] = v_resize(x[k+1], _nk[k+1]);
    u[k] = v_resize(u[k], _mk[k]);
    p[k] = v_resize(p[k], _nk[k+1]);

    //   u = -Rux*x -Ru 
    v1 = mv_mlt(Rux[k], x[k], v1);
    v1 = v_add(v1, Ru[k], v1);
    u[k] = sv_mlt(-1.0, v1, u[k]);

    //   x[k+1] = fx*x+fu*u+f
    k1 = k+1;
    knm = _nmk[k]; 
    kn = _nks[k];  
    if ( _a_sparse )
      v1 = bspv_mlt(_A_ori, kn, knm, _nk[k+1], x[k], v1);
    else
      v1 = mv_mlt(fx[k], x[k], v1);
    x[k1] = v_add(f[k], v1, x[k1]);
    if ( _a_sparse )
      v1 = bspv_mlt(_A_ori, kn, knm+_nk[k], _nk[k+1], u[k], v1);
    else
      v1 = mv_mlt(fu[k], u[k], v1);
    x[k1] = v_add(x[k1], v1, x[k1]);

    //   [y; yb] = -(Ryx*x + Ry)
    if ( ( Ry[k]->dim+yb->dim > 0 ) || ( ( k == 0 ) && ( _fixed_x0 ) ) ) {
      v1 = mv_mlt(Ryx[k], x[k], v1);
      v1 = v_add(v1, Ry[k], v1);
      v1 = sv_mlt(-1.0, v1, v1);
      v2 = v_concat(v1, yb, v2);
      v1 = mv_mlt(PQ[k], v2, v1);
      if ( ( k == 0 ) && ( _fixed_x0 ) ) {
	v2 = v_resize(v2, _nk[k]);
	v2 = vm_mltadd(Vx[0], yb, cbx[0], 1.0, v2);
	v2 = mv_mltadd(v2, x[0], Vxx[0], 1.0, v2);
	v2 = sv_mlt(-1.0, v2, v2);
	v2 = v_concat(v2, v1, v2);
	v1 = v_copy1(v2, v1);
      }
      y[k] = v_resize(y[k], a[k]->dim);
      y[k] = v_move(v1, 0, a[k]->dim, y[k], 0);
      yb = v_resize(yb, v1->dim-y[k]->dim);
      yb = v_move(v1, a[k]->dim, v1->dim-a[k]->dim, yb, 0);
    }

    //   p = Vxx*x + Vx + cbx*yb
    p[k] = mv_mlt(Vxx[k1], x[k1], p[k]);
    p[k] = v_add(p[k], Vx[k1], p[k]);
    if ( yb->dim > 0 ) {
      v1 = vm_mlt(cbx[k1], yb, v1);
      p[k] = v_add(p[k], v1, p[k]);
    }
  }

  //   y[_kmax] = yb
  y[_kmax] = v_resize(y[_kmax], yb->dim);
  y[_kmax] = v_copy1(yb, y[_kmax]);
}
