/*
 * Hqp_IpLQDOCP.h --
 *   - Solution of the Newton step equations for the Interior Point 
 *     algorithm.
 *   - Extended Riccati method for equality constrained linear-quadratic
 *     discrete-time optimal control problems.
 *
 * E. Arnold  09/18/96
 *            03/04/97  iterative improvement in base class
 */

/*
    Copyright (C) 1996--1998  Eckhard Arnold

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

#ifndef Hqp_IpLQDOCP_H
#define Hqp_IpLQDOCP_H

#include "Hqp_IpMatrix.h"
#include "Meschach.h"

#include "t_mesch.h"
#include "meschext_ea.h"

//----------   Hqp_IpLQDOCP   ----------

class Hqp_IpLQDOCP: public Hqp_IpMatrix {

 protected:

  int _kmax; 
  int _kmax1;
  Real _wz_tol;
  Real _ge_tol;
  int _a_sparse;
  int _fixed_x0;
  int _logging;

  TVECP x; 
  TVECP u; 
  TVECP p; 
  TVECP gx; 
  TVECP gu; 
  TVECP f; 
  TVECP Ru; 
  TVECP Vx;
  
  TMATP fx; 
  TMATP fu; 
  TMATP Rux; 
  TMATP Vxx; 
  TMATP CH_Guu; 
  TPERMP CH_Guu_p; 
  TPERMP CH_Guu_b; 
  IVECP CH_Guu_nr;
  TMATP Gxu; 
  TMATP Guu;
  
  MATP m1; 
  MATP m2; 
  MATP Gxx;
  VECP v1; 
  VECP v2; 
  VECP Gx; 
  VECP Gu;
  
  TVECP a; 
  TVECP Ry; 
  TVECP y; 
  TVECP cb;
  MATP ax;
  MATP au;
  TMATP PQ; 
  TMATP Ryx; 
  TMATP cbx; 
  TMATP ctx; 
  TMATP S;
  MATP Z;
  
  PERM *pivot;
  PERM *blocks;
  VECP ct; 
  VECP yb;

  IVECP _ra;
  IVECP _rc;
  IVECP _rak;
  IVECP _rck;
  IVECP _rcka;
  IVECP _rca;
  TIVECP _raki;
  TIVECP _rcki;
  TIVECP _rckai;

  SPMAT *_Q_ori;
  SPMAT *_A_ori;
  SPMAT *_C_ori;
  VECP  _Z_ori;
  VECP  _W_ori;
  VECP  _WZ_ori;
  SPMAT *_CT;
  SPMAT *_CTC;
  SPMAT *_QCTC;
  SPMAT *_AT;

  VECP yny;
  MATP db;
  TMATP d11;
  TMATP d12;
  MATP d22;
  TMATP Ruyb;
  TMATP Ryyb;
  
  TMATP CD;
  TPERMP CD_p;

  //   residuum
  VECP _res1;
  VECP _res2;
  VECP _res3;
  VECP _res4;

  //   different number of variables per stage  02/22/97
  IVECP _nk;
  IVECP _mk;
  IVECP _nmk;
  IVECP _nks;
  IVECP _indk;

  //   diagonal scaling
  VECP scx0;
  TVECP sc;

 public:
  Hqp_IpLQDOCP();
  ~Hqp_IpLQDOCP();
  
  int Get_Dim(const Hqp_Program *);
  int Check_Structure(const Hqp_Program *);
  void Get_Constr_Dim(const Hqp_Program *);
  void resize();
  void free();
  void dump(char *);
  void dump() { dump("dump.dat"); }
  void init(const Hqp_Program *);
  void update(const Hqp_Program *);
  void factor(const Hqp_Program *, const VEC *, const VEC *);
  void step(const Hqp_Program *, const VEC *, const VEC *,
            const VEC *, const VEC *, const VEC *, const VEC *,
            VEC *, VEC *, VEC *, VEC *);

  void Residuum(const VEC *, const VEC *, const VEC *, const VEC *,
		const VEC *, const VEC *, const VEC *, const VEC *, const int);

  void FormGxx(int);
  void FormGxxSp(int);
  void FormGx(int);
  void FormGxSp(int);
  void Formax(int);
  void ExRiccatiFactor(void);
  void ExRiccatiSolve(void);
  void ExRiccatiFactorSc(void);
  void ExRiccatiSolveSc(void);

  char	*name() {return "LQDOCP";}
};

#endif
