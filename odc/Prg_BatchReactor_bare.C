/*
 * Prg_BatchReactor_bare.C -- class definition
 *
 * rf, 13/12/16
 */

#include "Prg_BatchReactor_bare.h"

#include <If_Real.h>

IF_CLASS_DEFINE("BatchReactor_bare", Prg_BatchReactor_bare, Omu_Program);

//--------------------------------------------------------------------------
Prg_BatchReactor_bare::Prg_BatchReactor_bare()
{
  set_K(40); // can be modified via prg_K
  _kinf = 0.5;

  // interface elements
  _ifList.append(new If_Real("prg_kinf", &_kinf));
}

//--------------------------------------------------------------------------
void Prg_BatchReactor_bare::setup_stages(IVECP ks, VECP ts)
{
  stages_alloc(ks, ts, K(), 1, 0.0, 1.0);
}

//--------------------------------------------------------------------------
void Prg_BatchReactor_bare::setup_struct(int k,
					 const Omu_VariableVec &x,
					 const Omu_VariableVec &u,
					 Omu_DependentVec &xt, Omu_DependentVec &F,
					 Omu_DependentVec &f,
					 Omu_Dependent &f0, Omu_DependentVec &c)
{
  // consistic just takes states from optimizer
  m_ident(xt.Jx);
  m_zero(xt.Ju);
  xt.set_linear();

  // continuous provides explicit ODE
  if (k < _K) {
    m_zero(F.Jdx);
    for (unsigned int i = 0; i < F.Jdx->m; i++)
      F.Jdx[i][i] = -1.0;
  }
  F.set_linear(Omu_Dependent::WRT_dx);

  // update just takes final states from continuous
  m_zero(f.Ju);
  m_zero(f.Jx);
  m_ident(f.Jxf);
  f.set_linear();

  // linear objective
  v_zero(f0.gx);
  v_zero(f0.gxf);
  if (k < _K) {
    v_zero(f0.gu);
  }
  else
    f0.gx[1] = -1.0;
  f0.set_linear();
}

//--------------------------------------------------------------------------
void Prg_BatchReactor_bare::setup(int k,
				  Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  x.alloc(2);
  if (k < K())
    u.alloc(1);

  // initial control inputs and bounds
  if (k < K()) {
    u.initial[0] = 1.0;
    u.min[0] = 0.0;
    u.max[0] = 5.0;
  }
  // initial states and initial state constraints
  if (k == 0) {
    x.min[0] = x.max[0] = x.initial[0] = 1.0;
    x.min[1] = x.max[1] = x.initial[1] = 0.0;
  }
  // path constraint
  else {
    x.initial[0] = 0.5;
    x.initial[1] = 0.5;
    x.min[0] = x.min[1] = 0.0;
    x.max[1] = x.max[1] = 1.0;
  }
}

//--------------------------------------------------------------------------
void Prg_BatchReactor_bare::consistic(int kk, double t,
				      const Omu_StateVec &x, const Omu_Vec &u,
				      Omu_DependentVec &xt)
{
  v_copy(x, xt);
}

//--------------------------------------------------------------------------
void Prg_BatchReactor_bare::continuous(int kk, double t,
				       const Omu_StateVec &x, const Omu_Vec &u,
				       const Omu_StateVec &dx, Omu_DependentVec &F)
{
  F[0] = -(u[0] + _kinf*u[0]*u[0])*x[0] - dx[0];
  F[1] = u[0]*x[0] - dx[1];

  // provide Jacobians if required
  if (F.is_required_J()) {
    F.Jx[0][0] =  -u[0] - _kinf*u[0]*u[0];
    F.Jx[1][0] =  u[0];
    F.Jx[0][1] =  0.0;
    F.Jx[1][1] =  0.0;
    F.Ju[0][0] = -x[0] - 2.0*_kinf*u[0]*x[0];
    F.Ju[1][0] =  x[0];
  }
}

//--------------------------------------------------------------------------
void Prg_BatchReactor_bare::update(int kk,
				   const Omu_StateVec &x, const Omu_Vec &u,
				   const Omu_StateVec &xf,
				   Omu_DependentVec &f, Omu_Dependent &f0,
				   Omu_DependentVec &c)
{
  if (kk < _KK)
    v_copy(xf, f);
  else
    f0 = -x[1];
}

//==========================================================================
