/*
 * Prg_ASCEND.C -- class definition
 *
 * rf, 12/17/99
 */

/*
    Copyright (C) 1999  Ruediger Franke

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

#include "Prg_ASCEND.h"
#include "Hqp_Program.h"
#include "Hqp_SqpSolver.h"

#include <string.h>

#include <If_Class.h>
#include <If_Int.h>

IF_CLASS_DEFINE("ASCEND", Prg_ASCEND, Hqp_SqpProgram);

//-------------------------------------------------------------------------
Prg_ASCEND::Prg_ASCEND()
{
  _nvars = 0;
  _nrels = 0;
  _vars = NULL;
  _rels = NULL;
  _obj = NULL;
  _safe_calc = TRUE;
  _var_lb = v_resize(v_get(1), 0);
  _var_ub = v_resize(v_get(1), 0);
  _var_asc2hqp = iv_resize(iv_get(1), 0);
  _derivatives = v_resize(v_get(1), 0);
  _var_master_idxs = iv_resize(iv_get(1), 0);
  _var_solver_idxs = iv_resize(iv_get(1), 0);
  _Inf = 1e38;	// todo: should be initialized from ASCEND data
  _me_bounds = 0;
  _m_bounds = 0;
  memset(&_slv_status, 0, sizeof(slv_status_t));

  // create interface elements
  _ifList.append(new If_Int("prg_safe_calc", &_safe_calc));

  slv_register_client(Prg_ASCEND::slv_register, NULL, NULL);
}

//-------------------------------------------------------------------------
Prg_ASCEND::~Prg_ASCEND()
{
  iv_free(_var_solver_idxs);
  iv_free(_var_master_idxs);
  v_free(_derivatives);
  iv_free(_var_asc2hqp);
  v_free(_var_ub);
  v_free(_var_lb);
}

//-------------------------------------------------------------------------
int Prg_ASCEND::slv_register(SlvFunctionsT *sft)
{
  if (sft == NULL)  {
    FPRINTF(stderr,"Prg_ASCEND::slv_register called with NULL pointer\n");
    return 1;
  }
  fprintf(stderr, "Prg_ASCEND::slv_register called\n");

  sft->name = "HQP";
  sft->ccreate = Prg_ASCEND::slv_create;
  sft->cdestroy = Prg_ASCEND::slv_destroy;
  sft->getstatus = Prg_ASCEND::slv_get_status;
  sft->presolve = Prg_ASCEND::slv_presolve;
  //sft->solve = Prg_ASCEND::slv_solve;
  sft->iterate = Prg_ASCEND::slv_iterate;
  return 0;
}

//-------------------------------------------------------------------------
SlvClientToken Prg_ASCEND::slv_create(slv_system_t, int *status_index)
{
  extern Hqp_SqpSolver *theSqpSolver;  
  fprintf(stderr, "Prg_ASCEND::slv_create called\n");
  *status_index = 0;
  return (SlvClientToken)theSqpSolver;
}

//-------------------------------------------------------------------------
int Prg_ASCEND::slv_destroy(slv_system_t, SlvClientToken)
{
  fprintf(stderr, "Prg_ASCEND::slv_destroy called\n");
  return 0;
}

//-------------------------------------------------------------------------
void Prg_ASCEND::slv_get_status(slv_system_t, SlvClientToken clt,
				slv_status_t *status)
{
  //fprintf(stderr, "Prg_ASCEND::slv_get_status called\n");

  Hqp_SqpSolver *sqp = (Hqp_SqpSolver *)clt;
  Prg_ASCEND *prg = (Prg_ASCEND *)sqp->prg();
  // todo: assert prg != NULL && prg->name() == "ASCEND"

  *status = prg->_slv_status;
}

//-------------------------------------------------------------------------
void Prg_ASCEND::slv_presolve(slv_system_t system, SlvClientToken clt)
{
  fprintf(stderr, "Prg_ASCEND::slv_presolve called\n");

  Hqp_SqpSolver *sqp = (Hqp_SqpSolver *)clt;
  Prg_ASCEND *prg = (Prg_ASCEND *)sqp->prg();
  // todo: assert prg != NULL && prg->name() == "ASCEND"

  prg->_slv_system = system;
  prg->setup();
  prg->init_x();

  sqp->init();

  prg->_slv_status.calc_ok = TRUE;
  prg->_slv_status.ok = TRUE;
  prg->_slv_status.ready_to_solve = TRUE;
  prg->_slv_status.iteration = 0;
}

//-------------------------------------------------------------------------
void Prg_ASCEND::slv_iterate(slv_system_t, SlvClientToken clt)
{
  //fprintf(stderr, "Prg_ASCEND::slv_iterate called\n");

  Hqp_SqpSolver *sqp = (Hqp_SqpSolver *)clt;
  Prg_ASCEND *prg = (Prg_ASCEND *)sqp->prg();
  // todo: assert prg != NULL && prg->name() == "ASCEND"

  // currently call Tcl command hqp_solve
  // todo: replace with calls to sqp->{qp_update, qp_solve, step}()
  extern Tcl_Interp *theInterp;
  char *tcl_channel = "stdout";
  sqp->set_max_iters(sqp->iter() + 1);	// allow one iteration
  if (Tcl_VarEval(theInterp, "hqp_solve ", tcl_channel, NULL) == TCL_OK) {
    prg->_slv_status.converged = TRUE;
    prg->_slv_status.ready_to_solve = FALSE;
  }
  else if (strcmp(theInterp->result, "iters") != 0) {
    prg->_slv_status.ok = FALSE;
    prg->_slv_status.ready_to_solve = FALSE;
    Tcl_VarEval(theInterp, "puts ", tcl_channel,
		" \"HQP Error: ", theInterp->result, "\"", NULL);
  }

  prg->_slv_status.iteration = sqp->iter();
}

//-------------------------------------------------------------------------
void Prg_ASCEND::slv_solve(slv_system_t, SlvClientToken clt)
{
  fprintf(stderr, "Prg_ASCEND::slv_solve called\n");

  Hqp_SqpSolver *sqp = (Hqp_SqpSolver *)clt;
  Prg_ASCEND *prg = (Prg_ASCEND *)sqp->prg();
  // todo: assert prg != NULL && prg->name() == "ASCEND"

  // currently call Tcl command hqp_solve
  // todo: replace with call to sqp->solve()
  extern Tcl_Interp *theInterp;
  Tcl_Eval(theInterp, "hqp_solve");

  prg->_slv_status.ready_to_solve = FALSE;
  prg->_slv_status.iteration = sqp->iter();
}

//-------------------------------------------------------------------------
void Prg_ASCEND::setup()
{
  int n, me, m;
  int i, j, row_idx, idx;
  SPMAT *J;
  int nincidences;
  const struct var_variable **incidences;

  // obtain ASCEND system
  // todo: should check that system can be solved with HQP (e.g. no integers)
  _nvars = slv_get_num_solvers_vars(_slv_system);
  _vars = slv_get_solvers_var_list(_slv_system);
  _nrels = slv_get_num_solvers_rels(_slv_system);
  _rels = slv_get_solvers_rel_list(_slv_system);
  _obj = slv_get_obj_relation(_slv_system);

  // count number of optimization variables and bounds
  _var_lb = v_resize(_var_lb, _nvars);
  _var_ub = v_resize(_var_ub, _nvars);
  _var_asc2hqp = iv_resize(_var_asc2hqp, _nvars);
  _derivatives = v_resize(_derivatives, _nvars);
  _var_master_idxs = iv_resize(_var_master_idxs, _nvars);
  _var_solver_idxs = iv_resize(_var_solver_idxs, _nvars);
  n = 0;
  me = 0;
  m = 0;
  for (i = 0; i < _nvars; i++) {
    _var_lb[i] = var_lower_bound(_vars[i]); 
    _var_ub[i] = var_upper_bound(_vars[i]);
    /*
    var_write_name(_slv_system, _vars[i], stderr);
    fprintf(stderr, ":\t%i,\t%g,\t%g\n", var_fixed(_vars[i]),
            _var_lb[i], _var_ub[i]);
    */
    if (var_fixed(_vars[i])) {
      _var_asc2hqp[i] = -1;
    }
    else {
      _var_asc2hqp[i] = n++;
      if (_var_lb[i] == _var_ub[i])
	++me;
      else {
	if (_var_lb[i] > -_Inf)
	  ++m;
	if (_var_ub[i] < _Inf)
	  ++m;
      }
    }
  }

  // consider bounds as linear constraints (i.e. no Jacobian update)
  _me_bounds = me;
  _m_bounds = m;

  // count number of HQP constraints
  for (i = 0; i < _nrels; i++) {
    if (rel_equal(_rels[i]))
      ++me;	// equality constraint
    else
      ++m;	// inequality constraint
  }

  // allocate QP approximation and optimization variables vector
  _qp->resize(n, me, m);
  _x = v_resize(_x, n);

  // allocate sparse structure for bounds
  // (write constant elements in Jacobians)
  me = m = 0;
  for (i = 0; i < _nvars; i++) {
    idx = _var_asc2hqp[i];
    if (idx < 0)
      continue;
    if (_var_lb[i] == _var_ub[i]) {
      row_idx = me++;
      sp_set_val(_qp->A, row_idx, idx, 1.0);
    }
    else {
      if (_var_lb[i] > -_Inf) {
	row_idx = m++;
	sp_set_val(_qp->C, row_idx, idx, 1.0);
      }
      if (_var_ub[i] < _Inf) {
	row_idx = m++;
	sp_set_val(_qp->C, row_idx, idx, -1.0);
      }
    }
  }
  
  // allocate sparse structure for general constraints
  // (just insert dummy values; actual values are set in update method)
  for (i = 0; i < _nrels; i++) {
    if (rel_equal(_rels[i])) {
      row_idx = me++;
      J = _qp->A;
    }
    else {
      row_idx = m++;
      J = _qp->C;
    }
    nincidences = rel_n_incidences(_rels[i]);
    incidences = rel_incidence_list(_rels[i]);
    for (j = 0; j < nincidences; j++) {
      idx = _var_asc2hqp[var_sindex(incidences[j])];
      if (idx >= 0)
	sp_set_val(J, row_idx, idx, 1.0);
    }      
  }

  // todo: setup sparse structure of Hessian
  // for now initialize something resulting in dense BFGS update
  for (j = 0; j < n-1; j++) {
    sp_set_val(_qp->Q, j, j, 0.0);
    sp_set_val(_qp->Q, j, j+1, 0.0);
  }
  sp_set_val(_qp->Q, j, j, 0.0);
}

//-------------------------------------------------------------------------
void Prg_ASCEND::init_x()
{
  int i, idx;
  for (i = 0; i < _nvars; i++) {
    idx = _var_asc2hqp[i];
    if (idx >= 0)
      _x[idx] = var_value(_vars[i]);
  }
}

//-------------------------------------------------------------------------
void Prg_ASCEND::update_bounds()
{
  int j;
  int m, me, idx, row_idx;

  m = me = 0;
  for (j = 0; j < _nvars; j++) {
    idx = _var_asc2hqp[j];
    if (idx < 0)
      continue;
    if (_var_lb[j] == _var_ub[j]) {
      row_idx = me++;
      _qp->b[row_idx] = _x[idx] - _var_lb[j];
    }
    else {
      if (_var_lb[j] > -_Inf) {
	row_idx = m++;
	_qp->d[row_idx] = _x[idx] - _var_lb[j];
      }
      if (_var_ub[j] < _Inf) {
	row_idx = m++;
	_qp->d[row_idx] = _var_ub[j] - _x[idx];
      }
    }
  }
}

//-------------------------------------------------------------------------
void Prg_ASCEND::update_fbd()
{
  int i, idx, row_idx;
  int m, me;
  int32 calc_ok;
  double residual, mult;

  // write current optimization variables to ASCEND
  for (i = 0; i < _nvars; i++) {
    idx = _var_asc2hqp[i];
    if (idx >= 0)
      var_set_value(_vars[i], _x[idx]);
  }

  // update variable bounds
  update_bounds();

  // update objective
  _f = relman_eval(_obj, &calc_ok, _safe_calc);
  _slv_status.calc_ok &= calc_ok;

  // update general constraints
  me = _me_bounds;
  m = _m_bounds;
  for (i = 0; i < _nrels; i++) {
    residual = relman_eval(_rels[i], &calc_ok, _safe_calc);
    _slv_status.calc_ok &= calc_ok;
    if (rel_equal(_rels[i])) {
      row_idx = me++;
      _qp->b[row_idx] = residual;
    }
    else {
      row_idx = m++;
      mult = rel_less(_rels[i])? -1.0: 1.0;
      _qp->d[row_idx] = mult * residual;
    }
  }
}

//-------------------------------------------------------------------------
void Prg_ASCEND::update(const VECP y, const VECP z)
{
  var_filter_t vfilter;
  int i, j, row_idx, idx;
  int m, me;
  int32 count;
  int status, newel;
  real64 residual;
  double mult;
  SPMAT *J;

  // write current optimization variables to ASCEND
  for (j = 0; j < _nvars; j++) {
    idx = _var_asc2hqp[j];
    if (idx >= 0)
      var_set_value(_vars[j], _x[idx]);
  }

  // update variable bounds
  update_bounds();

  // initialize vfilter so that all optimization variables are considered
  vfilter.matchbits = (VAR_ACTIVE | VAR_INCIDENT | VAR_SVAR | VAR_FIXED);
  vfilter.matchvalue = (VAR_ACTIVE | VAR_INCIDENT | VAR_SVAR);

  // update objective and derivatives
  status = relman_diff_grad(_obj, &vfilter, _derivatives->ve,
			    _var_master_idxs->ive, _var_solver_idxs->ive,
			    &count, &residual, _safe_calc);
  if (status == 0) // todo: this contradicts docu in solver/relman.h!?
    _slv_status.calc_ok &= FALSE;
  _f = residual;
  for (j = 0; j < count; j++) {
    idx = _var_asc2hqp[_var_solver_idxs[j]];
    if (idx >= 0)
      _qp->c[idx] = _derivatives[j];
  }

  // update general constraints
  me = _me_bounds;
  m = _m_bounds;
  newel = 0;
  for (i = 0; i < _nrels; i++) {
    status = relman_diff_grad(_rels[i], &vfilter, _derivatives->ve,
			      _var_master_idxs->ive, _var_solver_idxs->ive,
			      &count, &residual, _safe_calc);
    if (status == 0) // todo: this contradicts docu in solver/relman.h!?
      _slv_status.calc_ok &= FALSE;
    if (rel_equal(_rels[i])) {
      row_idx = me++;
      mult = 1.0;
      _qp->b[row_idx] = residual;
      J = _qp->A;
    }
    else {
      row_idx = m++;
      mult = rel_less(_rels[i])? -1.0: 1.0;
      _qp->d[row_idx] = mult * residual;
      J = _qp->C;
    }
    sprow_zero(J->row + row_idx);
    for (j = 0; j < count; j++) {
      idx = _var_asc2hqp[_var_solver_idxs[j]];
      if (idx >= 0)
	newel |= sp_update_val(J, row_idx, idx, mult * _derivatives[j]);
    }      
  }
  if (newel) {
    // the sparsity structure has changed
    _slv_status.ok = FALSE;
    _slv_status.inconsistent = TRUE;
    _slv_status.ready_to_solve = FALSE;
  }
}

//=========================================================================
