/*
 * Hqp_LPSolve.C -- class definition
 *
 * rf, 4/11/09
 *
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


#include "Hqp_LPSolve.h"
#include "Hqp_SqpProgram.h"

// clear conflicting Meschach defines 
#undef REAL
#undef FLOAT
#undef DOUBLE
#undef set_row

namespace LPSOLVE {
#include <lp_lib.h>
}

#include <If_Int.h>
#include <If_Real.h>
#include <If_String.h>
#include <If_Method.h>

typedef If_Method<Hqp_LPSolve> If_Cmd;

#define GET_SET_CB(vartype, prefix, name) \
  GET_CB(vartype, prefix, name), \
  IF_SET_CB(vartype, Hqp_LPSolve, set_##name)

#define GET_CB(vartype, prefix, name) \
  prefix#name, \
  IF_GET_CB(vartype, Hqp_LPSolve, name)

IF_CLASS_DEFINE("LPSolve", Hqp_LPSolve, Hqp_MipSolver);

//--------------------------------------------------------------------------
void __WINAPI logfunction(lprec *lp, void *, char *buf)
{
  Tcl_VarEval(theInterp, "puts -nonewline {", buf, "}", NULL);
  Tcl_Eval(theInterp, "update idletasks");
}

//--------------------------------------------------------------------------
Hqp_LPSolve::Hqp_LPSolve()
{
  _lp = LPSOLVE::make_lp(0, 0);
  _timeout = LPSOLVE::get_timeout(_lp);
  _gap = LPSOLVE::get_mip_gap(_lp, TRUE);
  _logging = LPSOLVE::get_verbose(_lp);
  _dump_format = strdup("lp");

  _ifList.append(new If_Real(GET_SET_CB(double, "mip_", timeout)));
  _ifList.append(new If_Real(GET_SET_CB(double, "mip_", gap)));
  _ifList.append(new If_Int(GET_SET_CB(int, "mip_", logging)));
  _ifList.append(new If_String(GET_SET_CB(const char *, "mip_", dump_format)));
  _ifList.append(new If_Cmd("mip_dump", &Hqp_LPSolve::dump, this));
}

//--------------------------------------------------------------------------
Hqp_LPSolve::~Hqp_LPSolve()
{
  free(_dump_format);
  LPSOLVE::delete_lp(_lp);
}

//--------------------------------------------------------------------------
int Hqp_LPSolve::init(int, char *[], char **retString)
{
  if (!_prg) {
    m_error(E_NULL, "Hqp_LPSolve::init");
  }
  Hqp_Program *qp = _prg->qp();
  const VECP x = _prg->x();

  int ncols = qp->x->dim;
  int nrows = qp->b->dim + qp->d->dim;
  int *colno = NULL, i, j, jdx;
  REAL *row = NULL;
  int ret = 0;

  SPROW *sprow;
  row_elt *elt;
  double val;

  // update the qp
  VECP y = v_zero(v_get(qp->b->dim));
  VECP z = v_zero(v_get(qp->d->dim));
  _prg->update(y, z);
  v_free(z);
  v_free(y);

  // create new lp_solve model
  LPSOLVE::delete_lp(_lp);
  _lp = LPSOLVE::make_lp(0, ncols);
  if(_lp == NULL) {
    ret = 1;
    if (retString)
      *retString = "Failed to construct a linear program";
  }
  // initialize solver options
  if (ret == 0) {
    LPSOLVE::set_timeout(_lp, _timeout);
    LPSOLVE::set_mip_gap(_lp, TRUE, _gap);
    LPSOLVE::set_verbose(_lp, _logging);
  }
  LPSOLVE::put_logfunc(_lp, logfunction, NULL);

  // alloc memory for setup of lp
  if (ret == 0) {
    // allocate memory large enough for one row
    colno = (int *) malloc(ncols * sizeof(*colno));
    row = (REAL *) malloc(ncols * sizeof(*row));
    if((colno == NULL) || (row == NULL)) {
      ret = 2;
      if (retString)
        *retString = "Failed to allocate memory";
    }
  }

  // initialize variables
  //   - give names
  //   - remove default lower bounds of zero
  //   - specify integer variables, incl. weights for branch&bound
  for (i = 0; i < ncols; i++) {
    //LPSOLVE::set_col_name(_lp, i + 1, "x");
    LPSOLVE::set_unbounded(_lp, i + 1);
    if (qp->x_int[i] > 0)
      LPSOLVE::set_int(_lp, i + 1, TRUE);
    row[i] = (REAL)qp->x_int[i];
  }
  LPSOLVE::set_var_weights(_lp, row);

  // setup optimization objective and constraints

  // cx --> min
  if (ret == 0) {
    // set the objective function
    j = 0;
    for (i = 0; i < ncols; i++) {
      if (qp->c[i] != 0.0) {
        colno[j] = i + 1;
        row[j++] = qp->c[i];
      }
    }
    if(!LPSOLVE::set_obj_fnex(_lp, j, row, colno)) {
      ret = 4;
      if (retString)
        *retString = "Failed to set the objective function";
    }
  }

  if(ret == 0) {
    // specify minimization problem
    LPSOLVE::set_minim(_lp);
  }

  // add rows for constraints
  LPSOLVE::set_add_rowmode(_lp, TRUE);

  // Ax + b == 0
  if (ret == 0) {
    for (i = 0; i < (int)qp->b->dim; i++) {
      sprow = qp->A->row + i;
      val = -qp->b[i];
      jdx = 0;
      for (j = 0; j < (int)sprow->len; j++) {
        elt = sprow->elt + j;
        if (elt->val != 0.0) {
          colno[jdx] = elt->col + 1;
          row[jdx++] = elt->val;
          val += elt->val * x[elt->col];
        }
      }
      // add the row to lpsolve
      if(!LPSOLVE::add_constraintex(_lp, jdx, row, colno, EQ, val)) {
        ret = 3;
        if (retString)
          *retString = "Failed to add an equality constraint row";
        break;
      }
    }
  }

  // Cx + d >= 0
  // x_min <= x <= x_max
  if (ret == 0) {
    for (i = 0; i < (int)qp->d->dim; i++) {
      sprow = qp->C->row + i;
      val = -qp->d[i];
      jdx = 0;
      for (j = 0; j < (int)sprow->len; j++) {
        elt = sprow->elt + j;
        if (elt->val != 0.0) {
          colno[jdx] = elt->col + 1;
          row[jdx++] = elt->val;
          val += elt->val * x[elt->col];
        }
      }
      // add the row to lpsolve
      if (jdx == 1) {
        // just a bound
        if (row[0] > 0.0) {
          LPSOLVE::set_lowbo(_lp, colno[0],
                             max(val/row[0], LPSOLVE::get_lowbo(_lp, colno[0])));
        }
        else {
          LPSOLVE::set_upbo(_lp, colno[0],
                            min(val/row[0], LPSOLVE::get_upbo(_lp, colno[0])));
        }
      }
      else {
        // general constraint
        if(!LPSOLVE::add_constraintex(_lp, jdx, row, colno, GE, val)) {
          ret = 3;
          if (retString)
            *retString = "Failed to add an inequality constraint row";
          break;
        }
      }
    }
  }

  // stop adding rows
  LPSOLVE::set_add_rowmode(_lp, FALSE);

  // free allocated memory
  if(row != NULL)
    free(row);
  if(colno != NULL)
    free(colno);

  if (ret != 0)
    return IF_ERROR;

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_LPSolve::solve(int argc, char *argv[], char **retString)
{
  int ret;

  // init
  ret = init(argc, argv, retString);
  if (ret != IF_OK)
    return ret;

  // solve
  ret = LPSOLVE::solve(_lp);
  if(ret != OPTIMAL && ret != SUBOPTIMAL) {
    if (retString)
      *retString = "No solution found!";
    return IF_ERROR;
  }

  // read result
  LPSOLVE::get_variables(_lp, _prg->qp()->x->ve);
  _prg->set_x(_prg->qp()->x);

  if (retString) {
    if (ret == OPTIMAL) 
      *retString = "optimal";
    else
      *retString = "suboptimal";
  }
  return IF_OK;
}

//--------------------------------------------------------------------------
void Hqp_LPSolve::set_dump_format(const char *value)
{
  if (strcmp(value, "lp") != 0 && strcmp(value, "mps") != 0)
    m_error(E_FORMAT, "Hqp_LPSolve::set_dump_format(\"lp\" or \"mps\")");

  free(_dump_format);
  _dump_format = strdup(value);
}

//-------------------------------------------------------------------------
void Hqp_LPSolve::dump()
{
  // dump generated lp
  if (strcmp(_dump_format, "lp") == 0)
    LPSOLVE::write_lp(_lp, "mip_dump.lp");
  else if (strcmp(_dump_format, "mps") == 0)
    LPSOLVE::write_mps(_lp, "mip_dump.mps");
  else
    m_error(E_FORMAT, "Hqp_LPSolve::dump");
}

//==========================================================================
