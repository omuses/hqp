/*
 * Omu_IntOdeTs.h --
 *   -- class for integrating an autonomous Ode over a stage
 *   -- using Taylor series expansion of Ode
 *   -- derived from adolc driver routine forodec
 *
 * hl, 28/01/98
 */

/*
    Copyright (C) 1997--1999  Hartmut Linke

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

#ifndef Omu_IntOdeTs_H
#define Omu_IntOdeTs_H

#include "Omu_Integrator.h"

/**
 *  Solving autonomous ordinary differential equations using 
 *  Taylor series expansion. The implementation is adapted 
 *  from the driver routine forodec of the adol-c 
 *  automatic differentiation code.
 */

//--------------------------------------------------------------------------
class Omu_IntOdeTs: public Omu_Integrator {

 public:

  Omu_IntOdeTs(int ); 	///< constructor
  Omu_IntOdeTs(); 	///< constructor
  ~Omu_IntOdeTs(); 	///< destructor

  /**
   * @name Implementation of predefined methods.
   * @see Omu_Integrator
   */

  //@{

  char *name() {return "OdeTs";}

  void init_stage(int k, const Omu_VariableVec &x, 
		  const Omu_VariableVec &u,
		  const Omu_DependentVec &Fc, bool sa);

  void solve(int kk, double tstart, double tend,
	     const Omu_VariableVec &x, const Omu_VariableVec &u,
	     Omu_Program *sys, Omu_DependentVec &cF, Omu_StateVec &cx);

  //@}

 protected:

  /**
   * maximum degree of the taylor series expansion of the ode
   */

  int           _max_deg;

  /**
   * time scaling of the taylor series
   */
  double	_tau;

 private:

  bool          _adtaylor;
  bool          _multiple_record;
  bool          _check_autonomous;

  int           _output;

  int           _max_deg0;

  double        _rho;

  VECP          _a;

  Omu_SVec      _cxp;

  MATP          _Sxd;
  MATP          _Sxd0;
  MATP          _Sxx;
  MATP          _Sxx0;
  MATP          _Sxu;
  MATP          _Sxu0;
  MATP          _Sh;

  int           _nnz;
  short         **_nz;         // sparsity pattern

  double        *_aa1;
  double        *_ab1;
  double        **_aa;
  double        **_ab;
  double        **_ai;
  double        ***_aA;
  double        ***_aB;

  void		resize(int );
  void          free_adtaylor();

};

#endif




