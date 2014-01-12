/**
 * @file Omu_IntIMP.h --
 *    integrate an ODE over a stage using the implicit midpoint rule
 *
 * E. Arnold   1999-04-12
 *             2000-05-12 _rtol, _atol -> Omu_Integrator 
 *             2000-05-30 step size control
 *             2001-08-16 prg_int_nsteps --> _stepsize
 *             2003-01-02 modified Newton's method from Hairer/Wanner
 *             2003-08-25 Omu_Integrator
 *             2003-09-06 step size control by Richardson extrapolation
 *
 */

/*
    Copyright (C) 1999--2014  Eckhard Arnold

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

#ifndef Omu_IntIMP_H
#define Omu_IntIMP_H

#include "Omu_Integrator.h"

/**
 * Integrator for (stiff) ordinary differential equations (ODEs). The solver
 * is employed internally to numerically integrate continuous-time model
 * equations over sample periods and discrete-time stages.
 * The implicit midpoint rule is implemented using modified Newton's method 
 * described by Hairer/Wanner.
 */  
class Omu_IntIMP: public Omu_Integrator {

 public:

    Omu_IntIMP();     ///< constructor
    ~Omu_IntIMP();    ///< destructor

    const char *name() {return "IMP";}

    // interface routines
    virtual void init(int k,
		      const Omu_StateVec &xc, const Omu_Vec &q,
		      const Omu_DependentVec &Fc, bool sa);

    virtual void solve(int kk, double tstart, double tend,
		       Omu_StateVec &xc, Omu_StateVec &dxc, Omu_Vec &q,
		       Omu_DependentVec &Fc);

 private:

    void resize();
    int  step(double tstep, double dt, VECP y, MATP yq, double tol);
    void fixedstep(double tstep, double dt, VECP y, MATP yq);

    VECP	  _x;
    VECP	  _y;
    VECP	  _k1;
    VECP	  _y0;
    VECP	  _res;
    VECP	  _fh;
    VECP	  _z;
    VECP	  _zp;
    VECP	  _z_old;
    MATP	  _xs;
    MATP	  _yy;
    MATP	  _yyn;
    MATP	  _yq;
    MATP	  _yq1;
    MATP	  _yq2;
    MATP	  _yqq;
    PERM          *_ppivot;
    VECP	  _y1;
    VECP	  _y2;

    int           _maxiters;
    double        _hinit;
    double        _dt;
    double        _eta;
    double        _kappa;
    int           _ixgz;
    int           _max_modnewtonsteps;
    int           _modnewtonsteps;
    bool          _correct_der; 
    int           _nsing;
    int           _max_sing;


    void sys(double t, VECP x, VECP xp);
    void sys_jac(double t, VECP x, VECP xp, MATP fx);
    void sys_jac(double t, VECP x, VECP xp, MATP fx, MATP fu);
    int lufac_jac(double gamma, double delta, MATP fx);

    int              _kk;
    Omu_Vec	     *_q_ptr;
    Omu_SVec	     _dxc;
    Omu_StateVec     *_xc_ptr;
    Omu_DependentVec *_Fc_ptr;


};  

#endif
