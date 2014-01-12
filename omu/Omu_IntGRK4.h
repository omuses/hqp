/*
 * @file Omu_IntGRK4.h --
 *    integrate a (stiff) ODE over a stage using a linear-implicit RK method
 *
 * Reference: FORTRAN code ros4.f:
 *
 *    NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC)
 *    SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS  MY'=F(X,Y).
 *    THIS IS AN EMBEDDED ROSENBROCK METHOD OF ORDER (3)4  
 *    (WITH STEP SIZE CONTROL).
 *    C.F. SECTION IV.7
 *
 *    AUTHORS: E. HAIRER AND G. WANNER
 *             UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
 *             CH-1211 GENEVE 24, SWITZERLAND 
 *             E-MAIL:  HAIRER@CGEUGE51.BITNET,  WANNER@CGEUGE51.BITNET
 *    
 *    THIS CODE IS PART OF THE BOOK:
 *        E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
 *        EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
 *        SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
 *        SPRINGER-VERLAG (1990)               
 *     
 *    VERSION OF NOVEMBER 17, 1992
 *
 * E. Arnold   2000-05-25 C++ version
 *             2003-01-03
 *             2003-09-07 Omu_Integrator
 *
 */

/*
    Copyright (C) 2000-2014   Eckhard Arnold

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

#ifndef Omu_IntGRK4_H
#define Omu_IntGRK4_H

#include "Omu_Integrator.h"

/**
 * Integrator for (stiff) ordinary differential equations (ODEs). The solver
 * is employed internally to numerically integrate continuous-time model
 * equations over sample periods and discrete-time stages.
 * The linear-implicit RK method is derived from FORTRAN code ros4.f
 * described by Hairer/Wanner.
 */  
class Omu_IntGRK4: public Omu_Integrator {

 public:

    Omu_IntGRK4();     ///< constructor
    ~Omu_IntGRK4();    ///< destructor

    const char *name() {return "GRK4";}

    // interface routines
    virtual void init(int k,
		      const Omu_StateVec &xc, const Omu_Vec &q,
		      const Omu_DependentVec &Fc, bool sa);

    virtual void solve(int kk, double tstart, double tend,
		       Omu_StateVec &xc, Omu_StateVec &dxc, Omu_Vec &q,
		       Omu_DependentVec &Fc);

 private:

    double A21, A31, A32, C21, C31, C32, C41, C42, C43, B1, B2, B3, B4;
    double E1, E2, E3, E4, GAMMA, C2, C3, D1, D2, D3, D4;
    VECP _y, _ynew, _fx, _dy, k1, k1_ori, k2, k3, k4, y1;
    MATP _yy;
    MATP _yyn;
    MATP _yq;
    MATP _yq1;
    MATP _yqn;
    MATP _ys;
    PERM  *_ppivot;

    double _uround, _safe, _beta, _fac1, _fac2, _facc1, _facc2;
    long int _nmax, _nstiff, _nfcn, _nstep, _naccpt, _nrejct, _nsing;
    int _npar, _coeffs, _max_sing;
    bool _sensrk4;
    
    double _xold, _posneg, _h, _x, _xend, _hmax;
    double _hinit, _hmaxinit;

    int              _kk;
    Omu_Vec	     *_q_ptr;
    Omu_SVec	     _dxc;
    Omu_StateVec     *_xc_ptr;
    Omu_DependentVec *_Fc_ptr;

 private:
    void resize();
    void simulation();
    int coeffs();
    void sys(double t, const VECP x, VECP xp);
    void sys_jac(double t, const VECP x, VECP xp, MATP fx);
    void sys_jac(double t, const VECP x, VECP xp, MATP fx, MATP fu);
    int  lufac_jac(double gamma, double delta, MATP fx);
    void lusolve_jac(MATP A, VECP b, VECP x);
    void update_sens(const MATP fx, const VECP s, double fac, 
		     const VECP ds, const MATP fu, VECP sp);

};  

#endif
