/*
 * Omu_IntDopri5.h --
 *   -- integrate ODE over a stage using Dormand-Prince
 *
 * Reference: FORTRAN code dopri5.f:
 *
 *    NUMERICAL SOLUTION OF A SYSTEM OF FIRST 0RDER
 *    ORDINARY DIFFERENTIAL EQUATIONS  Y'=F(X,Y).
 *    THIS IS AN EXPLICIT RUNGE-KUTTA METHOD OF ORDER (4)5  
 *    DUE TO DORMAND & PRINCE (WITH STEPSIZE CONTROL AND
 *    DENSE OUTPUT).
 *
 *    AUTHORS: E. HAIRER AND G. WANNER
 *             UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
 *             CH-1211 GENEVE 24, SWITZERLAND 
 *             E-MAIL:  HAIRER@UNI2A.UNIGE.CH,  WANNER@UNI2A.UNIGE.CH
 *    
 *    THIS CODE IS DESCRIBED IN:
 *        E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
 *        DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
 *        SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
 *        SPRINGER-VERLAG (1993)               
 *     
 *    VERSION OF MARCH 24, 1993
 *
 *
 * E. Arnold   1999-01-23 C++ version
 *                        continuous output not used!
 *             1999-04-23 egcs, g77, libg2c *
 *             2000-03-29 Real -> double
 *                        error estimation using first _n components only
 *                        _rtol, _atol -> Omu_Integrator
 *
 */

/*
    Copyright (C) 1997--2000   Eckhard Arnold

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

#ifndef Omu_IntDopri5_H
#define Omu_IntDopri5_H

#include "Omu_IntODE.h"

//--------------------------------------------------------------------------
class Omu_IntDopri5: public Omu_IntODE {

 public:

  Omu_IntDopri5();
  ~Omu_IntDopri5();

  char *name() {return "Dopri5";}

  // interface routine
  void ode_solve(double tstart, VECP y, const VECP u, double tend);

 private:

   /*  RUNGE-KUTTA coefficients of DORMAND and PRINCE (1980) */
  static const double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0/9.0,
    a21 = 0.2, a31 = 3.0/40.0, a32 = 9.0/40.0,
    a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0,
    a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0,
    a54 = -212.0/729.0, a61 = 9017.0/3168.0, a62 = -355.0/33.0,
    a63 = 46732.0/5247.0, a64 = 49.0/176.0, a65 = -5103.0/18656.0,
    a71 = 35.0/384.0, a73 = 500.0/1113.0, a74 = 125.0/192.0,
    a75 = -2187.0/6784.0, a76 = 11.0/84.0, e1 = 71.0/57600.0,
    e3 = -71.0/16695.0, e4 = 71.0/1920.0, e5 = -17253.0/339200.0,
    e6 = 22.0/525.0, e7 = -1.0/40.0,  
  /* dense output of SHAMPINE (1986) */
    d1 = -12715105075.0/11282082432.0,
    d3 = 87487479700.0/32700410799.0,
    d4 = -10690763975.0/1880347072.0,
    d5 = 701980252875.0/199316789632.0,
    d6 = -1453857185.0/822651844.0,
    d7 = 69997945.0/29380423.0;
  static const int iord = 5;

  VECP _y, k1, k2, k3, k4, k5, k6, y1, ysti;
  VECP cont, cont1, cont2, cont3, cont4, cont5;
  VECP _u;
  PERM *_icont;

  double _uround, _safe, _beta, _fac1, _fac2, _facc1, _facc2;
  long int _nmax, _nstiff, _nfcn, _nstep, _naccpt, _nrejct;
 
  double _xold, _posneg, _h, _x, _xend, _hmax;
  /*  double _rtol, _atol; */
  double _hinit, _hmaxinit;

 public:
  VECP contd5(const double x, VECP cont);
  void f(const double x, const VECP y, VECP dy) { syseq(x, y, _u, dy); };
  int fout(const int nr, const double x, VECP y) { return 0;}
  void set_cont(PERM *icont);
  void set_initstep(double h, double hmax = 0.0); 
  double hinit();
  int simulation();
};  

#endif
