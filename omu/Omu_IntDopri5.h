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
 *             2002-04-09 remove static const from class declaration
 *
 */

/*
    Copyright (C) 1997--2014   Eckhard Arnold

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

  const char *name() {return "Dopri5";}

  // interface routine
  void ode_solve(double tstart, VECP y, const VECP u, double tend);

 private:
  static const double c2, c3, c4, c5,
      a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, 
      a61, a62, a63, a64, a65, a71, a73, a74, a75, a76, 
      e1, e3, e4, e5, e6, e7,  
      d1, d3, d4, d5, d6, d7;
  static const int iord;

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
