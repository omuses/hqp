/*
 * Omu_Program.h --
 *   -- declare an interface for a multistage optimization problem
 *      described by differential and algebraic equations
 *
 * rf, 16/1/97
 */

/*
    Copyright (C) 1997--2001  Ruediger Franke

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

#ifndef Omu_Program_H
#define Omu_Program_H

#include <adouble.h>
#include <Meschach.h>
#include <Omu_Vector.h>

#include <If_Class.h>

IF_BASE_DECLARE(Omu_Program);

//--------------------------------------------------------------------------
class Omu_Program {

 public:

  Omu_Program();
  virtual ~Omu_Program();

  //
  // method interface
  //

  /**
   * Setup stages and time sampling.
   */
  virtual void setup_stages();

  /**
   * Routine called by default implementation of
   * public setup_stages()
   */
  virtual void setup_stages(IVECP ks, VECP ts);

  virtual void setup(int k,
		     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c) = 0;

  virtual void init_simulation(int k,
			       Omu_Vector &x, Omu_Vector &u);

  virtual void update(int kk, 
		      const adoublev &x, const adoublev &u,
		      adoublev &f, adouble &f0, adoublev &c) = 0;

  /**
   * High-level consistic (consistent initial conditions) routine
   * for the initialization of continuous-time states from
   * optimization variables x and u.
   * The routine is called via the default low-level version
   * before the first call to continuous in each sample period.
   * The default implementation copies x to xt.
   */
  virtual void consistic(int kk, double t,
			 const adoublev &x, const adoublev &u,
			 adoublev &xt);

  /**
   * Low-level consistic (consistent initial conditions) routine
   * for the user specification of derivative information.
   * The default implementation calls the high-level version
   * and uses ADOL-C for the determination of xtx and xtu.
   * If a calling routine passes no xtx and xtu (or null pointers),
   * then no Jacobians are calculated at all.
   * If a calling routine passes xtx and xtu, then the matrices are
   * initialized with zeros.
   */
  virtual void consistic(int kk, double t,
			 const Omu_Vector &x, const Omu_Vector &u,
			 VECP xt, MATP xtx = MNULL, MATP xtu = MNULL);

  /**
   * High-level continuous routine for specifying differential equations
   * of the form F(x,u,\dot{x})=0.
   * The default implementation is empty to indicate that no differential
   * equations are used.
   */
  virtual void continuous(int kk, double t,
			  const adoublev &x, const adoublev &u,
			  const adoublev &xp, adoublev &F);

  /**
   * Low-level continuous routine,
   * which may be implemented in addition to high-level continuous
   * by a derived class to provide Jacobians for differential equations.
   * The default implementation calls the high-level continuous()
   * and employs ADOL-C for automatic Jacobian calculation.
   * If a calling routine passes no Fxp (or a null pointer),
   * then only Fx and Fu are calculated.
   * If a calling routine passes no Fx, Fu, and Fxp (or null pointers),
   * then no Jacobians are calculated at all.
   * If a calling routine passes Fx, Fu, or Fxp, then the matrices are
   * initialized with zeros.
   */
  virtual void continuous(int kk, double t,
			  const VECP x, const VECP u,
			  const VECP xp, VECP F,
			  MATP Fx = MNULL, MATP Fu = MNULL,
			  MATP Fxp = MNULL);

  /**
   * This routine is intended for a calling module,
   * which might want to avoid calling the low-level continuous
   * if it is not overridden (e.g. as it has then a faster sensitivity
   * calculation).
   * @return true if low-level continuous is overloaded
   * (currently low-level continuous must have been called once before 
   *  false may be returned)
   */
  bool has_low_level_continuous();

  //
  // member access routines
  //

  int		K() {return _K;}
  void		set_K(int K) {_K = K;}
  int		KK() {return _KK;}
  void		set_KK(int KK) {_KK = KK;}
  const IVECP 	ks() {return _ks;}
  const VECP 	ts() {return _ts;}
  int	 	ks(int k) {return _ks[k];}
  Real 		ts(int kk) {return _ts[kk];}

  //
  // textual name of a problem definition
  //

  virtual char *name() = 0;

 protected:

  //
  // provide a container to keep track of interface elements
  //

  If_List	_ifList;

  //
  // predefined basic variables
  //

  int 	_K;	// number of stages
  int 	_KK;	// number of sample periods over all stages

  //		     
  // service routines
  //

  void 	stages_alloc(IVECP ks, VECP ts, int K, int sps,
		     double t0 = 0.0, double tf = 1.0);

private:
  IVECP _ks;	// starting indices of stages
  VECP 	_ts;	// time steps

  void		ad_alloc(int ndep, int nindep);
  void		ad_realloc(int ndep, int nindep);
  void		ad_free();

  int		_max_ndep;
  int		_max_nindep;
  MATP 		_U;
  MATP		_Z;

  bool 		_has_low_level_continuous;
};  

#endif
