/*
 * Prg_DID.h -- 
 *   - double integrator with state constraint
 *
 * rf, 1/29/95
 *
 * rf, 2/10/98
 *   - updated for HQP 1.3
 *
 * rf, 7/11/00
 *   - updated for HQP 1.7
 */

#ifndef Prg_DID_H
#define Prg_DID_H

#include <Hqp_Docp.h>

class Prg_DID: public Hqp_Docp {
 protected:
  /**
   * @name implement interface of Hqp_Docp
   */
  //@{
  virtual void setup_horizon(int &k0, int &kf);

  virtual void setup_vars(int k,
			  VECP x, VECP x_min, VECP x_max, IVECP x_int,
			  VECP u, VECP u_min, VECP u_max, IVECP u_int,
			  VECP c, VECP c_min, VECP c_max);

  virtual void setup_struct(int k, const VECP x, const VECP u,
			    MATP fx, MATP fu, IVECP f_lin,
			    VECP f0x, VECP f0u, int &f0_lin,
			    MATP cx, MATP cu, IVECP c_lin,
			    MATP Lxx, MATP Luu, MATP Lxu);

  virtual void update_vals(int k, const VECP x, const VECP u,
			   VECP f, Real &f0, VECP c);

  virtual void update_stage(int k, const VECP x, const VECP u,
			    VECP f, Real &f0, VECP c,
			    MATP fx, MATP fu, VECP f0x, VECP f0u,
			    MATP cx, MATP cu,
			    const VECP rf, const VECP rc,
			    MATP Lxx, MATP Luu, MATP Lxu);

  //@}

  int  _kmax; 		///< number of sampling intervals
  bool _with_cns; 	///< treat overshoot with additional constraint

 public:

  Prg_DID();

  const char *name() {return "DID";}
};


#endif
