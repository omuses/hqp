/*
 * Prg_Bio.h -- 
 *   - Hqp/Omuses example: bio-technological process
 *
 * Adapted from:
 *    @PhdThesis{Pfaff:91,
 *       author = {Pfaff, M.},
 *       title =  {Entwurf optimaler {S}teuerstrategien ausgew{\"a}hlter 
 *                 biotechnologischer {P}rozesse anhand aggregierter 
 *                 kinetischer {M}odelle},
 *       school = {Technische Hochschule Ilmenau},
 *       year =   1991}
 *
 * E. Arnold 10/20/96
 *           2003-02-16 odc example
 */

/*
    Copyright (C) 1996--2003  Eckhard Arnold

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

#ifndef Prg_Bio_H
#define Prg_Bio_H

#include <Omu_Program.h>

//--------------------------------------------------------------------------
class Prg_Bio: public Omu_Program {
    
 protected:

    double _fscale;
    VECP   _y;
    bool   _controller;
    double _uinit;
  
    // model parameters
    double pimax, ks, kis, kip, kd, yps, kappa, cdos, kp, kap, kos;
    double cs0, p0, x0, v0, Fs, Fsmin, Fsmax, tf;

    void setup_stages(IVECP ks, VECP ts);

    void setup(int k,
	       Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

    void init_simulation(int k,
			 Omu_Vector &x, Omu_Vector &u);

    void update(int kk,
		const adoublev &x, const adoublev &u,
		adoublev &f, adouble &f0, adoublev &c);

    void continuous(int kk, double t,
		    const adoublev &x, const adoublev &u, const adoublev &xp,
		    adoublev &F);

    // nonlinear state controller for generation of initial control trajectory
    double Controller(const VECP, const double, const double);

 public:

    Prg_Bio();
    ~Prg_Bio();

    char *name() {return "Bio";}
};  

#endif
