/*
 * cg_sfun.h: code generation include file for HQP.
 * Currently macros for checking optional S-function methods are defined.
 * This allows to write the same code for an inlined S-function and for 
 * an external MEX S-function
 *
 * (Simulink is a registered trademark of The MathWorks, Inc.)
 *
 * rf, 05/05/2001
 */

/*
    Copyright (C) 1994--2001  Ruediger Franke

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

#if defined(MDL_START)
#define ssGetmdlStart(S) mdlStart
#else
#define ssGetmdlStart(S) NULL
#endif

#if defined(MDL_INITIALIZE_CONDITIONS)
#define ssGetmdlInitializeConditions(S) mdlInitializeConditions
#else
#define ssGetmdlInitializeConditions(S) NULL
#endif

#if defined(MDL_UPDATE)
#define ssGetmdlUpdate(S) mdlUpdate
#else
#define ssGetmdlUpdate(S) NULL
#endif

#if defined(MDL_DERIVATIVES)
#define ssGetmdlDerivatives(S) mdlDerivatives
#else
#define ssGetmdlDerivatives(S) NULL
#endif

