/* 
 * iftest.C --
 *    test application for interface elements
 *
 *  rf, 2/6/97
 */

/*
    Copyright (C) 1994--2002  Ruediger Franke

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

#include <tcl.h>
#include <string.h>
#include <stdio.h>

#include "If_Int.h"
#include "If_Bool.h"
#include "If_Real.h"
#include "If_RealVec.h"
#include "If_RealMat.h"
#include "If_IntVec.h"
#include "If_String.h"
#include "If_StdString.h"
#include "If_Procedure.h"
#include "If_Method.h"

/*
 * the main function
 */

int
main(int argc, char **argv)
{
  Tcl_Main(argc, argv, Tcl_AppInit);
  return 0;
}

/*
 * Tcl_AppInit function
 */

int if_int = 1234;
double if_real = 5.6;
bool if_bool = true;
VECP if_realVec = VNULL;
MATP if_realMat = MNULL;
IVECP if_intVec = IVNULL;
const char *if_string = "initial if_string value";

void if_procedure()
{
  printf("if_procedure called\n");
}

class A {
 public:
  int   _i;
  IVECP _intVecp;
  char *_string;
  std::string _stdString;

  A(): _i(11),
       _stdString("initial stdString value") {
    _intVecp = iv_get(2);
    _string = strdup("initial string value");
  }

  ~A() {
    free(_string);
    iv_free(_intVecp);
  }

  int i() const {
    return _i;
  }

  const IVECP intVecp() const {
    return _intVecp;
  }
  void set_intVecp(const IVECP value) {
    iv_copy(value, _intVecp);
  }

  const char *string() const {return _string;}
  void set_string(const char *value) {
    free(_string);
    _string = strdup(value);
  }

  const std::string &stdString() const {return _stdString;}
  void set_stdString(const std::string &value) {_stdString = value;}

  void if_method()
  {
    m_error(E_INPUT, "A::if_method");
  }
};

A a;

int
Tcl_AppInit(Tcl_Interp *interp)
{
  if (Tcl_InitStubs(interp, "8.1", 0) == NULL) {
    return TCL_ERROR;
  }

  if (Tcl_Init(interp) == TCL_ERROR) {
    return TCL_ERROR;
  }

  theInterp = interp;

  new If_Int("if_int", &if_int);
  new If_Int("if_intCb",
	     new If_GetCb<int, A>(&A::i, &a));

  new If_Bool("if_bool", &if_bool);
  new If_Real("if_real", &if_real);

  if_realVec = v_get(5);
  new If_RealVec("if_realVec", &if_realVec);

  if_realMat = m_get(2,3);
  new If_RealMat("if_realMat", &if_realMat);

  if_intVec = iv_get(3);
  new If_IntVec("if_intVec", &if_intVec);

  new If_IntVec("if_intVecp",
                new If_GetCb<const IVECP, A>(&A::intVecp, &a),
                new If_SetCb<const IVECP, A>(&A::set_intVecp, &a));

  new If_String("if_string", &if_string);
		
  new If_String("if_stringCb",
		new If_GetCb<const char *, A>(&A::string, &a),
		new If_SetCb<const char *, A>(&A::set_string, &a));

  new If_StdString("if_stdString",
		   new If_GetCb<const std::string&, A>(&A::stdString, &a),
		   new If_SetCb<const std::string&, A>(&A::set_stdString, &a));

  new If_Procedure("if_procedure", &if_procedure);

  new If_Method<A>("if_method", &A::if_method, &a);

  Tcl_SetVar(interp, "tcl_rcFileName", "~/.tclshrc", TCL_GLOBAL_ONLY);
  return TCL_OK;
}
