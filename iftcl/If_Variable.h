/**
 *  @file If_Variable.h
 *     base classes for interface variables
 *
 *  rf, 6/22/94
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

#ifndef If_Variable_H
#define If_Variable_H

#include "If_Element.h"

/** Construct a read callback */
#define IF_GET_CB(vartype, classtype, name) \
  new If_GetCb<vartype, classtype>(&classtype::name, this)

/** Construct a write callback */
#define IF_SET_CB(vartype, classtype, name) \
  new If_SetCb<vartype, classtype>(&classtype::name, this)

/** Interface for a read callback. */
template <class VarType>
class If_GetIf {
 public:
  virtual ~If_GetIf() {} 	   ///< destructor
  virtual VarType get() const = 0; ///< get method
};


/** Read callback for specific class type. */
template <class VarType, class ClassType>
class If_GetCb: public If_GetIf<VarType> {

 protected:
  VarType 	(ClassType::*_get)() const; ///< pointer to callback method
  ClassType	*_object; 		    ///< pointer to object

 public:
  /** constructor */
  If_GetCb(VarType (ClassType::*n_get)() const, ClassType *n_object) {
    assert(n_object && n_get);
    _get = n_get;
    _object = n_object;
  }
  /** get method */
  VarType get() const {
    return (_object->*_get)();
  }
};

/** Read callback directly accessing through a pointer. */
template <class VarType>
class If_GetPt: public If_GetIf<VarType> {
 protected:
  VarType *_varPtr; 	///< pointer to accessed variable

 public:
  /** constructor */
  If_GetPt(VarType *varPtr): _varPtr(varPtr) {}
  /** get method */
  VarType get() const {return *_varPtr;}
};

/** Read callback directly accessing through a pointer of different type. */
template <class VarType, class InternalType>
class If_GetPt2: public If_GetIf<VarType> {
 protected:
  InternalType *_varPtr; 	///< pointer to accessed variable

 public:
  /** constructor */
  If_GetPt2(InternalType *varPtr): _varPtr(varPtr) {}
  /** get method assuming a valid conversion from InternalType to VarType */
  VarType get() const {return *_varPtr;}
};


/** Interface for a write callback. */
template <class VarType>
class If_SetIf {
 public:
  virtual ~If_SetIf() {} 		///< destructor
  virtual void set(VarType newVal) = 0; ///< set method
};

/** Write callback for specific class type. */
template <class VarType, class ClassType>
class If_SetCb: public If_SetIf<VarType> {
 protected:
  void		(ClassType::*_set)(VarType); 	///< pointer to callback
  ClassType	*_object; 			///< pointer to object

 public:
  /** constructor */
  If_SetCb(void (ClassType::*n_set)(VarType), ClassType *n_object) {
    assert(n_object && n_set);
    _set = n_set;
    _object = n_object;
  }
  /** set method */
  void set(VarType newVal) {
    (_object->*_set)(newVal);
  }
};

/** Write callback directly accessing through a pointer. */
template <class VarType>
class If_SetPt: public If_SetIf<VarType> {
 protected:
  VarType *_varPtr; 	///< pointer to accessed variable

 public:
  /** constructor */
  If_SetPt(VarType *varPtr): _varPtr(varPtr) {}
  /** set method */
  void set(VarType value) {*_varPtr = value;}
};


/**
 *  Interface variable of abstract type VarType managing callbacks for
 *  reading and writing. A derived class should pass dynamically created 
 *  callback objects to the constructor. In order to invoke the callbacks,
 *  the methods get and set should be called from within the
 *  implementations for getTclObj and setTclObj, respectively
 *  (as defined by If_VariableCmd). 
 *  The destructor of this base class will delete the callback objects.
 */
template <class VarType>
class If_Variable: public If_Element {

 protected:
  If_GetIf<VarType>	*_getCb;  ///< callback for reading
  If_SetIf<VarType>	*_setCb;  ///< callback for writing

  /** Constructor taking callback methods as arguments. */
  If_Variable(const char *ifName, If_GetIf<VarType> *getCb,
	      If_SetIf<VarType> *setCb = NULL)
    :If_Element(ifName) {
    _getCb = getCb;
    _setCb = setCb;
  }

  /** Process command invokation */
  int invoke(Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[]) {
    switch (objc) {
    case 1:
      if (_getCb)
	return getTclObj(interp);
      else {
	Tcl_AppendResult(interp, "read permission denied for ",
			 ifName(), NULL);
	return TCL_ERROR;
      }

    case 2:
      if (_setCb)
	return setTclObj(interp, objv[1]);
      else {
	Tcl_AppendResult(interp, "write permission denied for ",
			 ifName(), NULL);
	return TCL_ERROR;
      }

    default:
      Tcl_AppendResult(interp, "wrong # args, should be: ",
		       ifName(), " [new value]", NULL);
      return TCL_ERROR;
    }
    return TCL_ERROR;
  }

  /** Get the value of the interface variable
      and create a Tcl object containing it. */
  virtual int	getTclObj(Tcl_Interp *) = 0;

  /** Convert Tcl object to internal type and
      set the value of the interface variable. */
  virtual int	setTclObj(Tcl_Interp *, Tcl_Obj *CONST objPtr) = 0;

 public:
  /** Destructor */
  virtual ~If_Variable() {
    delete _setCb;
    delete _getCb;
  }

  /** Get value of interface variable */
  virtual VarType get() const {
    assert(_getCb);
    return _getCb->get();
  }

  /** Set value of interface variable */
  virtual void set(VarType value) {
    assert(_setCb);
    _setCb->set(value);
  }
};


#endif
