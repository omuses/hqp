##############################################################################
## Makefile.am -- Process this file with automake to produce Makefile.in
## Revision: $Id: Makefile.am,v 1.1 2004/10/13 15:05:42 e_arnold Exp $
##
## Copyright (C) 2004
##               Technical University Dresden
##               Department of Mathematics
##               Institute of Scientific Computing
##############################################################################

CFLAGS               = @ac_adolc_cflags@
CXXFLAGS             = @ac_adolc_cxxflags@
LDADD                = -ladolc
AM_LDFLAGS           = -L../../../adolc

INCLUDES             = -I$(top_builddir)

MAINTAINERCLEANFILES = Makefile.in *~

if ADDEXA
noinst_PROGRAMS      = odexam od2exam od3exam

odexam_SOURCES       = vfunc_robertson.cpp vgenodemain.cpp ../clock/myclock.cpp

od2exam_SOURCES      = vfunc_eutroph.cpp vgenodemain.cpp ../clock/myclock.cpp

od3exam_SOURCES      = vfunc_fic.cpp vgenodemain.cpp ../clock/myclock.cpp
endif