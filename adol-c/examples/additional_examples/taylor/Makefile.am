##############################################################################
## Makefile.am -- Process this file with automake to produce Makefile.in
## Revision: $Id: Makefile.am,v 1.2 2004/10/14 13:29:52 e_arnold Exp $
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
noinst_PROGRAMS      = taylorexam inversexam accessexam trigger coordinates

taylorexam_SOURCES   = taylorexam.cpp

inversexam_SOURCES   = inversexam.cpp

accessexam_SOURCES   = accessexam.cpp

trigger_SOURCES      = trigger.cpp

coordinates_SOURCES  = coordinates.cpp
endif