##############################################################################
## Makefile.am -- Process this file with automake to produce Makefile.in
## subdirectory drivers
## Revision: $Id: Makefile.am,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
##
## Copyright (C) 2004
##               Technical University Dresden
##               Department of Mathematics
##               Institute of Scientific Computing
##############################################################################

CFLAGS                    = @ac_adolc_cflags@
CXXFLAGS                  = @ac_adolc_cxxflags@
AM_LDFLAGS                = -no-undefined

MAINTAINERCLEANFILES      = Makefile.in *~

INCLUDES                  = -I$(top_builddir)

noinst_LTLIBRARIES        = libdrivers.la

libdriversincludedir      = $(pkgincludedir)/drivers

libdriversinclude_HEADERS = drivers.h odedrivers.h taylor.h

libdrivers_la_SOURCES     = drivers.c driversf.c odedrivers.c odedriversf.c \
                            taylor.c
