##############################################################################
## Makefile.am -- Process this file with automake to produce Makefile.in
## subdirectory tapedoc
## Revision: $Id: Makefile.am,v 1.2 2004/10/14 13:29:48 e_arnold Exp $
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

noinst_LTLIBRARIES        = libtapedoc.la

libtapedocincludedir      = $(pkgincludedir)/tapedoc

libtapedocinclude_HEADERS = tapedoc.h

libtapedoc_la_SOURCES     = tapedoc.c
