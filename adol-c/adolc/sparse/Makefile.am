##############################################################################
## Makefile.am -- Process this file with automake to produce Makefile.in
## subdirectory sparse
## Revision: $Id: Makefile.am,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
##
## Copyright (C) 2004
##               Technical University Dresden
##               Department of Mathematics
##               Institute of Scientific Computing
##############################################################################

CFLAGS                   = @ac_adolc_cflags@
CXXFLAGS                 = @ac_adolc_cxxflags@
AM_LDFLAGS               = -no-undefined

INCLUDES                 = -I$(top_builddir)

MAINTAINERCLEANFILES     = Makefile.in *~

EXTRA_DIST               = int_for.c int_rev.c

noinst_LTLIBRARIES       = libsparse.la

libsparseincludedir      = $(pkgincludedir)/sparse

libsparseinclude_HEADERS = jacutils.h sparse.h

libsparse_la_SOURCES     = int_forward_s.c int_forward_t.c \
                           int_reverse_s.c int_reverse_t.c \
                           jacutils.c sparse.cpp sparsedrivers.c
