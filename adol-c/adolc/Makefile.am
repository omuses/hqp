##############################################################################
## Makefile.am -- Process this file with automake to produce Makefile.in
## Revision: $Id: Makefile.am,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
##
## Copyright (C) 2004
##               Technical University Dresden
##               Department of Mathematics
##               Institute of Scientific Computing
##############################################################################

CFLAGS               = @ac_adolc_cflags@
CXXFLAGS             = @ac_adolc_cxxflags@
AM_LDFLAGS           = -no-undefined

INCLUDES             = -I$(top_builddir)

SUBDIRS              = drivers sparse tapedoc

MAINTAINERCLEANFILES = Makefile.in *~

EXTRA_DIST           = uni5_for.c fo_rev.c ho_rev.c \
                       taputil_p.h tayutil_p.h

pkginclude_HEADERS   = adolc.h adalloc.h adouble.h adutils.h adutilsc.h \
                       avector.h common.h convolut.h dvlparms.h fortutils.h \
                       interfaces.h oplate.h taputil.h tayutil.h usrparms.h

lib_LTLIBRARIES      = libadolc.la

libadolc_la_LDFLAGS  = -version-info 0:0:0

libadolc_la_SOURCES  = adalloc.c adouble.cpp avector.cpp \
                       convolut.c fortutils.c \
                       interfaces.cpp interfacesf.c \
                       taputil.c taputilc.cpp \
                       tayutil.c tayutilc.cpp \
                       zos_forward.c fos_forward.c fov_forward.c \
                       hos_forward.c hov_forward.c hov_wk_forward.c \
                       fos_reverse.c fov_reverse.c \
                       hos_reverse.c hos_ov_reverse.c hov_reverse.c \
                       forward_partx.c

libadolc_la_LIBADD   = drivers/libdrivers.la sparse/libsparse.la \
                       tapedoc/libtapedoc.la
