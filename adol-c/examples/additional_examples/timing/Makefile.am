##############################################################################
## Makefile.am -- Process this file with automake to produce Makefile.in
## Revision: $Id: Makefile.am,v 1.2 2004/10/14 13:29:52 e_arnold Exp $
##
## Copyright (C) 2004
##               Technical University Dresden
##               Department of Mathematics
##               Institute of Scientific Computing
##############################################################################

CFLAGS                  = @ac_adolc_cflags@
CXXFLAGS                = @ac_adolc_cxxflags@
LDADD                   = -ladolc
AM_LDFLAGS              = -L../../../adolc

INCLUDES                = -I$(top_builddir)

EXTRA_DIST              = rotations.h

MAINTAINERCLEANFILES    = Makefile.in *~

if ADDEXA
noinst_PROGRAMS         = detexam speelpenning griewankexam michalewitzexam \
                          rosenbrockexam powexam helmholtzexam shuttlexam \
                          gearexam pargearexam simplevec eutrophexam \
                          robertsonexam ficexam experimental

detexam_SOURCES         = sfunc_determinant.cpp sgenmain.cpp \
                          ../clock/myclock.cpp

speelpenning_SOURCES    = sfunc_speelpenning.cpp sgenmain.cpp \
                          ../clock/myclock.cpp

griewankexam_SOURCES    = sfunc_griewank.cpp sgenmain.cpp \
                          ../clock/myclock.cpp

michalewitzexam_SOURCES = sfunc_michalewitz.cpp sgenmain.cpp \
                          ../clock/myclock.cpp

rosenbrockexam_SOURCES  = sfunc_rosenbrock.cpp sgenmain.cpp \
                          ../clock/myclock.cpp

powexam_SOURCES         = sfunc_power.cpp sgenmain.cpp \
                          ../clock/myclock.cpp

helmholtzexam_SOURCES   = sfunc_helmholtz.cpp sgenmain.cpp \
                          ../clock/myclock.cpp

shuttlexam_SOURCES      = vfunc_shuttle.cpp vgenmain.cpp \
                          ../clock/myclock.cpp

gearexam_SOURCES        = vfunc_gear.cpp vgenmain.cpp \
                          ../clock/myclock.cpp

pargearexam_SOURCES     = vfunc_pargear.cpp rotations.cpp vgenmain.cpp \
                          ../clock/myclock.cpp

simplevec_SOURCES       = vfunc_simplevec.cpp vgenmain.cpp \
                          ../clock/myclock.cpp

eutrophexam_SOURCES     = ../ode/vfunc_eutroph.cpp vgenmain.cpp \
                          ../clock/myclock.cpp

robertsonexam_SOURCES   = ../ode/vfunc_robertson.cpp vgenmain.cpp \
                          ../clock/myclock.cpp

ficexam_SOURCES         = ../ode/vfunc_fic.cpp vgenmain.cpp \
                          ../clock/myclock.cpp

experimental_SOURCES    = sfunc_experimental.cpp sgenmain.cpp \
                          ../clock/myclock.cpp
endif