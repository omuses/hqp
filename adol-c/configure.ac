##############################################################################
# configure.ac -- Process this file with autoconf to produce configure
# Revision: $Id: configure.ac,v 1.1 2004/10/13 15:05:26 e_arnold Exp $
#
# Copyright (C) 2003
#               Technical University Dresden
#               Department of Mathematics
#               Institute of Scientific Computing
##############################################################################
adolc_version    = 1
adolc_subversion = 9
adolc_patchlevel = 0

AC_PREREQ(2.57)
AC_INIT(adolc,
        1.9.0,
        [adol-c@math.tu-dresden.de])
AC_REVISION($Revision: 1.1 $)
AC_CONFIG_SRCDIR([adolc/adolc.h])
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_HEADERS([adolc/config.h])
AC_PREFIX_DEFAULT(${HOME}/adolc_base)
AM_INIT_AUTOMAKE

# checks for programs
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL

# Checks for libraries
AC_CHECK_LIB([m], [pow])

# checks for header files
AC_HEADER_STDC
AC_HEADER_TIME
AC_CHECK_HEADERS([malloc.h stdlib.h string.h])

# checks for types
AC_C_CONST
AC_C_INLINE
AC_STRUCT_TM

# Checks for library functions
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt])

# substitutions
AC_SUBST(ac_aux_dir)

# ADOL-C configuration
AC_DEFINE(ADOLC_VERSION,1,[ADOL-C Version])
AC_DEFINE(ADOLC_SUBVERSION,9,[ADOL-C Subversion])
AC_DEFINE(ADOLC_PATCHLEVEL,0,[ADOL-C Patchlevel])

AC_MSG_CHECKING(whether to use ADOL-C debug mode)
AC_ARG_ENABLE(debug,
              [
AC_HELP_STRING([--enable-debug],[enable ADOL-C debug mode [default=no]])],
  [AC_DEFINE(ADOLC_DEBUG,1,[ADOL-C Debug mode])
   AC_MSG_RESULT(yes)
   adolc_debug="yes"],
  [AC_MSG_RESULT(no)
   adolc_debug="no"])

# ADOL-C Examples
AC_MSG_CHECKING(whether to build documented examples)
AC_ARG_WITH(docexa,
            [
AC_HELP_STRING([--with-docexa],[build documented examples [default=no]])],
            [docexa="yes"
             AC_MSG_RESULT(yes)],
            [docexa="no"
             AC_MSG_RESULT(no)])
AM_CONDITIONAL(DOCEXA,[test x${docexa} = xyes])

AC_MSG_CHECKING(whether to build additional examples)
AC_ARG_WITH(addexa,
            [AC_HELP_STRING([--with-addexa],
                            [build additional examples [default=no]])],
            [addexa="yes"
             AC_MSG_RESULT(yes)],
            [addexa="no"
             AC_MSG_RESULT(no)])
AM_CONDITIONAL(ADDEXA,[test x${addexa} = xyes])

# ADOL-C Compilation flags
AC_MSG_CHECKING(which CFLAGS to use)
AC_ARG_WITH(cflags,
            [
Compilation Flags:]
AC_HELP_STRING([--with-cflags=FLAGS],
               [use CFLAGS=FLAGS (default: -O3 -Wall -ansi)]),
            [ac_adolc_cflags=$withval
             AC_SUBST(ac_adolc_cflags)
             AC_MSG_RESULT($ac_adolc_cflags)],
            [ac_adolc_cflags="-O3 -Wall -ansi"
             AC_SUBST(ac_adolc_cflags)
             AC_MSG_RESULT($ac_adolc_cflags)])

AC_MSG_CHECKING(which CXXFLAGS to use)
AC_ARG_WITH(cxxflags,
            [AC_HELP_STRING([--with-cxxflags=FLAGS],
                            [use CXXFLAGS=FLAGS (default: -O3 -Wall)])],
            [ac_adolc_cxxflags=$withval
             AC_SUBST(ac_adolc_cxxflags)
             AC_MSG_RESULT($ac_adolc_cxxflags)],
            [ac_adolc_cxxflags="-O3 -Wall"
             AC_SUBST(ac_adolc_cxxflags)
             AC_MSG_RESULT($ac_adolc_cxxflags)])

# output
AC_CONFIG_FILES(Makefile
                adolc/Makefile
                adolc/drivers/Makefile 
                adolc/sparse/Makefile
                adolc/tapedoc/Makefile
                doc/Makefile
                examples/Makefile
                examples/additional_examples/Makefile
                examples/additional_examples/clock/Makefile
                examples/additional_examples/gauss/Makefile
                examples/additional_examples/hessmat/Makefile
                examples/additional_examples/lufact/Makefile
                examples/additional_examples/scal/Makefile
                examples/additional_examples/speelpenning/Makefile
                examples/additional_examples/taylor/Makefile
                examples/additional_examples/detexam/Makefile
                examples/additional_examples/helm/Makefile
                examples/additional_examples/lighthouse/Makefile
                examples/additional_examples/ode/Makefile
                examples/additional_examples/sparse/Makefile
                examples/additional_examples/tapesave/Makefile
                examples/additional_examples/pow/Makefile
                examples/additional_examples/timing/Makefile)

AC_OUTPUT

# echo configuration
echo \
"-----------------------------------------------------------------------------
Configuration:

  Source code location:       ${srcdir}
  C compiler:                 ${CC}
  C++ compiler:               ${CXX}
  Install path:               ${prefix}
  
  CFLAGS:                     ${ac_adolc_cflags}
  CXXFLAGS:                   ${ac_adolc_cxxflags}

  Use ADOL-C debug mode:      ${adolc_debug}

  Build documented examples:  ${docexa}
  Build additional examples:  ${addexa}

  See adolc/config.h for further configuration information.
-----------------------------------------------------------------------------"

# check for make install problems resulting from missing write permission
if test x${prefix} = x${HOME}/adolc_base ; then
echo \
"
*****************************************************************************

  Installation (make install) will copy all files to:
     ${HOME}/adolc_base

  To successfully compile and run programs using the shared adolc library do
  the following things:
     compiling:
        add \"-I\${HOME}/adolc_base/include\" to your compiler call
     linking:
        add \"-L\${HOME}/adolc_base/lib -ladolc\" to your linker call
     executing (do one of the following things):
        * add \${HOME}/adolc_base to your LD_LIBRARY_PATH variable
        * ask your system administrator for adding \${HOME}/adolc_base to the
          global file containing library search paths (/etc/ld.so.conf)

     (or use the static library by replacing
        \"-L\${HOME}/adolc_base/lib -ladolc\" with
        \"\${HOME}/adolc_base/lib/libadolc.a\" for linking)

  See README for instructions on how to change to a nonlocal installation!

*****************************************************************************
"
fi
