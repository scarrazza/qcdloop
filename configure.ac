#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
#define([gitversion], [esyscmd([sh -c "git rev-parse --short HEAD | tr -d '\n'"])])
#AC_INIT(qcdloop, [gitversion], stefano.carrazza@cern.ch)
AC_INIT(qcdloop, 2.0.2, stefano.carrazza@cern.ch)
AC_CONFIG_MACRO_DIR([m4])
AM_MAINTAINER_MODE
AC_CONFIG_AUX_DIR([autoconfig])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

LT_INIT()

AC_CONFIG_HEADERS([src/qcdloop/config.h])
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([bin/Makefile]) 
AC_CONFIG_FILES([bin/qcdloop-config])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([examples/Makefile])
AM_INIT_AUTOMAKE([subdir-objects])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Check for libraries.
AX_CXX_COMPILE_STDCXX_11([noext])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h sys/time.h quadmath.h])

# Check for compile flag
AX_CHECK_COMPILE_FLAG([-fext-numeric-literals],[CPPFLAGS="$CPPFLAGS -fext-numeric-literals"])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday memset mkdir pow sqrt strstr])

AC_DEFINE(QL_NAMES,[],['Disable ql names conflict'])

AC_OUTPUT
