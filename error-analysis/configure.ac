#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([eax-homework], [1.0], [idomagalski@berkeley.edu])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/gauss_stats.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h fitsio.h math.h time.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_OUTPUT(Makefile src/Makefile scripts/Makefile)
