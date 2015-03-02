#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([muon-lifetime], [1.0], [idomagalski@berkeley.edu])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/mu2fits.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O

# Checks for libraries.
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([gslcblas], [cblas_dgemm])
AC_CHECK_LIB([gsl], [gsl_blas_dgemm])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h\
                  string.h\
                  unistd.h\
                  libgen.h\
                  ctype.h\
                  math.h\
                  time.h\
                  fitsio.h\
                  dirent.h\
                  sys/types.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([strtol])

AC_OUTPUT(Makefile src/Makefile scripts/Makefile)
