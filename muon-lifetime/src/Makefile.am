# Compiler options
AM_CFLAGS = -Wall -I ../incl -I../../lib/incl
AM_LDFLAGS = -lm -lcfitsio -lgsl -lgslcblas
AUTOMAKE_OPTIONS = subdir-objects

# Binaries
bin_PROGRAMS = mu2fits catmufits timeres eff_dig dcalib plinear getlifetime
mu2fits_SOURCES 	= mu2fits.c \
					  ../../lib/src/misc.c
catmufits_SOURCES 	= catmufits.c \
					  muo_common.c \
					  ../../lib/src/misc.c
timeres_SOURCES 	= timeres.c \
					  muo_common.c \
					  ../../lib/src/fit_gaus.c \
					  ../../lib/src/misc.c
eff_dig_SOURCES 	= eff_dig.c \
					  muo_common.c \
					  ../../lib/src/map.c \
					  ../../lib/src/fit_gaus.c \
					  ../../lib/src/misc.c
dcalib_SOURCES		= dcalib.c \
					  muo_common.c \
					  ../../lib/src/misc.c
plinear_SOURCES 	= plinear.c \
					  muo_common.c \
					  ../../lib/src/misc.c
getlifetime_SOURCES = getlifetime.c \
					  muo_common.c \
					  ../../lib/src/fit_exp.c \
					  ../../lib/src/map.c \
					  ../../lib/src/filter.c \
					  ../../lib/src/misc.c
