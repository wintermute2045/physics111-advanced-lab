#!/usr/bin/env python2

################################################################################
## This code is for the atomic physics lab from the UC Berkeley Physics Advanced
## Lab class
## Copyright (C) 2013  Rachel Domagalski: idomagalski@berkeley.edu
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import sys as _sys
from libc.stdlib cimport malloc, free
from libc.math   cimport sqrt

# GSL libraries I want to wrap
cdef extern from "gsl/gsl_fit.h":
    int gsl_fit_linear (const double *x,
                        const size_t xstride,
                        const double *y,
                        const size_t ystride,
                        size_t n,
                        double *c0,
                        double *c1,
                        double *cov00,
                        double *cov01,
                        double *cov11,
                        double *sumsq)
    int gsl_fit_wlinear (const double *x,
                         const size_t xstride,
                         const double *w,
                         const size_t wstride,
                         const double *y,
                         const size_t ystride,
                         size_t n,
                         double *c0,
                         double *c1,
                         double *cov00,
                         double *cov01,
                         double *cov11,
                         double *chisq)

def linear_fit(xdata, ydata, weights=None):
    """
    Wrapper function to the GSL fitting functions.
    """
    cdef int xlen = len(xdata)
    cdef int ylen = len(ydata)
    cdef int wlen
    if xlen != ylen:
        print 'xdata and ydata must have equal lengths.'
        _sys.exit(1)

    cdef double *x, *y, *w
    cdef double c0 = 0, c1 = 0, cov00 = 0, cov01 = 0, cov11 = 0, chisq = 0

    # Give data to the C arrays
    x = <double *> malloc(xlen * sizeof(double))
    y = <double *> malloc(ylen * sizeof(double))
    for i in range(xlen):
        x[i] = xdata[i]
        y[i] = ydata[i]

    # Compute the linear least squares fit.
    if weights is None:
        gsl_fit_linear(x, 1, y, 1, xlen,
                       &c0, &c1, &cov00, &cov01, &cov11, &chisq)
    else:
        wlen = len(weights)
        if wlen != xlen:
            print 'Weight length does not match data.'
            _sys.exit(1)

        # Assign the weights
        w = <double *> malloc(wlen * sizeof(double))
        for i in range(wlen):
            w[i] = weights[i]

        # compute the fit.
        gsl_fit_wlinear(x, 1, y, 1, w, 1, xlen,
                        &c0, &c1, &cov00, &cov01, &cov11, &chisq)
        free(w)

    # Free memory and return the results.
    free(x)
    free(y)
    slope = c1
    intercept = c0
    err_slope = sqrt(cov11)
    err_intercept = sqrt(cov00)
    return (slope, err_slope, intercept, err_intercept, chisq)
