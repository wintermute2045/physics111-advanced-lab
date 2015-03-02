#!/usr/bin/env python2.7

################################################################################
## This code is part of the analysis for the muon lifetime lab.
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

import sys
import pyfits
from math import sqrt

# Run everything
if __name__ == '__main__':
    print 'Combining results from the following files:'
    for i in sys.argv[1:]:
        print i
    print

    # Open the inputted FITS files, get lifetimes and errors.
    fitslist = map(pyfits.open, sys.argv[1:])
    lbda     = map(lambda x: x[1].header['EXPLBDA'], fitslist)
    lbdaerr  = map(lambda x: x[1].header['ERRLBDA'], fitslist)
    thresh   = map(lambda x: x[1].header['CHTH0'],   fitslist)
    tau      = map(lambda x: 1 / x, lbda)
    errtau   = map(lambda x: x[0]**2 * x[1], zip(tau, lbdaerr))
    weights  = map(lambda x: 1 / x**2, errtau)

    # Combine results.
    errcomb  = sum(weights)
    taucomb  = sum(map(lambda x, y: x*y, tau, weights)) / errcomb
    errcomb  = 1 / sqrt(errcomb)

    # Display results
    for i in range(len(sys.argv[1:])):
        print 'Channel 0 Threshold:', thresh[i], 'V'
        print 'Mean muon lifetime:', tau[i], 'microseconds'
        print 'Error on mean muon lifetime', errtau[i], 'microseconds'
        print

    print 'Combining all lifetimes:'
    print 'Mean muon lifetime:', taucomb, 'microseconds'
    print 'Error on mean muon lifetime', errcomb, 'microseconds'
