#!/usr/bin/env python2.7

################################################################################
## This code is part of the analysis for the optical pumping lab
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

import os                   as _os
import sys                  as _sys
import scipy                as _sp
import matplotlib.pyplot    as _plt
from opt_pump import *

datadir = _os.path.dirname(_os.path.dirname(_os.path.abspath(_sys.argv[0])))
datadir += '/data'

def main():
    # Open the data files
    outward = datadir + '/85Outward.csv'
    inward  = datadir + '/85Inward.csv'

    outdata = _sp.genfromtxt(outward, delimiter=',', skip_header=1)
    indata  = _sp.genfromtxt(inward,  delimiter=',', skip_header=1)

    # Create plots
    pcoeffs, ncoeffs = mkplot(outdata, 'Current adjusted increasing', 'b.')
    mkplot(indata,  'Current adjusted decreasing', 'g.')

    # Set display options
    _plt.legend()
    _plt.xlabel('Current (A)')
    _plt.ylabel('Resonance (MHz)')
    _plt.title(r'$^{85}Rb$')

    _plt.show()

    # Calculate nuclear spin
    fit_b, fit_m, freq_err, err_b, err_m = pcoeffs
    measured_spin, err_spin = nuc_spin(fit_m, err_m, 0.29, 0.005, 135)
    spin = get_spin(measured_spin)

    # Calculate the earth's field
    outfield, sigoutfield = earth_field(outdata, spin)
    infield, siginfield   = earth_field(indata,  spin)

    sigfield = 1.0/sigoutfield**2 + 1.0/siginfield**2
    sigfield = 1.0/_sp.sqrt(sigfield)
    field    = outfield/sigoutfield**2 + infield/siginfield**2
    field   *= sigfield**2

    return (spin, measured_spin, field, sigfield, freq_err, ncoeffs[2])

if __name__ == '__main__':
    main()
