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

# This is just a testing plot done in class with very limited data. Accurate
# answers are not to be expected with the data passed into here.

import sys
import scipy as sp
import scipy.stats as sps
import matplotlib.pyplot as plt

def mkplot(current, frequency, isotope):
    """
    Make a plot of current vs frequency for an isotope of rubidium.
    """
    frequency = map(lambda x: 1e-6 * x, frequency)
    errors    = map(lambda x: 0.01, range(len(current)))

    # Make the plot
    plt.figure()
    plt.plot(frequency, current, 'ko')
    plt.errorbar(frequency, current, yerr=errors, fmt=None, ecolor='#555753')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Current (A)')
    plt.title(r'$^{' + str(isotope) + r'}Rb$')

    # Fit the curve
    slope, intercept, r, p, std_err = sps.linregress(current, frequency)
    fity = [min(current), max(current)]
    fitx = map(lambda x: slope*x + intercept, fity)
    plt.plot(fitx, fity, 'r')

    # Get nuclear spin, 135 turns in helmholtz coil of radius 27.5 cm
    nuc_spin = 2.799*0.9e-2 * 135 / (slope * 0.275)
    nuc_spin = 0.5 * (nuc_spin - 1)
    print 'Nuclear spin:', nuc_spin

    # Get earth's B field
    if isotope == 85:
        spin = 3.0 / 2
    else:
        spin = 5.0 / 2
    earth_field = (2*spin + 1) / 2.799 * intercept

    print 'Earth\'s magnetic field:', earth_field, 'gauss\n'

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Invalid arguments.'
        sys.exit(1)

    # Open the file
    infile = sys.argv[1]
    data = sp.genfromtxt(infile, skip_header=1, delimiter=',')

    frequency     = map(lambda x: x[0],      data)
    lower_current = map(lambda x: -0.1*x[2], data)
    upper_current = map(lambda x: -0.1*x[4], data)

    mkplot(lower_current, frequency, 85)
    mkplot(upper_current, frequency, 87)
    plt.show()
