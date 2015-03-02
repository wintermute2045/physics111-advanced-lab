#!/usr/bin/env python2.7

################################################################################
## This code is for problem 3 of the EAX homework for Physics 111 Advanced Lab
## at UC Berkeley.
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
import getopt
import pyfits
import matplotlib.pyplot as plt

def mkplot(fitshdu, title, filename):
    """
    Make a histogram from a given fits hdu.
    """
    print 'Saving plot to', filename
    plt.figure()
    xdata = map(lambda x: x[0], fitshdu.data)
    ydata = map(lambda x: x[1], fitshdu.data)
    plt.plot(xdata, ydata, 'wo')

    # plot fit
    m = fitshdu.header['M']
    b = fitshdu.header['B']
    fity = map(lambda x: m*x + b, xdata)
    plt.plot(xdata, fity, 'r')

    # plot error bars if they exit
    if len(fitshdu.data[0]) == 3:
        errors = map(lambda x: x[2], fitshdu.data)
        plt.errorbar(xdata, ydata, yerr=errors, fmt=None, ecolor='#000000')

    # Add a title and some crap like that
    plt.xlabel('Current (A)')
    plt.ylabel('Frequency (MHz)')
    plt.title(title)
    plt.savefig(filename, bbox_inches='tight')

if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:p')
    except getopt.GetoptError as err:
        print str(err)
        sys.exit(1)
    basename = ''
    display_plots = False
    for opt, arg in opts:
        if opt == '-i':
            basename = arg
        elif opt == '-p':
            display_plots = True
    if basename == '':
        print 'Script cannot produce output without filename base.'
        sys.exit(1)

    plot_titles = ['',
            'Unweighted fit',
            r'Fit with $\sigma = 0.01$ MHz',
            r'Fit with $\sigma = 1$ MHz',
            r'Fit with $\sigma = 0.3 + 0.3f$ MHz']

    # Open file and make plots
    hdulist = pyfits.open(basename + '.fits')
    for i in range(1, len(hdulist)):
        mkplot(hdulist[i], plot_titles[i], basename + '_' + str(i) + '.pdf')
    if display_plots:
        plt.show()
