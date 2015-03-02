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
import getopt
import pyfits
import matplotlib.pyplot as plt

# Run everything
if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:o:q')
    except getopt.GetoptError as err:
        print str(err)
        sys.exit(1)
    infile = ''
    outfile = ''
    display_plots = True
    for opt, arg in opts:
        if opt == '-i':
            infile = arg
        elif opt == '-o':
            outfile = arg
        elif opt == '-q':
            display_plots = False
    if infile == '' or outfile == '':
        print 'Script cannot produce output without filename base.'
        sys.exit(1)

    # Open file and get header info
    print 'Generating plot from', infile + '.'
    hdulist = pyfits.open(infile)
    nch0    = hdulist[1].header['NCH0']
    cslope  = hdulist[1].header['FITM']
    cinter  = hdulist[1].header['FITB']
    intime  = map(lambda x: x[0], hdulist[1].data)
    outtime = map(lambda x: x[1], hdulist[1].data)
    errtime = map(lambda x: x[2], hdulist[1].data)
    hdulist.close()

    if nch0:
        plot_title = 'Clock calibration using Channel 0'
    else:
        plot_title = 'Clock calibration using Channel 1'

    # Make plots
    plt.plot(intime, outtime, 'wo')
    plt.errorbar(intime, outtime, yerr=errtime, fmt=None, ecolor='#555753')

    # Add fitted curve to the plot
    fitx = [0, int(max(intime))+1]
    fity = map(lambda x: cslope * x + cinter, fitx)
    plt.plot(fitx, fity, 'r')
    plt.axis([0, int(fity[1])+1, 0, int(fity[1])+1])

    # Add labels and stuff
    plt.xlabel(r'Input time delay ($\mu s$)')
    plt.ylabel(r'Measured time delay ($\mu s$)')
    plt.title(plot_title)

    print 'Saving figure to', outfile + '.'
    plt.savefig(outfile, bbox_inches='tight', format='pdf')
    if display_plots:
        plt.show()
