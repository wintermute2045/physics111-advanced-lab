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
import numpy as np
import matplotlib.pyplot as plt

def fitcurve(pars):
    """
    Fits the log of the bins.
    """
    def __fc__(arg):
        val = pars[0] + pars[1] * arg
        return 10**val
    return __fc__

def subtract_data(nrows, col1, col2, data):
    """
    Subtract data and print progress of all of the subtraction.
    """
    def __subtract__(i):
        progress = 100 * float(i+1) / nrows
        print 'Percent complete:', str(round(progress, 2)) + '%\r',
        difference = data[i][col2] - data[i][col1]
        return 1e6 * difference
    return __subtract__


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
    print 'This may take a while.\n'
    hdulist  = pyfits.open(infile)
    hrange   = (0.0, 30.0)
    nbins    = hdulist[1].header['NBINS']
    nch0     = hdulist[1].header['NCH0']
    mean     = hdulist[1].header['HMEAN']
    sigma    = hdulist[1].header['HSIGMA']
    fitb     = hdulist[1].header['FITB']
    fitm     = hdulist[1].header['FITM']

    if nch0:
        col1 = 0
        col2 = 3
        plot_title = 'Channel 0 Efficiency'
    else:
        col1 = 6
        col2 = 9
        plot_title = 'Channel 1 Efficiency'

    data  = filter(lambda x: abs(x[col2+1] - mean) < sigma, hdulist[1].data)
    nrows = len(data)
    data  = map(subtract_data(nrows, col1, col2, hdulist[1].data), range(nrows))
    hdulist.close()

    # Create the histogram
    n, bins, patches = plt.hist(data, nbins, range=hrange,
            histtype='step',color='k')
    offset    = (bins[nbins] - bins[0]) / (2*nbins)
    bins      = map(lambda x: x + offset, bins[:-1])
    error     = map(np.sqrt, n)
    errorbars = plt.errorbar(bins, n, yerr=error, fmt=None, ecolor='#555753')

    # Add the fitted data to the plot.
    fpars = [fitb, fitm]
    xfit    = map(lambda x: x*hrange[1] / 1000.0, range(1001))
    fitdata = map(fitcurve(fpars), xfit)
    plt.plot(xfit, fitdata, 'r')

    # Add labels and stuff
    plt.xlabel(r'$\Delta t$ ($\mu s$)')
    plt.ylabel('Counts')
    plt.title(plot_title)
    plt.yscale('log')

    print 'Saving figure to', outfile + '.'
    plt.savefig(outfile, bbox_inches='tight', format='pdf')
    if display_plots:
        plt.show()
