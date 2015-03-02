#!/usr/bin/env python2.7

################################################################################
## This code is for problem 4 of the EAX homework for Physics 111 Advanced Lab
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
import scipy as sp
import scipy.special as ss
import scipy.optimize as spo
import matplotlib.pyplot as plt

def find_true(items, condition):
    """
    Finds the first element of a list satisfying some condition. Returns -1 if
    no items satisfy the condition

    (a -> Bool) -> [a] -> Int
    """
    if len(items) == 0:
        return -1
    elif condition(items[0]):
        return 0
    else:
        return 1 + find_true(items[1:], condition)

def gaussian(x, A, mean, sigma):
    """
    Gaussian distribution for fitting purposes.
    """
    return A * sp.exp(-(x - mean)**2 / (2. * sigma**2))

def eval_gaus(par):
    """
    For evaluating gaussians in maps.
    """
    def __gaus__(x):
        return gaussian(x, par[0], par[1], par[2])
    return __gaus__

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

    # Create histogram from data
    print 'Creating histogram.'
    hdulist = pyfits.open(basename + '.fits')
    data = map(lambda x: x[0], hdulist[1].data)
    nbins = hdulist[1].header['HBINS']
    bin_range = (hdulist[1].header['HMIN'], hdulist[1].header['HMAX'])
    offset = (bin_range[1] - bin_range[0]) / (2*nbins)
    n, bins, patches = plt.hist(data,bins=nbins,range=bin_range,histtype='step')
    bins = map(lambda x: x + offset, bins[:-1])
    error = map(sp.sqrt, n)
    errorbars = plt.errorbar(bins, n, yerr=error, fmt=None)

    # Only use the histogram bins with data in plotting the fitted data
    start = find_true(n, lambda x: x > 0)
    end = 1 + reduce(lambda i,j: j if n[j] else i, range(len(n)))

    # Plot gaussian with fitted data
    magnitude = hdulist[1].header['HMAG']
    mean = hdulist[1].header['HMEAN']
    sigma = hdulist[1].header['SIGMA']
    gpars = [magnitude, mean, sigma]
    gauss_x = map(lambda x: bins[start] + x*(bins[end]-bins[start])/1000.0,
            range(1001))
    gauss_y = map(eval_gaus(gpars), gauss_x)
    gauss_plot = plt.plot(gauss_x, gauss_y, 'k--')

    # Plot titles
    title = 'Energies from gamma ray experiment'
    plt.xlabel('Energy')
    plt.ylabel('Counts')
    plt.title(title)
    plt.legend((patches[0], errorbars, gauss_plot[0]),
            ('Data', 'Errors', 'Gaussian fit'))

    print 'Saving plot to', basename + '.pdf'
    plt.savefig(basename + '.pdf', bbox_inches='tight')
    if display_plots:
        plt.show()

