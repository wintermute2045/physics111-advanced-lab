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

def mkhist(fitshdu, colnum, nbins, filename):
    """
    Make a histogram from a given fits hdu.
    """
    print 'Saving plot to', filename
    plt.figure()
    data = map(lambda x: x[colnum], fitshdu.data)
    n, bins, patches = plt.hist(data, bins=nbins, range=(-1,1), histtype='step')
    nsamples = int(fitshdu.header['NSAMPLES'])
    title = 'Means from sampling from a gaussian ' + str(nsamples) + ' times'
    plt.xlabel('Samples')
    plt.ylabel('Counts')
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

    # Open file and make plots
    hdulist = pyfits.open(basename + '.fits')
    for i in range(1, len(hdulist)):
        mkhist(hdulist[i], 0, 200, basename + '_' + str(i) + '.pdf')
    if display_plots:
        plt.show()
