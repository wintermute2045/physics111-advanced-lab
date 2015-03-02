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

def correct_time(nrows, data, pars):
    """
    Correct the time delays according to the calibration.
    """
    def __ct__(i):
        progress = 100 * float(i+1) / nrows
        print 'Time correction:', str(round(progress, 2)) + '%\r',
        return (data[i] - pars[0]) / pars[1]
    return __ct__

def expcurve(pars):
    """
    Exponential curve for easy use in maps.
    """
    def __ecurve__(arg):
        return pars[0] * np.exp(-pars[1]*arg) + pars[2]
    return __ecurve__


def threshold_plot(hdu, xtitle, plot_title, filename):
    """
    Make a plot of the lifetime as a function of the minimum threshold.
    """
    thresh   = map(lambda x: x[0],           hdu.data)
    lifetime = map(lambda x: 1 / x[4],       hdu.data)
    errlft   = map(lambda x: x[7] / x[4]**2, hdu.data)

    plt.figure()
    plt.plot(thresh, lifetime, 'b.')
    plt.errorbar(thresh, lifetime, yerr=errlft, fmt=None)
    plt.xlabel(xtitle)
    plt.ylabel(r'$\mu$ lifetime ($\mu s$)')
    plt.title(plot_title)
    print 'Saving figure to', filename + '.'
    plt.savefig(filename, bbox_inches='tight', format='pdf')


def subtract_data(nrows, col1, col2, data):
    """
    Subtract data and print progress of all of the subtraction.
    """
    def __subtract__(i):
        progress = 100 * float(i+1) / nrows
        print 'Subtraction:', str(round(progress, 2)) + '%\r',
        difference = data[i][col2] - data[i][col1]
        return 1e6 * difference
    return __subtract__


# Run everything
if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'c:i:o:q')
    except getopt.GetoptError as err:
        print str(err)
        sys.exit(1)
    infile    = ''
    outfile   = ''
    calibfile = ''
    display_plots = True
    for opt, arg in opts:
        if   opt == '-c':
            calibfile = arg
        elif opt == '-i':
            infile    = arg
        elif opt == '-o':
            outfile   = arg
        elif opt == '-q':
            display_plots = False
    if infile == '' or outfile == '':
        print 'Script cannot produce output without filename base.'
        sys.exit(1)
    if calibfile == '':
        print 'No calibration file specified.'
        sys.exit(1)

    # Get the calibration infomation
    hdulist    = pyfits.open(calibfile)
    cslope     = hdulist[1].header['FITM']
    cinter     = hdulist[1].header['FITB']
    calib_pars = [cinter, cslope]
    hdulist.close()

    # Open file and get header info
    print 'Generating plots from', infile + '.'
    print 'This may take a while.\n'
    hdulist  = pyfits.open(infile)
    hrange   = (0, hdulist[1].header['RANGE'])
    cutstart = hdulist[1].header['CUTSTART']
    cutend   = hdulist[1].header['CUTEND']
    nbins    = hdulist[1].header['NBINS']
    expmag   = hdulist[1].header['EXPMAG']
    explbda  = hdulist[1].header['EXPLBDA']
    explim   = hdulist[1].header['EXPLIM']
    thresh   = hdulist[1].header['CHTH0']
    data     = filter(lambda x: x[1] > 0.05 and x[4] > 0.05, hdulist[1].data)
    data     = filter(lambda x: x[2] < 1e-7 and x[5] < 1e-7, data)
    nrows    = len(data)
    data     = map(subtract_data(nrows, 0, 3, data),      range(nrows))
    data     = map(correct_time(nrows, data, calib_pars), range(nrows))

    # Create histogrm and make plots
    heights, bins = np.histogram(data, nbins, range=hrange)
    bin_width = (bins[-1] - bins[0]) / float(nbins)
    bins = bins[:-1]
    bins = map(lambda x: x + bin_width/2.0, bins)
    error = map(np.sqrt, heights)

    plot_title = 'Channel 0 Threshold = ' + str(thresh) + ' V'
    labels = ['Fit region', 'Background', 'Exponential fit']

    plt.figure()
    plt.plot(bins[cutstart:cutend], heights[cutstart:cutend],
            'k', label=labels[0])
    plt.plot(bins[:cutstart], heights[:cutstart], '#787878',)
    plt.plot(bins[cutend:], heights[cutend:], '#787878', label=labels[1])

    # Add fit to the plot
    fitpars = [expmag, explbda, explim]
    fitdata = map(expcurve(fitpars), bins)
    plt.plot(bins, fitdata, 'r', label=labels[2])

    # Titles
    plt.xlabel(r'Decay time ($\mu s$)')
    plt.ylabel('Counts')
    plt.title(plot_title)
    plt.yscale('log')
    plt.legend()

    # Save the plot
    plotfile = outfile + '_all.pdf'
    print 'Saving figure to', plotfile + '.'
    plt.savefig(plotfile, bbox_inches='tight', format='pdf')

    # Plot the fit region only
    plt.figure()
    plt.plot(bins[cutstart:cutend], heights[cutstart:cutend],
            'k.', label='Data')
    plt.errorbar(bins[cutstart:cutend], heights[cutstart:cutend],
            yerr=error[cutstart:cutend], fmt=None, ecolor='#555753')
    plt.plot(bins[cutstart:cutend], fitdata[cutstart:cutend],
            'r', label=labels[2])

    # Titles
    plt.xlabel(r'Decay time ($\mu s$)')
    plt.ylabel('Counts')
    plt.title(plot_title)
    plt.legend()

    # Save the plot
    plotfile = outfile + '_cut.pdf'
    print 'Saving figure to', plotfile + '.'
    plt.savefig(plotfile, bbox_inches='tight', format='pdf')

    # Get rid of peak in the beginning
    plt.figure()
    plt.plot(bins[cutstart:cutend], heights[cutstart:cutend],
            'k', label=labels[0])
    plt.plot(bins[cutend:], heights[cutend:], '#787878', label=labels[1])
    plt.plot(bins[cutstart:], fitdata[cutstart:], 'r', label=labels[2])

    # Titles
    plt.xlabel(r'Decay time ($\mu s$)')
    plt.ylabel('Counts')
    plt.title(plot_title)
    plt.legend()

    # Save the plot
    plotfile = outfile + '_cut_bgnd.pdf'
    print 'Saving figure to', plotfile + '.\n'
    plt.savefig(plotfile, bbox_inches='tight', format='pdf')

    # Create the plot of lifetime as a function of threshold.
    print 'Plotting threshold tests.\n'
    xtitle   = ['Threshold for pulse 1 (V)',
                'Threshold for pulse 2 (V)',
                'Threshold for both pulses (V)']
    filename = [outfile + '_thresh1.pdf',
                outfile + '_thresh2.pdf',
                outfile + '_boththresh.pdf']
    for i in range(3):
        threshold_plot(hdulist[i+2], xtitle[i], plot_title, filename[i])

    # Close fits file and optionally display plots.
    hdulist.close()
    if display_plots:
        plt.show()

