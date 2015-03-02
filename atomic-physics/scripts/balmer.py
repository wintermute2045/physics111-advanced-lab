#!/usr/bin/env python2.7

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

import os
import sys
import getopt
import scipy as sp
import scipy.optimize as spo
import matplotlib.pyplot as plt

# Cython shared object
import gsl_wlsq as lsq

# Simple linear function used for curve fitting
linear = lambda x, m, b: m*x + b

def usage(code):
    """
    Explain the inputs of the program.
    """
    print 'python balmer.py -H hydrogen.txt -M mercury.txt [-q]'
    print 'Flags:'
    print '    -h: Show this help message and exit.'
    print '    -H: Hydrogen data file.'
    print '    -m: Mercury data file.'
    print '    -q: Suppress plots.'
    sys.exit(1)

def chisq(observed, expected, sigma=None):
    """
    I have no idea what the hell is going on in the chi^2 test that
    scipy ships with...

    Hey look! I have everything that I need my least squares
    test to do in only 4 lines of code! IT JUST WERKS!
    """
    diffsq = (observed - expected)**2
    if sigma != None:
        diffsq /= sigma**2
    return sum(diffsq)

def correct_lbda(slope, intercept):
    """
    Correct the wavelength from the plot of measured lbda vs true
    lbda. This function was written for maps.
    """
    def _correct_lbda(lbda_obs):
        lbda = (lbda_obs - intercept) / slope
        return lbda
    return _correct_lbda

def get_err(xdata, ydata, coefs):
    """
    Get the error on the ydata assuming that the fit gives a chi^2/ndf
    of 1. This function also returns the errors on the fit
    coefficients. It shall be assumed that there are no errors on the
    x data, since that's the only situation where this function will be
    actually used. It shall also be assumed that the coefs variable is
    in the form (slope, intercept). I'm also assuming that the xdata
    and ydata lists have the same length.
    """
    # Convert x and y lists to arrays for easier manipulation
    xlen   = len(xdata)
    ndf    = xlen - 2
    xdata  = sp.array(xdata)
    ydata  = sp.array(ydata)
    xmean  = sp.mean(xdata)
    xsqsum = sum(xdata**2)

    # Get the average of the ydata
    m, b  = tuple(coefs)
    sigma = sum((ydata - m*xdata - b)**2)
    sigma = sp.sqrt(sigma / ndf)

    # Get the average of the coefficients.
    sigma_m = sp.sqrt(sigma**2 / xsqsum)
    sigma_b = sp.sqrt(sigma**2 / xlen + (xmean * sigma_m)**2)

    return (sigma, sigma_m, sigma_b)

def get_vac(lbda, sigma):
    """
    Get the vacuum wavelegnth for an air wavelength.
    """
    # Index of refraction
    n = 1 + 6432.8e-8
    n += 2929810. / (15768. - sigma**2)
    n += 25540. / (4228. - sigma**2)
    lbda_vac = n * lbda
    return lbda_vac

def max_snd(tup_list):
    """
    Get the max tuple by comparing the second element.
    """
    if len(tup_list) == 1:
        return tup_list[0]
    else:
        head = tup_list[0]
        max_tail = max_snd(tup_list[1:])
        if head[1] > max_tail[1]:
            return head
        else:
            return max_tail

def qs_fst(tuplist):
    """
    Run a quicksort on a list of tuples based on the first item.
    """
    if tuplist == []:
        return []
    else:
        fst = tuplist[0]
        smaller = filter(lambda t: t[0] <= fst[0], tuplist[1:])
        larger  = filter(lambda t: t[0] >  fst[0], tuplist[1:])
        return smaller + [fst] + larger

# Wavelength tables
lbda_mercury  = { 'red'    : 6149.50,
                  'yellow1': 5790.66,
                  'yellow2': 5769.60,
                  'green'  : 5460.74,
                  'blue'   : 4358.33,
                  'indigo' : 4046.56
                }
# Wavenumber tables
sigma_mercury = { 'red'    : 16256.986,
                  'yellow1': 17264.401,
                  'yellow2': 17327.439,
                  'green'  : 18397.479,
                  'blue'   : 22938.156,
                  'indigo' : 24705.376
                }
sigma_hydrogen = {'alpha'  : 15233.377,
                  'beta'   : 20564.758,
                  'gamma'  : 23032.505,
                  'delta'  : 23023.505,
                  'epsilon': 25181.336
                 }
# Quantum_numbers
quantum_num    = {'alpha'  : 3,
                  'beta'   : 4,
                  'gamma'  : 5,
                  'delta'  : 6,
                  'epsilon': 7,
                  'zeta'   : 8,
                  'eta'    : 9
                 }

if __name__ == '__main__':
    # Set up the environment.
    figdir = os.path.dirname(os.path.abspath(sys.argv[0]))
    figdir = os.path.dirname(figdir) + '/plots'
    os.system('mkdir -p ' + figdir)

    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hH:m:q')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)
    hydrogen = None
    mercury  = None
    display_plots = True
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-H':
            hydrogen = arg
        elif opt == '-m':
            mercury = arg
        elif opt == '-q':
            display_plots = False
    if hydrogen == None:
        print 'Hydrogen file missing!'
        usage(1)
    if mercury  == None:
        print 'Mercury file missing!'
        usage(1)

    # Open data files
    mercury  = map(tuple, sp.genfromtxt(mercury))
    hydrogen = map(tuple, sp.genfromtxt(hydrogen))
    m_lbda   = [merc[0] for merc in mercury]
    m_volt   = [merc[1] for merc in mercury]
    h_lbda   = [hydr[0] for hydr in hydrogen]
    h_volt   = [hydr[1] for hydr in hydrogen]

    # Get maximums
    merc_max = {
            'red'    : filter(lambda m: m[0] > 6000 and m[0] < 6400, mercury),
            'yellow2': filter(lambda m: m[0] > 5750 and m[0] < 6000, mercury),
            'green'  : filter(lambda m: m[0] > 5400 and m[0] < 5750, mercury),
            'blue'   : filter(lambda m: m[0] > 4350 and m[0] < 4750, mercury),
            'indigo' : filter(lambda m: m[0] > 4000 and m[0] < 4350, mercury)
            }
    merc_lbda = []
    merc_volt = []
    merc_true = []
    for key in merc_max.keys():
        if merc_max[key] != None:
            max_lbda = list(max_snd(merc_max[key]))
            max_lbda.append(lbda_mercury[key])
            merc_max[key] = tuple(max_lbda)
            merc_lbda.append(merc_max[key][0])
            merc_volt.append(merc_max[key][1])
            merc_true.append(merc_max[key][2])

    # Fit the curve and do a wavelength calibration on hydrogen
    # TODO include errors from line widths.
    hg_popt, hg_pcov = spo.curve_fit(linear, merc_true, merc_lbda)
    hg_m, hg_b = tuple(hg_popt)
    h_lbda   = map(correct_lbda(hg_m, hg_b), h_lbda)
    hydrogen = zip(h_lbda, h_volt)

    # Get peaks for hydrogen lines
    hydr_max = {
            'alpha'  : filter(lambda h: h[0] > 6500 and h[0] < 7000, hydrogen),
            'beta'   : filter(lambda h: h[0] > 4600 and h[0] < 5100, hydrogen),
            'gamma'  : filter(lambda h: h[0] > 4250 and h[0] < 4500, hydrogen),
            'delta'  : filter(lambda h: h[0] > 4000 and h[0] < 4250, hydrogen),
            'epsilon': filter(lambda h: h[0] > 3900 and h[0] < 4050, hydrogen),
            'zeta'   : filter(lambda h: h[0] > 3880 and h[0] < 3920, hydrogen)
            #'eta'    : filter(lambda h: h[0] > 3790 and h[0] < 3850, hydrogen)
            }
    hydr_lbda = []
    hydr_volt = []
    hydr_qntm  = []
    for key in hydr_max.keys():
        if hydr_max[key] != None:
            hydr_max[key] = max_snd(hydr_max[key])
            hydr_lbda.append(hydr_max[key][0])
            hydr_volt.append(hydr_max[key][1])
            hydr_qntm.append(quantum_num[key])
    hydr_lbda  = 1.0  * sp.array(hydr_lbda)
    hydr_qntm  = 1.0  / sp.array(hydr_qntm)**2
    invr_lbda  = 1e10 / hydr_lbda

    # TODO add sigma to the fitter
    hydr_lsq = lsq.linear_fit(hydr_qntm, invr_lbda)
    rydberg  = hydr_lsq[0]
    offset   = hydr_lsq[2]

    # Calculate the weighted average of the rydberg constant
    coefs = (rydberg, offset)
    err_inv_lbda, err_m, err_b = get_err(hydr_qntm, invr_lbda, coefs)
    rydberg1 = abs(rydberg)
    rydberg2 = abs(4 * offset)
    err_b = 4 * err_b # Simply scale the error because the value is scaled
    weight1, weight2 = 1.0 / err_m**2, 1.0 / err_b**2
    sum_weights= weight1 + weight2
    weighted_sum = weight1 * rydberg1 + weight2 * rydberg2
    rydberg_ave = weighted_sum / sum_weights
    rydberg_err = sp.sqrt(1.0 / sum_weights)

    # Plot some mercury results
    plt.figure()
    plt.plot(m_lbda, m_volt, label='Measured spectrum')
    plt.plot(merc_lbda, merc_volt, 'ro', label='Measured peaks')
    plt.plot(merc_true, merc_volt, 'go', label='Actual peaks')

    plt.legend(loc=(0.02,0.02))
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Voltage (V)')
    plt.title('Spectrum calibration with Mercury')
    plt.yscale('log')
    plt.savefig(figdir + '/merc_spec.pdf', bbox_inches='tight')

    # Wavelength calibration plot
    plt.figure()
    plt.plot(merc_true, merc_lbda, 'ko', label='Data')

    # Create fit lines
    fit_hg_x = [min(merc_true), max(merc_true)]
    fit_hg_y = [hg_m * x + hg_b for x in fit_hg_x]
    plt.plot(fit_hg_x, fit_hg_y, 'r', label='Calibration fit')

    plt.legend(loc=(0.02,0.85))
    plt.xlabel(r'$\lambda_{true}$')
    plt.ylabel(r'$\lambda_{obs}$')
    plt.title('Spectrum calibration with Mercury')
    plt.savefig(figdir + '/merc_calib.pdf', bbox_inches='tight')

    # Ballmer series plots
    plt.figure()
    plt.plot(h_lbda, h_volt, label='Spectrum')
    plt.plot(hydr_lbda, hydr_volt, 'r.', label='Balmer peaks')

    plt.legend(loc=(0.02,0.85))
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Voltage (V)')
    plt.title('Balmer series')
    plt.yscale('log')
    plt.savefig(figdir + '/hydrogen_spec.pdf', bbox_inches='tight')

    # Let's see if quantum mechanics really works
    plt.figure()
    plt.plot(hydr_qntm, invr_lbda, 'ko', label='Data')

    # Create fit lines
    fitx = [min(hydr_qntm), max(hydr_qntm)]
    fity = [rydberg * x + offset for x in fitx]
    plt.plot(fitx, fity, 'r', label='Rydberg fit')

    plt.legend()
    plt.xlabel(r'$1 / n^2$')
    plt.ylabel(r'$1 / \lambda\ (m^{-1})$')
    plt.title('Balmer series')
    plt.savefig(figdir + '/rydberg.pdf', bbox_inches='tight')

    # Finally! I never thought we would make it!
    print 'Rydberg constant 1:', rydberg_ave, '+=', rydberg_err, 'm^{-1}'
    if display_plots:
        plt.show()
