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
import scipy.stats          as _sps
import matplotlib.pyplot    as _plt

# Quickly get elements out of tuples.
def fst(tup): return tup[0]
def snd(tup): return tup[1]
def trd(tup): return tup[2]

def earth_field(data, spin):
    """
    Calculate the Earth's magnetic field from the data.
    """
    current   = map(fst, data)
    frequency = map(trd, data)
    positive  = filter(lambda t: t[0] > 0, zip(current, frequency))
    negative  = filter(lambda t: t[0] < 0, zip(current, frequency))

    # Calculate the field
    pairs = []
    for curr, freq in positive:
        neg_freq = filter(lambda t: abs(t[0]) == curr, negative)[0][1]
        pairs.append((freq, neg_freq))
    fields      = map(lambda t: t[0]-t[1], pairs)
    fields      = map(lambda x: 0.5 * x * (2*spin+1) / 2.799, fields)
    mean_field  = _sp.mean(fields)
    sigma_field = _sp.std(fields, ddof=1) / _sp.sqrt(len(fields))
    return (mean_field, sigma_field)

def get_err(xdata, ydata, coefs):
    """
    This function computes the error in the data necessry for chi^2 to
    be optimal. The format of the coefficients tuple is (b, m).
    """
    if len(coefs) != 2:
        print 'Invalid fit coefficients.'
        _sys,exit(1)
    if len(xdata) != len(ydata):
        print 'X and Y datasets must have the same legnth.'
        _sys.exit(1)

    # Get the coefficients
    fitb, fitm = tuple(coefs)
    ndf = len(xdata) - 2

    # Compute the error
    sigma = 0
    for i in range(ndf+2):
        sigma += (ydata[i] - fitm*xdata[i] - fitb) ** 2
    sigma /= float(ndf)

    return _sp.sqrt(sigma)

def get_err_coeff(xdata, ydata, erry):
    """
    Get the error on the linear fit coefficients.
    """
    meanx  = sum(_sp.array(xdata) / erry**2) / len(xdata)
    xprime = _sp.array(xdata) - meanx

    err_m = erry / _sp.sqrt(sum(xprime**2))
    err_b = _sp.sqrt(erry**2 / len(xdata) + (meanx * err_m)**2)

    return (err_b, err_m)

def get_spin(measured):
    """
    Determine the actual spin based off of the measured value.
    """
    measured *= 2
    spin1 = int(measured)
    spin2 = spin1 + 1
    if abs(measured-spin1) > abs(measured-spin2):
        return float(spin2) / 2
    else:
        return float(spin1) / 2

def mkplot(raw_data, plabel, color):
    """
    This function creates a plot from raw data collected in the optical
    pumping experiment.
    """
    current   = _sp.array(map(lambda l: -0.1*l[1], raw_data))
    frequency = _sp.array(map(lambda l:      l[2], raw_data))

    # Plot current and frequency
    if len(plabel):
        _plt.plot(current, frequency, color, label=plabel)
    else:
        _plt.plot(current, frequency, color)

    # Separate data and perform fits
    positive = filter(lambda t: t[0] >= 0, zip(current, frequency))
    pos_curr = map(fst, positive)
    pos_freq = map(snd, positive)
    negative = filter(lambda t: t[0] <  0, zip(current, frequency))
    neg_curr = map(fst, negative)
    neg_freq = map(snd, negative)

    # Perform fits and get the error on RF.
    pos_m, pos_b, r, p, std_err = _sps.linregress(pos_curr, pos_freq)
    neg_m, neg_b, r, p, std_err = _sps.linregress(neg_curr, neg_freq)
    pos_err = get_err(pos_curr, pos_freq, (pos_b, pos_m))
    neg_err = get_err(neg_curr, neg_freq, (neg_b, neg_m))
    perr_b, perr_m = get_err_coeff(pos_curr, pos_freq, pos_err)
    nerr_b, nerr_m = get_err_coeff(neg_curr, neg_freq, neg_err)

    # Plot the errors and the fit curve.
    p_err_array = [pos_err for p in pos_freq]
    n_err_array = [neg_err for n in neg_freq]
    _plt.errorbar(pos_curr, pos_freq, p_err_array, fmt=None, ecolor='#555753')
    _plt.errorbar(neg_curr, neg_freq, n_err_array, fmt=None, ecolor='#555753')

    # Plot the fit curves
    fit_pcurr = _sp.array([min(pos_curr), max(pos_curr)])
    fit_ncurr = _sp.array([min(neg_curr), max(neg_curr)])
    fit_pfreq = _sp.array(map(lambda x: pos_m*x + pos_b, fit_pcurr))
    fit_nfreq = _sp.array(map(lambda x: neg_m*x + neg_b, fit_ncurr))
    _plt.plot(fit_pcurr, fit_pfreq, color[0])
    _plt.plot(fit_ncurr, fit_nfreq, color[0])

    pos_coeffs = (pos_b, pos_m, pos_err, perr_b, perr_m)
    neg_coeffs = (neg_b, neg_m, neg_err, nerr_b, nerr_m)
    return (pos_coeffs, neg_coeffs)

def nuc_spin(slope, err_slope, radius, err_radius, nturns):
    """
    Calculate the nuclear spin and errors.
    """
    spin = 2.799 * 0.9e-2 * nturns / abs(slope * radius)
    spin = 0.5*(spin - 1)

    err_spin = (err_slope / slope)**2 + (err_radius / radius)**2
    err_spin = spin * _sp.sqrt(err_spin)

    return (spin, err_spin)
