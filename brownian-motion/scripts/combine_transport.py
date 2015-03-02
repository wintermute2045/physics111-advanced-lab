#!/usr/bin/env python2.7

################################################################################
## This code is for the brownian motion lab from the UC Berkeley Physics
## Advanced Lab class.
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
import cPickle
import numpy as np
import scipy.constants as spc
import brownian_motion_cells as bmc

def usage(code):
    """
    Display the allowed options for the program.
    """
    print 'python combine_transport.py [options] pickle_files'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    sys.exit(code)

def get_values(info_dicts, key):
    """
    Extract values out of the dictionaries.
    """
    err_key = 'err_' + key
    transport = np.array(map(lambda d: (d[key], d[err_key]), info_dicts))
    return transport.transpose()

# Run the simulation.
if __name__ == '__main__':
    # Parse arguments
    if len(sys.argv) == 1:
        usage(1)
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Default conditions for the simulation
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
    if len(args) == 0:
        print 'No input pickles supplied.'
        usage(1)

    # Open pickle files and compute the weigted means.
    info_dicts = []
    for arg in args:
        with open(os.path.abspath(arg)) as pkl:
            info_dicts.append(cPickle.load(pkl))

    # Compute the mean RMS velocity
    vrms      = get_values(info_dicts, 'vel_rms')
    work      = get_values(info_dicts, 'mean_work')
    work_rate = get_values(info_dicts, 'work_rate')

    ave_vrms, err_vrms           = bmc.mean_with_err(*vrms)
    ave_work, err_work           = bmc.mean_with_err(*work)
    ave_work_rate, err_work_rate = bmc.mean_with_err(*work_rate)

    atp_energy = 0.4 * 20.5 / spc.N_A
    num_myosin = ave_work / atp_energy
    err_myosin = err_work / atp_energy

    print 'RMS speed:', ave_vrms, 'm/s'
    print 'Uncertainty:', err_vrms, 'm/s'
    print
    print 'Average work:', ave_work, 'J'
    print 'Uncertainty:', err_work, 'J'
    print
    print 'Average number of myosin motors:', num_myosin
    print 'Uncertainty:', err_myosin
    print
    print 'Average work rate:', ave_work_rate, 'J/m'
    print 'Uncertainty:', err_work_rate, 'J/m'
