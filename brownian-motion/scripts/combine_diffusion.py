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
import brownian_motion_cells as bmc

def usage(code):
    """
    Display the allowed options for the program.
    """
    print 'python combine_diffusion.py [options] pickle_files'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -D: Diameter of each particle.'
    print '    -o: Output file to save the diffusion to.'
    print '    -T: Temperature of the system.'
    print '    -v: Viscosity of the solution.'
    sys.exit(code)

def get_diffusion(info_dicts, diff_type):
    """
    Get the diffusions from the list of info dictionaries.
    """
    # The pickled dictionaries have keys that begin with fit and ave.
    if diff_type != 'ave' and diff_type != 'fit':
        raise ValueError('Invalid diffusion type.')

    diff_key = diff_type + '_diffusion'
    err_key  = diff_key  + '_err'

    diffusion = np.array(map(lambda d: (d[diff_key], d[err_key]), info_dicts))
    return diffusion.transpose()

# Run the simulation.
if __name__ == '__main__':
    # Parse arguments
    if len(sys.argv) == 1:
        usage(1)
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hD:o:T:v:')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Default conditions for the simulation
    radius = None
    outfile = None
    temperature = 293.0
    viscosity = None
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-D':
            radius = float(arg) / 2
        elif opt == '-o':
            outfile = os.path.abspath(arg)
        elif opt == '-T':
            temperature = float(arg)
        elif opt == '-v':
            viscosity = float(arg)
    if len(args) == 0:
        print 'No input pickles supplied.'
        usage(1)

    # Compute the theoretical diffusion value if all required parameters have
    # been supplied.
    has_theory = False
    if radius != None and viscosity != None:
        diffusion_theory = bmc.theoretical(radius, viscosity, temperature)
        has_theory = True

    # Open pickle files and compute the weigted means.
    info_dicts = []
    for arg in args:
        with open(os.path.abspath(arg)) as pkl:
            info_dicts.append(cPickle.load(pkl))

    ave_diffusion = get_diffusion(info_dicts, 'ave')
    fit_diffusion = get_diffusion(info_dicts, 'fit')

    ave_diffusion, err_ave_diffusion = bmc.mean_with_err(*ave_diffusion)
    fit_diffusion, err_fit_diffusion = bmc.mean_with_err(*fit_diffusion)

    if outfile != None:
        diffusion = dict()
        diffusion['diffusion']   = fit_diffusion
        diffusion['uncertainty'] = err_fit_diffusion
        print 'Saving diffusion coefficient to', outfile + '.'
        with open(outfile, 'w') as pkl:
            cPickle.dump(diffusion, pkl)

    # Print some results.
    print
    if has_theory:
        print 'Particle radius:', radius, 'm'
        print 'Viscosity: %g Ps' % viscosity
        print 'Temperature:', temperature, 'K'
        print 'Theoretical diffusion:', diffusion_theory
        print
    # Assume that if the theoretical value is not present, then this script is
    # being used to analyze vesicles in the red onion cells.
    else:
        radius = 2 * 0.000000222
        # When solving for viscosity in terms of diffusion, eta and D switch
        # places in the eqation, so the bmc.theoretical function can be used to
        # find the viscosity.
        viscosity = bmc.theoretical(radius, fit_diffusion, temperature)
        print 'Particle radius:', radius, 'm'
        print 'Viscosity: %g cP' % (viscosity * 1e3)
        print 'Temperature:', temperature, 'K'
        print


    print 'Averaging displacements:'
    print 'Diffusion:', ave_diffusion
    print 'Uncertainty:', err_ave_diffusion
    print
    print 'Curve fitting:'
    print 'Diffusion:', fit_diffusion
    print 'Uncertainty:', err_fit_diffusion

