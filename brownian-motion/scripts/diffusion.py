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
import itertools as it
import numpy.random as npr
import scipy.optimize as spo
import scipy.constants as spc
import matplotlib.pyplot as plt
import brownian_motion_cells as bmc

def usage(code):
    """
    Display the allowed options for the program.
    """
    print 'python diffusion.py -d data.txt [options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -d: Data file of the particle tracks.'
    print '    -f: File to save a plot to.'
    print '    -D: Diameter of each particle.'
    print '    -p: Pickle file to save to.'
    print '    -q: Do not display plots, but print messages.'
    print '    -Q: Do not print messages and do not display plots.'
    print '    -T: Temperature of the system.'
    print '    -v: Viscosity of the solution.'
    sys.exit(code)

def get_diffusion_ave(time_step):
    """
    Get the diffusion from a track in the data. This is designed only
    for 2D data, since that is what gets recorded by the camera.

    Elements of the track arrays:
    x y time dx dy dt dr^2 TotalDisplacementSquared
    0 1 2    3  4  5  6    7
    """
    def __diffusion__(track):
        # The first line in the data file is the origin of the particle. It
        # needs to be discarded when computing the diffusion coefficient.
        disp_sq = track[6][1:]

        # The squared displacements were recorded for each step.
        nsteps = len(disp_sq)
        dispsq_mean = np.mean(disp_sq)
        dispsq_std  = np.std(disp_sq, ddof=1) / np.sqrt(nsteps)

        ndim = 2
        diffusion = dispsq_mean / (2 * ndim * time_step)
        err_diffusion = dispsq_std / (2 * ndim * time_step)

        # Return the results as a tuple.
        return (diffusion, err_diffusion)
    return __diffusion__

def linear(x, m, b):
    """
    Linear function for fitting.
    """
    return m*x + b

# Run the simulation.
if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hd:D:f:p:qQT:v:')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Default conditions for the simulation
    data_file = None
    radius = None
    fig_file = None
    pkl_file = None
    temperature = 293.0
    viscosity = None
    quiet = False
    extra_quiet = False
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-d':
            data_file = os.path.abspath(arg)
        elif opt == '-D':
            radius = float(arg) / 2
        elif opt == '-f':
            fig_file = os.path.abspath(arg)
        elif opt == '-p':
            pkl_file = os.path.abspath(arg)
        elif opt == '-q':
            quiet = True
        elif opt == '-Q':
            extra_quiet = True
        elif opt == '-T':
            temperature = float(arg)
        elif opt == '-v':
            viscosity = float(arg)
    if data_file == None:
        print 'No data supplied. Cannot compute diffusion coefficient.'
        sys.exit(1)
    else:
        print 'Reading data from', data_file + '.'

    # Compute the theoretical diffusion value if
    # all required parameters have been supplied.
    has_theory = False
    if radius != None and viscosity != None:
        diffusion_theory = bmc.theoretical(radius, viscosity, temperature)
        has_theory = True

    # Open and read the data file.
    all_tracks = bmc.read_data(data_file)
    # Only accept tracks where there is displacement in each step.
    all_tracks = filter(lambda t: all(t[6][1:]), all_tracks)
    nparticles = len(all_tracks)

    # Compute the diffusions.
    time_step = 0.0330033 # The time step isn't recorded uniformly.
    all_diffusion = map(get_diffusion_ave(time_step), all_tracks)
    diff_list     = [d[0] for d in all_diffusion]
    diff_err_list = [d[1] for d in all_diffusion]
    diffusion, err_diffusion = bmc.mean_with_err(diff_list, diff_err_list)

    # Create average displacements and fit a curve
    ave_dispsq = bmc.average_dispsq(all_tracks, 50)
    num_avg = len(ave_dispsq[0])
    time = np.linspace(1, num_avg, num_avg) * time_step
    popt, pcov = spo.curve_fit(linear, time, ave_dispsq[0], sigma=ave_dispsq[1])
    fit_data = popt[0] * time + popt[1]
    fit_diff = popt[0] / 4.0
    fit_err  = np.sqrt(pcov[0][0]) / 4.0
    axis_bounds  = [min(time), max(time)]
    axis_bounds += [0, max(max(fit_data), max(ave_dispsq[0]))]

    # Create a dictionary to save to and save it to a pickle.
    if pkl_file != None:
        brownian_info = dict()
        brownian_info['ave_diffusion'] = diffusion
        brownian_info['ave_diffusion_err'] = err_diffusion
        brownian_info['num_particles'] = nparticles
        brownian_info['fit_diffusion'] = fit_diff
        brownian_info['fit_diffusion_err'] = fit_err
        print
        print 'Writing results to', pkl_file
        with open(pkl_file, 'w') as pkl:
            cPickle.dump(brownian_info, pkl)

    # Print some results.
    if not extra_quiet:
        print
        if has_theory:
            print 'Particle radius:', radius, 'm'
            print 'Viscosity: %g Ps' % viscosity
            print 'Temperature:', temperature, 'K'
            print 'Theoretical diffusion:', diffusion_theory
            print

        print 'Averaging displacements:'
        print 'Diffusion:', diffusion
        print 'Uncertainty:', err_diffusion
        print 'Number of particles tracks:', nparticles
        print
        print 'Curve fitting:'
        print 'Diffusion:', fit_diff
        print 'Uncertainty:', fit_err

    if fig_file == None and (quiet or extra_quiet):
        sys.exit()

    # Create plots.
    colors = map(tuple, npr.random((nparticles, 3)))

    # Plot the square displacements.
    plt.figure()
    for i in range(nparticles):
        plt.plot(all_tracks[i][2][1:], all_tracks[i][7][1:], color=colors[i])

    plt.plot(time, ave_dispsq[0], '#555753', linewidth=4, label='Average')
    plt.plot(time, fit_data, 'k', linewidth=4, label='Linear fit')
    plt.axis(axis_bounds)

    ax = plt.gca()
    yticks = [1e12 * tick for tick in ax.get_yticks()]
    ax.set_yticklabels(map(str, yticks))
    plt.xlabel('Time (s)')
    plt.ylabel(r'Displacement Squared ($\mu m^2$)')
    plt.title('Combined Particle Displacements')
    plt.legend(loc='upper left')
    plt.tight_layout()

    if fig_file == None:
        plt.show()
    else:
        print
        print 'Saving plot to', fig_file + '.'
        plt.savefig(fig_file)
