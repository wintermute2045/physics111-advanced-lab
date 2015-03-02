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
import numpy.random as npr
import scipy.optimize as spo
import matplotlib.pyplot as plt
import brownian_motion_cells as bmc

def usage(code):
    """
    Display the allowed options for the program.
    """
    print 'python cell_bulk_flow.py -d data.txt [options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -d: Data file of the particle tracks.'
    print '    -f: File to save a plot to.'
    print '    -p: Pickle file to save to.'
    print '    -q: Do not display plots, but print messages.'
    print '    -Q: Do not print messages and do not display plots.'
    sys.exit(code)

def quadratic(x, a, b, c):
    """
    Quadratic function for fitting.
    """
    return a*x**2 + b*x + c

# Run the simulation.
if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hd:f:p:qQ')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Default conditions for the simulation
    data_file = None
    fig_file = None
    pkl_file = None
    quiet = False
    extra_quiet = False
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-d':
            data_file = os.path.abspath(arg)
        elif opt == '-f':
            fig_file = os.path.abspath(arg)
        elif opt == '-p':
            pkl_file = os.path.abspath(arg)
        elif opt == '-q':
            quiet = True
        elif opt == '-Q':
            extra_quiet = True
    if data_file == None:
        print 'No data supplied. Cannot compute flow.'
        sys.exit(1)
    else:
        print 'Reading data from', data_file + '.'

    # Open and read the data file.
    all_tracks = bmc.read_data(data_file)
    # Only accept tracks where there is displacement in each step.
    all_tracks = filter(lambda t: all(t[6][1:]), all_tracks)
    nparticles = len(all_tracks)

    # Compute the diffusions.
    time_step = 0.0330033 # The time step isn't recorded uniformly.

    # Create average displacements and fit a curve
    ave_dispsq = bmc.average_dispsq(all_tracks, 50)
    num_avg = len(ave_dispsq[0])
    time = np.linspace(1, num_avg, num_avg) * time_step
    popt, pcov = spo.curve_fit(quadratic, time, ave_dispsq[0],
            sigma=ave_dispsq[1])
    fit_data = popt[0]*time**2 + popt[1]*time + popt[2]
    axis_bounds  = [min(time), max(time)]
    axis_bounds += [0, max(max(fit_data), max(ave_dispsq[0]))]

    vel_sq = popt[0]
    err_vel_sq = np.sqrt(pcov[0][0])
    speed = np.sqrt(vel_sq)
    err_speed = err_vel_sq / (2 * speed)

    # Create a dictionary to save to and save it to a pickle.
    if pkl_file != None:
        bulk_flow_info = dict()
        bulk_flow_info['speed_sq'] = vel_sq
        bulk_flow_info['err_speed_sq'] = err_vel_sq
        bulk_flow_info['speed'] = speed
        bulk_flow_info['err_speed'] = err_speed
        print
        print 'Writing results to', pkl_file
        with open(pkl_file, 'w') as pkl:
            cPickle.dump(bulk_flow_info, pkl)

    # Print some results.
    if not extra_quiet:
        print
        print 'Number of particles:', nparticles
        print
        print 'Speed of the flow:', speed, 'm/s'
        print 'Uncertainty:', err_speed, 'm/s'
        print
        print 'Speed squared:', vel_sq, '(m/s)^2'
        print 'Uncertainty:', err_vel_sq, '(m/s)^2'

    if fig_file == None and (quiet or extra_quiet):
        sys.exit()

    # Create plots.
    colors = map(tuple, npr.random((nparticles, 3)))

    # Plot the square displacements.
    plt.figure()
    for i in range(nparticles):
        plt.plot(all_tracks[i][2][1:], all_tracks[i][7][1:], color=colors[i])

    plt.plot(time, ave_dispsq[0], '#555753', linewidth=4, label='Average')
    plt.plot(time, fit_data, 'k', linewidth=4, label='Quadratic fit')
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
