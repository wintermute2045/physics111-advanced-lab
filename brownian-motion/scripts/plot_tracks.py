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
import scipy.constants as spc
import matplotlib.pyplot as plt
import brownian_motion_cells as bmc

def usage(code):
    """
    Display the allowed options for the program.
    """
    print 'python active_transport.py -d data.txt -D diffusion.pkl [options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -d: Data file of the particle tracks.'
    print '    -f: File to save a plot to.'
    print '    -m: Minimum number of steps in a track.'
    sys.exit(code)

# Run the simulation.
if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hd:f:m:')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Default conditions for the simulation
    data_file = None
    fig_file  = None
    min_steps = None
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-d':
            data_file = os.path.abspath(arg)
        elif opt == '-f':
            fig_file  = os.path.abspath(arg)
        elif opt == '-m':
            min_steps = int(arg)
    if data_file == None:
        print 'No data supplied.'
        usage(1)
    else:
        print 'Reading data from', data_file + '.'

    # Open and read the data file.
    all_tracks = bmc.read_data(data_file)
    # Only accept tracks where there is displacement in each step.
    all_tracks = filter(lambda t: all(t[6][1:]), all_tracks)
    if min_steps != None:
        all_tracks = filter(lambda t: len(t[0])-1 >= min_steps, all_tracks)
    nparticles = len(all_tracks)

    colors = map(tuple, npr.random((nparticles, 3)))

    # Plot the particle trajectories.
    plt.figure()
    for i in range(nparticles):
        plt.plot(*all_tracks[i][:2], color=colors[i])

    ax = plt.gca()
    xticks = [1e6 * tick for tick in ax.get_xticks()]
    yticks = [1e6 * tick for tick in ax.get_yticks()]
    ax.set_xticklabels(map(str, xticks))
    ax.set_yticklabels(map(str, yticks))
    plt.xlabel(r'X position ($\mu m$)')
    plt.ylabel(r'Y position ($\mu m$)')
    plt.title('Combined Particle Tracks')
    plt.tight_layout()

    if fig_file == None:
        plt.show()
    else:
        print
        print 'Saving plot to', fig_file + '.'
        plt.savefig(fig_file)
