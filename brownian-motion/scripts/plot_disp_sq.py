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
    print 'python plot_disp_sq.py -d data.txt [options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -d: Data file of the particle tracks.'
    print '    -f: File to save a plot to.'
    print '    -x: Maximum x value.'
    print '    -y: Maximum y value.'
    sys.exit(code)

# Run the simulation.
if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hd:f:p:qQx:y:')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Default conditions for the simulation
    data_file = None
    fig_file = None
    pkl_file = None
    quiet = False
    extra_quiet = False
    max_x = None
    max_y = None
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
        elif opt == '-x':
            max_x = float(arg)
        elif opt == '-y':
            max_y = float(arg)
    if data_file == None:
        print 'No data supplied. Cannot compute diffusion coefficient.'
        usage(1)
    else:
        print 'Reading data from', data_file + '.'
    if max_x != None and max_y == None:
        print 'Cannot set coordinates without a max y value.'
        usage(1)

    # Open and read the data file.
    all_tracks = bmc.read_data(data_file)
    # Only accept tracks where there is displacement in each step.
    all_tracks = filter(lambda t: all(t[6][1:]), all_tracks)
    nparticles = len(all_tracks)

    # Create plots.
    colors = map(tuple, npr.random((nparticles, 3)))

    # Plot the square displacements.
    plt.figure()
    for i in range(nparticles):
        plt.plot(all_tracks[i][2][1:], all_tracks[i][7][1:], color=colors[i])

    if max_x != None and max_y != None:
        plt.axis([0, max_x, 0, max_y])
    elif max_y != None:
        time  = map(lambda l: l[2], all_tracks)
        max_x = max(map(max, time))
        plt.axis([0, max_x, 0, max_y])

    ax = plt.gca()
    yticks = [1e12 * tick for tick in ax.get_yticks()]
    ax.set_yticklabels(map(str, yticks))
    plt.xlabel('Time (s)')
    plt.ylabel(r'Displacement Squared ($\mu m^2$)')
    plt.title('Combined Particle Displacements')
    plt.tight_layout()

    if fig_file == None:
        plt.show()
    else:
        print
        print 'Saving plot to', fig_file + '.'
        plt.savefig(fig_file)
