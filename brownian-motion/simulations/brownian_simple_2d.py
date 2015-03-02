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

import sys
import getopt
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

def usage(code):
    """
    Display the allowed options for the program.
    """
    print 'python brownian.py [options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -n: Number of steps in the random walk.'
    sys.exit(code)

class Position():
    """
    Generic struct-like class for the positions of particles.
    """
    def __init__(self, nsteps):
        self.xpts   = np.cumsum(npr.normal(0, 1, nsteps))
        self.ypts   = np.cumsum(npr.normal(0, 1, nsteps))
        self.distsq = self.xpts**2 + self.ypts**2

if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hn:')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Default the number of steps.
    nsteps = 1000
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-n':
            nsteps = int(arg)

    # Get a track.
    position = Position(nsteps)

    # Plot of position
    plt.figure()
    plt.plot(position.xpts, position.ypts)
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    plt.title('position versus time in 2D')

    # Plot this distances.
    plt.figure()
    plt.plot(position.distsq)

    plt.show()
