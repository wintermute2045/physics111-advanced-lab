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
import numpy as np
import matplotlib.pyplot as plt
import brownian_motion_cells as bmc

def usage(code):
    """
    Display the allowed options for the program.
    """
    print 'python num_tracks.py -d data.txt [options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -d: Data file of the particle tracks.'
    print '    -q: Do not display plots.'
    sys.exit(code)

# Run the simulation.
if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hd:')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Default conditions for the simulation
    data_file = None
    quiet = False
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-d':
            data_file = os.path.abspath(arg)
        elif opt == '-q':
            quiet = True
    if data_file == None:
        print 'No data supplied.'
        sys.exit(1)

    # Open and read the data file.
    all_tracks = bmc.read_data(data_file)

    lengths = map(lambda l: len(l[0])-1, all_tracks)
    max_len = max(lengths)

    print 'Total number of tracks:', len(lengths)

    n, bins, patches = plt.hist(lengths, max_len, (0, max_len))

    plt.show()
