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

import numpy as _np
import itertools as _it
import scipy.constants as _spc

def average_dispsq(all_tracks, min_tracks=10):
    """
    Create the average dispacement squared for all tracks.
    """
    disp_sq = map(lambda l: l[7][1:], all_tracks)
    longest_track = max(map(len, disp_sq))

    # Add a list of Nones to the end of the length of each track so that zip
    # will not truncate any track.
    for i in range(len(disp_sq)):
        track_len   = len(disp_sq[i])
        disp_sq[i]  = list(disp_sq[i])
        disp_sq[i] += [None for i in range(longest_track - track_len)]

    # Create a list of displacement squares for each time step.
    disp_list = map(list, map(lambda l: filter(lambda x: x, l), zip(*disp_sq)))
    disp_list = list(_it.takewhile(lambda x: len(x) >= min_tracks, disp_list))

    ave_list = _np.array(map(lambda l: mean_with_err(l), disp_list))
    return ave_list.transpose()

def mean_with_err(data, err=None):
    """
    Computes the mean, plus error on the mean for data with errors.
    """
    data = _np.array(data)
    if err == None:
        mean = _np.mean(data)
        err_mean = _np.std(data, ddof=1) / _np.sqrt(len(data))
    else:
        err  = _np.array(err)
        weights = 1.0 / err**2
        mean, sumweights = _np.average(data, weights=weights, returned=True)
        err_mean = 1 / _np.sqrt(sumweights)

    return(mean, err_mean)

def read_data(filename):
    """
    Read data from a file produced by the Brownian motion software.

    Format of the input file:
    x y time dx dy dt dr^2 TotalDisplacementSquared
    """
    # Open the data file.
    try:
        infile = open(filename, 'r')
    except IOError:
        raise IOError('Input file ' + filename + ' does not exist!')

    # The first line just contains a header. It can be thrown away.
    infile.readline()

    # Main loop that records the data.
    all_tracks = []
    num_steps = infile.readline()
    while num_steps != '':
        # DOS line endings were a stupid idea.
        # What the fuck was Microsoft even thinking?
        if num_steps[-2:] == '\r\n':
            num_steps = int(num_steps[:-2])
        else:
            num_steps = int(num_steps[:-1])

        # Collect steps in the particle track.
        particle = []
        for i in range(num_steps):
            data_line = infile.readline()

            # Detect DOS or UNIX line endings...
            if data_line[-2:] == '\r\n':
                data_line = data_line[:-2]
            else:
                data_line = data_line[:-1]

            particle.append(map(float, data_line.split()))

        # Save data and then get the next line of steps from the file.
        all_tracks.append(_np.array(particle).transpose())
        num_steps = infile.readline()

    # Close file and return data.
    infile.close()
    return all_tracks

def theoretical(radius, eta, temp):
    """
    Theoretical diffusion value.
    """
    diffusion = _spc.k * temp / (3 * _spc.pi * eta * radius)
    return diffusion
