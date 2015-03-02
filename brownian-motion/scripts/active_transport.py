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
    print '    -D: Input file containing the diffusion coefficient.'
    print '    -f: File to save a plot to.'
    print '    -m: Minimum number of steps in a track.'
    print '    -n: Number of bins in the histograms.'
    print '    -q: Do not display plots, but print messages.'
    print '    -Q: Do not print messages and do not display plots.'
    print '    -T: Temperature of the system.'
    sys.exit(code)

def gaussian(x, a, mean, sigma):
    """
    Gaussian function for fitting.
    """
    return a * np.exp(-(x - mean)**2 / (2.0 * sigma**2))

def get_path_length(track):
    """
    Get the line integral of a particle's trajectory.
    """
    disp_sq = np.array(track[6])
    return np.sum(np.sqrt(disp_sq))

def get_particle_speedsq(track):
    """
    Get the speed that a particle is moving, plus the uncertainty.
    """
    time_step = 0.0330033
    disp_sq = track[6][1:]
    num_steps = len(disp_sq)
    time = np.linspace(1, num_steps, num_steps) * time_step
    popt, pcov = spo.curve_fit(quadratic, time, disp_sq)
    if popt[0] <= 0:
        return 0
    else:
        return popt[0]

def quadratic(x, a, b, c):
    """
    Quadratic function for fitting.
    """
    return a*x**2 + b*x + c

def stokes_work(temperature, diffusion):
    """
    Compute the amount of work required to transport a particle.
    """
    def __sw__(track):
        time_step = 0.0330033
        dispsq = np.array(track[6][1:])
        speed = dispsq / time_step
        work_per_step = 2 * spc.k * temperature * speed / diffusion
        return np.sum(work_per_step)
    return __sw__

# Run the simulation.
if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hd:D:f:m:n:p:qQT:')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Default conditions for the simulation
    data_file = None
    diff_file = None
    fig_file  = None
    pkl_file  = None
    min_steps = None
    nbins = 50
    temp = 293.0
    quiet = False
    extra_quiet = False
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-d':
            data_file = os.path.abspath(arg)
        elif opt == '-D':
            diff_file = os.path.abspath(arg)
        elif opt == '-f':
            fig_file  = os.path.abspath(arg)
        elif opt == '-m':
            min_steps = int(arg)
        elif opt == '-n':
            nbins = int(arg)
        elif opt == '-p':
            pkl_file  = os.path.abspath(arg)
        elif opt == '-q':
            quiet = True
        elif opt == '-Q':
            extra_quiet = True
        elif opt == '-T':
            temp = float(arg)
    if data_file == None:
        print 'No data supplied.'
        usage(1)
    else:
        print 'Reading data from', data_file + '.'

    if diff_file == None:
        print 'No diffusion information supplied.'
        usage(1)
    else:
        print 'Reading diffusion information from', diff_file + '.'

    # Open and read the data file.
    all_tracks = bmc.read_data(data_file)
    # Only accept tracks where there is displacement in each step.
    all_tracks = filter(lambda t: all(t[6][1:]), all_tracks)
    if min_steps != None:
        all_tracks = filter(lambda t: len(t[0])-1 >= min_steps, all_tracks)
    nparticles = len(all_tracks)

    # Read the diffusion and compute the stokes work for each particle.
    with open(diff_file) as pkl:
        diffusion_dict = cPickle.load(pkl)
    diffusion = diffusion_dict['diffusion']
    err_diffusion = diffusion_dict['uncertainty']

    work = np.array(map(stokes_work(temp, diffusion), all_tracks))
    mean_work = np.mean(work)
    err_mean_work = np.std(work, ddof=1) / np.sqrt(len(work))
    max_work  = np.max(work)

    path_lengths = np.array(map(get_path_length, all_tracks))
    longest_path = np.max(path_lengths)
    energy_per_atp = 20.5e3 / spc.N_A

    # Get particles speeds
    particle_speedsq = np.array(map(get_particle_speedsq, all_tracks))
    speedsq = np.array(filter(lambda x: x, particle_speedsq))

    # Compute the work with the mean speed.
    mean_speedsq = np.mean(speedsq)
    vrms = np.sqrt(mean_speedsq)
    err_vrms = np.std(speedsq, ddof=1) / (2 * np.sqrt(len(speedsq)) * vrms)

    # Print some results.
    if not extra_quiet:
        print
        print 'RMS velocity:', vrms, 'm/s'
        print 'Work required to transport particles:', mean_work, 'J'
        print 'Maximum amount of work for a particle:', max_work, 'J'
        print 'Longest path length:', longest_path, 'm'

    if fig_file == None and (quiet or extra_quiet):
        sys.exit()

    # Get the work per path length and fit to a gaussian
    plt.figure()
    work_rate = work / path_lengths
    n, bins, patches = plt.hist(work_rate, nbins, (0, np.max(work_rate)/2.0),
            histtype='step')
    plt.xlabel('Work / path length (J/m)')
    plt.ylabel('Counts')
    plt.title('Average transport rate of the myosin motors')

    # Gaussian fit
    n = np.array(n)
    bins = np.array(bins[:-1])
    bin_width = (bins[1] - bins[0]) / 2.0
    bins += bin_width
    err_n = np.sqrt(n)
    plt.errorbar(bins, n, yerr=err_n, fmt=None)

    hist = zip(bins, n, err_n)
    hist = np.array(filter(lambda t: t[1] > 1, hist)).transpose()
    popt, pcov = spo.curve_fit(gaussian, hist[0], hist[1],
            p0=[1.0, 1e-13, 1e-13], sigma=hist[2])
    fit_xdata = np.linspace(bins[0], bins[-1], 1000)
    fit_ydata = gaussian(fit_xdata, *popt)
    gaus_plot = plt.plot(fit_xdata, fit_ydata, 'k')
    plt.legend((patches[0], gaus_plot[0]), ('Data', 'Gaussian fit'))
    plt.tight_layout()

    mean_work_rate     = popt[1]
    err_mean_work_rate = np.sqrt(pcov[1][1])

    if fig_file:
        filename = fig_file + '_workrate.pdf'
        print
        print 'Saving plot to', filename + '.'
        plt.savefig(filename)

    if not extra_quiet:
        print
        print 'Average work per path length:', mean_work_rate, 'J/m'
        print 'Uncertainty:', err_mean_work_rate

    # Create a dictionary to save to and save it to a pickle.
    if pkl_file != None:
        transport_info = dict()
        transport_info['vel_rms'] = vrms
        transport_info['err_vel_rms'] = err_vrms
        transport_info['mean_work'] = mean_work
        transport_info['err_mean_work'] = err_mean_work
        transport_info['work_rate'] = mean_work_rate
        transport_info['err_work_rate'] = err_mean_work_rate
        print
        print 'Writing results to', pkl_file
        with open(pkl_file, 'w') as pkl:
            cPickle.dump(transport_info, pkl)

    # Plot the speed squared.
    plt.figure()
    speeds = np.sqrt(speedsq)
    n, bins, patches = plt.hist(speeds, nbins, (0, np.max(speeds)/2.0))
    plt.xlabel(r'Speed ($\mu m / s$)')
    plt.ylabel('Counts')
    plt.title('Speed distribution of vesicles.')
    ax = plt.gca()
    xticks = [1e6 * tick for tick in ax.get_xticks()]
    ax.set_xticklabels(map(str, xticks))
    plt.tight_layout()

    if fig_file:
        filename = fig_file + '_speeds.pdf'
        print
        print 'Saving plot to', filename + '.'
        plt.savefig(filename)

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

    if fig_file:
        filename = fig_file + '_tracks.pdf'
        print
        print 'Saving plot to', filename + '.'
        plt.savefig(filename)

    # Plot the square displacements.
    plt.figure()
    for i in range(nparticles):
        plt.plot(all_tracks[i][2][1:], all_tracks[i][7][1:], color=colors[i])

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
        filename = fig_file + '_dispsq.pdf'
        print
        print 'Saving plot to', filename + '.'
        plt.savefig(filename)
