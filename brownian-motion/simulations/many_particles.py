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
import scipy.optimize as spo
import scipy.constants as spc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from Particle import *

def usage(code):
    """
    Display the allowed options for the program.
    """
    print 'python many_particles.py [options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -d: Number of dimensions in the simulation.'
    print '    -D: Diameter of each particle.'
    print '    -n: Number of steps in the random walk.'
    print '    -N: Number of particles in the simulation.'
    print '    -q: Do not display plots.'
    print '    -t: Time difference of each step.'
    print '    -T: Temperature of the system.'
    print '    -v: Viscosity of the solution.'
    sys.exit(code)

def average_displacement(displacements):
    """
    Get the average of each particle's displacement at each step.
    """
    average = np.array(map(lambda l: mean_with_err(l), zip(*displacements)))
    return average.transpose()

# Run the simulation.
if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'd:D:hn:N:qt:T:v:')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Default conditions for the simulation
    ndim = 2
    nsteps = 50
    nparticles = 10
    radius = 1e-6
    time_step = 0.1
    temperature = 293.0
    viscosity = 1e-3
    quiet = False
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-d':
            ndim = int(arg)
        elif opt == '-D':
            radius = float(arg) / 2
        elif opt == '-n':
            nsteps = int(arg)
        elif opt == '-N':
            nparticles = int(arg)
        elif opt == '-q':
            quiet = True
        elif opt == '-t':
            time_step = float(arg)
        elif opt == '-T':
            temperature = float(arg)
        elif opt == '-v':
            viscosity = float(arg)
    if ndim < 1 or ndim > 3:
        print 'Invalid number of dimensions.'
        sys.exit(1)

    # Run the simulation on a list of particles.
    particles = []
    for i in range(nparticles):
        # Create a particle.
        particles.append(Particle(ndim, radius, viscosity, temperature))
        particles[-1].generate_motion(nsteps, time_step)
        particles[-1].get_diffusion()

    # Get the mean diffusion value from the simulation
    diffusion     = np.array([p.diffusion.data for p in particles])
    err_diffusion = np.array([p.diffusion.std_err for p in particles])
    mean_diff, err_mean_diff = mean_with_err(diffusion)

    err_theory  =  particles[0].diffusion.theory * np.sqrt(2)
    err_theory /= np.sqrt(ndim*nsteps*nparticles)

    # Curve fitting
    time_steps = particles[0].time
    displacement = [p.get_pos_sq() for p in particles]
    ave_displacements = average_displacement(displacement)
    scale_factor = 2 * particles[0].length_scale**2 / particles[0].time_step
    theory_data = time_steps * scale_factor
    #popt, pcov = spo.curve_fit(linear, time_steps, ave_displacements[0],
    #        sigma=ave_displacements[1])

    print 'Number of particles:', nparticles
    print 'Number of dimensions:', ndim
    print 'Number of steps:', nsteps
    print 'Radius of particles:', radius, 'um'
    print 'Viscosity of solution:', viscosity, 'Ps'
    print 'Temperature of solution:', temperature, 'K'
    print 'Time step:', time_step, 's'
    print
    print 'Diffusion:', mean_diff
    print 'Uncertainty:', err_mean_diff
    print 'Theoretical diffusion:', particles[0].diffusion.theory
    print 'Theoretical error:', err_theory

    if quiet:
        sys.exit()

    # Plot colors
    colors = map(tuple, npr.random((nparticles, 3)))

    # Plot the random walks
    fig = plt.figure()
    # Enable 3D plotting.
    if ndim == 3:
        plt3d = fig.gca(projection='3d')
    else:
        ax = plt.gca()

    for i in range(nparticles):
        if ndim == 3:
            plt3d.plot(*particles[i].position, color=colors[i])
        else:
            plt.plot(*particles[i].position, color=colors[i])
    # Format the titles
    if ndim == 1:
        plt.xlabel('Step number')
        plt.ylabel(r'X position ($\mu m$)')
        yticks = [1e6 * tick for tick in ax.get_yticks()]
        ax.set_yticklabels(map(str, yticks))
    elif ndim == 2:
        plt.xlabel(r'X position ($\mu m$)')
        plt.ylabel(r'Y position ($\mu m$)')
        xticks = [1e6 * tick for tick in ax.get_xticks()]
        yticks = [1e6 * tick for tick in ax.get_yticks()]
        ax.set_xticklabels(map(str, xticks))
        ax.set_yticklabels(map(str, yticks))
    if ndim == 3:
        plt3d.set_xlabel(r'X position ($\mu m$)')
        plt3d.set_ylabel(r'Y position ($\mu m$)')
        plt3d.set_zlabel(r'Z position ($\mu m$)')
        xticks = [1e6 * tick for tick in plt3d.get_xticks()]
        yticks = [1e6 * tick for tick in plt3d.get_yticks()]
        zticks = [1e6 * tick for tick in plt3d.get_zticks()]
        plt3d.set_xticklabels(map(str, xticks))
        plt3d.set_yticklabels(map(str, yticks))
        plt3d.set_zticklabels(map(str, zticks))
    plt.title('Combined Particle Tracks')
    plt.tight_layout()

    # Plot the displacements.
    plt.figure()
    ax = plt.gca()
    for i in range(nparticles):
        plt.plot(time_steps, displacement[i], color=colors[i])
    plt.plot(time_steps, theory_data,
            color='#555753', linewidth=4, label='Theoretical')
    plt.plot(time_steps, ave_displacements[0],
            color='#000000', linewidth=4, label='Average')
    plt.legend(loc='upper left')
    plt.xlabel('Time (s)')
    plt.ylabel(r'Displacement Squared ($\mu m^2$)')
    plt.title('Combined Particle Displacements')
    yticks = [1e12 * tick for tick in ax.get_yticks()]
    ax.set_yticklabels(map(str, yticks))
    plt.tight_layout()

    # Plot the diffusion calculations
    plt.figure()
    ax = plt.gca()
    particle_num  = np.arange(1, nparticles+1)
    plt.plot(particle_num, diffusion, 'ro')
    plt.errorbar(particle_num, diffusion, err_diffusion, ecolor='r', fmt=None)

    # Make the mean into a repetition of nparticle items
    mean_diff = np.zeros(nparticles) + mean_diff
    err_diff_up   = mean_diff + err_mean_diff
    err_diff_down = mean_diff - err_mean_diff

    # Add the mean and its error to the plot.
    plt.plot(particle_num, mean_diff, 'b',
            linewidth=4, label='Average value of D')
    plt.plot(particle_num, err_diff_up,   'g')
    plt.plot(particle_num, err_diff_down, 'g')
    plt.legend(loc='upper left')
    plt.xlabel('Particle')
    plt.ylabel(r'Diffusion ($\mu m^2/s$)')
    plt.title('Diffusion coefficients for many particles')
    yticks = [1e12 * tick for tick in ax.get_yticks()]
    ax.set_yticklabels(map(str, yticks))
    plt.tight_layout()

    plt.show()
