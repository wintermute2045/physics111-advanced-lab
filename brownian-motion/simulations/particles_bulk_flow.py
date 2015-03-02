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
import scipy.constants as spc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from Particle import *

def usage(code):
    """
    Display the allowed options for the program.
    """
    print 'python particles_bulk_flow.py [options]'
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
    print '    -x: Flow in the x direction.'
    print '    -y: Flow in the y direction.'
    print '    -z: Flow in the z direction.'
    sys.exit(code)

def average_displacement(displacements):
    """
    Get the average of each particle's displacement at each step.
    """
    average  = np.zeros(displacements[0].shape)
    average  = reduce(lambda x,y: x+y, displacements, average)
    average /= len(displacements)
    return average

# Run the simulation.
if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'd:D:hn:N:qt:T:v:x:y:z:')
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
    flow_coefs = [0.0, 0.0, 0.0]
    #flow_coefs = [0.2, 0.05, 0.0]
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
        elif opt == '-x':
            flow_coefs[0] = float(arg)
        elif opt == '-y':
            flow_coefs[1] = float(arg)
        elif opt == '-z':
            flow_coefs[2] = float(arg)
    if ndim < 1 or ndim > 3:
        print 'Invalid number of dimensions.'
        sys.exit(1)

    # Run the simulation on a list of particles.
    particles = []
    for i in range(nparticles):
        # Create a particle.
        particles.append(Particle(ndim, radius, viscosity, temperature))
        particles[-1].generate_motion(nsteps, time_step)
        particles[-1].add_bulk_flow(flow_coefs[:ndim], time_step, nsteps)
        particles[-1].get_diffusion(flow=True)

    # Get the mean diffusion value from the simulation
    diffusion     = np.array([p.diffusion.data for p in particles])
    err_diffusion = np.array([p.diffusion.std_err for p in particles])
    vel_sq = np.array([p.flow.vel_sq for p in particles])
    err_vel_sq = np.array([p.flow.err_vel_sq for p in particles])

    mean_diff, err_mean_diff = mean_with_err(diffusion)
    mean_vel_sq, err_vel_sq = mean_with_err(vel_sq)

    # Get the average flow components.
    mean_flow = []
    err_flow  = []
    mean_vel  = []
    for i in range(ndim):
        flow = map(lambda p: p.flow.data[i], particles)
        velocity = map(lambda p: p.flow.velocity[i], particles)
        mean_flow.append(np.mean(flow))
        err_flow.append(np.std(flow, ddof=1) / np.sqrt(len(flow)))
        mean_vel.append(np.mean(velocity))

    diff_corr = mean_diff - time_step * mean_vel_sq / (2.0 * ndim)
    err_diff_corr  = err_mean_diff ** 2
    err_diff_corr += (time_step * err_vel_sq / (2.0 * ndim)) ** 2
    err_diff_corr  = np.sqrt(err_diff_corr)

    print 'Number of particles:', nparticles
    print 'Number of dimensions:', ndim
    print 'Number of steps:', nsteps
    print 'Radius of particles:', radius, 'um'
    print 'Viscosity of solution:', viscosity, 'Ps'
    print 'Temperature of solution:', temperature, 'K'
    print 'Time step:', time_step, 's'
    print
    print 'Uncorrected diffusion:', mean_diff
    print 'Uncertainty:', err_mean_diff
    print 'Corrected diffusion:', diff_corr
    print 'Uncertainty:', err_diff_corr
    print 'Theoretical diffusion:', particles[0].diffusion.theory
    #print
    #print 'WARNING: Corrected diffusion might not be',
    #print '(probably isn\'t) accurate.'

    print
    print 'Measured speed of flow:', np.sqrt(particles[0].flow.vel_sq), 'm/s'
    print '| Exact flow (m) |',
    print 'Measured flow (m) |',
    print ' Error (m)  |',
    print 'Velocity (m/s) |'
    for i in range(ndim):
        coef = particles[0].flow.actual[i]
        flow = mean_flow[i]
        error_flow = err_flow[i]
        velocity = mean_vel[i]
        params = (coef, flow, error_flow, velocity)
        print '|   %#6g  |    %#6g    | %#6g |   %#6g  |' % params

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
    time_steps = particles[0].time
    displacement = [p.get_pos_sq() for p in particles]
    ave_displacements = average_displacement(displacement)

    vel_sq = particles[0].flow.vel_sq_th
    scale_factor = 2 * particles[0].length_scale**2 / particles[0].time_step
    disp_theory = scale_factor * time_steps + vel_sq * time_steps**2

    plt.figure()
    ax = plt.gca()
    for i in range(nparticles):
        plt.plot(time_steps, displacement[i], color=colors[i])
    plt.plot(time_steps, disp_theory,
            color='#555753', linewidth=4, label='Theoretical')
    plt.plot(time_steps, ave_displacements,
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
