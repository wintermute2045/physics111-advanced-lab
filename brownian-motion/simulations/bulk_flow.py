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
    print 'python bulk_flow.py [options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -d: Number of dimensions in the simulation.'
    print '    -D: Diameter of each particle.'
    print '    -n: Number of steps in the random walk.'
    print '    -q: Do not display plots.'
    print '    -t: Time difference of each step.'
    print '    -T: Temperature of the system.'
    print '    -v: Viscosity of the solution.'
    print '    -x: Flow in the x direction.'
    print '    -y: Flow in the y direction.'
    print '    -z: Flow in the z direction.'
    sys.exit(code)

# Run the simulation.
if __name__ == '__main__':
    # Parse arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'd:D:hn:qt:T:v:x:y:z:')
    except getopt.GetoptError as err:
        print str(err)
        usage(1)

    # Default conditions for the simulation
    ndim = 2
    nsteps = 1000
    radius = 1e-6
    time_step = 0.1
    temperature = 293.0
    viscosity = 1e-3
    flow_coefs = [0.2, 0.05, 0.0]
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

    # Get a track.
    small_sphere = Particle(ndim, radius, viscosity, temperature)
    small_sphere.add_bulk_flow(flow_coefs[:ndim], time_step, nsteps)

    # The flow is computed in get_diffusion when correcting for it.
    # No need to call get_bulk_flow.
    small_sphere.get_diffusion(flow=True)

    # Display some results:
    print 'Number of dimensions:', ndim
    print 'Number of steps:', nsteps
    print 'Radius of particles:', radius, 'um'
    print 'Viscosity of solution:', viscosity, 'Ps'
    print 'Temperature of solution:', temperature, 'K'
    print 'Time step:', time_step, 's'
    print
    print 'Theoretical diffusion:', small_sphere.diffusion.theory
    print 'Diffusion (flow-corrected):', small_sphere.diffusion.flow_corr
    print 'Diffusion (uncorrected):', small_sphere.diffusion.data
    print 'Standard error:', small_sphere.diffusion.std_err
    print 'Absolute error:', small_sphere.diffusion.abs_err
    print 'Absolute error from corrected diffusion:',
    print small_sphere.diffusion.theory - small_sphere.diffusion.flow_corr
    print
    print 'Measured speed of flow:', np.sqrt(small_sphere.flow.vel_sq), 'm/s'
    print '| Exact flow (m) |',
    print 'Measured flow (m) |',
    print ' Error (m)  |',
    print 'Velocity (m/s) |'
    for i in range(ndim):
        coef = small_sphere.flow.actual[i]
        flow = small_sphere.flow.data[i]
        err_flow = small_sphere.flow.error[i]
        velocity = small_sphere.flow.velocity[i]
        params = (coef, flow, err_flow, velocity)
        print '|   %#6g  |    %#6g    | %#6g |   %#6g  |' % params

    if quiet:
        sys.exit()

    # Plot the random walks
    fig = plt.figure()
    # Enable 3D plotting.
    if ndim == 3:
        plt3d = fig.gca(projection='3d')
        plt3d.plot(*small_sphere.position)
    else:
        plt.plot(*small_sphere.position)
        ax = plt.gca()

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
    plt.title('Single Particle Track')
    plt.tight_layout()

    # Plot displacements
    plt.figure()
    plt.plot(small_sphere.time, small_sphere.get_pos_sq())

    # Compare observed to the theoretical predictions.
    vel_sq = small_sphere.flow.vel_sq_th
    scale_factor = 2 * small_sphere.length_scale**2 / small_sphere.time_step
    disp_theory  = scale_factor * small_sphere.time
    disp_theory += vel_sq * small_sphere.time**2

    plt.plot(small_sphere.time, disp_theory)
    ax = plt.gca()
    plt.xlabel('Time (s)')
    plt.ylabel(r'Displacement Squared ($\mu m^2$)')
    plt.title('Single Particle Displacement Under Bulk Flow')
    yticks = [1e12 * tick for tick in ax.get_yticks()]
    ax.set_yticklabels(map(str, yticks))
    plt.tight_layout()

    plt.show()
