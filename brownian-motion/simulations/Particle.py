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
import numpy.random as _npr
import scipy.constants as _spc

# Generic data class to store variables.
class _Struct:
    pass

class EmptyVariable(Exception):
    pass

class Particle:
    def __init__(self, ndim, radius, viscosity, temperature):
        # Get the diffusion coefficients and the time
        self.ndim = ndim
        diffusion  = _spc.k * temperature / (3 * _spc.pi * viscosity * radius)
        self.diffusion = _Struct()
        self.diffusion.theory = diffusion

        # Create empty variables for data that will be computed later
        self.flow = None
        self.time = None
        self.position = None
        self.time_step = None
        self.displacement = None
        self.length_scale = None
        self.diffusion.data = None
        self.diffusion.std_err = None
        self.diffusion.abs_err = None

    def add_bulk_flow(self, flow, time_step, nsteps=None):
        """
        Add a bulk flow or create a trajectory with a bulk flow.
        """
        if len(flow) != self.ndim:
            raise ValueError('Flow parameters do not match the dimensions.')

        if self.displacement == None:
            # Check for errors.
            if nsteps == None:
                raise ValueError('Must iterate if displacement list is empty.')

            # Compute displacements
            self.generate_motion(nsteps, time_step)

        # Add flow to the displacement
        flow = _np.array([self.length_scale * f for f in flow])
        velocity = flow / time_step
        vel_sq = _np.sum(velocity**2)

        self.flow = _Struct()
        self.flow.actual = flow
        self.flow.vel_sq_th = vel_sq
        self.displacement = map(lambda x,y: x+y, self.displacement, flow)
        self.position = map(_np.cumsum, self.displacement)

    def generate_motion(self, nsteps, time_step):
        """
        Generate a random walk for a given number of steps.
        """
        # Theoretical mean displacement of particles. This sets the scale
        length_scale = _np.sqrt(2.0 * self.diffusion.theory * time_step)

        self.length_scale = length_scale
        self.time_step = time_step

        # Compute the displacements and positions of the particle.
        self.time = time_step * _np.r_[1:nsteps+1]
        self.displacement = [_npr.normal(0, length_scale, nsteps)
                for i in range(self.ndim)]
        self.position = map(_np.cumsum, self.displacement)

    def get_bulk_flow(self):
        """
        Estimate the flow coefficients from the data.
        """
        # This will only work after the generate_motion command has been run...
        if self.displacement == None:
            raise EmptyVariable('No displacements available!')

        nsteps = self.time.shape[0]
        time_step = self.time_step

        flow = _np.array([_np.mean(d) for d in self.displacement])
        err_flow = _np.array([_np.std(d, ddof=1) for d in self.displacement])
        err_flow /= _np.sqrt(nsteps)

        velocity = flow / time_step
        err_velocity = err_flow / time_step

        vel_sq = _np.sum(velocity**2)
        err_vel_sq = 2 * _np.sqrt(_np.sum(velocity**2 * err_velocity**2))

        self.flow.data = flow
        self.flow.error = err_flow
        self.flow.velocity = velocity
        self.flow.err_vel = err_velocity
        self.flow.vel_sq = vel_sq
        self.flow.err_vel_sq = err_vel_sq

    def get_diffusion(self, flow=False):
        """
        Get the diffusion from the sampled distances.
        """
        # This will only work after the generate_motion command has been run...
        if self.displacement == None:
            raise EmptyVariable('No displacements available!')

        scale = 2.0 * self.ndim * self.time_step
        disp_sq = self.get_disp_sq()
        diffusion = _np.mean(disp_sq) / scale
        self.diffusion.data = diffusion

        # Get the error on the (uncorrected) diffusion
        std_err = _np.std(disp_sq, ddof=1) / scale / _np.sqrt(disp_sq.shape[0])
        abs_err = self.diffusion.theory - diffusion
        self.diffusion.std_err = std_err
        self.diffusion.abs_err = abs_err

        # Testing the corrected diffusion formula I derived.
        if flow:
            self.get_bulk_flow()
            vel_sq = self.flow.vel_sq
            offset = self.time_step * vel_sq / (2 * self.ndim)
            flow_corr = diffusion - offset
            self.diffusion.flow_corr = flow_corr

            # Compute the error. Only use this when dealing with one track.
            ndim = float(self.ndim)
            time_step = self.time_step
            err_vel_sq = self.flow.err_vel_sq
            err_flow_corr  = std_err ** 2
            err_flow_corr += (time_step * err_vel_sq / (2 * ndim)) ** 2
            err_flow_corr  = _np.sqrt(err_flow_corr)
            self.diffusion.err_flow_corr = err_flow_corr

    def get_disp_sq(self):
        """
        Generate a list of the squares of each displacment in the
        random walk.
        """
        # This will only work after the generate_motion command has been run...
        if self.displacement == None:
            raise EmptyVariable('No displacements available!')

        zero_list = _np.zeros(self.displacement[0].shape)
        return reduce(lambda x,y: x + y**2, self.displacement, zero_list)

    def get_pos_sq(self):
        """
        Generate a list of the squares of the displacements from the
        origin.
        """
        # This will only work after the generate_motion command has been run...
        if self.position == None:
            raise EmptyVariable('No positions available!')

        zero_list = _np.zeros(self.position[0].shape)
        return reduce(lambda x,y: x + y**2, self.position, zero_list)

def linear(x, m, b):
    """
    Function for linear fitting.
    """
    return m*x + b

# Averaging function.
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
