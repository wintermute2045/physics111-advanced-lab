#!/usr/bin/env python2.7

################################################################################
## This code is part of the analysis for the optical pumping lab
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

import plot85
import plot87
import scipy as sp

if __name__ == '__main__':
    sp85, msp85, field85, errfield85, pfreq_err85, nfreq_err85 = plot85.main()
    sp87, msp87, field87, errfield87, pfreq_err87, nfreq_err87 = plot87.main()

    # Combine the field measurements
    error  = 1.0/errfield85**2 + 1.0/errfield87**2
    error  = 1.0/sp.sqrt(error)
    field  = field85/errfield85**2 + field87/errfield87**2
    field *= error**2
    field  = round(field, 3)
    error  = round(error, 3)

    print 'Nuclear spin of ^{85}Rb:', sp85, msp85
    print 'Nuclear spin of ^{87}Rb:', sp87, msp87
    print 'Earth\'s magnetic field:', field, 'Gauss \pm', error, 'Gauss'

    print
    print 'Uncertainties of RF for ^{85}Rb:', pfreq_err85, nfreq_err85
    print 'Uncertainties of RF for ^{87}Rb:', pfreq_err87, nfreq_err87
