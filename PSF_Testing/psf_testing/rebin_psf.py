""" @file rebin_psf.py

    Created 21 Sep 2015

    Function to rebin a subsampled PSF

    ---------------------------------------------------------------------

    Copyright (C) 2015 Bryan R. Gillis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import division

import numpy as np
from scipy.signal import convolve

from psf_testing import magic_values as mv

def rebin(subsampled_image,
          cdf,
          x_shift=0,
          y_shift=0,
          subsampling_factor = mv.default_subsampling_factor):
    
    # Get the size of the final image
    ss_nx, ss_ny = np.shape(subsampled_image)
    
    rb_nx = 2 * (((ss_nx-np.abs(x_shift))//subsampling_factor - 1) // 2) + 1 
    rb_ny = 2 * (((ss_ny-np.abs(y_shift))//subsampling_factor - 1) // 2) + 1 
    
    x_offset = (ss_nx - rb_nx*subsampling_factor)//2 - x_shift
    y_offset = (ss_ny - rb_ny*subsampling_factor)//2 - y_shift
    
    # Shift is in star - model
    
    # Make the rebinned array
    rebinned_array = np.zeros((rb_nx,rb_ny))
    for xi in xrange(rb_nx):
        for yi in xrange(rb_ny):
            xm = x_offset + subsampling_factor*(xi+1) - int(subsampling_factor/2)
            ym = y_offset + subsampling_factor*(yi+1) - int(subsampling_factor/2)
            rebinned_array[xi,yi] = subsampled_image[xm:xm+subsampling_factor,
                                                     ym:ym+subsampling_factor ].sum()
    
    # Convolve it with the charge diffusion kernel
    rebinned_diffused_array = convolve(rebinned_array, cdf, mode="same")
    
    return rebinned_diffused_array