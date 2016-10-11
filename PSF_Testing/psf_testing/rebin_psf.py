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

from copy import deepcopy
import numpy as np
from scipy.signal import convolve
import galsim

from psf_testing import magic_values as mv

try:
    import cIceBRGpy
except ImportError as _e:
    from Release import cIceBRGpy
    
def galsim_rebin(a,
                 x_shift = 0,
                 y_shift = 0,
                 subsampling_factor = mv.default_subsampling_factor,
                 guiding_error_mag1 = 0.,
                 guiding_error_mag2 = 0.,
                 guiding_error_angle = 0.,):
    
    # Interpret the subsampled image as an interpolated image
    ss_prof = galsim.InterpolatedImage(galsim.Image(a,scale=1./subsampling_factor))
    
    # Convolve by guiding error if appropriate
    if guiding_error_mag1 != 0. or guiding_error_mag2:
        guiding_error_prof = galsim.Box(guiding_error_mag1,guiding_error_mag2)
        guiding_error_prof = guiding_error_prof.rotate(guiding_error_angle*galsim.degrees)
        
        ss_prof = galsim.Convolve([ss_prof,guiding_error_prof])
    
    # Get the rebinned shape
    ss_nx, ss_ny = np.shape(a)
    rb_nx = 2*((ss_nx//subsampling_factor-1)//2) + 1
    rb_ny = 2*((ss_ny//subsampling_factor-1)//2) + 1
    
    rb_array = np.zeros((rb_nx,rb_ny))
    rb_image = galsim.Image(rb_array,scale=1.)
    
    ss_prof.drawImage(rb_image,offset=(float(y_shift)/subsampling_factor,float(x_shift)/subsampling_factor))
    
    return rb_array

def rebin(a,
          cdf,
          x_shift=0,
          y_shift=0,
          subsampling_factor = mv.default_subsampling_factor,
          conserve=False,
          use_galsim=True,
          guiding_error_mag1=0.,
          guiding_error_mag2=0.,
          guiding_error_angle=0.,):
    
    if (guiding_error_mag1 != 0. or guiding_error_mag2 != 0.) and not use_galsim:
        raise Exception("Error: Galsim must be used for rebinning if guiding error is nonzero.")

    if use_galsim:
        rebinned_array = galsim_rebin(a,x_shift,y_shift,subsampling_factor,guiding_error_mag1,
                                      guiding_error_mag2,guiding_error_angle)
    else:
        # If we want to conserve, do so by operating on a copy of the array
        if(conserve):
            a = deepcopy(a)
        else:
            # Ensure it's contiguous
            a = np.ascontiguousarray(a)
    
        # Use the proper function for the data type
        if a.dtype == 'float32':
            f = cIceBRGpy.rebin_float
        elif a.dtype == 'float64':
            f = cIceBRGpy.rebin_double
        elif a.dtype == 'int32':
            f = cIceBRGpy.rebin_int
        elif a.dtype == 'int64':
            f = cIceBRGpy.rebin_long
        elif a.dtype == 'uint32':
            f = cIceBRGpy.rebin_uint
        elif a.dtype == 'uint64':
            f = cIceBRGpy.rebin_ulong
        else:
            a = np.asarray(a,dtype='float32')
            f = cIceBRGpy.rebin_float
    
        new_shape_c = f(a, y_shift, x_shift, subsampling_factor) # Note swap of indices here
        new_shape = new_shape_c[1], new_shape_c[0]
    
        # Resort the new array into the proper shape
        new_size = np.product(new_shape)
    
        rebinned_array = np.reshape(np.ravel(a)[0:new_size], new_shape)
    
    if np.any(np.isnan(rebinned_array)):
        pass
    assert(not np.any(np.isnan(rebinned_array)))
    
    # Convolve it with the charge diffusion kernel
    rebinned_diffused_array = convolve(rebinned_array, cdf, mode="same")
    
    return rebinned_diffused_array