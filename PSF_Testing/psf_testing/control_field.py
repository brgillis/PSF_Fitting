""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing/control_field.py

    Created 16 Sep 2016

    @TODO: File docstring

    ---------------------------------------------------------------------

    Copyright (C) 2016 Bryan R. Gillis

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

import numpy as np
import galsim

from psf_testing import magic_values as mv
from psf_testing.get_model_psf import get_model_psf
from psf_testing.psf_model_scheme import psf_model_scheme

spec_types = np.linspace(0,18,19,endpoint=True)
spec_types_pdf = spec_types**3
spec_types_pdf[0] = spec_types_pdf[-1] = 0
spec_types_pdf /= spec_types_pdf.sum()
    

def make_control_field(image_filename,
                       
                       random_seed = 0,
                       
                       image_shape=mv.default_image_shape,
                       
                       num_stars=1000,
                       min_star_mag=mv.default_min_star_mag,
                       max_star_mag=mv.default_max_star_mag,
                       
                       exp_time=1298.,
                       scale=mv.pixel_scale,
                       gain=mv.gain,
                       zeropoint=mv.zeropoint,
                       read_noise=mv.read_noise,
                       sky_level=mv.sky_level,
                       suppress_noise=False,
                       num_grid_points=(0,0),
                       
                       focus=mv.default_focus,
                       tinytim_params=None,
                       subsampling_factor=mv.default_subsampling_factor,
                       randomize_spectral_type=True,
                       use_cache=True):
    
    if tinytim_params is None:
        tinytim_params = mv.default_tinytim_params
    tinytim_params['subsampling_factor'] = subsampling_factor
        
    scheme = psf_model_scheme(focus = focus,
                              num_grid_points = num_grid_points,
                              image_shape = image_shape)
        
    image = galsim.Image(image_shape[0],image_shape[1],scale=scale)
    
    im_zeropoint = zeropoint + 2.5*np.log10(exp_time)
    
    u = galsim.random.UniformDeviate(galsim.BaseDeviate(random_seed))
    spec_f = galsim.DistDeviate(galsim.BaseDeviate(random_seed+1), function=galsim.LookupTable(spec_types,spec_types_pdf))
    
    for _i in range(num_stars):
        
        # Choose a random position
        xp = u()*image_shape[0]
        yp = u()*image_shape[1]
        
        # Choose a random magnitude
        mag = min_star_mag + u()*(max_star_mag-min_star_mag)
        
        # Choose a random spectral type if desired
        if randomize_spectral_type:
            spec_type = (1,int(round(spec_f())))
        else:
            spec_type = mv.default_model_psf_spec_type
            _ = u() # Advance the deviate to keep seeding the same
        
        # Calculate the flux for this star given its magnitude
        flux = 10.0**(0.4*(im_zeropoint-mag))
        
        psf_data = get_model_psf(xp,yp,
                                 scheme=scheme,
                                 tinytim_params=tinytim_params,
                                 use_cache=use_cache,
                                 spec_type=spec_type)
        
        psf_data *= flux/psf_data.sum()
        
        psf_ny, psf_nx = np.shape(psf_data)
        psf_x_arm = int((psf_nx - 1.)/2)
        psf_y_arm = int((psf_ny - 1.)/2)
        
        image_window = image.array
        psf_window = psf_data
        
        xp = int(xp)
        yp = int(yp)
        
        if yp < psf_y_arm:
            psf_window = psf_window[psf_y_arm-yp:,:]
            image_window = image_window[:yp+psf_y_arm+1,:]
        elif yp > image_shape[1]-psf_y_arm-1:
            psf_window = psf_window[:image_shape[1]-psf_y_arm-1-yp,:]
            image_window = image_window[yp-psf_y_arm:,:]
        else:
            image_window = image_window[yp-psf_y_arm:yp+psf_y_arm+1,:]
        
        if xp < psf_x_arm:
            psf_window = psf_window[:,psf_x_arm-xp:]
            image_window = image_window[:,:xp+psf_x_arm+1]
        elif xp > image_shape[0]-psf_x_arm-1:
            psf_window = psf_window[:,:image_shape[0]-psf_x_arm-1-xp]
            image_window = image_window[:,xp-psf_x_arm:]
        else:
            image_window = image_window[:,xp-psf_x_arm:xp+psf_x_arm+1]
            
        if not np.shape(image_window)==np.shape(psf_window):
            pass
        image_window += psf_window
        
        pass
        
    # Set up the image's header
    image.header = {}
    image.header[mv.header_chip_keyword] = tinytim_params["chip"]
    image.header[mv.header_exp_time_keyword] = exp_time
    image.header[mv.header_gain_keyword] = gain
    image.header[mv.header_obs_date_keyword] = 0.
    image.header[mv.header_obs_time_keyword] = 0.
    image.header[mv.header_ra_keyword] = 0.
    image.header[mv.header_dec_keyword] = 0.
    
    # Add noise to the image
    if not suppress_noise:
        image.addNoise(galsim.CCDNoise(rng, gain=gain, read_noise=read_noise, sky_level=sky_level))

    # Output the image
    image.write(image_filename,clobber=True)    
    
    return