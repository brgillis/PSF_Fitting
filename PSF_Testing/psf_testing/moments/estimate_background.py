""" @file estimate_background.py

    Created 18 Sep 2015

    Functions to get the background level and noise of an image.

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

import numpy as np
from psf_testing import magic_values as mv
from psf_testing.remove_outliers import remove_outliers
from psf_testing.moments.centre_image import centre_image

def get_background_level_and_noise(image,
                                   rmin=mv.min_stamp_size,
                                   rmax=mv.min_stamp_size+2,
                                   weight_func=mv.default_prim_weight_func):
    """ Get the background level and noise of an image by sampling pixels sufficiently distant from the centre
        and removing outliers
    
        Requires: image <2d array>
        
        Returns: background_level <float>, background_noise <float> (as standard deviation)
    """
    
    x_array, y_array = centre_image(image,weight_func)[2:4]
    
    r_array = np.sqrt(x_array**2+y_array**2)
    
    # Get pixels around the edge of the image
    edge_pixels = np.logical_and(r_array.ravel() > rmin, r_array.ravel() <= rmax)
    edge_pixel_values = image.ravel()[edge_pixels]
    
    masked_edge_pixels = remove_outliers(edge_pixel_values)
    
    good_pixels = masked_edge_pixels[~masked_edge_pixels.mask]
    
    return np.mean(good_pixels), np.std(good_pixels)

def get_background_level(image,*args,**kwargs):
    """ Get the background level of an image by sampling pixels around the edge
        and removing outliers.
    
        Requires: image <2d array>
        
        Returns: background_level <float>
    """
    
    return get_background_level_and_noise(image,*args,**kwargs)[0]

def get_background_noise(image,*args,**kwargs):
    """ Get the background noise of an image by sampling pixels around the edge
        and removing outliers.
    
        Requires: image <2d array>
        
        Returns: background_noise <float> (as standard deviation)
    """
    
    return get_background_level_and_noise(image,*args,**kwargs)[1]