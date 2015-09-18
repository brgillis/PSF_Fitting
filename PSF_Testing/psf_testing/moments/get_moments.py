""" @file get_psf_testing.moments.py

    Created 17 Sep 2015

    Function to get the moments of an image.

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

from psf_testing.moments.centre_image import centre_image
from psf_testing.moments.coords import get_coords_of_array
from psf_testing.moments.estimate_background import get_background_noise
from psf_testing.moments.make_weight_mask import make_weight_mask

def get_moments_and_variances(image,
                weight_func = lambda x,y : 1.,
                weight_mask = None,
                xc = None,
                yc = None,
                background_noise = None,
                gain = mv.gain):
    
    if((xc is None) or (yc is None)):
        # Don't overwrite weight_mask if it's given; trust the user
        if(weight_mask is None):
            xc, yc, x_array, y_array, weight_mask = centre_image(image, weight_func)
        else:
            xc, yc, x_array, y_array, _ = centre_image(image, weight_func)
        x2_array = np.square(x_array)
        y2_array = np.square(y_array)
        xy_array = x_array * y_array
    else:
        nx, ny = np.shape(image)
        x_array, y_array, x2_array, y2_array, xy_array = get_coords_of_array(nx=nx,
                                               ny=ny,
                                               xc=xc,
                                               yc=yc)
        if(weight_mask is None):
            weight_mask = make_weight_mask(weight_func=weight_func,
                                       nx=nx,
                                       ny=ny,
                                       xc=xc,
                                       yc=yc,
                                       x_array=x_array,
                                       y_array=y_array)
    
    # Get the unnormalized moments
    
    m0 = (image*weight_mask).sum()
    
    mx = (x_array*image*weight_mask).sum()
    my = (y_array*image*weight_mask).sum()
    
    mxx = (x2_array*image*weight_mask).sum()
    myy = (y2_array*image*weight_mask).sum()
    mxy = (xy_array*image*weight_mask).sum()
    
    # Normalize the moments and store these values in a tuple of tuples
    
    Mx = mx/m0
    My = my/m0
    
    Mxx = mxx/m0
    Myy = myy/m0
    Mxy = mxy/m0
    
    moments = ((m0,),
               (Mx,My),
               (Mxx,Myy,Mxy))
        
    # Now calculate the errors
    
    if(background_noise is None):
        background_noise = get_background_noise(image)
    
    # Store some new arrays we'll use first
    square_weight_mask = np.square(weight_mask)
    image_var = image/gain + np.square(background_noise)
    x4_array = np.square(x2_array)
    y4_array = np.square(y2_array)
    x2y2_array = np.square(xy_array)
    
    var_m0 = (image_var * square_weight_mask).sum()
    
    var_mx = (x2_array * image_var * square_weight_mask).sum()
    var_my = (y2_array * image_var * square_weight_mask).sum()
    
    var_mxx = (x4_array * image_var * square_weight_mask).sum()
    var_myy = (y4_array * image_var * square_weight_mask).sum()
    var_mxy = (x2y2_array * image_var * square_weight_mask).sum()
    
    # And now get the variances of the normalized moments
    square_m0 = np.square(m0)
    quart_m0 = np.square(square_m0)
    
    var_Mx = var_mx/square_m0 + np.square(mx)*var_m0/quart_m0
    var_My = var_my/square_m0 + np.square(my)*var_m0/quart_m0
    
    var_Mxx = var_mxx/square_m0 + np.square(mxx)*var_m0/quart_m0
    var_Myy = var_myy/square_m0 + np.square(myy)*var_m0/quart_m0
    var_Mxy = var_mxy/square_m0 + np.square(mxy)*var_m0/quart_m0
    
    variances = ((var_m0,),
                 (var_Mx,var_My),
                 (var_Mxx,var_Myy,var_Mxy))
    
    return moments, variances