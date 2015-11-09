""" @file get_centroid_of_image.py

    Created 17 Sep 2015

    Function to get the centroid of an image, where the dipole moments are zero.

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

from psf_testing.moments.make_weight_mask import make_weight_mask
from psf_testing.moments.coords import get_x_and_y_of_array
from psf_testing import magic_values as mv

def centre_image(image,
                 weight_func = lambda x,y: np.ones_like(x),
                 precision = 0.01):
    """ Get the centroid of an image for a given weight function.
    
        Requires: image <ndarray>
                  weight_func <function> (must be able to take ndarrays and work element-wise)
                  precision <float> (to what fraction of a pixel should we get the centroid?)
                  
        Returns: xc <float>,
                 yc <float>,
                 x_array <ndarray> (x coords relative to centroid),
                 y_array <ndarray> (y coords relative to centroid),
                 weight_mask <ndarray>
                 m0 <float> (monopole moment of image)
    """
    # Get the shape of the image
    nx, ny = np.shape(image)
    
    # Start with the center
    xc = (nx-1.)/2.
    yc = (ny-1.)/2.
    
    for _ in range(mv.max_loop_iterations):
        
        # Get coordinate arrays
        x_array, y_array = get_x_and_y_of_array(nx, ny, xc, yc)
        
        if not (xc>=0 or xc<0):
            pass
        
        # Make the weight mask
        weight_mask = make_weight_mask(weight_func, nx, ny, xc, yc, x_array, y_array)
        
        # Get the dipole moments
        m0 = (image*weight_mask).sum()
        if not (m0>0):
            pass
        dx = (x_array*image*weight_mask).sum()/m0
        dy = (y_array*image*weight_mask).sum()/m0
        
        # Get the larger distance we'll move
        d = max((np.abs(dx),np.abs(dy)))
        
        # If we're under the precision, break
        if(d<precision):
            break
        
        # Move the centroid
        xc += dx
        yc += dy
        
    # Return the results 
    return (xc,
            yc,
            x_array,
            y_array,
            weight_mask,
            m0)