""" @file make_weight_mask.py

    Created 17 Sep 2015

    Function to make a weight mask given a function and centre position

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

from psf_testing import magic_values as mv
from moments.coords import get_x_and_y_of_array

def make_weight_mask(weight_func,
                     nx = 2*mv.min_stamp_size+1,
                     ny = 2*mv.min_stamp_size+1,
                     xc = None,
                     yc = None,
                     x_array = None,
                     y_array = None):
    """ Make a weight mask from a weight function and an optional centre point.
    
        Requires: weight_func <binary function> (Function of x,y. Must be able to handle numpy
                                                 ndarrays element-wise.)
        Optional: nx <int> (result array x size)
                  ny <int> (result array y size)
                  xc <int> (x centre position. If not provided, will use centre of array)
                  yc <int> (y centre position. If not provided, will use centre of array)
                  x_array <ndarray> (array of x coordinates. Will be generated if not provided
                                     so that it can be reused)
                  y_array <ndarray> (array of x coordinates. Will be generated if not provided
                                     so that it can be reused)
                                    
        Returns: <ndarray> (Weight mask normalized so sum is 1)
    """
    
    if(xc is None):
        xc = (nx-1.)/2.
    if(yc is None):
        yc = (ny-1.)/2.
    
    # Generate x and y arrays if needed
    if((x_array is None) or (y_array is None)):
        x_array, y_array = get_x_and_y_of_array(nx, ny, xc, yc)
    
    weight_mask = weight_func(x_array,y_array)
    
    # Normalize it
    weight_mask = weight_mask/weight_mask.sum()
    
    return weight_mask