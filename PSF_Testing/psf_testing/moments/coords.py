""" @file coords.py

    Created 17 Sep 2015

    Methods to get coordinate arrays.

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

def get_x_and_y_of_array(nx, ny, xc=None, yc=None):
    """ Gets arrays for x and y coordinates relative to the centre of an array.
    
        Requires: nx <int> (x size of desired array)
                  ny <int> (y size of desired array)
        Optional: xc <float> (x centre value)
                  yc <float> (y centre value)
                  
        Returns <ndarray>, <ndarray> (x, y)
    """
    if(xc is None):
        xc = (nx-1.)/2.
    if(yc is None):
        yc = (ny-1.)/2.
    
    indices = np.indices((nx, ny), dtype=int)
    x_array = indices[0] - xc
    y_array = indices[1] - yc
    
    return x_array, y_array

def get_coords_of_array(nx, ny, xc=None, yc=None):
    """ Gets arrays for coordinates relative to the centre of an array, plus 2nd order values.
    
        Requires: nx <int> (x size of desired array)
                  ny <int> (y size of desired array)
        Optional: xc <float> (x centre value)
                  yc <float> (y centre value)
                  
        Returns <ndarray> (x),
                <ndarray> (y),
                <ndarray> (x^2),
                <ndarray> (y^2),
                <ndarray> (xy)
    """
    
    x_array, y_array = get_x_and_y_of_array(nx, ny, xc, yc)
    
    x2_array = np.square(x_array)
    y2_array = np.square(y_array)
    xy_array = x_array * y_array
    
    return (x_array, 
            y_array,
            x2_array,
            y2_array,
            xy_array)