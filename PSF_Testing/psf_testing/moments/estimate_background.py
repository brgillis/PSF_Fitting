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
from psf_testing.remove_outliers import remove_outliers

def get_background_level_and_noise(image):
    """ Get the background level and noise of an image by sampling pixels around the edge
        and removing outliers.
    
        Requires: image <2d array>
        
        Returns: background_level <float>, background_noise <float> (as standard deviation)
    """
    
    # Get pixels around the edge of the image
    edge_pixels = np.concatenate((image[0,:],image[-1,:],image[1:-2,0],image[1:-2,-1]))
    
    masked_edge_pixels = remove_outliers(edge_pixels)
    
    good_pixels = masked_edge_pixels[~masked_edge_pixels.mask]
    
    return np.mean(good_pixels), np.std(good_pixels)

def get_background_level(image):
    """ Get the background level of an image by sampling pixels around the edge
        and removing outliers.
    
        Requires: image <2d array>
        
        Returns: background_level <float>
    """
    
    return get_background_level_and_noise(image)[0]

def get_background_noise(image):
    """ Get the background noise of an image by sampling pixels around the edge
        and removing outliers.
    
        Requires: image <2d array>
        
        Returns: background_noise <float> (as standard deviation)
    """
    
    return get_background_level_and_noise(image)[1]