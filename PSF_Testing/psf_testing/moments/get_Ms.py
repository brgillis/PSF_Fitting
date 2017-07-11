""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/moments/get_Ms.py

    Created 17 Sep 2015

    Function to get the quantities we wish to compare for an image (actual
    star or model PSF)

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
from psf_testing.moments.Msize import get_Msize
from psf_testing.moments.centre_image import centre_image
from psf_testing.moments.coords import get_x_and_y_of_array
from psf_testing.moments.get_moments import get_moments
from psf_testing.moments.make_weight_mask import make_weight_mask


def get_m0_and_Ms(image,
                  prim_weight_func=mv.default_prim_weight_func,
                  sec_weight_func=mv.default_sec_weight_func,
                  xc=None,
                  yc=None):
        
    nx, ny = np.shape(image)

    if (xc is None) or (yc is None):
        xc, yc, x_array, y_array, prim_weight_mask, _ = centre_image(image, prim_weight_func)
    else:
        x_array, y_array = get_x_and_y_of_array(nx=nx, ny=ny, xc=xc, yc=yc)
        prim_weight_mask = make_weight_mask(weight_func=prim_weight_func,
                                       nx=nx,
                                       ny=ny,
                                       xc=xc,
                                       yc=yc,
                                       x_array=x_array,
                                       y_array=y_array)
    sec_weight_mask = make_weight_mask(weight_func=sec_weight_func,
                                   nx=nx,
                                   ny=ny,
                                   xc=xc,
                                   yc=yc,
                                   x_array=x_array,
                                   y_array=y_array)

    # Get the moments and variances
    ((m0,), (Mx, My), (_Mxx, _Myy, Mxy, Mplus)) = get_moments(image=image,
                                                      prim_weight_func=prim_weight_func,
                                                      prim_weight_mask=prim_weight_mask,
                                                      sec_weight_func=sec_weight_func,
                                                      sec_weight_mask=sec_weight_mask,
                                                      xc=xc,
                                                      yc=yc)

    # Get the M values from the moments, plus errors from the variances

    scale = mv.pixel_scale
    square_scale = np.square(scale)

    Mx *= scale

    My *= scale

    Mplus *= square_scale

    Mcross = 2. * Mxy * square_scale

    # Get Qsize and its error now
    Msize = get_Msize(image=image,
                      prim_weight_func=prim_weight_func,
                      sec_weight_func=sec_weight_func,
                      xc=xc,
                      yc=yc)

    # Put the Q values into numpy arrays
    Mxy = np.array([Mx, My])
    
    Mpcs = np.array([Mplus, Mcross, Msize])

    return m0, Mxy, Mpcs
