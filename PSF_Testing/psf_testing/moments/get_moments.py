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
from psf_testing.moments.make_weight_mask import make_weight_mask

def get_moments(image,
                              prim_weight_func=mv.default_prim_weight_func,
                              sec_weight_func=mv.default_sec_weight_func,
                              prim_weight_mask=None,
                              sec_weight_mask=None,
                              xc=None,
                              yc=None):

    nx, ny = np.shape(image)

    if (xc is None) or (yc is None):
        # Don't overwrite weight_mask if it's given; trust the user
        if prim_weight_mask is None:
            xc, yc, x_array, y_array, prim_weight_mask, _ = centre_image(image, prim_weight_func)
        else:
            xc, yc, x_array, y_array, _, _ = centre_image(image, prim_weight_func)

        x2_array = np.square(x_array)
        y2_array = np.square(y_array)
        xy_array = x_array * y_array
    else:
        x_array, y_array, x2_array, y2_array, xy_array = get_coords_of_array(nx=nx,
                                               ny=ny,
                                               xc=xc,
                                               yc=yc)
        if prim_weight_mask is None:
            prim_weight_mask = make_weight_mask(weight_func=prim_weight_func,
                                           nx=nx,
                                           ny=ny,
                                           xc=xc,
                                           yc=yc,
                                           x_array=x_array,
                                           y_array=y_array)
    if sec_weight_mask is None:
        sec_weight_mask = make_weight_mask(weight_func=sec_weight_func,
                                       nx=nx,
                                       ny=ny,
                                       xc=xc,
                                       yc=yc,
                                       x_array=x_array,
                                       y_array=y_array)
    plus_array = x2_array - y2_array

    m0 = np.zeros(2)
    mx = np.zeros(2)
    my = np.zeros(2)
    mxx = np.zeros(2)
    myy = np.zeros(2)
    mxy = np.zeros(2)
    mplus = np.zeros(2)

    # Get the unnormalized moments

    for i, weight_mask in zip(range(2), (prim_weight_mask, sec_weight_mask)):

        m0[i] = (image * weight_mask).sum()

        mx[i] = (x_array * image * weight_mask).sum()
        my[i] = (y_array * image * weight_mask).sum()

        mxx[i] = (x2_array * image * weight_mask).sum()
        myy[i] = (y2_array * image * weight_mask).sum()
        mxy[i] = (xy_array * image * weight_mask).sum()
        mplus[i] = (plus_array * image * weight_mask).sum()

    # Normalize the moments and store these values in a tuple of tuples

    Mx = mx / m0
    My = my / m0

    Mxx = mxx / m0
    Myy = myy / m0
    Mxy = mxy / m0
    Mplus = mplus / m0

    moments = ((m0,),
               (Mx, My),
               (Mxx, Myy, Mxy, Mplus))

    return moments
