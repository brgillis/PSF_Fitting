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
                              prim_weight_func=mv.default_prim_weight_func,
                              sec_weight_func=mv.default_sec_weight_func,
                              prim_weight_mask=None,
                              sec_weight_mask=None,
                              xc=None,
                              yc=None,
                              background_noise=None,
                              gain=mv.gain):

    nx, ny = np.shape(image)

    if (xc is None) or (yc is None):
        # Don't overwrite weight_mask if it's given; trust the user
        if prim_weight_mask is None:
            xc, yc, x_array, y_array, prim_weight_mask, m0 = centre_image(image, prim_weight_func)
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

    # Now calculate the errors

    if background_noise is None:
        background_noise = get_background_noise(image)

    # Store some new arrays we'll use first
    weight_mask = (prim_weight_mask, sec_weight_mask)

    image_var = np.abs(image) / gain + np.square(background_noise)
    x4_array = np.square(x2_array)
    y4_array = np.square(y2_array)
    x2y2_array = np.square(xy_array)
    plus2_array = np.square(plus_array)

    var_m0 = np.zeros((2, 2))
    var_mx = np.zeros((2, 2))
    var_my = np.zeros((2, 2))
    var_mxx = np.zeros((2, 2))
    var_myy = np.zeros((2, 2))
    var_mxy = np.zeros((2, 2))
    var_mplus = np.zeros((2, 2))

    for i in range(2):
        for j in range(2):

            weighted_image_var = image_var * weight_mask[i] * weight_mask[j]

            var_m0[i, j] = weighted_image_var.sum()

            var_mx[i, j] = (x2_array * weighted_image_var).sum()
            var_my[i, j] = (y2_array * weighted_image_var).sum()

            var_mxx[i, j] = (x4_array * weighted_image_var).sum()
            var_myy[i, j] = (y4_array * weighted_image_var).sum()
            var_mxy[i, j] = (x2y2_array * weighted_image_var).sum()
            var_mplus[i, j] = (plus2_array * weighted_image_var).sum()

            square_m0 = m0[i] * m0[j]

    # And now get the variances of the normalized moments
    quart_m0 = np.square(square_m0)

    var_Mx = var_mx / square_m0 + np.square(mx) * var_m0 / quart_m0
    var_My = var_my / square_m0 + np.square(my) * var_m0 / quart_m0

    var_Mxx = var_mxx / square_m0 + np.square(mxx) * var_m0 / quart_m0
    var_Myy = var_myy / square_m0 + np.square(myy) * var_m0 / quart_m0
    var_Mxy = var_mxy / square_m0 + np.square(mxy) * var_m0 / quart_m0
    var_Mplus = var_mplus / square_m0 + np.square(mplus) * var_m0 / quart_m0

    variances = ((var_m0,),
                 (var_Mx, var_My),
                 (var_Mxx, var_Myy, var_Mxy, var_Mplus))

    return moments, variances
