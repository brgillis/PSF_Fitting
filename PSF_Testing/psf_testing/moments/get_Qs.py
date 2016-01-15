""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/moments/get_Qs.py

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
from psf_testing.moments.Qsize import get_Qsize_and_var
from psf_testing.moments.centre_image import centre_image
from psf_testing.moments.coords import get_x_and_y_of_array
from psf_testing.moments.estimate_background import get_background_noise
from psf_testing.moments.get_moments import get_moments_and_variances
from psf_testing.moments.make_weight_mask import make_weight_mask


def get_m0_and_Qs(image,
                  prim_weight_func=mv.default_prim_weight_func,
                  sec_weight_func=mv.default_sec_weight_func,
                  xc=None,
                  yc=None,
                  background_noise=None,
                  gain=mv.gain):

    if background_noise is None:
        background_noise = get_background_noise(image)
        
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
    moments_and_variances = get_moments_and_variances(image=image,
                                                      prim_weight_func=prim_weight_func,
                                                      prim_weight_mask=prim_weight_mask,
                                                      sec_weight_func=sec_weight_func,
                                                      sec_weight_mask=sec_weight_mask,
                                                      xc=xc,
                                                      yc=yc,
                                                      background_noise=background_noise,
                                                      gain=gain)

    ((m0,), (Mx, My), (_Mxx, _Myy, Mxy, Mplus)) = \
        moments_and_variances[0]
    ((covar_m0,), (covar_Mx, covar_My),
     (_covar_Mxx, _covar_Myy, covar_Mxy, covar_Mplus)) = moments_and_variances[1]

    # Get the Q values from the moments, plus errors from the variances

    scale = mv.pixel_scale
    square_scale = np.square(scale)
    quart_scale = np.square(square_scale)
    
    err_m0 = np.sqrt(np.diag(covar_m0))

    Qx = Mx * scale
    covar_Qx = covar_Mx * square_scale
    err_Qx = np.sqrt(np.diag(covar_Qx))

    Qy = My * scale
    covar_Qy = covar_My * square_scale
    err_Qy = np.sqrt(np.diag(covar_Qy))

    Qplus = Mplus * square_scale
    covar_Qplus = covar_Mplus * quart_scale
    err_Qplus = np.sqrt(np.diag(covar_Qplus))

    Qcross = 2. * Mxy * square_scale
    covar_Qcross = 4. * covar_Mxy * quart_scale
    err_Qcross = np.sqrt(np.diag(covar_Qcross))

    # Get Qsize and its error now
    Qsize, err_Qsize, covar_Qsize = get_Qsize_and_var(image=image,
                                         prim_weight_func=prim_weight_func,
                                         sec_weight_func=sec_weight_func,
                                         xc=xc,
                                         yc=yc,
                                         background_noise=background_noise,
                                         gain=gain)

    # Put the Q values into numpy arrays
    Qxy = np.array([Qx, Qy])
    err_Qxy = np.array([err_Qx, err_Qy])
    covar_Qxy = np.array([covar_Qx, covar_Qy])
    
    Qpcs = np.array([Qplus, Qcross, Qsize])
    err_Qpcs = np.array([err_Qplus, err_Qcross, err_Qsize])
    covar_Qpcs = np.array([covar_Qplus, covar_Qcross, covar_Qsize])

    return m0, err_m0, covar_m0, Qxy, err_Qxy, covar_Qxy, Qpcs, err_Qpcs, covar_Qpcs
