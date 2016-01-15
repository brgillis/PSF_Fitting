""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing/moments/get_Qsize.py

    Created 18 Sep 2015

    Functions to get the "Qsize" parameter of an image and its error.

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

def get_light_distribution(image,
                           prim_weight_func=mv.default_prim_weight_func,
                           xc=None,
                           yc=None,
                           x_array=None,
                           y_array=None,
                           background_noise=None,
                           gain=mv.gain):
    
    dmax = np.max(np.shape(image)) // 2
    
    if x_array is None or y_array is None:
        if(xc is None) or (yc is None):
            xc, yc, _, _, _, _ = centre_image(image, prim_weight_func)
        
        nx, ny = np.shape(image)
        x_array, y_array, _, _, _ = get_coords_of_array(nx=nx, ny=ny, xc=xc, yc=yc)
        
    if background_noise is None:
        background_noise = get_background_noise(image)
    
    # Initialize things so we can loop through radially

    # Get the radial distances of all pixels from the centre
    r_array = np.sqrt(np.square(x_array) + np.square(y_array))

    # Get raveled arrays for the image and radius
    flattened_image = np.ravel(image)
    flattened_rs = np.ravel(r_array)
    
    # Initialize arrays and quantities we'll sum
    I = []
    N = []
    N_lt = []
    I_mean = np.zeros(dmax)
    I_lt_mean = np.zeros(dmax)
    var_I_mean = np.zeros(dmax)
    W = np.zeros(dmax)
    var_W = np.zeros(dmax)
    covar_W = np.zeros((dmax, dmax))
    for ri in xrange(dmax):

        # Start by binning pixels by radial distances
        I.append(flattened_image[ np.logical_and(flattened_rs <= ri + 0.5,
                                                  flattened_rs > ri - 0.5) ])
        N.append(np.size(I[ri]))

        # Get the number of contained pixels
        N_lt.append(0)
        for rj in xrange(ri):
            N_lt[ri] += N[rj]

        if N[-1] > 0:
            # If we have any pixels here, get the mean, contained mean, and variances
            I_mean[ri] = np.mean(I[ri])

            if N_lt[ri] == 0:
                I_lt_mean[ri] = 0
            else:
                for rj in xrange(ri):
                    I_lt_mean[ri] += (I_mean[rj] * N[rj]) / N_lt[ri]


            var_I_mean[ri] = np.sum(np.abs(I[ri]) / gain + np.square(background_noise)) \
                                / (np.square(N[ri]))

        # Get W for this bin
        if N_lt[ri] > 0:
            W[ri] = I_lt_mean[ri] - I_mean[ri]

        # Get the covariance of this bin with itself and every smaller bin

        # For the bin with itself, we want the first term to be the variance of the contained mean
        # For this, we need a sum of a function on all lesser bins.
        if N_lt[ri] == 0:
            covar_W[ri, ri] = 0.
        else:
            covar_W[ri, ri] = np.sum(np.square(N[0:ri]) * var_I_mean[0:ri]) \
                              / np.square(N_lt[ri]) + var_I_mean[ri]

        var_W[ri] = covar_W[ri, ri]

        # Now, the covariance with all smaller bins
        for rj in xrange(ri):
            if (N_lt[ri] == 0) or (N_lt[rj] == 0):
                covar_W[ri, rj] = 0
            else:
                covar_W[ri, rj] = np.sum(np.square(N[0:rj]) * var_I_mean[0:rj]) \
                                    / (N_lt[ri] * N_lt[rj]) - \
                                 N[rj] * var_I_mean[rj] \
                                    / N_lt[ri]

            covar_W[rj, ri] = covar_W[ri, rj]

    return N, I_mean, W, covar_W
    

def get_Qsize_and_var(image,
                      prim_weight_func=mv.default_prim_weight_func,
                      sec_weight_func=mv.default_sec_weight_func,
                      xc=None,
                      yc=None,
                      background_noise=None,
                      gain=mv.gain):

    nx, ny = np.shape(image)
    dmax = np.max((nx, ny)) // 2

    if(xc is None) or (yc is None):
        xc, yc, _, _, _, _ = centre_image(image, prim_weight_func)
        
    xc = np.round(xc,0)
    yc = np.round(yc,0)
        
    x_array, y_array, _, _, _ = get_coords_of_array(nx=nx, ny=ny, xc=xc, yc=yc)

    if background_noise is None:
        background_noise = get_background_noise(image)

    # Get a 1-d version of the weight function
    def prim_radial_weight_func(r):
        return prim_weight_func(r, np.zeros_like(r)) # Assuming it's circular
    def sec_radial_weight_func(r):
        return sec_weight_func(r, np.zeros_like(r)) # Assuming it's circular
    
    N, _I_mean, W, covar_W = get_light_distribution(image=image,
                                prim_weight_func=prim_weight_func,
                                x_array=x_array,
                                y_array=y_array,
                                background_noise=background_noise,
                                gain=gain)

    # Put this into the calculation for Qsize
    ri_array = np.linspace(start=0., stop=dmax, num=dmax, endpoint=False)

    Qsize = np.zeros(2)
    covar_Qsize = np.zeros((2, 2))

    for i, weight_func in zip(range(2), (prim_radial_weight_func, sec_radial_weight_func)):

        w_ri_array = weight_func(ri_array)*N

        weighted_W = W * w_ri_array

        Qsize_numerator = (weighted_W * ri_array).sum()
        square_Qsize_numerator = np.square(Qsize_numerator)

        Qsize_denominator = weighted_W.sum()
        square_Qsize_denominator = np.square(Qsize_denominator)
        cube_Qsize_denominator = Qsize_denominator * square_Qsize_denominator
        quart_Qsize_denominator = Qsize_denominator * cube_Qsize_denominator

        # Get Qsize from this
        Qsize[i] = Qsize_numerator / Qsize_denominator

        # Now get the variance on Qsize

        w_ri_ri_array = ri_array * w_ri_array

        for j, weight_func in zip(range(2), (prim_radial_weight_func, sec_radial_weight_func)):

            other_w_ri_array = weight_func(ri_array)*N
            other_w_ri_ri_array = ri_array * other_w_ri_array

            var_Qsize_numerator = (np.outer(w_ri_ri_array, other_w_ri_ri_array) \
                                     * covar_W).sum()
            var_Qsize_denominator = (np.outer(w_ri_array, other_w_ri_array) \
                                     * covar_W).sum()

            covar_Qsize_num_denom = (np.outer(w_ri_ri_array, other_w_ri_array) \
                                     * covar_W).sum()

            covar_Qsize[i, j] = var_Qsize_numerator / square_Qsize_denominator \
                         + square_Qsize_numerator * var_Qsize_denominator / quart_Qsize_denominator \
                         - 2. * Qsize_numerator * covar_Qsize_num_denom / cube_Qsize_denominator
                         
    err_Qsize = np.sqrt(np.diag(covar_Qsize))

    return (Qsize * mv.pixel_scale,
            err_Qsize * mv.pixel_scale,
            covar_Qsize * np.square(mv.pixel_scale))
