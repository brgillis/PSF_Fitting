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

def get_light_distribution(image,
                           prim_weight_func=mv.default_prim_weight_func,
                           xc=None,
                           yc=None,
                           x_array=None,
                           y_array=None):
    
    dmax = np.max(np.shape(image)) // 2
    
    if x_array is None or y_array is None:
        if(xc is None) or (yc is None):
            xc, yc, _, _, _, _ = centre_image(image, prim_weight_func)
        
        nx, ny = np.shape(image)
        x_array, y_array, _, _, _ = get_coords_of_array(nx=nx, ny=ny, xc=xc, yc=yc)
    
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
    W = np.zeros(dmax)
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


        # Get W for this bin
        if N_lt[ri] > 0:
            W[ri] = I_lt_mean[ri] - I_mean[ri]

    return N, I_mean, W
    

def get_Msize(image,
                      prim_weight_func=mv.default_prim_weight_func,
                      sec_weight_func=mv.default_sec_weight_func,
                      xc=None,
                      yc=None):

    nx, ny = np.shape(image)
    dmax = np.max((nx, ny)) // 2

    if(xc is None) or (yc is None):
        xc, yc, _, _, _, _ = centre_image(image, prim_weight_func)
        
    try:
        xc = int(np.round(xc,0))
    except ValueError:
        pass
    yc = int(np.round(yc,0))
        
    x_array, y_array, _, _, _ = get_coords_of_array(nx=nx, ny=ny, xc=xc, yc=yc)

    # Get a 1-d version of the weight function
    def prim_radial_weight_func(r):
        return prim_weight_func(r, np.zeros_like(r)) # Assuming it's circular
    def sec_radial_weight_func(r):
        return sec_weight_func(r, np.zeros_like(r)) # Assuming it's circular
    
    N, _I_mean, W = get_light_distribution(image=image,
                                prim_weight_func=prim_weight_func,
                                x_array=x_array,
                                y_array=y_array)
    
    # Get similar stats for a point source image
    ps_image = np.zeros_like(image)
    ps_image[xc,yc] = 1
    
    N_ps, _I_mean_ps, W_ps = get_light_distribution(image=ps_image,
                                                    prim_weight_func=prim_weight_func,
                                                    x_array=x_array,
                                                    y_array=y_array)

    # Put this into the calculation for Qsize
    ri_array = np.linspace(start=0., stop=dmax, num=dmax, endpoint=False)

    msize = np.zeros(2)
    msize_ps = np.zeros(2)

    for i, weight_func in zip(range(2), (prim_radial_weight_func, sec_radial_weight_func)):

        w_ri_array = weight_func(ri_array)*N
        w_ri_ps_array = weight_func(ri_array)*N_ps

        weighted_W = W * w_ri_array
        weighted_W_ps = W_ps * w_ri_ps_array

        # Get Qsize from this
        msize[i] = (weighted_W * ri_array).sum() / weighted_W.sum()
        msize_ps[i] = (weighted_W_ps * ri_array).sum() / weighted_W_ps.sum()
        
    # Subtract the point source result to normalize
    msize -= msize_ps

    return (msize * mv.pixel_scale)**2
