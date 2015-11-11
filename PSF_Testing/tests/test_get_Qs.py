""" @file test_get_Qs.py

    Created 18 Sep 2015

    Unit tests for functions to calculate Q values

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

import pytest

import numpy as np
from scipy import signal

from psf_testing.moments.get_Qs import get_m0_and_Qs

nx = ny = 41
image_sigma = 3.0
weight_sigma = 5.0
gain = 2.0
flux = 10000.0
background_noise = 10.0

test_shear = 0.2

seed = 151357

@pytest.fixture(scope="module")
def noiseless_circ_Gaussian():
    unnormed_Gaussian = np.outer(signal.gaussian(nx, image_sigma), signal.gaussian(ny, image_sigma))
    return unnormed_Gaussian/unnormed_Gaussian.sum()

@pytest.fixture(scope="module")
def noiseless_2x_circ_Gaussian():
    unnormed_Gaussian = np.outer(signal.gaussian(nx, 2*image_sigma), signal.gaussian(ny, image_sigma))
    return unnormed_Gaussian/unnormed_Gaussian.sum()

@pytest.fixture(scope="module")
def noiseless_4x_circ_Gaussian():
    unnormed_Gaussian = np.outer(signal.gaussian(nx, 4*image_sigma), signal.gaussian(ny, image_sigma))
    return unnormed_Gaussian/unnormed_Gaussian.sum()

@pytest.fixture(scope="module")
def noiseless_horiz_Gaussian():
    unnormed_Gaussian = np.outer(signal.gaussian(nx, image_sigma*(1+test_shear)),
                                 signal.gaussian(ny, image_sigma*(1-test_shear)))
    return unnormed_Gaussian/unnormed_Gaussian.sum()

@pytest.fixture(scope="module")
def circ_Gaussian_image(noiseless_circ_Gaussian):
    scaled_Gaussian = flux*noiseless_circ_Gaussian
    
    noise_per_pixel = np.sqrt(scaled_Gaussian/gain + np.square(background_noise))
    
    np.random.seed(seed)
    Gaussian_image = scaled_Gaussian + noise_per_pixel * np.random.randn(nx,ny)
    
    return Gaussian_image

@pytest.fixture(scope="module")
def horiz_Gaussian_image(noiseless_horiz_Gaussian):
    scaled_Gaussian = flux*noiseless_horiz_Gaussian
    
    noise_per_pixel = np.sqrt(scaled_Gaussian/gain + np.square(background_noise))
    
    np.random.seed(seed)
    Gaussian_image = scaled_Gaussian + noise_per_pixel * np.random.randn(nx,ny)
    
    return Gaussian_image

def test_circular_Gaussian(noiseless_circ_Gaussian,circ_Gaussian_image):
    
    tolerance = 2.0
    
    _m0, _err_m0, model_Qs, _err_model_Qs = get_m0_and_Qs(noiseless_circ_Gaussian)
    _m0, _err_m0, image_Qs, err_image_Qs = get_m0_and_Qs(circ_Gaussian_image)
    
    Zs = (image_Qs-model_Qs)/np.array((err_image_Qs[:,0,0],err_image_Qs[:,1,1])).transpose()
    
    assert(np.all(np.abs(Zs[0]) < tolerance))
    assert(np.all(np.abs(Zs[1]) < tolerance))
    assert(np.all(np.abs(Zs[2]) < tolerance))
    assert(np.all(np.abs(Zs[3]) < tolerance))
    assert(np.all(np.abs(Zs[4]) < tolerance))

def test_horiz_Gaussian(noiseless_horiz_Gaussian,horiz_Gaussian_image):
    
    tolerance = 2.0
    
    _m0, _err_m0, model_Qs, _err_model_Qs = get_m0_and_Qs(noiseless_horiz_Gaussian)
    _m0, _err_m0, image_Qs, err_image_Qs = get_m0_and_Qs(horiz_Gaussian_image,)
    
    Zs = (image_Qs-model_Qs)/np.array((err_image_Qs[:,0,0],err_image_Qs[:,1,1])).transpose()
    
    assert(np.all(np.abs(Zs[0]) < tolerance))
    assert(np.all(np.abs(Zs[1]) < tolerance))
    assert(np.all(np.abs(Zs[2]) < tolerance))
    assert(np.all(np.abs(Zs[3]) < tolerance))
    assert(np.all(np.abs(Zs[4]) < tolerance))
    
def test_err_calcs(noiseless_circ_Gaussian):
    
    tolerance = 0.2
    
    num_test_Gaussians = 100
    
    all_Qs = []
    
    np.random.seed(seed)
    
    scaled_Gaussian = flux*noiseless_circ_Gaussian
        
    noise_per_pixel = np.sqrt(scaled_Gaussian/gain + np.square(background_noise))
    
    for _i in xrange(num_test_Gaussians):
        
        Gaussian_image = scaled_Gaussian + noise_per_pixel * np.random.randn(nx,ny)
        _m0, _err_m0, Qs, _err_Qs = get_m0_and_Qs(Gaussian_image,
                                                 xc = (nx-1)/2.,
                                                 yc = (ny-1)/2.,
                                                 background_noise=background_noise,
                                                 gain=gain)
        
        all_Qs.append(Qs)
        
    _m0, _err_m0, _Qs, mean_err_Qs = get_m0_and_Qs(scaled_Gaussian,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.,
                                             background_noise=background_noise,
                                             gain=gain)
    
    mean_err_Qs = np.array((mean_err_Qs[:,0,0],mean_err_Qs[:,1,1])).transpose()
        
    std_Qs = np.std(all_Qs,axis=0)
    
    assert(np.all(np.abs(mean_err_Qs[0]-std_Qs[0])/std_Qs[0] < tolerance))
    assert(np.all(np.abs(mean_err_Qs[1]-std_Qs[1])/std_Qs[1] < tolerance))
    assert(np.all(np.abs(mean_err_Qs[2]-std_Qs[2])/std_Qs[2] < tolerance))
    assert(np.all(np.abs(mean_err_Qs[3]-std_Qs[3])/std_Qs[3] < tolerance))
    assert(np.all(np.abs(mean_err_Qs[4]-std_Qs[4])/std_Qs[4] < tolerance))
    
    pass

def test_Qsize_scaling(noiseless_circ_Gaussian,noiseless_2x_circ_Gaussian,
                       noiseless_4x_circ_Gaussian):
    
    # Get size measurement for 1x size Gaussian
    _m0, _err_m0, Qs_1x, _err_Qs = get_m0_and_Qs(flux*noiseless_circ_Gaussian,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.,
                                             background_noise=background_noise,
                                             gain=gain)
    Qsize_1x = Qs_1x[4]
    
    # Get size measurement for 2x size Gaussian
    _m0, _err_m0, Qs_2x, _err_Qs = get_m0_and_Qs(flux*noiseless_2x_circ_Gaussian,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.,
                                             background_noise=background_noise,
                                             gain=gain)
    Qsize_2x = Qs_2x[4]
    
    # Get size measurement for 4x size Gaussian
    _m0, _err_m0, Qs_4x, _err_Qs = get_m0_and_Qs(flux*noiseless_4x_circ_Gaussian,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.,
                                             background_noise=background_noise,
                                             gain=gain)
    Qsize_4x = Qs_4x[4]
    
    assert(np.all(Qsize_4x > Qsize_2x) and np.all(Qsize_2x>Qsize_1x), "Qsize doesn't scale with size of " +
           "input Gaussian.")
    
    pass