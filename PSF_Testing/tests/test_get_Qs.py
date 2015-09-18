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
background_noise = 1.0

test_shear = 0.2

seed = 151357

@pytest.fixture(scope="module")
def test_weight_func():
    def func(x,y):
        r = np.sqrt(np.square(x)+np.square(y))
        return np.exp(-np.square(r) / (2. * np.square(weight_sigma)))
    return func

@pytest.fixture(scope="module")
def noiseless_circ_Gaussian():
    unnormed_Gaussian = np.outer(signal.gaussian(nx, image_sigma), signal.gaussian(ny, image_sigma))
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

def test_circular_Gaussian(noiseless_circ_Gaussian,circ_Gaussian_image,test_weight_func):
    
    tolerance = 2.0
    
    _m0, _err_m0, model_Qs, _err_model_Qs = get_m0_and_Qs(noiseless_circ_Gaussian,
                                                          weight_func=test_weight_func)
    _m0, _err_m0, image_Qs, err_image_Qs = get_m0_and_Qs(circ_Gaussian_image,
                                                         weight_func=test_weight_func)
    
    Zs = (image_Qs-model_Qs)/err_image_Qs
    
    assert(np.abs(Zs[0]) < tolerance)
    assert(np.abs(Zs[1]) < tolerance)
    assert(np.abs(Zs[2]) < tolerance)
    assert(np.abs(Zs[3]) < tolerance)
    assert(np.abs(Zs[4]) < tolerance)

def test_horiz_Gaussian(noiseless_horiz_Gaussian,horiz_Gaussian_image,test_weight_func):
    
    tolerance = 2.0
    
    _m0, _err_m0, model_Qs, _err_model_Qs = get_m0_and_Qs(noiseless_horiz_Gaussian,
                                                          weight_func=test_weight_func)
    _m0, _err_m0, image_Qs, err_image_Qs = get_m0_and_Qs(horiz_Gaussian_image,
                                                         weight_func=test_weight_func)
    
    Zs = (image_Qs-model_Qs)/err_image_Qs
    
    assert(np.abs(Zs[0]) < tolerance)
    assert(np.abs(Zs[1]) < tolerance)
    assert(np.abs(Zs[2]) < tolerance)
    assert(np.abs(Zs[3]) < tolerance)
    assert(np.abs(Zs[4]) < tolerance)