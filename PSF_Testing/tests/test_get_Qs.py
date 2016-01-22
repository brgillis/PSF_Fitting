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
bg_level = 1

test_shear = 0.2

seed = 151357

def check_close(a,b,m_tol,a_tol=0):
    assert (a >= b*(1-m_tol)-a_tol) and (a <= b*(1+m_tol)+a_tol)

def check_all_close(la,lb,m_tol,a_tol=0):
    for a, b in zip(la,lb):
        check_close(a,b,m_tol,a_tol)

@pytest.fixture(scope="module")
def noiseless_circ_Gaussian():
    unnormed_Gaussian = np.outer(signal.gaussian(nx, image_sigma), signal.gaussian(ny, image_sigma))
    return unnormed_Gaussian/unnormed_Gaussian.sum()

@pytest.fixture(scope="module")
def noiseless_circ_2x_Gaussian():
    unnormed_Gaussian = np.outer(signal.gaussian(nx, 2*image_sigma), signal.gaussian(ny, image_sigma))
    return unnormed_Gaussian/unnormed_Gaussian.sum()

@pytest.fixture(scope="module")
def noiseless_circ_4x_Gaussian():
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
def circ_2x_Gaussian_image(noiseless_circ_2x_Gaussian):
    scaled_Gaussian = flux*noiseless_circ_2x_Gaussian
    
    noise_per_pixel = np.sqrt(scaled_Gaussian/gain + np.square(background_noise))
    
    np.random.seed(seed)
    Gaussian_image = scaled_Gaussian + noise_per_pixel * np.random.randn(nx,ny)
    
    return Gaussian_image

@pytest.fixture(scope="module")
def circ_4x_Gaussian_image(noiseless_circ_4x_Gaussian):
    scaled_Gaussian = flux*noiseless_circ_4x_Gaussian
    
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

def test_circular_Gaussian(noiseless_circ_Gaussian, circ_Gaussian_image,
                              noiseless_circ_2x_Gaussian, circ_2x_Gaussian_image,
                              noiseless_circ_4x_Gaussian, circ_4x_Gaussian_image,):
    
    m_tol = 0.01
    a_tol = 0.01
    
    _model_m0, model_Qxy, model_Qpcs = get_m0_and_Qs(noiseless_circ_Gaussian)
    _image_m0, image_Qxy, image_Qpcs = get_m0_and_Qs(circ_Gaussian_image)
    
    for image_Q, model_Q in zip(image_Qxy, model_Qxy):
        check_all_close(image_Q,model_Q,m_tol,a_tol)
    for image_Q, model_Q in zip(image_Qpcs, model_Qpcs):
        check_all_close(image_Q,model_Q,m_tol,a_tol)

def test_horiz_Gaussian(noiseless_horiz_Gaussian,horiz_Gaussian_image):
    
    m_tol = 0.01
    a_tol = 0.01
    
    _model_m0, model_Qxy, model_Qpcs = get_m0_and_Qs(noiseless_horiz_Gaussian)
    _image_m0, image_Qxy, image_Qpcs = get_m0_and_Qs(horiz_Gaussian_image)
    
    for image_Q, model_Q in zip(image_Qxy, model_Qxy):
        check_all_close(image_Q,model_Q,m_tol,a_tol)
    for image_Q, model_Q in zip(image_Qpcs, model_Qpcs):
        check_all_close(image_Q,model_Q,m_tol,a_tol)

def test_Qsize_scaling(noiseless_circ_Gaussian, circ_Gaussian_image,
                              noiseless_circ_2x_Gaussian, circ_2x_Gaussian_image,
                              noiseless_circ_4x_Gaussian, circ_4x_Gaussian_image,):
    
    # Get size measurement for 1x size Gaussian
    _m0_1x, _Qxy_1x, Qpcs_1x = get_m0_and_Qs(flux*noiseless_circ_Gaussian,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.)
    Qsize_1x = Qpcs_1x[2]
    
    # Get size measurement for 2x size Gaussian
    _m0_1x, _Qxy_2x, Qpcs_2x = get_m0_and_Qs(flux*noiseless_circ_2x_Gaussian,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.)
    Qsize_2x = Qpcs_2x[2]
    
    # Get size measurement for 4x size Gaussian
    _m0_1x, _Qxy_4x, Qpcs_4x = get_m0_and_Qs(flux*noiseless_circ_4x_Gaussian,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.)
    Qsize_4x = Qpcs_4x[2]
    
    assert(np.all(Qsize_4x > Qsize_2x) and np.all(Qsize_2x>Qsize_1x), "Qsize doesn't scale with size of " +
           "input Gaussian.")
    
    # Get size measurement for 1x size noisy Gaussian
    _m0_1x, _Qxy_1x, Qpcs_1x = get_m0_and_Qs(flux*circ_Gaussian_image,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.)
    Qsize_1x = Qpcs_1x[2]
    
    # Get size measurement for 2x size noisy Gaussian
    _m0_1x, _Qxy_2x, Qpcs_2x = get_m0_and_Qs(flux*circ_2x_Gaussian_image,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.)
    Qsize_2x = Qpcs_2x[2]
    
    # Get size measurement for 4x size noisy Gaussian
    _m0_1x, _Qxy_4x, Qpcs_4x = get_m0_and_Qs(flux*circ_4x_Gaussian_image,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.)
    Qsize_4x = Qpcs_4x[2]
    
    assert(np.all(Qsize_4x > Qsize_2x) and np.all(Qsize_2x>Qsize_1x), "Qsize doesn't scale with size of " +
           "input Gaussian.")
    
    pass

def test_Qsize_independence(noiseless_circ_Gaussian, circ_Gaussian_image):
    
    m_tol = 0.000001
    
    # Get size measurement for 1x size Gaussian
    _m0_1x, _Qxy_1x, Qpcs_1x = get_m0_and_Qs(flux*noiseless_circ_Gaussian,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.)
    Qsize_1x = Qpcs_1x[2]
    
    # Get size measurement for 1x size, higher flux Gaussian
    _m0_1x, _Qxy_1x, Qpcs_1x = get_m0_and_Qs(2*flux*noiseless_circ_Gaussian,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.)
    Qsize_1x_double_flux = Qpcs_1x[2]
    
    check_all_close(Qsize_1x_double_flux, Qsize_1x, m_tol)
    
    # Get size measurement for 1x size, bg-added Gaussian
    _m0_1x, _Qxy_1x, Qpcs_1x = get_m0_and_Qs(flux*noiseless_circ_Gaussian+bg_level,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.)
    Qsize_1x_bg_added = Qpcs_1x[2]
    
    check_all_close(Qsize_1x_bg_added, Qsize_1x, m_tol)
    
    # Get size measurement for 1x size, double flux, bg-added Gaussian
    _m0_1x, _Qxy_1x, Qpcs_1x = get_m0_and_Qs(2*flux*noiseless_circ_Gaussian+bg_level,
                                             xc = (nx-1)/2.,
                                             yc = (ny-1)/2.)
    Qsize_1x_double_flux_bg_added = Qpcs_1x[2]
    
    check_all_close(Qsize_1x_double_flux_bg_added, Qsize_1x, m_tol)
    
    pass