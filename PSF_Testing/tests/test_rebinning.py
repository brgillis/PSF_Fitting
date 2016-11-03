""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/tests/test_rebinning.py

    Created 3 Oct 2016

    Tests to make sure rebinning is behaving as expected.

    ---------------------------------------------------------------------

    Copyright (C) 2016 Bryan R. Gillis

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

from psf_testing.rebin_psf import rebin

ss_factor = 5

@pytest.fixture()
def rebinning_test_array():
    a = np.zeros((25,25),dtype=float)
    
    # Have it step by 1 in the first index in the centre 5x5 coords
    for i in range(5):
        a[10+i,10:15] += i+1
    
    # Have it step by 10 in the second index in the centre 5x5 coords
    for i in range(5):
        a[10:15,10+i] += 10*(i+1)
        
    # Centre should now look like:
    # [[ 11.,  21.,  31.,  41.,  51.],
    #  [ 12.,  22.,  32.,  42.,  52.],
    #  [ 13.,  23.,  33.,  43.,  53.],
    #  [ 14.,  24.,  34.,  44.,  54.],
    #  [ 15.,  25.,  35.,  45.,  55.]]

    return a

@pytest.fixture()
def rebinning_unity_kernel():
    return np.array([[1.]])

def test_rebin_size(rebinning_test_array,rebinning_unity_kernel):
    
    ss_nx, ss_ny = np.shape(rebinning_test_array)
    
    # Test with no shift
    
    rb_array_0 = rebin(rebinning_test_array,
                     rebinning_unity_kernel,
                     subsampling_factor=ss_factor,
                     conserve=True)
    
    assert np.shape(rb_array_0) == (ss_nx/ss_factor,ss_ny/ss_factor)
    
    # Test that a small shift will decrease the size as expected
    
    rb_array_xm1 = rebin(rebinning_test_array,
                     rebinning_unity_kernel,
                     subsampling_factor=ss_factor,
                     x_shift=-1,
                     conserve=True)
    
    assert np.shape(rb_array_xm1) == (ss_nx/ss_factor-2,ss_ny/ss_factor)
    
    rb_array_xp1 = rebin(rebinning_test_array,
                     rebinning_unity_kernel,
                     subsampling_factor=ss_factor,
                     x_shift=1,
                     conserve=True)
    
    assert np.shape(rb_array_xp1) == (ss_nx/ss_factor-2,ss_ny/ss_factor)
    
    rb_array_ym1 = rebin(rebinning_test_array,
                     rebinning_unity_kernel,
                     subsampling_factor=ss_factor,
                     y_shift=-1,
                     conserve=True)
    
    assert np.shape(rb_array_ym1) == (ss_nx/ss_factor,ss_ny/ss_factor-2)
    
    rb_array_yp1 = rebin(rebinning_test_array,
                     rebinning_unity_kernel,
                     subsampling_factor=ss_factor,
                     y_shift=1,
                     conserve=True)
    
    assert np.shape(rb_array_yp1) == (ss_nx/ss_factor,ss_ny/ss_factor-2)
    
    rb_array_xp1_yp1 = rebin(rebinning_test_array,
                     rebinning_unity_kernel,
                     subsampling_factor=ss_factor,
                     x_shift=1,
                     y_shift=1,
                     conserve=True)
    
    assert np.shape(rb_array_xp1_yp1) == (ss_nx/ss_factor-2,ss_ny/ss_factor-2)
    
    # Test a larger shift now
    
    rb_array_xmss = rebin(rebinning_test_array,
                     rebinning_unity_kernel,
                     subsampling_factor=ss_factor,
                     x_shift=-ss_factor,
                     conserve=True)
    
    assert np.shape(rb_array_xmss) == (ss_nx/ss_factor-2,ss_ny/ss_factor)
    
    rb_array_xmss1 = rebin(rebinning_test_array,
                     rebinning_unity_kernel,
                     subsampling_factor=ss_factor,
                     x_shift=-ss_factor-1,
                     conserve=True)
    
    assert np.shape(rb_array_xmss1) == (ss_nx/ss_factor-4,ss_ny/ss_factor)

def test_rebin_values(rebinning_test_array,rebinning_unity_kernel):
    
    ss_nx, ss_ny = np.shape(rebinning_test_array)
    
    # Test with no shift
    
    rb_array_0 = rebin(rebinning_test_array,
                     rebinning_unity_kernel,
                     subsampling_factor=ss_factor,
                     conserve=True)
    
    assert np.sum(rb_array_0)==np.sum(rebinning_test_array)
    assert rb_array_0[2,2]==825
    
    # Test with a small shift
    
    rb_array_xm1 = rebin(rebinning_test_array,
                         rebinning_unity_kernel,
                         subsampling_factor=ss_factor,
                         x_shift=-1,
                         conserve=True)
    
    assert np.sum(rb_array_xm1)==np.sum(rebinning_test_array)
    assert rb_array_xm1[1,2]==670
    assert rb_array_xm1[0,2]==155
    
    rb_array_xp1 = rebin(rebinning_test_array,
                         rebinning_unity_kernel,
                         subsampling_factor=ss_factor,
                         x_shift=1,
                         conserve=True)
    
    assert np.sum(rb_array_xp1)==np.sum(rebinning_test_array)
    assert rb_array_xp1[1,2]==650
    assert rb_array_xp1[2,2]==175
    
    rb_array_ym1 = rebin(rebinning_test_array,
                         rebinning_unity_kernel,
                         subsampling_factor=ss_factor,
                         y_shift=-1,
                         conserve=True)
    
    assert np.sum(rb_array_ym1)==np.sum(rebinning_test_array)
    assert rb_array_ym1[2,1]==760
    assert rb_array_ym1[2,0]==65
    
    rb_array_yp1 = rebin(rebinning_test_array,
                         rebinning_unity_kernel,
                         subsampling_factor=ss_factor,
                         y_shift=1,
                         conserve=True)
    
    assert np.sum(rb_array_yp1)==np.sum(rebinning_test_array)
    assert rb_array_yp1[2,1]==560
    assert rb_array_yp1[2,2]==265
    
    rb_array_xm1ym1 = rebin(rebinning_test_array,
                            rebinning_unity_kernel,
                            subsampling_factor=ss_factor,
                            x_shift=-1,
                            y_shift=-1,
                            conserve=True)
    
    assert np.sum(rb_array_xm1ym1)==np.sum(rebinning_test_array)
    assert rb_array_xm1ym1[1,1]==616
    assert rb_array_xm1ym1[1,0]==54
    assert rb_array_xm1ym1[0,1]==144
    assert rb_array_xm1ym1[0,0]==11
    
    