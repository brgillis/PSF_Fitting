""" @file remove_outliers_test.py

    Created 18 Sep 2015

    @TODO: File docstring

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

from psf_testing.remove_outliers import remove_outliers

@pytest.fixture()
def array_with_outliers():
    
    return np.array([0.9,0.95,0.95,1.,1.,1.05,1.05,1.1,
                     10,
                     1000])

def test_my_test(array_with_outliers):
    
    trimmed_array = remove_outliers(array_with_outliers)
    
    assert(np.size(trimmed_array[~trimmed_array.mask]) == np.size(array_with_outliers)-2 )
    assert(np.mean(trimmed_array[~trimmed_array.mask]) == 1.)