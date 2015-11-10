""" @file brute_force_minimize.py

    Created 5 Nov 2015

    Function to minimize a one-parameter function through a brute force approach
    followed by a precision fitting.

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

def bf_minimize(func,
                args=[],
                kwargs={},
                min_input=mv.default_min_test_focus,
                max_input=mv.default_max_test_focus,
                test_points=mv.default_focus_samples,
                precision=mv.default_focus_precision,
                max_iter=100):

    assert precision > 0
    assert max_input > min_input
    assert test_points > 0

    # Set up test points
    test_points = np.linspace(start=min_input,
                              stop=max_input,
                              num=test_points)

    step_size = test_points[1] - test_points[0]

    # Define a function with the args and kwargs passed to it
    def afunc(x):
        return func(x,*args,**kwargs)
    
    # Define a vectorized version of the function
    vfunc = np.vectorize(afunc, otypes=[float])

    # Get the best point in the initial set
    test_outs = vfunc(test_points)

    best_i = np.argmin(test_outs)
    best_out = test_outs[best_i]
    best_point = test_points[best_i]

    # Now zoom in, looking on either side of this best point and repeating

    test_outs = np.zeros(3)
    i = 0

    while(np.abs(step_size) > precision and i<max_iter):
        step_size /= 2
        i += 1

        test_points = np.linspace(start=best_point - step_size,
                                  stop=best_point + step_size,
                                  num=3)

        test_outs[0] = vfunc(test_points[0])
        test_outs[1] = best_out
        test_outs[2] = vfunc(test_points[2])

        best_i = np.argmin(test_outs)
        best_out = test_outs[best_i]
        best_point = test_points[best_i]

    return best_point
