#!/usr/bin/env python

""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/test_grid_convergence.py

    Created 26 Sep 2016

    Runs a series of focus fittings to see how subsampling converges.

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

from copy import deepcopy
from multiprocessing import Pool, cpu_count
import os
import sys

import numpy as np
from psf_testing import magic_values as mv
from psf_testing.smart_logging import get_default_logger
from psf_testing.test_psf import test_psf

image_size = (4096, 2048)

class test_psf_caller(object):
    def __init__(self, image_filename, *args, **kwargs):
        self.image_filename = image_filename
        self.args = args
        self.kwargs = kwargs
    def __call__(self, x):
        try:
            if x == 0:
                num_grid_points = (0, 0)
            else:
                num_grid_points = tuple(np.divide(image_size, x))
            test_psf(self.image_filename, num_grid_points=num_grid_points, results_tag="gp" + str(x), *self.args, **self.kwargs)
        except Exception as e:
            logger = get_default_logger()
            logger.error("Exception when trying grid points " + str(x) + ": " + str(e))
            raise

image_dir = "/disk2/brg/Data/HST_Fields/"
image_filename = "control_image_n_rs.fits"
results_dir = "/disk2/brg/Data/HST_Fields/grid_convergence_testing"
subsampling_factor = 13
parallelize = True

def main(argv):
    """ @TODO main docstring
    """

    grid_sizes = (2048, 1024, 512, 256, 128, 64, 32, 1)

    # Check if we're debugging
    try:
        import pydevd as unused_import
        debugging = True
    except ImportError:
        debugging = False

    caller = test_psf_caller(os.path.join(image_dir, image_filename),
                             subsampling_factor=subsampling_factor,
                             min_lowest_separation=1.0,
                             min_class_star=0.01,
                             min_star_mag=22.,
                             max_star_mag=25.,
                             focus=mv.default_focus,
                             norm_errors=True,
                             min_star_snr=50.,
                             results_dir=results_dir,
                             parallelize=True)

    for grid_size in grid_sizes:
        caller(grid_size)


if __name__ == "__main__":
    main(sys.argv)
