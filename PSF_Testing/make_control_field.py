#!/usr/bin/env python

""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/make_control_field.py

    Created 19 Sep 2016

    @TODO: File docstring

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

import sys
from psf_testing.control_field import make_control_field

def main(argv):
    """ @TODO main docstring
    """
    make_control_field("/disk2/brg/Data/HST_Fields/control_image_n_b3.fits",
                       random_seed=3,
                       num_grid_points=(0, 0),
                       sky_level=0.,
                       read_noise=51., # Approximately right for background, but doesn't include unresolved sources
                       suppress_noise=True,
                       num_stars=1000,
                       binary_fraction=0.3,
                       subsampling_factor=20,
                       randomize_spectral_type=False,
                       use_cache=True)
    pass

if __name__ == "__main__":
    main(sys.argv)
