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
    # Check if we're debugging
    try:
        import pydevd as _
        debugging = True
    except ImportError:
        debugging = False

    num_fields = 100
    i_start = 0
    mode = "none"

    if len(argv)>1:
        num_fields = int(argv[1])
    if len(argv)>2:
        i_start = int(argv[2])

    if len(argv)>3:
        mode = argv[3]

    tag = ""
    base_image_root = None
    binary_fraction = 0.0
    binary_r_max = 1.
    guiding_error_mag1=0.
    guiding_error_mag2=0.
    randomize_spectral_type=False
    suppress_noise = False

    if mode!="none":
        tag = "_" + mode

    if mode == "gb" or mode == "full":
        base_image_root = "/home/brg/Data/HST_Fields/mock_galaxies_"
    if mode == "b3" or mode == "full":
        binary_fraction = 0.3
    if mode == "b32":
        binary_fraction = 0.3
        binary_r_max = 2.
    if mode == "blur":
        guiding_error_mag1 = 0.010/0.05
    if mode == "blur2" or mode == "full":
        guiding_error_mag1 = 0.010/0.05
        guiding_error_mag2 = 0.010/0.05
    if mode == "rs" or mode == "full":
        randomize_spectral_type = True
    if mode == "nf":
        suppress_noise = True
        


    def run_for_i(i):
        focus = -6.0 + (i%100)/10.

        print("Generating image " + str(i) + " with focus " + str(focus))

        if base_image_root is not None:
            base_image=base_image_root+str(i)+".fits"
        else:
            base_image = None

        make_control_field("/home/brg/Data/HST_Fields/control_image_"+str(i)+tag+".fits",
                   random_seed=10*(i+1),
                   num_grid_points=(0, 0),
                   focus=focus,
                   subtracted_sky_level=400., # Approximately right for background, but doesn't include unresolved sources
                   unsubtracted_sky_level=5.,
                   read_noise=20.8, 
                   base_image=base_image,
                   suppress_noise=suppress_noise,
                   num_stars=10,
                   binary_fraction=binary_fraction,
                   binary_r_max=binary_r_max,
                   subsampling_factor=10,
                   guiding_error_mag1=guiding_error_mag1,
                   guiding_error_mag2=guiding_error_mag2,
                   guiding_error_angle=0.,
                   randomize_spectral_type=randomize_spectral_type,
                   use_cache=True,
                   parallelize=True)
        
    for i in range(num_fields):
        run_for_i(i+i_start)

if __name__ == "__main__":
    main(sys.argv)
