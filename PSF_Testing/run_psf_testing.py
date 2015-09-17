#!/usr/bin/env python

""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/Run_PSF_Testing.py

    Created 16 Sep 2015

    Main module for running PSF testing.

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

import sys
import argparse

from psf_testing import magic_values as mv
from psf_testing.test_psf import test_psf

def main(argv):
    """ @TODO main docstring
    """
    parser = argparse.ArgumentParser()
    
    # Image filename
    parser.add_argument("image_filename",type=str,
                        help="The filename of the image to test the PSF model on.")
    
    # Star selection
    parser.add_argument("--min_class_star", type=float, default=mv.default_min_class_star,
                        help="The minimum class_star value for an object to be categorized as a star.")
    parser.add_argument("--min_mag", type=float, default=mv.default_min_star_mag,
                        help="The minimum magnitude for stars to be tested.")
    parser.add_argument("--max_mag", type=float, default=mv.default_max_star_mag,
                        help="The minimum magnitude for stars to be tested.")
    parser.add_argument("--min_lowest_separation", type=float, default=None,
                        help="The minimum required separation of a star from any other object in arcsec.")
    
    # Focus fitting
    parser.add_argument("--focus", type=float, default=None,
                        help="If given, will only test this focus value. Otherwise will fit best focus.")
    parser.add_argument("--min_focus", type=float, default=mv.default_min_test_focus,
                        help="Minimum focus value to test.")
    parser.add_argument("--max_focus", type=float, default=mv.default_max_test_focus,
                        help="Maximum focus value to test.")
    parser.add_argument("--focus_samples", type=float, default=mv.default_focus_samples,
                        help="Initial number of focus values to test between min and max before trying " +
                             "to find the precise value.")
    parser.add_argument("--focus_precision", type=float, default=mv.default_focus_precision,
                        help="Desired precision for the fit focus value. eg. if 0.1, then the focus " +
                             "will be fit to +/- 0.1.")
    
    # SExtractor data/files
    parser.add_argument("--sex_data_path",type=str, default=mv.default_sex_data_path,
                        help="Path of data required by sextractor (ie. template .cfg files etc.).")
    parser.add_argument("--cleanup_sex_files", action="store_true",
                        help="Cleanup generated intermediate sextractor files after execution.")
    
    # TinyTim data/files
    parser.add_argument("--tinytim_data_path",type=str, default=mv.default_tinytim_data_path,
                        help="Path where PSFs generated by TinyTim will be stored.")
    parser.add_argument("--cleanup_tinytim_files", action="store_true",
                        help="Cleanup generated TinyTim PSFs after execution. Note that enabling this" +
                             " can greatly slow down repeated runs.")
    parser.add_argument("--force_tinytim_update", action="store_true",
                        help="Force update of TinyTim PSFs.")
    
    # Execute command-line parsing
    parser.parse_args(args=argv)
    
    if(parser.focus is not None):
        test_single_focus = True
    else:
        test_single_focus = False
    
    # Pass the cline-args to the test_psf function, which carries out the testing
    test_psf(image_filename = parser.image_filename,
             
             min_class_star = parser.min_class_star,
             min_star_mag = parser.min_star_mag,
             max_star_mag = parser.max_star_mag,
             min_lowest_separation = parser.min_lowest_separation,
             
             test_single_focus = test_single_focus,
             test_focus = parser.focus,
             min_test_focus = parser.min_focus,
             max_test_focus = parser.max_focus,
             test_focus_samples = parser.focus_samples,
             test_focus_precision = parser.test_focus_precision,
             
             sex_data_path = parser.sex_data_path,
             cleanup_sex_files = parser.cleanup_sex_files,
             
             tinytim_data_path = parser.tinytim_data_path,
             cleanup_tinytim_files = parser.cleanup_tinytim_files,
             force_tinytim_update = parser.force_tinytim_update)

if __name__ == "__main__":
    main(sys.argv)