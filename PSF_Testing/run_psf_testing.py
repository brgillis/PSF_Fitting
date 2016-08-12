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
from os.path import join
from multiprocessing import Pool, cpu_count
from copy import deepcopy

from psf_testing import magic_values as mv
from psf_testing.test_psf import test_psf
from psf_testing.smart_logging import get_default_logger

class test_psf_caller(object):
    def __init__(self,*args,**kwargs):
        self.args = args,
        self.kwargs = kwargs
    def __call__(self,x):
        try:
            test_psf(x,*self.args,**self.kwargs)
        except Exception as e:
            logger = get_default_logger()
            logger.error("Exception when processing file " + str(x) + ": " + str(e))
            if "usable stars" not in str(e):
                raise
            
def main(argv):
    """ @TODO main docstring
    """
    parser = argparse.ArgumentParser()
    
    # Image filename
    parser.add_argument("--image_filename",type=str, default=None,
                        help="The filename of the image to test the PSF model on.")
    parser.add_argument("--image_list_filename",type=str, default=None,
                        help="The filename which contains a list of images to test.")
    parser.add_argument("--image_dir",type=str, default=None,
                        help="The directory where images can be found.")
    
    # Star selection
    parser.add_argument("--min_class_star", type=float, default=mv.default_min_class_star,
                        help="The minimum class_star value for an object to be categorized as a star.")
    parser.add_argument("--min_star_mag", type=float, default=mv.default_min_star_mag,
                        help="The minimum magnitude for stars to be tested.")
    parser.add_argument("--max_star_mag", type=float, default=mv.default_max_star_mag,
                        help="The minimum magnitude for stars to be tested.")
    parser.add_argument("--min_star_size", type=float, default=mv.default_min_star_size,
                        help="The minimum size for stars to be tested (FWHM in arcsec).")
    parser.add_argument("--max_star_size", type=float, default=mv.default_max_star_size,
                        help="The minimum size for stars to be tested (FWHM in arcsec).")
    parser.add_argument("--min_lowest_separation", type=float, default=None,
                        help="The minimum required separation of a star from any other object in arcsec.")
    parser.add_argument("--min_star_snr", type=float, default=mv.default_min_star_snr,
                        help="The minimum signal-to-noise ratio for this star to be tested")
    
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
    parser.add_argument("--focus_sample_x_points", type=int, default=mv.default_num_grid_points[0],
                        help="Number of sample points on the image to get model foci for " +
                             "along the x axis.")
    parser.add_argument("--focus_sample_y_points", type=int, default=mv.default_num_grid_points[1],
                        help="Number of sample points on the image to get model foci for " +
                             "along the y axis.")
    
    # Parameter fitting
    parser.add_argument("--fit_all_params", action="store_true",
                        help="If true, will also fit astigmatism parameters. This will greatly increase runtime.")
    parser.add_argument("--focus_penalty_sigma", type=float, default=mv.default_focus_penalty_sigma,
                        help="How lenient to be in letting focus vary from default values. If 0, will impose no penalty.")
    parser.add_argument("--penalty_sigma", type=float, default=mv.default_penalty_sigma,
                        help="How lenient to be in letting params vary from default values. If 0, will impose no penalty.")
    
    # SExtractor data/files
    parser.add_argument("--sex_data_path",type=str, default=mv.default_sex_data_path,
                        help="Path of data required by sextractor (ie. template .cfg files etc.).")
    parser.add_argument("--cleanup_sex_files", action="store_true",
                        help="Cleanup generated intermediate sextractor files after execution.")
    
    # TinyTim data/files
    parser.add_argument("--tinytim_path",type=str, default=mv.default_tinytim_path,
                        help="Path to the TinyTim executable.")
    parser.add_argument("--tinytim_data_path",type=str, default=mv.default_tinytim_data_path,
                        help="Path where PSFs generated by TinyTim will be stored.")
    parser.add_argument("--cleanup_tinytim_files", action="store_true",
                        help="Cleanup generated TinyTim PSFs after execution. Note that enabling this" +
                             " can greatly slow down repeated runs.")
    
    parser.add_argument("--force_update", action="store_true",
                        help="Force update of Sextractor catalogues and TinyTim PSFs.")
    
    parser.add_argument("--logging_level", type=str, default=mv.default_logging_level,
                        help="Level of logging info to display. Default 'info'.")
    
    # Other options
    parser.add_argument("--seed", type=int, default=None,
                        help="Base seed for random number generation.")
    parser.add_argument("--refresh_only", action="store_true",
                        help="Only run for images if the results file doesn't already exist.")
    parser.add_argument("--norm_errors", action="store_true",
                        help="Normalize errors - X^2 will no longer weight differently depending on scatter.")
    
    # Execute command-line parsing
    args = parser.parse_args()
    
    if(args.focus is not None):
        test_single_focus = True
    else:
        test_single_focus = False
        
    logger = get_default_logger()
    logger.setLevel(args.logging_level.upper())
    
    # Check if we're debugging
    try:
        import pydevd as unused_import
        debugging = True
    except ImportError:
        debugging = False
        
    # Set up kwargs used for all versions
    kwargs = deepcopy(vars(args))
    
    del (kwargs["focus_sample_x_points"],kwargs["focus_sample_y_points"])
    kwargs["num_grid_points"] = (args.focus_sample_x_points, args.focus_sample_y_points)
    
    del (kwargs["image_filename"], kwargs["image_list_filename"], kwargs["image_dir"], kwargs["logging_level"])
    
    # Pass the cline-args to the test_psf function, which carries out the testing
    
    if args.image_filename is not None:
        test_psf_caller(
                 # parallelize = (not debugging),
                 parallelize = (not debugging) and (not args.fit_all_params),
                 **kwargs )(args.image_filename)
        image_filenames = [args.image_filename]
    else:
        # We'll test a list of images
        image_filenames = []
        with open(args.image_list_filename) as fi:
            for line in fi:
                for word in line.split():
                    image_filename = join(args.image_dir, word)
                    image_filenames.append(image_filename)
        
        if args.refresh_only:
            nproc = 1
            parallelize_each_image = not debugging
        else:
            nproc = max((cpu_count()-1,1))
            parallelize_each_image = False
        
        if nproc == 1:
            for image_filename in image_filenames:
                test_psf_caller(parallelize = parallelize_each_image,
                                **kwargs )(image_filename)
        else:
        
            pool = Pool(processes=nproc,maxtasksperchild=1)
            _ = pool.map(test_psf_caller(parallelize = False,
                                         **kwargs ),
                     image_filenames,chunksize=1)
            pool.close()
            pool.join()
    
    logger.info("Execution complete.")

if __name__ == "__main__":
    main(sys.argv)
