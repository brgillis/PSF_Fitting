""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing/arg_parser.py

    Created 5 Oct 2016

    Provides a function to get a standard argument parser for test_psf.

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

import argparse
from psf_testing import magic_values as mv
from copy import deepcopy

def parse_and_get_kwargs(parser,special_keys=None):
    
    args = parser.parse_args()
    
    kwargs = deepcopy(vars(args))
    
    del (kwargs["focus_sample_x_points"],kwargs["focus_sample_y_points"])
    kwargs["num_grid_points"] = (args.focus_sample_x_points, args.focus_sample_y_points)
    
    # Split out kwargs we'll be handling differently
    special_kwargs = {}
    if special_keys is not None:
        for key in special_keys:
            special_kwargs[key] = kwargs[key]
            del kwargs[key]
            
    optical_params = {}
    for key in mv.default_params:
        optical_params[key] = kwargs[key]
        del kwargs[key]
            
    optical_param_slopes = {}
    for key in mv.default_params:
        slope_key = key+"_slope"
        optical_param_slopes[slope_key] = kwargs[slope_key]
        del kwargs[slope_key]
            
    kwargs["optical_params"] = optical_params
    kwargs["optical_param_slopes"] = optical_param_slopes
            
    return kwargs, special_kwargs

def get_arg_parser():
    
    parser = argparse.ArgumentParser()
    
    # Image filename
    parser.add_argument("--image_filename",type=str, default=None,
                        help="The filename of the image to test the PSF model on.")
    parser.add_argument("--image_list_filename",type=str, default=None,
                        help="The filename which contains a list of images to test.")
    parser.add_argument("--image_dir",type=str, default=None,
                        help="The directory where images can be found.")
    
    # What to do with results
    parser.add_argument("--results_dir",type=str, default=None,
                        help="The directory to store results data in.")
    parser.add_argument("--results_tag",type=str, default=None,
                        help="A tag to add to results filenames.")
    
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
    parser.add_argument("--min_focus", type=float, default=mv.default_min_focus,
                        help="Minimum focus value to test.")
    parser.add_argument("--max_focus", type=float, default=mv.default_max_focus,
                        help="Maximum focus value to test.")
    parser.add_argument("--focus_samples", type=int, default=mv.default_focus_samples,
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
    parser.add_argument("--subsampling_factor",type=int,default=mv.default_subsampling_factor,
                        help="The subsampling factor to use for TinyTim PSFs.")
    parser.add_argument("--galsim_rebin", action="store_true",
                        help="Use galsim to perform rebinning via interpolation. Will be slower but more accurate.")
    
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
    parser.add_argument("--disable_parallelization", action="store_true",
                        help="If set, will not use multiprocessing at any point.")
    
    # Optical parameters
    for optical_param in mv.default_params:
        parser.add_argument("--"+optical_param,type=float, default=mv.default_params[optical_param],
                            help="Optical parameter; will be ignored if fitting parameters.")
        parser.add_argument("--"+optical_param+"_slope",type=float, default=0,
                            help="Optical parameter slope against focus; will be ignored if fitting parameters.")
    
    return parser