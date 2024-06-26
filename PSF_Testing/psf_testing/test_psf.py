""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing/test_psf.py

    Created 16 Sep 2015

    Primary function for testing psfs.

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

    You should have received a copy of the GNU General Pmvublic License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os

import numpy as np
from astropy.io import fits

from psf_testing import magic_values as mv
from utility.check_updates import make_update_marker
from psf_testing.extract_stamp import extract_stamp_for_star
from psf_testing.fit_focus import fit_best_focus_and_test_psf
from psf_testing.fit_params import fit_best_params_and_test_psf
from utility.smart_logging import get_default_logger
from psf_testing.moments.centre_image import centre_image
from psf_testing.moments.estimate_background import get_background_level_and_noise
from psf_testing.moments.get_Ms import get_m0_and_Ms
from psf_testing.report_results import report_results
from psf_testing.stacking import make_and_save_stacks
from psf_testing.star_selection.image_info import (get_chip, get_exp_time, get_ra, get_dec,
                                                   get_gain, get_obs_time)
from psf_testing.star_selection.sextractor_utility import get_stars_in_image
from psf_testing.test_psf_for_params import test_psf_for_params

def test_psf(image_filename,

             results_filename=None,
             results_dir=None,
             results_tag=None,

             min_class_star=mv.default_min_class_star,
             min_star_mag=mv.default_min_star_mag,
             max_star_mag=mv.default_max_star_mag,
             min_star_size = mv.default_min_star_size,
             max_star_size = mv.default_max_star_size,
             min_lowest_separation=mv.default_min_lowest_separation,
             min_star_snr=mv.default_min_star_snr,

             fit_all_params=False,
             focus_penalty_sigma=mv.default_focus_penalty_sigma,
             penalty_sigma=mv.default_penalty_sigma,
             
             focus=None,
             optical_params=mv.default_params,
             optical_param_slopes=mv.default_slopes,
             
             min_focus=mv.default_min_focus,
             max_focus=mv.default_max_focus,
             focus_samples=mv.default_focus_samples,
             focus_precision=mv.default_focus_precision,
             num_grid_points=mv.default_num_grid_points,

             sex_data_path=mv.default_sex_data_path,
             cleanup_sex_files=True,

             tinytim_path=mv.default_tinytim_path,
             tinytim_data_path=mv.default_tinytim_data_path,
             cleanup_tinytim_files=False,
             force_update=False,
             save_stacks=True,
             refresh_only=False,
             subsampling_factor=mv.default_subsampling_factor,
             galsim_rebin=mv.default_galsim_rebin,
             
             parallelize=False,
             
             norm_errors=False,
             
             seed=None):
        
    logger = get_default_logger()

    # Mark that we need an update if we're forcing an update
    if force_update:
        make_update_marker()

    filename_root = image_filename.replace(mv.image_extension, "")

    # Set up the filename root for the results
    
    if results_dir is None:
        results_root = filename_root
    else:
        results_root = os.path.join(results_dir,os.path.split(filename_root)[-1])
    
    if results_tag is not None:
        results_root += "_" + results_tag
        
    if seed is not None:
        results_root += "_" + str(seed)
        
    # If refreshing only, check if we can skip this
    if refresh_only:
        if os.path.isfile(results_root+mv.results_tail):
            logger.info(results_root+mv.results_tail + " already exists, so skipping.")
            return
    
    logger.info("Testing " + image_filename + ".")

    # Start by inspecting the image and getting needed details about it
    image = fits.open(image_filename)[0]

    chip = get_chip(image)
    exp_time = get_exp_time(image)
    gain = get_gain(image)
    obs_time = get_obs_time(image)
    
    # Keep track of a list of files we'll want to cleanup when done
    files_to_cleanup = []

    # Get a list of the isolated stars in the image by running SExtractor on it
    stars = get_stars_in_image(image_filename=image_filename,
                               exp_time=exp_time,
                               min_class_star=min_class_star,
                               min_star_mag=min_star_mag,
                               max_star_mag=max_star_mag,
                               min_star_size=min_star_size,
                               max_star_size=max_star_size,
                               min_lowest_separation=min_lowest_separation,
                               min_star_snr=min_star_snr,
                               sex_data_path=sex_data_path,
                               files_to_cleanup=files_to_cleanup,
                               cleanup_sex_files=cleanup_sex_files)

    # Set up the weight functions we'll use
    prim_weight_func = mv.default_prim_weight_func
    sec_weight_func = mv.default_sec_weight_func

    # Get general data on all stars
    for star in stars:

        # Extract the star's postage stamp if possible
        try:
            star.stamp = extract_stamp_for_star(star=star,
                                            image=image)
        except AssertionError as _e:
            star.valid = False
            continue

        background_level, star.background_noise = get_background_level_and_noise(star.stamp)

        star.stamp -= background_level

        try:
            star.xc, star.yc, star.x_array, star.y_array, star.prim_weight_mask, star.m0 = \
                centre_image(image=star.stamp, weight_func=prim_weight_func)
            
            # Check the star's centring is good enough
            x_shift = np.abs( star.xc - (np.shape(star.stamp)[0]-1)/2 )
            y_shift = np.abs( star.yc - (np.shape(star.stamp)[1]-1)/2 )
            max_shift = np.max((np.abs(x_shift),np.abs(y_shift)))
            if max_shift > 4:
                star.valid = False
                continue
        except AssertionError as _e:
            star.valid = False
            continue
        except Exception as e:
            if "cannot be centred" in str(e):
                star.valid = False
                continue
            else:
                raise
                

        (star.m0, star.Mxy, star.Mpcs) = \
            get_m0_and_Ms(image=star.stamp,
                          prim_weight_func=prim_weight_func,
                          sec_weight_func=sec_weight_func,
                          xc=star.xc,
                          yc=star.yc)

        star.chip = chip

        # If the monopole term is negative or zero, mark it as invalid
        if star.m0[0] <= 0:
            star.valid = False
            continue
        
    # Set up kwargs to pass along regardless of function called
    tinytim_params = {"tinytim_path":tinytim_path,
                      "tinytim_data_path":tinytim_data_path,
                      "subsampling_factor":subsampling_factor,
                      "chip":chip,
                      "galsim_rebin":galsim_rebin}
    
    kwargs = {"stars":stars,

              "image_filename":image_filename,
              "image":image,

              "num_grid_points":num_grid_points,

              "prim_weight_func":prim_weight_func,
              "sec_weight_func":sec_weight_func,

              "tinytim_params":tinytim_params,

              "gain":gain,
              "save_models":save_stacks,
              "files_to_cleanup":files_to_cleanup,
              
              "parallelize":parallelize,
              
              "norm_errors":norm_errors,
              "seed":seed}

    kwargs_and_params = kwargs.copy()
    # Are we fitting all params?
    kwargs_and_params.update(optical_params)
    
    if fit_all_params:
        # Start with a focus fit to get a good starting position
        test_results, _ = fit_best_focus_and_test_psf(focus_penalty_sigma=focus_penalty_sigma,

                                                      min_focus=min_focus,
                                                      max_focus=max_focus,
                                                      focus_samples=focus_samples,
                                                      focus_precision=focus_precision,
                
                                                      **kwargs_and_params)
        
        kwargs_and_params["focus"]=test_results[0]["focus"]
        
        # Proceed to the full fit
        test_results, fitting_record = fit_best_params_and_test_psf(focus_penalty_sigma=focus_penalty_sigma,
                                                                    penalty_sigma=penalty_sigma,
                                                                    **kwargs_and_params)
    # Are we testing a specific focus value?
    elif focus is not None:
        kwargs_and_params.update(optical_param_slopes)
        test_results = test_psf_for_params(focus=focus,
                                           **kwargs_and_params)
        fitting_record = None
    # Otherwise, call the focus fitting function
    else:
        kwargs_and_params.update(optical_param_slopes)
        test_results, fitting_record = fit_best_focus_and_test_psf(focus_penalty_sigma=focus_penalty_sigma,

                                                                  min_focus=min_focus,
                                                                  max_focus=max_focus,
                                                                  focus_samples=focus_samples,
                                                                  focus_precision=focus_precision,
                
                                                                  **kwargs_and_params)
    
    # Report the results
    report_results(test_results=test_results,
                   fitting_record=fitting_record,
                   chip=chip,
                   obs_time=obs_time,
                   exp_time=exp_time,
                   ra=get_ra(image),
                   dec=get_dec(image),
                   filename_root=results_root)

    # Save stacks if desired
    if save_stacks:
        make_and_save_stacks(stars=stars,
                             filename_root=results_root,
                             stack_size=(2 * mv.default_weight_rmax + 1),
                             header={"NSTAR":test_results[1][1],
                                     "FOCUS":test_results[0]["focus"]})

    # Remove all files in the cleanup list
    if cleanup_tinytim_files:
        for filename in files_to_cleanup:
            os.remove(filename)
    
    # Cleanup unneeded objects
    del stars, image
    import gc; gc.collect(2)
        
    logger.info("Finished analysing " + image_filename + ".")

    return
