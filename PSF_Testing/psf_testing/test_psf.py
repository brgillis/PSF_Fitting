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

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from astropy.io import fits
import os

from psf_testing import magic_values as mv
from psf_testing.image_info import get_chip, get_exp_time
from psf_testing.sextractor_utility import get_stars_in_image

def test_psf(image_filename,
             
             min_class_star = mv.default_min_class_star,
             min_star_mag = mv.default_min_star_mag,
             max_star_mag = mv.default_max_star_mag,
             min_lowest_separation = mv.default_min_lowest_separation,
             
             test_single_focus = False,
             test_focus = None,
             min_test_focus = mv.default_min_test_focus,
             max_test_focus = mv.default_max_test_focus,
             test_focus_samples = mv.default_focus_samples,
             test_focus_precision = mv.default_focus_precision,
             
             sex_data_path = mv.default_sex_data_path,
             cleanup_sex_files = True,
             
             tinytim_data_path = mv.default_tinytim_data_path,
             cleanup_tinytim_files = False,
             force_tinytim_update = False,
             **kwargs):
    
    # Start by inspecting the image and getting needed details about it
    image = fits.open(image_filename)[0]
    
    chip = get_chip(image)
    exp_time = get_exp_time(image)
    
    # Keep track of a list of files we'll want to cleanup when done
    files_to_cleanup = []
    
    # Get a list of the isolated stars in the image by running SExtractor on it
    stars = get_stars_in_image(image_filename = image_filename,
                               exp_time = exp_time,
                               min_class_star = min_class_star,
                               min_star_mag = min_star_mag,
                               max_star_mag = max_star_mag,
                               min_lowest_separation = min_lowest_separation,
                               sex_data_path = sex_data_path,
                               files_to_cleanup = files_to_cleanup,
                               cleanup_sex_files = cleanup_sex_files,
                               **kwargs)
    
    # If we're testing a single focus value, do that now
    if(test_single_focus):
        test_results = test_psf_for_focus(stars = stars,
                                          
                                          test_focus = test_focus,
                                          
                                          tinytim_data_path = tinytim_data_path,
                                          cleanup_tinytim_files = cleanup_tinytim_files,
                                          force_tinytim_update = force_tinytim_update,
                                          
                                          files_to_cleanup = files_to_cleanup,
                                          
                                          **kwargs)
    # Otherwise, call the fitting function
    else:
        test_results = fit_best_focus_and_test_psf(stars=stars,
                                                   
                                                   min_test_focus=min_test_focus,
                                                   max_test_focus=max_test_focus,
                                                   test_focus_samples=test_focus_samples,
                                                   test_focus_precision=test_focus_precision,
                                          
                                                   tinytim_data_path = tinytim_data_path,
                                                   cleanup_tinytim_files = cleanup_tinytim_files,
                                                   force_tinytim_update = force_tinytim_update,
                                                  
                                                   files_to_cleanup = files_to_cleanup,
                                          
                                                   **kwargs)
        
    # Do something with the results
    report_results(test_results,**kwargs)
    
    # Remove all files in the cleanup list
    for filename in files_to_cleanup:
        os.remove(filename)
            
    return