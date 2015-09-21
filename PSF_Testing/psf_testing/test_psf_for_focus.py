""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing/test_psf_for_focus.py

    Created 21 Sep 2015

    Function to test the PSF model for a single focus value.

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
from psf_testing.get_model_psf import get_model_psf_for_star
from psf_testing.psf_model_scheme import psf_model_scheme

def test_psf_for_focus(stars,
                       
                       image_filename,
                       chip,
                       image = None,
                       
                       test_focus = mv.default_test_focus,
                       num_grid_points = mv.default_num_grid_points,
                       
                       weight_func = mv.default_weight_func,
                                  
                       tinytim_path = mv.default_tinytim_path,
                       tinytim_data_path = mv.default_tinytim_data_path,
                       cleanup_tinytim_files = False,
                      
                       files_to_cleanup = None,
                      
                       **kwargs):

    if(image is None):
        from astropy.io import fits
        image = fits.open(image_filename)[0]

    if(files_to_cleanup is None):
        files_to_cleanup = []
        
    # Initialize arrays for results per star tested
    star_m0s = []
    star_m0_errs = []
    star_Qs = []
    star_Q_errs = []
    psf_m0s = []
    psf_m0_errs = []
    comp_Zs = []
    
    # Get the image shape, and reverse its ordering to x,y in fits ordering
    image_shape = np.shape(image)
    image_shape = image_shape[1], image_shape[0]
    
    # Set up the focus generating scheme
    model_scheme = psf_model_scheme(focus=test_focus,
                                    num_grid_points=num_grid_points,
                                    image_shape=image_shape)
    
    # Loop through stars and get results for each
    for star in stars:
        
        if not star.valid:
            continue
        
        model_psf = get_model_psf_for_star(star=star,
                                           scheme=model_scheme,
                                           weight_func=weight_func,
                                           tinytim_path=tinytim_path,
                                           tinytim_data_path=tinytim_data_path)