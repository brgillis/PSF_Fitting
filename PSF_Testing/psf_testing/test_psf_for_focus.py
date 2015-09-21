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

from psf_testing import magic_values as mv
from psf_testing.moments.centre_image import centre_image
from psf_testing.extract_stamp import extract_stamp_for_star

def test_psf_for_focus(stars,
                       
                       image_filename,
                       chip,
                       image = None,
                       
                       test_focus = mv.default_test_focus,
                                  
                       tinytim_data_path = mv.default_tinytim_data_path,
                       cleanup_tinytim_files = False,
                       force_tinytim_update = False,
                      
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
    
    # Set up the weight function we'll use
    weight_func = mv.default_weight_func
    
    # Loop through stars and get results for Each
    for star in stars:
        
        # Extract the star's postage stamp
        star.stamp = extract_stamp_for_star(star=star,
                                            image=image)
        
        star_xc, star_yc, star_x_array, star_y_array, star_weight_mask, star.m0 = \
            centre_image(image=image, weight_func=weight_func)
        
        psf_stamp = get_psf_for_star(star=star,
                                     chip=chip,
                                     star_xc=star_xc,
                                     star_yc=star_yc,
                                     star_m0=star_m0)