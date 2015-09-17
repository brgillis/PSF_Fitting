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

def test_psf(image_filename,cleanup=False,**kwargs):
    
    # Start by inspecting the image and getting needed details about it
    image = fits.open(image_filename)[0]
    
    chip = get_chip(image)
    exp_time = get_exp_time(image)
    
    # Keep track of a list of files we'll want to cleanup when done
    files_to_cleanup = []
    
    # Get a list of the isolated stars in the image by running SExtractor on it
    stars = get_stars_in_image(image_filename, exp_time, files_to_cleanup, **kwargs)
    
    # ...
    
    # If we're cleaning up, remove all files in the cleanup list
    if(cleanup):
        for filename in files_to_cleanup:
            os.remove(filename)
            
    return