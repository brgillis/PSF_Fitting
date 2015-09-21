""" @file get_model_psf.py

    Created 21 Sep 2015

    Function to generate a model psf for a given star.

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
from psf_testing.check_updates import file_needs_update

def get_model_psf_for_star(star,
                           scheme,
                           tinytim_data_path = mv.default_tinytim_data_path):
    """ Gets a model psf for a given star and the chip it was detected on
    
        Requires: star <star>
                  scheme <psf_model_scheme> (plan for generating model psfs)
        Optional: tinytim_data_path <string> (location where TinyTim data will be stored)
                  
        Returns: model_psf <star> (with m_err, Q values, and Q errors not yet determined)
    """
    
    # Get the position we'll generate the model PSF for
    psf_position = scheme.get_position_to_use( star.sky_object.x_pix, star.sky_object.y_pix )
    
    # Determine the name for the subsampled model PSF file
    subsampled_name = tinytim_data_path + "subsampled_psf_x-" + psf_position[0] + \
                        "_y-" + psf_position[1] + "_f-" + scheme.focus + \
                        "_c-" + star.chip + mv.image_extension 
                        
    # Check if we need to update this file, or if we can reuse the existing version
    if(file_needs_update(subsampled_name)):
        
        # We'll need to update it, so we'll call TinyTim to generate a PSF model
        make_subsampled_psf_model(filename = subsampled_name,
                                  xp = psf_position[0],
                                  yp = psf_position[1],
                                  focus = scheme.focus,
                                  chip = star.chip)