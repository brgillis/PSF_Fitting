""" @file sextractor_utility.py

    Created 15 Sep 2015

    Various functions for setting up and using sextractor.

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

from math import log10
from os.path import isfile, join
import subprocess as sbp

from astropy.io import fits

from psf_testing import magic_values as mv
from psf_testing.io import replace_multiple_in_file
from psf_testing.star_selection.image_info import get_exp_time
from psf_testing.star_selection.star_selection import get_objects_from_cat, get_stars, get_isolated_stars
from psf_testing.check_updates import file_needs_update

def get_mag_zeropoint(exp_time,instrument_zeropoint=mv.zeropoint):
    """ Gets the magnitude zeropoint for an exposure from the exposure time and the instrument's
        zeropoint.
        
        Requires: exp_time <float> (Exposure time)
        Optional: instrument_zeropoint (The instrument's zeropoint. Defaults to the zeropoint of
                    HST's WFC (...) if not given)
        Returns: <float> (The magnitude zeropoint for this exposure)
    """
    
    return instrument_zeropoint + 2.5 * log10(exp_time)

def make_cfg_file(output_cfg_filename,
                  output_catalog_filename,
                  exp_time=None,
                  mag_zeropoint=None,
                  template_filename=mv.sex_field_template_cfg_filename,
                  data_path=mv.default_sex_data_path,
                  gain="2.0",
                  pixel_scale="0"):
    """ Makes a cfg file for SExtractor, using a given template file.
    
        Requires: output_cfg_filename <string>
                  output_catalog_filename <string> (The file SExtractor should output the cat to)
                  exp_time <float> (The exposure time for the image to be analyzed)
        Optional: template_filename <string> (The name of the template to use)
                  data_path <string> (The location where needed SExtractor files are.)
        
        Returns: None
    """
    
    if(mag_zeropoint is None):
        mag_zeropoint = get_mag_zeropoint(exp_time, mv.zeropoint)
    
    replace_multiple_in_file(input_filename=join(data_path,template_filename),
                             output_filename=output_cfg_filename,
                             input_strings=[mv.sex_template_cfg_output_tag,
                                            mv.sex_template_cfg_path_tag,
                                            mv.sex_template_cfg_zeropoint_tag,
                                            mv.sex_template_cfg_gain_tag,
                                            mv.sex_template_cfg_pixel_scale_tag],
                             output_strings=[output_catalog_filename,
                                             data_path,
                                             str(mag_zeropoint),
                                             str(gain),
                                             str(pixel_scale)])
    
    return

def run_sextractor(image_filename, sex_cfg_filename, sex_cat_name=None):
    """ Runs sextractor on a given image with a given config file.
    
        Requires: image_filename <string>
                  sex_cfg_filename <string>
                  
        Returns: None
        
        Side-effects: Runs SExtractor, which means the catalog in the config file
            will be overwritten.
    """
    
    if(sex_cat_name is not None):
        # Delete any old catalog if it exists
        cmd = "rm -f " + sex_cat_name
        sbp.call(cmd, shell=True)
    
    # Call SExtractor
    cmd = "sex " + image_filename + " -c " + sex_cfg_filename
    sbp.call(cmd, shell=True)
    
    if(sex_cat_name is not None):
        # Check that the catalog was successfully created
        assert(isfile(sex_cat_name))
    
    return

def get_stars_in_image(image_filename,
                       exp_time=None,
                       min_class_star=mv.default_min_class_star,
                       min_star_mag=mv.default_min_star_mag,
                       max_star_mag=mv.default_max_star_mag,
                       min_lowest_separation=None,
                       sex_data_path=mv.default_sex_data_path,
                       files_to_cleanup=None,
                       cleanup_sex_files=False):
    
    if(files_to_cleanup is None):
        files_to_cleanup = []
    
    # Get the root of the filename
    image_filename_root = image_filename.replace(mv.image_extension,"")
    
    # Get the desired name of the sex cfg file
    sex_cfg_filename = image_filename_root + mv.sex_cfg_tail
    
    # Get the desired name of the sex cat file
    sex_cat_filename = image_filename_root + mv.sex_cat_tail
    
    if(cleanup_sex_files):
        files_to_cleanup.append(sex_cfg_filename)
        files_to_cleanup.append(sex_cat_filename)
    
    # Get the exposure time from the image if necessary
    if exp_time is None:
        exp_time = get_exp_time(fits.open(image_filename))
    
    # Run SExtractor if needed
    if(file_needs_update(sex_cat_filename)):
        # Make the cfg file for running SExtractor
        make_cfg_file(output_cfg_filename=sex_cfg_filename,
                      output_catalog_filename=sex_cat_filename,
                      exp_time=exp_time,
                      template_filename=mv.sex_field_template_cfg_filename,
                      data_path=sex_data_path)
        
        # Call SExtracter
        run_sextractor(image_filename, sex_cfg_filename, sex_cat_name=sex_cat_filename)
    
    # Get all objects in the image
    objects_in_image = get_objects_from_cat(sex_cat_filename)
    
    # Get a list of those which are stars in the right magnitude range
    stars = get_stars(objects_in_image,min_class_star,min_star_mag,max_star_mag)
    
    # Get only isolated stars
    isolated_stars = get_isolated_stars(stars,objects_in_image,min_lowest_separation)
    
    return isolated_stars
    