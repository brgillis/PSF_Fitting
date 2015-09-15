""" @file sextractor_utility.py

    Created 15 Sep 2015

    @TODO: File docstring

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

from psf_testing import magic_values as mv
from psf_testing.io import replace_multiple_in_file

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
                  exp_time,
                  template_filename=mv.sex_field_template_cfg_filename,
                  data_path=mv.sex_data_path):
    """ Makes a cfg file for SExtractor, using a given template file.
    
        Requires: output_cfg_filename <string>
                  output_catalog_filename <string> (The file SExtractor should output the cat to)
                  exp_time <float> (The exposure time for the image to be analyzed)
        Optional: template_filename <string> (The name of the template to use)
                  data_path <string> (The location where needed SExtractor files are.)
        
        Returns: None
    """
    
    mag_zeropoint = str(get_mag_zeropoint(exp_time, mv.zeropoint))
    
    replace_multiple_in_file(input_filename=template_filename,
                             output_filename=output_cfg_filename,
                             input_strings=[mv.sex_template_cfg_output_tag,
                                            mv.sex_template_cfg_path_tag,
                                            mv.sex_template_cfg_zeropoint_tag],
                             output_strings=[output_catalog_filename,
                                             data_path,
                                             mag_zeropoint])
    
    return
        