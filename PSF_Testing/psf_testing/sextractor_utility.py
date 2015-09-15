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

import magic_values as mv

def get_mag_zeropoint(exp_time,instrument_zeropoint=mv.zeropoint):
    """ Gets the magnitude zeropoint for an exposure from the exposure time and the instrument's
        zeropoint.
        
        Requires: exp_time <float> (Exposure time)
        Optional: instrument_zeropoint (The instrument's zeropoint. Defaults to the zeropoint of
                    HST's WFC (...) if not given..)
        Returns: <float> The magnitude zeropoint for this exposure.
    """
    
    return instrument_zeropoint + 2.5 * log10(exp_time)

def make_cfg_file(output_filename,
                  exp_time,
                  template_filename=mv.sex_field_template_cfg_filename,
                  data_path=mv.sex_data_path):
    
    mag_zeropoint = str(get_mag_zeropoint(exp_time, mv.zeropoint))
    
    with open(template_filename, 'r') as template_file:
        with open(output_filename, 'w') as output_file:
            for line in template_file:
                output_file.write(line.replace('A', 'Orange'))
        