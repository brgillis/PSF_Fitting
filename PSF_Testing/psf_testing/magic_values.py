""" @file magic_values.py

    Created 15 Sep 2015

    Magic values for the PSF_Testing project.

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

# Instrument/exposure details
gain = 2.0
zeropoint = 26.08 # F606W filter in ABmags from www.stsci.edu/hst/wfc3/phot_zp_lbn

# Image header keywords
header_chip_keyword = "CCDCHIP"
header_exp_time_keyword = "EXPTIME"
header_obs_date_keyword = "DATE-OBS"
header_obs_time_keyword = "TIME-OBS"

# Sextractor config/parameter filenames
sex_field_par_filename = "hst_acs_wfc.par"
sex_field_template_cfg_filename = "hst_acs_wfc.cfg"
sex_psf_par_filename = "psf.par"
sex_psf_template_cfg_filename = "psf.cfg"

# Sextractor cfg template tags
sex_template_cfg_path_tag = "REPLACEME_PATH"
sex_template_cfg_output_tag = "REPLACEME_OUTPUT_CAT"
sex_template_cfg_zeropoint_tag = "REPLACEME_MAG_ZEROPOINT"

# Other Sextractor values
sex_data_path = "./data"
