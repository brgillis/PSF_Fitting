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

import numpy as np

# Misc runtime paramters
max_loop_iterations = 100

# Filename conventions
image_extension = ".fits"
sex_cfg_tail = "_sex.cfg"
sex_cat_tail = "_objects.cat"

# Specific filenames
update_marker_filename = ".PSF_Testing_update_marker"

# Instrument/exposure details
pixel_scale = 0.05
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
default_sex_data_path = "./data"

# Default values for selecting stars
default_min_class_star = 0.95
default_min_star_mag = 21.0
default_max_star_mag = 27.0
default_min_lowest_separation = None

# Values for making objects from a catalog line
min_stamp_size = 20
sex_cat_xp_col = 1
sex_cat_yp_col = 2
sex_cat_ra_col = 3
sex_cat_dec_col = 4
sex_cat_class_star_col = 16
sex_cat_mag_col = 17
sex_cat_flux_col = 14
sex_cat_xp_min_col = 5
sex_cat_xp_max_col = 7
sex_cat_yp_min_col = 6
sex_cat_yp_max_col = 8
sex_cat_num_cols = 24

# Focus fitting values
default_test_focus = 0.0
default_min_test_focus = -6.0
default_max_test_focus = 6.0
default_focus_samples = 7
default_focus_precision = 0.05

# Default weight function for measuring star/model moments
default_weight_sigma = 3.0 # Pixels
def default_weight_func(x,y):
    np.exp(-(np.square(x)+np.square(y)) / (2. * np.square(default_weight_sigma)))

# TinyTim values
default_tinytim_data_path = "/disk2/brg/Data/HST_Fields/PSF_models"
