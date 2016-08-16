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
table_extension = ".fits"
sex_cfg_tail = "_sex.cfg"
sex_cat_tail = "_objects.cat"

results_tail = "_results" + table_extension

# Specific filenames
update_marker_filename = ".PSF_Testing_update_marker"

# Instrument/exposure details
pixel_scale = 0.05
gain = 2.0
zeropoint = 26.08 # F606W filter in ABmags from www.stsci.edu/hst/wfc3/phot_zp_lbn

# Image header keywords
header_chip_keyword = "CCDCHIP"
header_exp_time_keyword = "EXPTIME"
header_gain_keyword = "CCDGAIN"
header_obs_date_keyword = "DATE-OBS"
header_obs_time_keyword = "TIME-OBS"
header_ra_keyword = "RA_TARG"
header_dec_keyword = "DEC_TARG"

# Sextractor config/parameter filenames
sex_cmd = "sextractor "
sex_field_par_filename = "hst_acs_wfc.par"
sex_field_template_cfg_filename = "hst_acs_wfc.cfg"
sex_psf_par_filename = "psf.par"
sex_psf_template_cfg_filename = "psf.cfg"

# Sextractor cfg template tags
sex_template_cfg_path_tag = "REPLACEME_PATH"
sex_template_cfg_output_tag = "REPLACEME_OUTPUT_CAT"
sex_template_cfg_zeropoint_tag = "REPLACEME_MAG_ZEROPOINT"
sex_template_cfg_gain_tag = "REPLACEME_GAIN"
sex_template_cfg_pixel_scale_tag = "REPLACEME_PIXEL_SCALE"

# Other Sextractor values
default_sex_data_path = "./data"

# Default values for selecting stars
default_min_class_star = 0.95
default_min_star_mag = 21.0
default_max_star_mag = 27.0
default_min_star_size = 0.8
default_max_star_size = 5.0
default_min_lowest_separation = 1.0
default_min_star_snr = 50.0

# Values for making objects from a catalog line
min_stamp_size = 10
sex_cat_xp_col = 1
sex_cat_yp_col = 2
sex_cat_ra_col = 3
sex_cat_dec_col = 4
sex_cat_class_star_col = 16
sex_cat_mag_col = 17
sex_cat_flux_col = 20
sex_cat_flux_err_col = 21
sex_cat_xp_min_col = 5
sex_cat_xp_max_col = 7
sex_cat_yp_min_col = 6
sex_cat_yp_max_col = 8
sex_cat_num_cols = 24

# Focus fitting values
default_focus = 0.0
default_init_focus = -1.0
default_min_focus = -6.0
default_max_focus = 6.0
default_focus_samples = 4
default_focus_precision = 0.1
default_num_grid_points = (8, 4) # (x,y) in fits ordering, (y,x) in C ordering

# Param fitting values
default_params = {"z2":0.,
                  "z3":0.,
                  "astigmatism_0":0.031,
                  "astigmatism_45":0.028,
                  "coma_x":0.003,
                  "coma_y":0.001,
                  "clover_x":0.008,
                  "clover_y":0.018,
                  "spherical_3rd":-0.025,
                  "z12":0.,
                  "z13":0.,
                  "z14":0.,
                  "z15":0.,
                  "z16":0.,
                  "z17":0.,
                  "z18":0.,
                  "z19":0.,
                  "z20":0.,
                  "z21":0.,
                  "spherical_5th":0.009,
                  "kernel_adjustment":1.,
                  }
default_penalty_sigma = 0.005
default_focus_penalty_sigma = 10

default_image_shape = (4096, 2048) # (x,y) in fits ordering, (y,x) in C ordering
default_logging_level = "info"

# Default weight function for measuring star/model moments
default_weight_sigma = 3.0 # Pixels
if default_min_lowest_separation is not None:
    default_weight_rmax = default_min_lowest_separation / (2*pixel_scale)
else:
    default_weight_rmax = 10.0
def default_prim_weight_func(x, y):
    r2 = np.square(x) + np.square(y)
    return np.where(r2 > np.square(default_weight_rmax),0.0,
                    np.exp(-r2 / (2. * np.square(default_weight_sigma))))
def default_sec_weight_func(x, y):
    r2 = np.square(x) + np.square(y)
    return np.where(r2 > np.square(default_weight_rmax),0.0,1.0)

if(min_stamp_size<default_weight_rmax):
    min_stamp_size = default_weight_rmax

# TinyTim values
default_tinytim_data_path = "/disk2/brg/Data/HST_Fields/PSF_models"
default_tinytim_path = "/home/brg/Program_Files/tinytim-7.5"
default_model_psf_width = 2.0
default_model_psf_spec_type = (1, 15) # K-type star 15
default_detector = 15 # WFC
default_filter = "f606w"
default_subsampling_factor = 6

undistorted_model_tail = "00_psf.fits"
subsampled_model_tail = "00.fits"

charge_diffusion_kernel = np.array([])

# Header values to store in the subsampled model psf's header
ss_model_xc_label = "XP_CEN"
ss_model_yc_label = "YP_CEN"
ss_model_m0_label = "M0"
