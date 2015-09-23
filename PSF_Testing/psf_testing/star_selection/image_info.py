""" image_info.py

    Created 15 Sep 2015

    Functions to get specific information from a fits image

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

def get_chip(image):
    return image.header[mv.header_chip_keyword]

def get_exp_time(image):
    return image.header[mv.header_exp_time_keyword]

def get_gain(image):
    return image.header[mv.header_gain_keyword]

def get_obs_date(image):
    return image.header[mv.header_obs_date_keyword]

def get_obs_time(image):
    return image.header[mv.header_obs_time_keyword]