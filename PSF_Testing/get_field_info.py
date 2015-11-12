#!/usr/bin/env python

""" @file /get_field_info.py

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

import click
from psf_testing.smart_logging import get_default_logger
from astropy.io import ascii
from astropy.io import fits
from os.path import join

from psf_testing import magic_values as mv

@click.command()
@click.option("--fields_list",default="data/HST_ACS_WFC_no_drizzle_fields.txt",
              help="List of fields to get details for.")
@click.option("--fields_path",default="/disk2/brg/Data/HST_Fields",
              help="Location of fields to get details for.")
@click.option("--output_filename",default="results/fields_info.dat")
def main(**kwargs):
    """ @TODO main docstring
    """
    
    logger = get_default_logger()
    
    logger.debug("Entering get_field_info main function.")
    
    # Get the fields we want to look at
    fields_list = []
    with open(kwargs['fields_list'], 'r') as fields_list_file:

        # Read in the file, except for comment lines
        for line in fields_list_file:
            fields_list.append(line.strip())
            
    # Initialize lists of data for each thing we're interested in
    obs_dates = [] # Observation dates
    obs_times = [] # Observation times
    exp_times = [] # Exposure times
    filter1s = [] # First filter
    filter2s = [] # Second filter
    chips = [] # Chip
    
    # Loop through fields and get data for each field
    for field_filename in fields_list:
        full_field_filename = join(kwargs['fields_path'],field_filename)
        
        # Get the FITS header
        field_header = fits.open(full_field_filename)[0].header
        
        obs_dates.append(field_header[mv.header_obs_date_keyword])
        obs_times.append(field_header[mv.header_obs_time_keyword])
        exp_times.append(field_header[mv.header_exp_time_keyword])
        filter1s.append(field_header['FILTER1'])
        filter2s.append(field_header['FILTER2'])
        chips.append(field_header[mv.header_chip_keyword])
        
    # Create a table of the data to be output
    ascii.write([obs_dates, obs_times, exp_times, filter1s, filter2s, chips], kwargs['output_filename'],
                names=['observation_date', 'observation_time', 'exposure_time', 'first_filter',
                       'second_filter', 'chip'],
                format="commented_header")
    
    logger.info("Exiting get_field_info main function.")
        

if __name__ == "__main__":
    main()
