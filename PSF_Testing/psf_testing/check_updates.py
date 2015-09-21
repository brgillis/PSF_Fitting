""" @file check_updates.py

    Created 15 Sep 2015

    Functions to check if certain things need to be updated or not, based on the last
    time they were changed.

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

from os.path import getmtime
from time import time
import os.error

from psf_testing import magic_values as mv

def make_update_marker():
    """ Make an update marker so that any files which need updates will have to be
        recreated.
        
        Requires: (nothing)
        
        Returns: None
        
        Side-effects: Overwrites $mv.update_marker_filename    
    """
    
    with open(mv.update_marker_filename, "w") as fout:
        fout.write(str(time()))

def file_needs_update(filename):
    """ Check if a file needs to be updated, based on whether or not is was last
        updated more recently than the update marker.
        
        Requires: filename <string>
        
        Returns: <bool> (True if it needs to be updated, False otherwise)
    """
    
    # Get the last time the update marker was updated
    try:
        last_update_time = getmtime(mv.update_marker_filename)
    except os.error:
        # The update marker doesn't exist, so create it stamped with the current time
        make_update_marker()
        return True
    
    # Get the last modification time of the file
    try:
        file_last_updated = getmtime(filename)
    except os.error:
        # The file doesn't exist, so we'll have to create it
        return True
    
    # If the file was last modified since the update marker was updated, it needs to be updated
    if(file_last_updated < last_update_time):
        return True
    else:
        return False