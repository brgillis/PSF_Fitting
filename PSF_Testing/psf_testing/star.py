""" star.py

    Created 16 Sep 2015

    A module containing a class which represents a star on an image.

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

class star(object):
    """
        @TODO classdocs
    """


    def __init__(self, sky_obj=None):
        """
            @TODO Constructor
        """
        
        self.sky_object = sky_obj
        
        self.stamp = None
        self.moments = None
        
    def get_position(self):
        return self.sky_object.get_position()