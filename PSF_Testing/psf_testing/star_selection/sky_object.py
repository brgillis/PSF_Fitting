""" @file sky_object.py

    Created 15 Sep 2015

    A class representing a star or galaxy identified in an image.

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
from scipy.spatial import cKDTree

from psf_testing import magic_values as mv

class sky_object(object):
    """
        @TODO: classdocs
    """

    def __init__(self, cat_line=None):
        """
            @TODO: Constructor docs
        """
        
        if(cat_line is None):
        
            self.x_pix = None
            self.y_pix = None
            
            self.ra = None
            self.dec = None
            
            self.class_star = None
            
            self.mag = None
            self.flux = None
            
            self.stamp_size = mv.min_stamp_size
            
        else:
            
            properties = cat_line.split()
        
            self.x_pix = float(properties[1])
            self.y_pix = float(properties[2])
            
            self.ra = float(properties[3])
            self.dec = float(properties[4])
            
            self.class_star = float(properties[16])
            
            self.mag = float(properties[17])
            self.flux = float(properties[14])
            
            xp_size = float(properties[7]) - float(properties[5])
            yp_size = float(properties[8]) - float(properties[6])
        
            if(xp_size > yp_size):
                test_stamp_size = int(np.ceil(xp_size / 2.) + 1)
            else:
                test_stamp_size = int(np.ceil(yp_size / 2.) + 1)
            
            self.stamp_size = max((test_stamp_size,mv.min_stamp_size))
            
        self.lowest_separation = None
        
    def get_position_tuple(self):
        return (self.x_pix, self.y_pix)
        
    def get_lowest_separation(self,all_objects=None, objects_tree=None):
        
        # If we haven't been passed an objects list, simply use the cached value
        if(all_objects is None):
            assert(self.lowest_separation is not None)
            return self.lowest_separation
        
        # If we haven't been passed an objects_tree, generate one
        if(objects_tree is None):
            
            # Set up list of position tuples
            positions_tuples = []
            for obj in all_objects:
                positions_tuples.append(obj.get_position_tuple())
                
            # Create the cKDTree
            objects_tree = cKDTree(positions_tuples)
            
        # Query the cKDTree for the nearest neighbor
        query_result = objects_tree.query((self.get_position_tuple()),k=2)
        self.lowest_separation = query_result[0][1] * mv.pixel_scale
        
        return (self.lowest_separation, objects_tree)