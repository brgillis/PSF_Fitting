""" @file psf_model_scheme.py

    Created 21 Sep 2015

    Class which stores data for how we go about generating model PSFs.

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

class psf_model_scheme(object):
    """
        Class which stores data for how we go about generating model PSFs.
    """


    def __init__(self, focus=mv.default_test_focus,
                       num_grid_points=mv.default_num_grid_points,
                       image_shape=mv.default_image_shape):
        """ Optional: focus <float> (focus value for the model)
                      num_grid_points <float, float> (number of grid points to use for focus
                                                      determinations. If None, will use stars' exact
                                                      positions for models.)
                      image_shape <float, float> (Shape of image, using C ordering)
        
        """
        
        self.focus = focus
        
        if(num_grid_points is None):
            self.use_grid = False
        else:
            self.use_grid = True
            
            grid_nx = num_grid_points[0]
            grid_ny = num_grid_points[1]
            image_nx = image_shape[0]
            image_ny = image_shape[1]
            
            self.grid_stepx = float(image_nx)/grid_nx
            self.grid_stepy = float(image_ny)/grid_ny
            
    def get_position_to_use(self, xp, yp):
        
        if not self.use_grid:
            # If we're not using a grid, round to nearest int
            return ( int(xp+0.5), int(yp+0.5) )
        
        # Otherwise, determine the nearest grid point
        gx = int(float(xp)/self.grid_stepx)
        gy = int(float(yp)/self.grid_stepy)
        
        xp = int((gx + 0.5) * self.grid_stepx)
        yp = int((gy + 0.5) * self.grid_stepy)
        
        return xp, yp
        