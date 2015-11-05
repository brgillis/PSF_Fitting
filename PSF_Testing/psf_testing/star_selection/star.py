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

        # Data on position, etc of the star
        self.sky_object = sky_obj

        # Extracted postage stamp and background noise estimate of this star
        self.stamp = None
        self.background_noise = None

        # Data determined from testing this star
        self.xc = None
        self.yc = None

        self.x_array = None
        self.y_array = None

        self.prim_weight_mask = None
        self.sec_weight_mask = None

        self.m0 = None
        self.m0_err = None

        self.Qs = None
        self.Q_vars = None

        self.valid = True
        self.outlier = False

        # Data for the model PSF corresponding to this star's position
        self.model_psf = None

        self.model_m0 = None
        self.model_m0_var = None

        self.model_Qs = None
        self.model_Q_vars = None

        # Data for the corresponding noisy model PSF
        self.noisy_model_psf = None

        self.noisy_m0 = None
        self.noisy_m0_var = None

        self.noisy_model_Qs = None
        self.noisy_model_Q_vars = None


    def get_position(self):
        return self.sky_object.get_position()
