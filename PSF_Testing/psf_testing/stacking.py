""" @file stacking.py

    Created 6 Nov 2015

    Functions used to stack star and model PSF images

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

from ast import literal_eval

from astropy.io import fits

import numpy as np
from psf_testing import magic_values as mv


def make_stacks(stars, stack_size=(2 * mv.default_weight_rmax + 1)):

    stacks = {}

    for stack_name, image_name in zip(("star", "model", "noisy_model"),
                                      ("image", "model_image", "noisy_model_image")):
        stack = np.zeros((stack_size, stack_size))
        num_stars = 0
        for star in stars:
            if not star.valid or star.outlier:
                continue

            image = literal_eval("star." + image_name)
            image_shape = np.shape(image)

            dx = image_shape[0] - stack_size / 2
            dy = image_shape[1] - stack_size / 2

            if (dx >= 0) and (dy >= 0):
                stack += image[dx:-dx, dy:-dy]
            elif (dx < 0) and (dy < 0):
                stack[-dx:dx, -dy:dy] += image
            elif (dx >= 0) and (dy < 0):
                stack[:, -dy:dy] += image[dx:-dx, :]
            elif (dx < 0) and (dy >= 0):
                stack[-dx:dx, :] += image[:, dy:-dy]

            num_stars += 1

        assert num_stars > 0
        stack /= num_stars

        stacks[stack_name] = stack

    return stacks

def save_stacks(stacks, filename_root):

    for stack_name in stacks:
        filename = filename_root + "_" + stack_name + "_stack" + mv.image_extension

        hdu = fits.PrimaryHDU(stacks[stack_name])

        hdu.writeto(filename)

    return

def make_and_save_stacks(stars, filename_root, stack_size=(2 * mv.default_weight_rmax + 1)):

    stacks = make_stacks(stars, stack_size)

    save_stacks(stacks, filename_root)

    return
