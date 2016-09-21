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

from astropy.io import fits

import numpy as np
from psf_testing import magic_values as mv
from psf_testing.moments.estimate_background import get_background_level


def make_stacks(stars, stack_size=int(2 * mv.default_weight_rmax + 1),
                weight_func=mv.default_prim_weight_func):
    
    stack_size = int(stack_size)

    stacks = {}

    for stack_name, image_name in zip(("star", "model", "noisy_model"),
                                      ("stamp", "model_psf", "noisy_model_psf")):
        stack = np.zeros((stack_size, stack_size))
        num_stars = 0
        for star in stars:
            if (not star.valid) or star.outlier:
                continue

            image = eval("star." + image_name)
            image_shape = np.shape(image)
            
            # Center it around the proper central pixel
            if "star" in stack_name:
                xc = int(star.xc+0.5)
                yc = int(star.yc+0.5)
            else:
                xc = int(eval("star."+stack_name+"_xc")+0.5)
                yc = int(eval("star."+stack_name+"_yc")+0.5)
            
            image_radius = np.min((xc,image_shape[0]-xc-1,yc,image_shape[1]-yc-1))
            image_view = image[xc-image_radius:xc+image_radius+1,yc-image_radius:yc+image_radius+1]
            
            image_view_shape = np.shape(image_view)

            dx = (image_view_shape[0] - stack_size) // 2
            dy = (image_view_shape[1] - stack_size) // 2
            
            if dx > 0:
                image_view = image_view[dx:-dx,:]
            if dy > 0:
                image_view = image_view[:,dy:-dy]
                
            image_view = image_view/star.m0[0]

            if (dx >= 0) and (dy >= 0):
                stack += image_view
            elif (dx < 0) and (dy < 0):
                stack[-dx:dx, -dy:dy] += image_view
            elif (dx >= 0) and (dy < 0):
                stack[:, -dy:dy] += image_view
            elif (dx < 0) and (dy >= 0):
                stack[-dx:dx, :] += image_view

            num_stars += 1

        assert num_stars > 0
        
        # Subtract off any lingering background
        stack -= get_background_level(stack)
        
        # Normalize the stack by its number of stars
        stack /= num_stars

        stacks[stack_name] = stack
        
    # Add residual and noisy residual stacks
    stacks["residual"] = stacks["star"] - stacks["model"]
    stacks["noisy_residual"] = stacks["star"] - stacks["noisy_model"]

    return stacks

def save_stacks(stacks, filename_root, header=None):

    for stack_name in stacks:
        filename = filename_root + "_" + stack_name + "_stack" + mv.image_extension

        hdu = fits.PrimaryHDU(stacks[stack_name])
        
        if header is not None:
            for key in header:
                hdu.header[key] = header[key]

        hdu.writeto(filename, clobber=True)

    return

def make_and_save_stacks(stars, filename_root, stack_size=(2 * mv.default_weight_rmax + 1), header=None):

    stacks = make_stacks(stars, stack_size)

    save_stacks(stacks, filename_root, header)

    return stacks
