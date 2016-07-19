""" @file get_model_psf.py

    Created 21 Sep 2015

    Function to generate a model psf for a given star.

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

import time
import os

from astropy.io import fits

import numpy as np
from psf_testing.function_cache import lru_cache
from psf_testing import magic_values as mv
from psf_testing.check_updates import file_needs_update
from psf_testing.io import replace_multiple_in_file
from psf_testing.moments.centre_image import centre_image
from psf_testing.rebin_psf import rebin
import subprocess as sbp


def make_subsampled_psf_model(filename,
                              chip=1,
                              psf_size=mv.default_model_psf_width,
                              xp=mv.default_image_shape[0]//2,
                              yp=mv.default_image_shape[1]//2,
                              focus=0.0,
                              detector=mv.default_detector,
                              filter=mv.default_filter,
                              spec_type=mv.default_model_psf_spec_type,
                              subsampling_factor=mv.default_subsampling_factor,
                              weight_func=mv.default_prim_weight_func,
                              tinytim_path=mv.default_tinytim_path,
                              astigmatism_0=None,
                              astigmatism_45=None,
                              coma_x=None,
                              coma_y=None,
                              clover_x=None,
                              clover_y=None,
                              spherical_3rd=None,
                              z12=None,
                              z13=None,
                              z14=None,
                              z15=None,
                              z16=None,
                              z17=None,
                              z18=None,
                              z19=None,
                              z20=None,
                              z21=None,
                              spherical_5th=None,
                              shape=None,
                              files_to_cleanup=None,
                              use_cache=True):
    """ Generates a subsampled psf model for a given image position, focus, and chip.

        Requires: filename <string>
                  xp <int>
                  yp <int>
                  focus <float>
                  chip <int>
                  weight_func <function of x,y> (Used for determining the centre of the
                                                 subsampled model)
                  tinytim_path <string> location of the TinyTim executable

        Returns: None

        Side-effects: Overwrites $filename if it exists with the new subsampled PSF model
    """

    filename_base = filename.replace(mv.image_extension, "")
    
    pid = str(os.getpid())
    
    process_filename_base = filename_base + pid
    process_filename = process_filename_base + mv.image_extension
    
    # Check for a lock file on this, to make sure its safe to create
    lock_filename = process_filename_base + ".lock"
    if os.path.isfile(lock_filename):
        time_start = time.time()
        while (not os.path.isfile(process_filename)) and (time.time() - time_start < 60):
            time.sleep(1)
        
        # Give the file a chance to be fully written
        time.sleep(1)
        
        if os.path.isfile(process_filename):
            if use_cache:
                try:
                    return fits.open(process_filename)[0]
                except IOError as _e:
                    open(lock_filename, 'a').close()
            else:
                open(lock_filename, 'a').close()
        else:
            # Looks like a rogue lock, so go ahead and use it ourselves
            open(lock_filename, 'a').close()
    else:
        open(lock_filename, 'a').close()
    
    par_file = process_filename_base + ".par"

    # Set up the command to call tiny1
    cmd = "export TINYTIM=" + tinytim_path + "\n" + \
          tinytim_path + "/tiny1 " + par_file + " << EOF \n" + \
          str(detector) + "\n" + \
          str(chip) + "\n" + \
          str(xp) + " " + str(yp) + "\n" + \
          str(filter) + "\n" + \
          str(spec_type[0]) + "\n" + \
          str(spec_type[1]) + "\n" + \
          str(psf_size) + "\n" + \
          str(focus) + "\n" + \
          process_filename_base + "\nEOF"
    # Run the command to call tiny1
    sbp.call(cmd, shell=True)
    
    # Edit the parameter file to adjust coma and astigmatism if necessary
    str_to_replace = []
    replacements = []
    if astigmatism_0 is not None:
        str_to_replace.append("0.031    # Z5 = 0 degree astigmatism")
        replacements.append(str(astigmatism_0) + "    # Z5 = 0 degree astigmatism")
    if astigmatism_45 is not None:
        str_to_replace.append("0.028    # Z6 = 45 degree astigmatism")
        replacements.append(str(astigmatism_45) + "    # Z6 = 45 degree astigmatism")
    if coma_x is not None:
        str_to_replace.append("0.003    # Z7 = X (V2) coma")
        replacements.append(str(coma_x) + "    # Z7 = X (V2) coma")
    if coma_y is not None:
        str_to_replace.append("0.001    # Z8 = Y (V3) coma")
        replacements.append(str(coma_y) + "    # Z8 = Y (V3) coma")
    if clover_x is not None:
        str_to_replace.append("0.008    # Z9 = X clover")
        replacements.append(str(clover_x) + "    # Z7 = X (V2) clover")
    if clover_y is not None:
        str_to_replace.append("0.018    # Z10 = Y clover")
        replacements.append(str(clover_y) + "    # Z8 = Y (V3) clover")
    if spherical_3rd is not None:
        str_to_replace.append("-0.025    # Z11 = 3rd order spherical")
        replacements.append(str(spherical_3rd) + "    # Z11 = 3rd order spherical")
    if z12 is not None:
        str_to_replace.append("0.       # Z12 = 0 degree Spherical astigmatism")
        replacements.append(str(z12) + "       # Z12 = 0 degree Spherical astigmatism")
    if z13 is not None:
        str_to_replace.append("0.       # Z13 = 45 degree Spherical astigmatism")
        replacements.append(str(z13) + "       # Z13 = 45 degree Spherical astigmatism")
    if z14 is not None:
        str_to_replace.append("0.       # Z14 = X Ashtray")
        replacements.append(str(z14) + "       # Z14 = X Ashtray")
    if z15 is not None:
        str_to_replace.append("0.       # Z15 = Y Ashtray")
        replacements.append(str(z15) + "       # Z15 = Y Ashtray")
    if z16 is not None:
        str_to_replace.append("0.       # Z16")
        replacements.append(str(z16) + "       # Z16")
    if z17 is not None:
        str_to_replace.append("0.       # Z17")
        replacements.append(str(z17) + "       # Z17")
    if z18 is not None:
        str_to_replace.append("0.       # Z18")
        replacements.append(str(z18) + "       # Z18")
    if z19 is not None:
        str_to_replace.append("0.       # Z19")
        replacements.append(str(z19) + "       # Z19")
    if z20 is not None:
        str_to_replace.append("0.       # Z20")
        replacements.append(str(z20) + "       # Z20")
    if z21 is not None:
        str_to_replace.append("0.       # Z21")
        replacements.append(str(z21) + "       # Z21")
    if spherical_5th is not None:
        str_to_replace.append("0.009    # Z22 = 5th order spherical")
        replacements.append(str(spherical_5th) + "    # Z22 = 5th order spherical")
        
    if len(str_to_replace) > 0:
        replace_multiple_in_file(par_file, par_file + ".new", str_to_replace, replacements)
        replace_multiple_in_file(par_file + ".new", par_file, [], [])

    # Set up the command to call tiny2
    cmd = "export TINYTIM=" + tinytim_path + "\n" + \
          tinytim_path + "/tiny2 " + par_file
    # Run the command to call tiny2
    sbp.call(cmd, shell=True)

    # Set up the command to call tiny3
    cmd = "export TINYTIM=" + tinytim_path + "\n" + \
          tinytim_path + "/tiny3 " + par_file + " SUB=" + \
        str(subsampling_factor)
    # Run the command to call tiny3
    sbp.call(cmd, shell=True)

    # Remove the unnecessary undistorted model
    try:
        os.remove(process_filename_base + mv.undistorted_model_tail)
    except OSError as _e:
        pass

    # Open the subsampled model image
    try:
        subsampled_image = fits.open(process_filename_base + mv.subsampled_model_tail)
    except IOError as _e:
        # This normally shouldn't happen, but might if two threads somehow end up using
        # the same lock. Best case is to wait and see if the other thread handles this
        # for us
        time_start = time.time()
        while ((not os.path.isfile(process_filename_base + mv.subsampled_model_tail))
            and (time.time() - time_start < 120)):
            time.sleep(1)
        
        # Give the file a chance to be fully written
        time.sleep(1)
        
        if not os.path.isfile(process_filename_base + mv.subsampled_model_tail):
            raise
        else:
            subsampled_image = fits.open(process_filename_base + mv.subsampled_model_tail)
            

    # Define a modified weight function to work on subsampled pixels
    def ss_weight_func(x, y):
        return weight_func(x/subsampling_factor, y/subsampling_factor)

    # Get the centre of this image
    ss_xc, ss_yc, _ss_x_array, _ss_y_array, _ss_weight_mask, ss_m0 = \
            centre_image(image=subsampled_image[0].data, weight_func=ss_weight_func)

    subsampled_image[0].header[mv.ss_model_xc_label] = ss_xc
    subsampled_image[0].header[mv.ss_model_yc_label] = ss_yc
    subsampled_image[0].header[mv.ss_model_m0_label] = ss_m0
    
    if shape is not None:
        full_shape = np.shape(subsampled_image[0].data)
        dx, dy = np.subtract(full_shape,shape) // 2
        
        if dx > 0 and dy > 0:
            subsampled_image[0].data = subsampled_image[0].data[dx:-dx,dy:-dy]
        elif dx > 0:
            subsampled_image[0].data = subsampled_image[0].data[dx:-dx,:]
        elif dy > 0:
            subsampled_image[0].data = subsampled_image[0].data[:,dy:-dy]
            
        full_shape = np.shape(subsampled_image[0].data)
        if full_shape[0] > shape[0]:
            subsampled_image[0].data = subsampled_image[0].data[1:,:]
        if full_shape[1] > shape[1]:
            subsampled_image[0].data = subsampled_image[0].data[:,1:]
            
    # Normalize the image
    subsampled_image[0].data /= subsampled_image[0].data.sum()

    # Write the image out to the proper filename
    if use_cache:
        subsampled_image.writeto(filename, clobber=True)

    if files_to_cleanup is not None:
        files_to_cleanup.append(filename)

    # Remove the old version, plus the .par and .tt3 files
    try:
        os.remove(process_filename_base + mv.subsampled_model_tail)
    except OSError as _e:
        pass
    try:
        os.remove(par_file)
    except OSError as _e:
        pass
    try:
        os.remove(par_file + ".new")
    except OSError as _e:
        pass
    try:
        os.remove(process_filename_base + ".tt3")
    except OSError as _e:
        pass
    
    try:
        os.remove(lock_filename)
    except OSError as _e:
        pass

    return subsampled_image[0]

@lru_cache(64)
def get_cached_subsampled_psf(tinytim_path,
                              tinytim_data_path,
                              weight_func,
                              psf_position,
                              chip,
                              focus,
                              use_cache,
                              astigmatism_0=None,
                              astigmatism_45=None,
                              coma_x=None,
                              coma_y=None,
                              clover_x=None,
                              clover_y=None,
                              spherical_3rd=None,
                              z12=None,
                              z13=None,
                              z14=None,
                              z15=None,
                              z16=None,
                              z17=None,
                              z18=None,
                              z19=None,
                              z20=None,
                              z21=None,
                              spherical_5th=None):

    # Determine the name for the subsampled model PSF file
    subsampled_name = os.path.join(tinytim_data_path, "subsampled_psf_x-" + str(psf_position[0]) + \
                        "_y-" + str(psf_position[1]) + "_f-" + str(focus) + \
                        "_c-" + str(chip))
    
    for (value, label) in ((astigmatism_0, "a0"),
                           (astigmatism_45, "a45"),
                           (coma_x, "cx"),
                           (coma_y, "cy"),
                           (clover_x, "clx"),
                           (clover_y, "cly"),
                           (spherical_3rd, "s3"),
                           (z12, "z12"),
                           (z13, "z13"),
                           (z14, "z14"),
                           (z15, "z15"),
                           (z16, "z16"),
                           (z17, "z17"),
                           (z18, "z18"),
                           (z19, "z19"),
                           (z20, "z20"),
                           (z21, "z21"),
                           (spherical_5th, "s5"),):
        if value is not None:
            subsampled_name += "_" + label + "-" + str(value)
    
    subsampled_name +=  mv.image_extension

    # Check if we need to update this file, or if we can reuse the existing version
    # if file_needs_update(subsampled_name) or len(params)>0:
    if file_needs_update(subsampled_name):

        # We'll need to update it, so we'll call TinyTim to generate a PSF model
        subsampled_model = make_subsampled_psf_model(filename=subsampled_name,
                                  xp=psf_position[0],
                                  yp=psf_position[1],
                                  focus=focus,
                                  use_cache=use_cache,
                                  chip=chip,
                                  weight_func=weight_func,
                                  tinytim_path=tinytim_path,
                                  astigmatism_0=astigmatism_0,
                                  astigmatism_45=astigmatism_45,
                                  coma_x=coma_x,
                                  coma_y=coma_y,
                                  clover_x=clover_x,
                                  clover_y=clover_y,
                                  spherical_3rd=spherical_3rd,
                                  z12=z12,
                                  z13=z13,
                                  z14=z14,
                                  z15=z15,
                                  z16=z16,
                                  z17=z17,
                                  z18=z18,
                                  z19=z19,
                                  z20=z20,
                                  z21=z21,
                                  spherical_5th=spherical_5th)

    else:

        # It doesn't need an update, so open up the old image
        try:
            subsampled_model = fits.open(subsampled_name)[0]
        except IOError as _e:
            # File is corrupt, so we'll regenerate it
            subsampled_model = make_subsampled_psf_model(filename=subsampled_name,
                                      xp=psf_position[0],
                                      yp=psf_position[1],
                                      focus=focus,
                                      chip=chip,
                                      weight_func=weight_func,
                                      tinytim_path=tinytim_path)
            
    return subsampled_model

def get_model_psf_for_star(star,
                           scheme,
                           weight_func=mv.default_prim_weight_func,
                           tinytim_path=mv.default_tinytim_path,
                           tinytim_data_path=mv.default_tinytim_data_path,
                           **params):
    """ Gets a model psf for a given star and the chip it was detected on

        Requires: star <star>
                  scheme <psf_model_scheme> (plan for generating model psfs)
        Optional: weight_func <function of x,y> (Used for determining the centre of the
                                                 subsampled model)
                  tinytim_path <string> (location of TinyTim executable)
                  tinytim_data_path <string> (location where TinyTim data will be stored)
                  **params <dict> (Extra parameters for describing the PSF)

        Returns: model_psf <star> (with m_err, Q values, and Q errors not yet determined)
    """
    
    use_cache = True
#     if len(params) > 0:
#         use_cache = False

    # Get the position we'll generate the model PSF for
    psf_position = scheme.get_position_to_use(star.x_pix, star.y_pix)
    
    focus = round(scheme.focus,5)
    
    rounded_params = {}
    
    for param in params:
        rounded_params[param] = round(params[param],5)

    subsampled_model = get_cached_subsampled_psf(tinytim_path,
                                                 tinytim_data_path,
                                                 weight_func,
                                                 psf_position,
                                                 star.chip,
                                                 focus,
                                                 use_cache=use_cache,
                                                 **rounded_params)
            

    # Get the charge diffusion kernel from the FITS comments

    if(mv.default_subsampling_factor > 1):
        fits_comments = subsampled_model.header['COMMENT']
    
        for test_i, s in enumerate(fits_comments):
            if 'following kernel' in s:
                i = test_i
                break
        if i == -1:
            raise Exception("Cannot find charge-diffusion kernel in fits file passed to " +
                            "read_kernel_from_fits")

        # kernel parameters are located in the three lines following to that index
        kernel = []
        for j in fits_comments[i + 1:i + 4]:
            kernel.append([float(x) for x in j.split()])
    
        # Convert to an ndarray
        kernel = np.asarray(kernel)
    else:
        kernel = np.array([[0.,0.,0.],[0.,1.,0.],[0.,0.,0.]])

    # Determine how far off the centre of the subsampled image is from the centre
    ss_ny, ss_nx = np.shape(subsampled_model.data)

    ss_model_d_xc = subsampled_model.header[mv.ss_model_xc_label] - (ss_nx - 1.) / 2
    ss_model_d_yc = subsampled_model.header[mv.ss_model_yc_label] - (ss_ny - 1.) / 2

    # Get the same for the star's postage stamp
    star_ny, star_nx = np.shape(star.stamp)

    star_d_xc = star.xc - (star_nx - 1.) / 2
    star_d_yc = star.yc - (star_ny - 1.) / 2

    # Determine how many subsampled pixels we'll have to shift the subsampled psf by
    x_shift = int(round(mv.default_subsampling_factor * (star_d_xc+0.5) - ss_model_d_xc - 0.5,0))
    y_shift = int(round(mv.default_subsampling_factor * (star_d_yc+0.5) - ss_model_d_yc - 0.5,0))
    
    # Check that the shifts are reasonable (within 2 non-subsampled pixels)
    max_shift = np.max((np.abs(x_shift),np.abs(y_shift)))/mv.default_subsampling_factor
    if max_shift > 2:
        raise Exception("Star's centring is too poor; requires too extreme of a shift.")

    # Get the rebinned PSF model
    rebinned_model = rebin(subsampled_model.data,
                           kernel,
                           x_shift=x_shift,
                           y_shift=y_shift,
                           subsampling_factor=mv.default_subsampling_factor)

    # Get the zeroth-order moment for the rebinned psf
    _xc, _yc, _, _, _, rb_model_m0 = centre_image(rebinned_model, weight_func=weight_func)
    

    scaled_model = rebinned_model * star.m0[0] / rb_model_m0

    return scaled_model
