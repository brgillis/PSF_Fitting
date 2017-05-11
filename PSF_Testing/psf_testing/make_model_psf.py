""" This module contains the needed functions to generate a Tiny Tim PSF, using either
    the best-fit optical parameters found by Gillis et al. (2017) or the linear
    relationship with the focus-secondary-mirror despace.
"""

import subprocess as sbp
import os
from psf_testing.magic_values import default_tinytim_path

# Default values
default_psf_position = (2048, 1024) # Center of the detector by default
default_focus = -1.0 # Approximately the middle of expected values
default_chip = 1 
default_spec_type = (1, 15) # Use spectrum for a K-type star by default
default_detector = 15 # ACS WFC
default_psf_size = 2.0
default_tinytim_path = "../../Program_Files/tinytim-7.5" # Adjust as needed for your own purposes
default_subsampling_factor = 7

--z2 -0.000800157 --z3 -0.00199172 --astigmatism_0 0.0445775 --astigmatism_45 0.0368361 --coma_x 0.011953 --coma_y 0.011953 --clover_x 0.00923067 --clover_y 0.0274982 --spherical_3rd -0.0290079 --z12 0.00169178 --z13 -0.00215789 --z14 0.00615519 --z15 0.00682615 --z16 -0.00226776 --z17 0.000538391 --z18 -0.00126454 --z19 -0.000948831 --z20 0.00168733 --z21 0.00207099 --spherical_5th 0.00984861 --kernel_adjustment 1.00105

# Default optical parameters
default_optical_params = {"z2":-0.000800157,
                          "z3":-0.00199172,
                          "astigmatism_0":0.0445775,
                          "astigmatism_45":0.0368361,
                          "coma_x":0.011953,
                          "coma_y":0.011953,
                          "clover_x":,
                          "clover_y":,
                          "spherical_3rd":,
                          "z12":,
                          "z13":,
                          "z14":,
                          "z15":,
                          "z16":,
                          "z17":,
                          "z18":,
                          "z19":,
                          "z20":,
                          "z21":,
                          "spherical_5th":,}

# Linear fit for optical parameters in (intercept, slope)

def make_subsampled_model_psf(filename,
                              psf_position,
                              focus = default_focus,
                              chip = default_chip,
                              spec_type = default_spec_type,
                              detector = default_detector,
                              psf_size=default_psf_size,
                              tinytim_path = default_tinytim_path,
                              linear_fit = False,
                              clobber = False,
                              **optical_params):
    """ @TODO function docstring

    """
    
    # If clobber is False, check if the desired file already exists
    if not clobber:
        if os.path.isfile(filename):
            raise IOError("File " + filename + " already exists. Set clobber=True if you wish to overwrite it.")
    
    # Create a directory to contain this project
    try:
        os.makedirs(os.path.split(filename)[0])
    except OSError as e:
        if not "[Errno 17] File exists:" in e:
            raise
        else:
            pass # No need to raise if the directory already exists

    filename_base = filename.replace(".fits", "")
    if not filename_base + ".fits" == filename:
        raise ValueError("Filename (" + filename + ") must end in '.fits'.")
    
    par_file = filename_base + ".par"

    # Set up the command to call tiny1 and execute it
    cmd = "export TINYTIM=" + tinytim_path + "\n" + \
          tinytim_path + "/tiny1 " + par_file + " << EOF \n" + \
          str(detector) + "\n" + \
          str(chip) + "\n" + \
          str(psf_position[0]) + " " + str(psf_position[1]) + "\n" + \
          str(filter) + "\n" + \
          str(spec_type[0]) + "\n" + \
          str(spec_type[1]) + "\n" + \
          str(psf_size) + "\n" + \
          str(focus) + "\n" + \
          filename_base + "\nEOF"
    sbp.call(cmd, shell=True)
    
    # Determine which optical parameters we'll be using
    
    # Edit the parameter file to adjust optical parameters as necessary
    strs_to_replace = []
    replacements = []
    if z2 is not None:
        strs_to_replace.append("0.       # Z2 = X (V2) tilt")
        replacements.append(str(z2) + "       # Z2 = X (V2) tilt")
    if z3 is not None:
        strs_to_replace.append("0.       # Z3 = Y (V3) tilt")
        replacements.append(str(z3) + "       # Z3 = Y (V3) tilt")
    if astigmatism_0 is not None:
        strs_to_replace.append("0.031    # Z5 = 0 degree astigmatism")
        replacements.append(str(astigmatism_0) + "    # Z5 = 0 degree astigmatism")
    if astigmatism_45 is not None:
        strs_to_replace.append("0.028    # Z6 = 45 degree astigmatism")
        replacements.append(str(astigmatism_45) + "    # Z6 = 45 degree astigmatism")
    if coma_x is not None:
        strs_to_replace.append("0.003    # Z7 = X (V2) coma")
        replacements.append(str(coma_x) + "    # Z7 = X (V2) coma")
    if coma_y is not None:
        strs_to_replace.append("0.001    # Z8 = Y (V3) coma")
        replacements.append(str(coma_y) + "    # Z8 = Y (V3) coma")
    if clover_x is not None:
        strs_to_replace.append("0.008    # Z9 = X clover")
        replacements.append(str(clover_x) + "    # Z7 = X (V2) clover")
    if clover_y is not None:
        strs_to_replace.append("0.018    # Z10 = Y clover")
        replacements.append(str(clover_y) + "    # Z8 = Y (V3) clover")
    if spherical_3rd is not None:
        strs_to_replace.append("-0.025    # Z11 = 3rd order spherical")
        replacements.append(str(spherical_3rd) + "    # Z11 = 3rd order spherical")
    if z12 is not None:
        strs_to_replace.append("0.       # Z12 = 0 degree Spherical astigmatism")
        replacements.append(str(z12) + "       # Z12 = 0 degree Spherical astigmatism")
    if z13 is not None:
        strs_to_replace.append("0.       # Z13 = 45 degree Spherical astigmatism")
        replacements.append(str(z13) + "       # Z13 = 45 degree Spherical astigmatism")
    if z14 is not None:
        strs_to_replace.append("0.       # Z14 = X Ashtray")
        replacements.append(str(z14) + "       # Z14 = X Ashtray")
    if z15 is not None:
        strs_to_replace.append("0.       # Z15 = Y Ashtray")
        replacements.append(str(z15) + "       # Z15 = Y Ashtray")
    if z16 is not None:
        strs_to_replace.append("0.       # Z16")
        replacements.append(str(z16) + "       # Z16")
    if z17 is not None:
        strs_to_replace.append("0.       # Z17")
        replacements.append(str(z17) + "       # Z17")
    if z18 is not None:
        strs_to_replace.append("0.       # Z18")
        replacements.append(str(z18) + "       # Z18")
    if z19 is not None:
        strs_to_replace.append("0.       # Z19")
        replacements.append(str(z19) + "       # Z19")
    if z20 is not None:
        strs_to_replace.append("0.       # Z20")
        replacements.append(str(z20) + "       # Z20")
    if z21 is not None:
        strs_to_replace.append("0.       # Z21")
        replacements.append(str(z21) + "       # Z21")
    if spherical_5th is not None:
        strs_to_replace.append("0.009    # Z22 = 5th order spherical")
        replacements.append(str(spherical_5th) + "    # Z22 = 5th order spherical")

    if len(strs_to_replace) > 0:
        replace_multiple_in_file(par_file, par_file + ".new", strs_to_replace, replacements)
        sbp.call("mv " + par_file + ".new " + par_file)

    return subsampled_model