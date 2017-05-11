""" This module contains the needed functions to generate a Tiny Tim PSF, using either
    the best-fit optical parameters found by Gillis et al. (2017) or the linear
    relationship with the focus-secondary-mirror despace.
"""

__all__ = ['make_subsampled_model_psf']

import subprocess as sbp
import os

# Default values
default_psf_position = (2048, 1024) # Center of the detector by default
default_focus = -1.0 # Approximately the middle of expected values
default_chip = 1 
default_spec_type = (1, 15) # Use spectrum for a K-type star by default
default_filter_name = 'f606w'
default_detector = 15 # ACS WFC
default_psf_size = 2.0
default_tinytim_path = "../../Program_Files/tinytim-7.5" # Adjust as needed for your own purposes
default_subsampling_factor = 7

# Default optical parameters
optical_params_means = {"z2":                -0.000800157,
                        "z3":                -0.00199172,
                        "astigmatism_0":      0.0445775,
                        "astigmatism_45":     0.0368361,
                        "coma_x":             0.011953,
                        "coma_y":             0.001, # Using default due to bug, will be fixed ASAP
                        "clover_x":           0.00923067,
                        "clover_y":           0.0274982,
                        "spherical_3rd":     -0.0290079,
                        "z12":                0.00169178,
                        "z13":               -0.00215789,
                        "z14":                0.00615519,
                        "z15":                0.00682615,
                        "z16":               -0.00226776,
                        "z17":                0.000538391,
                        "z18":               -0.00126454,
                        "z19":               -0.000948831,
                        "z20":                0.00168733,
                        "z21":                0.00207099,
                        "spherical_5th":      0.00984861,
                        "kernel_adjustment":  1.00105}

# Linear fit for optical parameters in (intercept, slope)
optical_params_int_and_slopes = {"z2":                (-0.002873648083,-0.00130053309303),
                                 "z3":                (-0.00181298435543,-0.00163556185266),
                                 "astigmatism_0":     (0.0406578112278,-0.00163787058868),
                                 "astigmatism_45":    (0.0325217048069,-0.00583735922479),
                                 "coma_x":            (0.0134420318745,0.00207411113338),
                                 "coma_y":            (0.001,0), # Using default due to bug, will be fixed ASAP
                                 "clover_x":          (0.0101213358836,0.00153365607983),
                                 "clover_y":          (0.0345197287386,0.00764902577666),
                                 "spherical_3rd":     (-0.029961124731,-0.00623843207348),
                                 "z12":               (0.00192029654149,0.000881703385712),
                                 "z13":               (-0.00261522941263,-0.0011394796251),
                                 "z14":               (0.00602192631498,0.00101694195336),
                                 "z15":               (0.00790605208103,0.000945861376067),
                                 "z16":               (-0.00335267943575,-0.00362187843696),
                                 "z17":               (0.00185168971574,-0.000174349121967),
                                 "z18":               (-0.000426898611297,0.000419751087548),
                                 "z19":               (-0.0022165718048,-0.0015019737541),
                                 "z20":               (0.00150504289848,0.000747263759456),
                                 "z21":               (-0.000996577577823,-0.00263556844535),
                                 "spherical_5th":     (0.0105835881518,0.00138130728956),
                                 "kernel_adjustment": (1.00265342235,0.000602841490909),}


def replace_multiple_in_file(input_filename, output_filename, input_strings, output_strings):
    """ Replaces every occurence of an input_string in input_filename with the corresponding
        output string and prints the results to $output_filename.
        
        Requires: input_filename <string>
                  output_filename <string>
                  input_strings <iterable of strings>
                  output_strings <iterable of strings>
                  
        Returns: None
        Side-effects: $output_filename is created/overwritten
    """
    
    with open(output_filename, "w") as fout:
        with open(input_filename, "r") as fin:
            for line in fin:
                new_line = line
                for input_string, output_string in zip(input_strings, output_strings):
                    if((input_string is None) or (output_string is None)):
                        continue
                    new_line = new_line.replace(input_string, output_string)
                fout.write(new_line)
                
    return

def make_subsampled_model_psf(filename,
                              psf_position = default_psf_position,
                              focus = default_focus,
                              chip = default_chip,
                              spec_type = default_spec_type,
                              detector = default_detector,
                              filter_name = default_filter_name,
                              psf_size=default_psf_size,
                              tinytim_path = default_tinytim_path,
                              subsampling_factor = default_subsampling_factor,
                              linear_fit = False,
                              clobber = True,
                              **optical_params):
    """ Generates a subsampled model PSF, using the desired (or default) optical parameters.
        For input parameters spec_type and detector, the allowed options can be seen through
        running tiny1
    
        @param filename <str> Desired filename of the generated PSF. If it already exists, the
                              'clobber' parameter will determine whether or not it will be
                              overwritten.
        @param psf_position <(float, float)> Position on the detector of the PSF in x, y
        @param focus <float> Focus-secondary-mirror despace of the PSF
        @param chip <int> Which chip the model PSF is for. Allowed values are 1 and 2
        @param spec_type <(int, *)> Spectral type of the PSF to generate. First value chooses type
                                    of spectrum, second chooses from options for this type
        @param detector <int> Index of detector to be used.
        @param filter_name <str> Name of the filter to use (eg. f606w)
        @param psf_size <float> Size of the PSF image in arcseconds
        @param tinytim_path <str> Location of the Tiny Tim executables
        @param subsampling_factor <int> Factor by which to subsample the PSF
        @param linear_fit <bool> If False, unspecified optical parameters will be given values
                                 based on the mean from Gillis et al. (2017)'s testing. If True,
                                 will use the linear fit from the analysis instead
        @param clobber <bool> Whether or not to overwrite the target file if it already exists.
        @param optical_params <dict> Optical parameters aside from focus for this PSF. If not
                                     specified here, defaults will be used based on the
                                     linear_fit parameter.
                                 

    """
    
    # If clobber is False, check if the desired file already exists
    if not clobber:
        if os.path.isfile(filename):
            raise IOError("File " + filename + " already exists. Set clobber=True if you wish to overwrite it.")
    
    # Create a directory to contain this project
    try:
        os.makedirs(os.path.split(filename)[0])
    except OSError as e:
        if not ("[Errno 17] File exists:" in str(e) or "[Errno 2] No such file or directory: ''" in str(e)):
            raise
        else:
            pass # No need to raise if the directory already exists

    filename_base = filename.replace(".fits", "")
    if not filename_base + ".fits" == filename:
        raise ValueError("Filename (" + filename + ") must end in '.fits'.")
    
    par_file = filename_base + ".par"
    tmp_par_file = filename_base + ".par.tmp"

    # Set up the command to call tiny1 and execute it
    cmd = "export TINYTIM=" + tinytim_path + "\n" + \
          tinytim_path + "/tiny1 " + tmp_par_file + " << EOF \n" + \
          str(detector) + "\n" + \
          str(chip) + "\n" + \
          str(psf_position[0]) + " " + str(psf_position[1]) + "\n" + \
          str(filter_name) + "\n" + \
          str(spec_type[0]) + "\n" + \
          str(spec_type[1]) + "\n" + \
          str(psf_size) + "\n" + \
          str(focus) + "\n" + \
          filename_base + "\nEOF"
    
    sbp.call(cmd, shell=True)
    
    # Determine which optical parameters we'll be using
    optical_params_to_use = {}
    for param in optical_params_means:
        if param in optical_params:
            optical_params_to_use[param] = optical_params[param]
        elif linear_fit:
            intercept, slope = optical_params_int_and_slopes[param]
            optical_params_to_use[param] = intercept + focus*slope
        else:
            optical_params_to_use[param] = optical_params_means[param]
    
    # Edit the parameter file to adjust optical parameters
    strs_to_replace = []
    replacements = []

    strs_to_replace.append("0.       # Z2 = X (V2) tilt")
    replacements.append(str(optical_params_to_use["z2"]) + "       # Z2 = X (V2) tilt")

    strs_to_replace.append("0.       # Z3 = Y (V3) tilt")
    replacements.append(str(optical_params_to_use["z3"]) + "       # Z3 = Y (V3) tilt")

    strs_to_replace.append("0.031    # Z5 = 0 degree astigmatism")
    replacements.append(str(optical_params_to_use["astigmatism_0"]) + "    # Z5 = 0 degree astigmatism")

    strs_to_replace.append("0.028    # Z6 = 45 degree astigmatism")
    replacements.append(str(optical_params_to_use["astigmatism_45"]) + "    # Z6 = 45 degree astigmatism")

    strs_to_replace.append("0.003    # Z7 = X (V2) coma")
    replacements.append(str(optical_params_to_use["coma_x"]) + "    # Z7 = X (V2) coma")

    strs_to_replace.append("0.001    # Z8 = Y (V3) coma")
    replacements.append(str(optical_params_to_use["coma_y"]) + "    # Z8 = Y (V3) coma")

    strs_to_replace.append("0.008    # Z9 = X clover")
    replacements.append(str(optical_params_to_use["clover_x"]) + "    # Z7 = X (V2) clover")

    strs_to_replace.append("0.018    # Z10 = Y clover")
    replacements.append(str(optical_params_to_use["clover_y"]) + "    # Z8 = Y (V3) clover")

    strs_to_replace.append("-0.025    # Z11 = 3rd order spherical")
    replacements.append(str(optical_params_to_use["spherical_3rd"]) + "    # Z11 = 3rd order spherical")

    strs_to_replace.append("0.       # Z12 = 0 degree Spherical astigmatism")
    replacements.append(str(optical_params_to_use["z12"]) + "       # Z12 = 0 degree Spherical astigmatism")

    strs_to_replace.append("0.       # Z13 = 45 degree Spherical astigmatism")
    replacements.append(str(optical_params_to_use["z13"]) + "       # Z13 = 45 degree Spherical astigmatism")

    strs_to_replace.append("0.       # Z14 = X Ashtray")
    replacements.append(str(optical_params_to_use["z14"]) + "       # Z14 = X Ashtray")

    strs_to_replace.append("0.       # Z15 = Y Ashtray")
    replacements.append(str(optical_params_to_use["z15"]) + "       # Z15 = Y Ashtray")

    strs_to_replace.append("0.       # Z16")
    replacements.append(str(optical_params_to_use["z16"]) + "       # Z16")

    strs_to_replace.append("0.       # Z17")
    replacements.append(str(optical_params_to_use["z17"]) + "       # Z17")

    strs_to_replace.append("0.       # Z18")
    replacements.append(str(optical_params_to_use["z18"]) + "       # Z18")

    strs_to_replace.append("0.       # Z19")
    replacements.append(str(optical_params_to_use["z19"]) + "       # Z19")

    strs_to_replace.append("0.       # Z20")
    replacements.append(str(optical_params_to_use["z20"]) + "       # Z20")

    strs_to_replace.append("0.       # Z21")
    replacements.append(str(optical_params_to_use["z21"]) + "       # Z21")

    strs_to_replace.append("0.009    # Z22 = 5th order spherical")
    replacements.append(str(optical_params_to_use["spherical_5th"]) + "    # Z22 = 5th order spherical")

    replace_multiple_in_file(tmp_par_file, par_file, strs_to_replace, replacements)

    # Set up the command to call tiny2
    cmd = "export TINYTIM=" + tinytim_path + "\n" + \
          tinytim_path + "/tiny2 " + par_file
    # Run the command to call tiny2
    sbp.call(cmd, shell=True)

    # Set up the command to call tiny3
    cmd = "export TINYTIM=" + tinytim_path + "\n" + \
          tinytim_path + "/tiny3 " + par_file + " SUB=" + \
          str(int(subsampling_factor))
    # Run the command to call tiny3
    sbp.call(cmd, shell=True)
    
    # PSF should be generated, now move it to the desired filename
    init_filename = filename_base + "00.fits"
    os.rename(init_filename,filename)

    # Clean up unneeded files. Silently suppress any errors here
    try:
        os.remove(filename_base + "00_psf.fits")
    except OSError as _e:
        pass
    try:
        os.remove(filename_base + ".tt3")
    except OSError as _e:
        pass
    try:
        os.remove(init_filename)
    except OSError as _e:
        pass
    try:
        os.remove(par_file)
    except OSError as _e:
        pass
    try:
        os.remove(tmp_par_file)
    except OSError as _e:
        pass

    return