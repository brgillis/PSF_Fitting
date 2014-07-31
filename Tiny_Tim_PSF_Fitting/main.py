#!/usr/bin/env python

import sys
oldpath = sys.path
sys.path = []
sys.path.append('/home/brg/lib/python2.7/site-packages')
sys.path.append('.')
sys.path = sys.path + oldpath
from os.path import isfile
import numpy as np
import subprocess as sbp
import multiprocessing as mtp

from rebin import rebin_and_shift_psf
from subtract_psf import subtract_psf
from cut_postage_stamp import cut_postage_stamp
from residual_chi2 import get_chi2
from whisker_shear import draw_whisker_shear
from call_tinytim import call_tinytim
from fits_functions import *

# Magic values
sextractor_cfg_name_end = "_sex.cfg"
sextractor_cat_name_end = "_sex_output.cat"
sextractor_cfg_template_name = "test_img.cfg"
sextractor_cat_template_name = "output.cat" # Make sure this matches the value in the config file sextractor_cfg_name
sextractor_psf_cfg_name_end = "_psf.cfg"
sextractor_psf_cat_name_end = "_psf_output.cat"
sextractor_psf_cfg_template_name = "psf.cfg"
sextractor_psf_cat_template_name = "psf_output.cat"

def main(argv):
    """ This is the main execution script for testing how well Tiny Tim PSFs fit the Hubble data.
    
        Version 1.1
       
        Requires: The full name of an image to process, given in the command line.
       
        Optional: Class_star threshold
                  Star_mag threshold
                  Number of grid psfs in x
                  Number of grid psfs in y (if these are omitted or set to -1,
                                            will calc psf for each star)
                  Whether or not to fit the best focus (True, False, or specific focus value in um)
                  Whether or not to cleanup files generated except for key output (True if so, blank or
                                                                            otherwise to not cleanup.)
                  Summary file to output fitted data to (Won't output if blank, will append one line to
                                                         it if entered)
                 
                  ex. python main.py file_name.fits 0.95 22 (No fit, use default focus of 0,
                                                                  calculate PSF for each star)
                  ex. python main.py file_name.fits 0.95 22 -1 -1 False (No fit, use default focus of 0,
                                                                  calculate PSF for each star
                                                                  (same as above))
                  ex. python main.py file_name.fits 0.95 22 -1 -1 3.0 (No fit, use focus offset of
                                                                       3.0 um, calculate PSF for each
                                                                       star)
                  ex. python main.py file_name.fits 0.95 22 8 4 True (fit focus, 8x4 grid)
                 
        Returns: (nothing)
       
        Side-effects: Creates number of .fits images for postage stamps, PSFs, and PSF residuals.
                      Creates two .dat ASCII tables listing the the multipole expansion moments
                          for the residuals,
                          along with size and shape measurements of PSFs and actual star images,
                          for both exponential weighting (*_moments.dat) and no weighting
                          (*_moments_wings.dat)
                      Creates a stacked residual image of all non-outlier stars minus their PSFs, divided by the
                           star fluxes.
                          
        Except: Will raise an exception if:
                    No image file name is passed in the command line
                    The image file to be used cannot be read
                    The various files to be created cannot be written
                  
    """

    # Magic numbers
    default_class_star_threshold = 0.95
    default_max_star_mag = 22
    default_focus = 0.0
    default_subsampling_factor= 5
    stamp_size_factor = 1 # How much bigger the stamps should be at minimum than the size of the detected star
    min_stamp_size = 20 # "radius"-like; actual minimum size will be twice this, plus 1
    
    # Initialize the params dictionary
    params = {}
    
    # First, check that we were passed a string in the command-line for an image file
    if(len(argv)) <= 1:
        raise Exception("Name of fits image file must be passed at command-line.\n" + \
                        "eg. python main.py HST_image.fits")
    
    # Check if a CLASS_STAR or star_mag threshold was passed in the command-line
    if(len(argv) <= 2):
        params['class_star_threshold'] = default_class_star_threshold
        params['max_star_mag'] = default_max_star_mag
    else:
        try:
            params['class_star_threshold'] = float(argv[2])
            if(params['class_star_threshold'] <= 0):
                params['class_star_threshold'] = default_class_star_threshold
        except:
            params['class_star_threshold'] = default_class_star_threshold
        if(len(argv) <= 3):
            params['max_star_mag'] = default_max_star_mag
        else:
            try:
                params['max_star_mag'] = float(argv[3])
                if(params['max_star_mag'] <= 0):
                    params['max_star_mag'] = default_max_star_mag
            except:
                params['max_star_mag'] = default_max_star_mag
                
    # Check if we're using the quick mode
    if(len(argv) <= 5):
        params['quick_mode'] = False
        params['x_psfs'] = -1
        params['y_psfs'] = -1
    else:
        try:
            params['quick_mode'] = True
            params['x_psfs'] = int(argv[4])
            params['y_psfs'] = int(argv[5])
        except:
            params['quick_mode'] = False
            params['x_psfs'] = -1
            params['y_psfs'] = -1
            
    if((params['x_psfs'] <= 0) or (params['y_psfs'] <= 0)):
        params['quick_mode'] = False
        params['x_psfs'] = -1
        params['y_psfs'] = -1
        print "Quick mode disabled. Calculating PSF for each star."
    else:
        params['quick_mode'] = True # Shouldn't be necessary, but just in case
        print "Quick mode enabled. Calculating a grid of " + str(params['x_psfs']) + "x" + \
              str(params['y_psfs']) + "psfs to use."

    # Check if we're fitting for best focus
    if(len(argv) <= 6):
        params['fit_focus'] = False
        params['focus'] = default_focus
    else:
        try:
            params['fit_focus'] = ( ( argv[6] == 'True' ) or (argv[6] == 'true') )
            try:
                params['focus'] = float(argv[6])
            except:
                params['focus'] = default_focus
        except:
            params['fit_focus'] = False
            params['focus'] = default_focus
            
    if(params['fit_focus']):
        print("Focus-fitting enabled.")
    else:
        print("Focus-fitting disabled. Using focus of " + str(params['focus']) + ".")
        
    # Check if we should clean up generated files when we're done
    if(len(argv) <= 7):
        params['cleanup'] = False
    else:
        try:
            params['cleanup'] = ( ( argv[7] == 'True' ) or (argv[7] == 'true') )
        except: params['cleanup'] = False
        
    if(len(argv) <= 8):
        params['summary_file_name'] = None
    else:
        try:
            params['summary_file_name'] = str(argv[8])
        except: params['summary_file_name'] = None
                
    params['image_file'] = argv[1]
    params['subsampling_factor'] = default_subsampling_factor
    
    # Does the name as given have the .fits extension?
    if(params['image_file'][-5:] == ".fits"):
        params['file_name_base'] = params['image_file'][0:-5]
    else:
        params['file_name_base'] = params['image_file']
        params['image_file'] = params['image_file'] + ".fits"
        
    # Get the size of the image and other info needed for tinytim
    get_image_info(params)

    # Now, we'll run sextractor on the passed image to identify stars
    
    # Run sextractor on the image
    sextractor_cat_name = params['file_name_base'] + sextractor_cat_name_end
    sextractor_cfg_name = params['file_name_base'] + sextractor_cfg_name_end
    cmd = "rm -f " + sextractor_cat_name
    sbp.call(cmd, shell=True)
    cmd = "awk '{sub(/" + sextractor_cat_template_name + "/,\"" + sextractor_cat_name + \
        "\")}; 1' " + sextractor_cfg_template_name + " > " + sextractor_cfg_name
    sbp.call(cmd, shell=True)
    cmd = "sex " + params['image_file'] + " -c " + sextractor_cfg_name
    sbp.call(cmd, shell=True)
    
    # Check that the output catalog exists now, to see if sextractor worked
    if(not isfile(sextractor_cat_name)):
        raise "Catalogue of stars could not be generated with sextractor.\n" + \
              "Check the input file and sextractor config files.\n"

    # Open the sextractor catalogue, and go through it to get a list of star positions
    object_lines = []
    with open(sextractor_cat_name, 'r') as cat_file:

        # Read in the file, except for comment lines
        for object_line in cat_file:
            object_line.strip()
            if((object_line[0] != '#') and (len(object_line) > 0) and (object_line[0] != '[')):
                object_lines.append(object_line.split())
                
    # Check which of these are likely to be stars
    stars = []
    i = 1
    for object_line in object_lines:
        class_star = float(object_line[16])
        star_mag = float(object_line[17])
        if((class_star > params['class_star_threshold']) and (star_mag < params['max_star_mag'])):
            xp = float(object_line[1])
            yp = float(object_line[2])
            
            dx = xp-np.round(xp)
            dy = yp-np.round(yp)
            xp = int(np.round(xp))
            yp = int(np.round(yp))
            
            xp_size = float(object_line[7]) - float(object_line[5])
            yp_size = float(object_line[8]) - float(object_line[6])
        
            if(xp_size > yp_size):
                stamp_size = np.ceil(xp_size * stamp_size_factor / 2) + 1
            else:
                stamp_size = np.ceil(yp_size * stamp_size_factor / 2) + 1
            
            if(stamp_size < min_stamp_size):
                stamp_size = min_stamp_size
                
            stamp_size = 2*stamp_size + 1
            
            # Check whether or not it's too close to an edge
            xl = int(xp) - int(stamp_size)
            yl = int(yp) - int(stamp_size)
            xm = int(xp) + int(stamp_size) + 1
            ym = int(yp) + int(stamp_size) + 1

            
            if((xl > 0) and (yl > 0) and (xm < params['x_size']) and (ym < params['y_size'])):
                a = float(object_line[21])
                b = float(object_line[22])
                theta = float(object_line[11])
                
                # Calculate ellipticity from a, b, and theta
                e = (a-b)/(a+b)
                e1 = e*np.cos(2*theta)
                e2 = e*np.sin(2*theta)
            
                flux_auto = float(object_line[19])
                flux_rad = float(object_line[14])
                kron_rad = float(object_line[13])*np.sqrt(a*b)
                
                n = flux_rad / kron_rad
                
                total_flux = flux_auto*flux_ratio(n)

                star = (i, xp, yp, class_star, star_mag, stamp_size, dx, dy, total_flux,
                        kron_rad, e1, e2)
                stars.append(star)
                i += 1
                
    # We now have our set of stars to work with
    
    # Now, either run once with the set focus or fit the focus
    try:
        if(params['fit_focus']):
            remove_psf_after_fit(stars, params)
        else:
            remove_psf_no_fit(stars, params)
    except:
        pass
    finally:
        if(params['cleanup']):
            cleanup(params['file_name_base'])
        
def get_image_info(params):
    """Gets the size (in pixels) of a .fits image file and other information needed by tinytim,
       using the image's header.
    
       Requires: params <dict> (Params dictionary, requiring the 'image_file' key)
       
       Returns: (nothing)
       
       Side-effects: params dictionary loaded with the following keys:
                         'x_size', 'y_size', 'chip', 'detector', 'filter'
                          
       Except: Unknown detector for the image
               Image cannot be accessed or is missing needed header values
    """
    
    fits_struct = pyf.open(params['image_file'])
    header = fits_struct[0].header
    params['x_size'] = header['NAXIS1']
    params['y_size'] = header['NAXIS2']
    params['chip'] = header['CCDCHIP']
    if(header['DETECTOR'].strip().lower()=="wfc"):
        params['detector'] = 15
    else:
        raise Exception("ERROR: Unknown detector for image: " + header['DETECTOR'] + "\n" +
                    "If tinytim can use this detector, please add its number to the get_image_info function.")
    params['filter'] = header['FILTER1'].strip().lower()
    
    
        
def cleanup(file_name_base):
    """ Cleans up intermediate image files generated by this script, including postage stamps, psfs,
        and residual images.
        
        NOTE: This may affect files generated from other runs. Only enable cleanup when you're sure
        it's what you want to do, as files deleted aren't usually recoverable. See below for list of
        file name patterns that will be deleted.
        
        Requires: file_name_base <string> (Root name for images to be cleaned up)
        
        Returns: (nothing)
        
        Side-effects: Deletes all files which match the following patterns:
                          file_name_base_stamp_*.fits
                          file_name_base_psf_*.fits
                          file_name_base_binned_*.fits
                          file_name_base_unbinned_*.fits
                          file_name_base_residual_*.fits (except file_name_base_residual_stack.fits)
                          file_name_base.par
                          file_name_base.tt3
                          file_name_base + sextractor_cfg_name_end
                          file_name_base + sextractor_cat_name_end
                          file_name_base + sextractor_psf_cfg_name_end
                          file_name_base + sextractor_psf_cat_name_end
    """
    
    # Move the residual stack file to preserve it if it exists
    cmd = "mv " + file_name_base + "_residual_stack.fits " + file_name_base + "_temp.fits"
    sbp.call(cmd, shell=True)
    
    # Delete files matching the right patterns
    cmd = "rm -f " + file_name_base + "_stamp_*.fits "+ \
                     file_name_base + "_psf_*.fits "+ \
                     file_name_base + "_binned_*.fits "+ \
                     file_name_base + "_unbinned_*.fits "+ \
                     file_name_base + "_residual_*.fits "+ \
                     file_name_base + sextractor_cfg_name_end + " " + \
                     file_name_base + sextractor_cat_name_end + " " + \
                     file_name_base + sextractor_psf_cfg_name_end + " " + \
                     file_name_base + sextractor_psf_cat_name_end + " " + \
                     file_name_base + ".par " + \
                     file_name_base + ".tt3 " + \
                     "check.fits"
    sbp.call(cmd,shell=True)
    
    # Move back the residual stack file
    cmd = "mv " + file_name_base + "_temp.fits " + file_name_base + "_residual_stack.fits"
    sbp.call(cmd, shell=True)

def get_image_size(image_file):
    """Gets the size (in pixels) of a .fits image file.
    
       Requires: image_file <string> (name of image)
       
       Returns: nx <int>, ny <int> (number of pixels in x- and y-dimensions, using PyFITS ordering,
                                    which is notably the reverse order of sextractor/tinytim.)
                          
       Side-effects: (none)
    """
    
    fits_struct = pyf.open(image_file)
    return np.shape(fits_struct[0].data)
    
def generate_psf(params, xp, yp, psf_file, subsampled_file):
    """Uses tinytim calling script to generate a psf at a given location in the image.
    
       Requires: params <dictionary> (parameter dictionary, properly set up with needed keys:
                                          'file_name_base', 'subsampling_factor', 'focus',
                                          'chip', 'detector', 'filter')
                 xp <int> (x-pixel coordinate (in sextractor/tinytim ordering))
                 yp <int> (y-pixel coordinate (in sextractor/tinytim ordering))
                 psf_file <string> (name of undistorted psf file to be generated)
                 subsampled_file <string> (name of subsampled, distorted psf file to be generated)
                 
       Returns: (nothing)
       
       Side-effects: Cleans and regenerates tinytim intermediate and output files
                     Creates and stores psf and subsampled files
                     Overwrites sextractor_psf_cat_name
                     
       Except: Will raise if:
                   params dictionary is missing needed keys
                      tinytim fails to generate the expected psf files for whatever reason.
                      sextractor fails to analyse the psf to find its barycentre
                      files that must be written to cannot be accessed or written
    """
    
    try:
        # Clean previous tinytim files if they exist
        cmd = "rm -f " + params['file_name_base'] + ".tt3 " + params['file_name_base'] + ".par " + \
              params['file_name_base'] + "00.fits" + params['file_name_base'] + "00_psf.fits"
        sbp.call(cmd, shell=True)
        
        # Invoke the calling function
        call_tinytim(params, xp, yp)
        
        # Check if the output from tinytim exists. If not, skip this star
        if(not isfile(params['file_name_base'] + "00_psf.fits")):
            raise Exception("Cannot create tinytime psfs.")
        if(not isfile(params['file_name_base'] + "00.fits")):
            raise Exception("Cannot create tinytime psfs.")
        
        # Move the psf file to store it
        cmd = "mv " + params['file_name_base'] + "00_psf.fits " + psf_file
        sbp.call(cmd, shell=True)
        
        # Move the unbinned distorted psf file to store it
        cmd = "mv " + params['file_name_base'] + "00.fits " + subsampled_file
        sbp.call(cmd, shell=True)
        
        # Now, use sextractor to determine the windowed barycentre of the unbinned PSF
        sextractor_psf_cat_name = params['file_name_base'] + sextractor_psf_cat_name_end
        sextractor_psf_cfg_name = params['file_name_base'] + sextractor_psf_cfg_name_end
        cmd = "rm -f " + sextractor_psf_cat_name
        sbp.call(cmd, shell=True)
        cmd = "awk '{sub(/" + sextractor_psf_cat_template_name + "/,\"" + sextractor_psf_cat_name + \
            "\")}; 1' " + sextractor_psf_cfg_template_name + " > " + sextractor_psf_cfg_name
        sbp.call(cmd, shell=True)
        cmd = "sex " + subsampled_file + " -c " + sextractor_psf_cfg_name
        sbp.call(cmd, shell=True)
        
        # Check that the output catalog exists now, to see if sextractor worked
        if(not isfile(sextractor_psf_cat_name)):
            raise Exception("PSF could not be analysed with sextractor.\n" + \
                  "Check the input file and sextractor config files.\n")
    
        # Open the sextractor catalogue, and go through it to get a list of star positions
        object_lines = []
        with open(sextractor_psf_cat_name, 'r') as cat_file:
    
            # Read in the file, except for comment lines
            for object_line in cat_file:
                object_line.strip()
                if((object_line[0] != '#') and (len(object_line) > 0) and (object_line[0] != '[')):
                    object_lines.append(object_line.split())
            
        # Assume first line is the one we want, since it should be the "brightest"
        xp = float(object_lines[0][0])
        yp = float(object_lines[0][1])
        a = float(object_lines[0][8])
        b = float(object_lines[0][9])
        theta = float(object_lines[0][10])
        kron_rad = float(object_lines[0][11]) * np.sqrt(a*b)/params['subsampling_factor']
        
        # Calculate ellipticity from a, b, and theta
        e = (a-b)/(a+b)
        e1 = e*np.cos(2*theta)
        e2 = e*np.sin(2*theta)
        
        # Add this info the psf's .fits header file
        # While we have it open, transpose the data to match the PyFITS ordering
        try:
            data = read_fits(subsampled_file)
            
            np.transpose(data[0].data)
            data[0].header['NAXIS1'],data[0].header['NAXIS2'] = \
                data[0].header['NAXIS2'],data[0].header['NAXIS1']
    
            # Note swap of xp and yp here, for the transposition
            data[0].header['XP_WIN'] = (yp,'x-pixel position of windowed barycentre')
            data[0].header['YP_WIN'] = (xp,'y-pixel position of windowed barycentre')
            data[0].header['E1_IMAGE'] = (e1,'e1 component of ellipticity in image frame')
            data[0].header['E2_IMAGE'] = (e2,'e2 component of ellipticity in image frame')
            data[0].header['KRON_RAD'] = (kron_rad,'Kron radius in non-subsampled pixel scale')
            
            add_comment(data, "Axes transposed for simplicity with PyFITS. Do not manipulate in " +
                              "conjunction with others without transposing axes back first.")
            write_fits(data, subsampled_file)
            
        except Exception, e:
            raise Exception("Could not write psf barycentre position to its .fits file:\n" +
              str(e))
        
    except KeyError, e:
        raise Exception("params dictionary passed to generate_psf is missing needed key:\n" +
              str(e))
    
    return
    
def load_psf_cache(params):
    """Generates a grid of PSFs at evenly spaced points (number determined by values within
       params dictionary).
       
       Requires: params <dictionary> (Must be properly loaded with needed keys:
                                          'x_size', 'y_size', 'x_psfs', 'y_psfs',
                                          'file_name_base', 'subsampling_factor', 'focus',
                                          'chip', 'detector', 'filter')
       
       Returns: (nothing)
       
       Side-effects: All side-effects of the generate_psf subroutine for each PSF generated
       
       Except: Will raise if params is not set up with all needed keys
    """
    
    try:
        for xi in xrange(0,params['x_psfs']):
            xp = (xi+0.5) * params['x_size'] / params['x_psfs']
            for yi in xrange(0,params['y_psfs']):
                yp = (yi+0.5) * params['y_size'] / params['y_psfs']
            
                psf_file = params['file_name_base'] + "_psf_" + str(xi) + "_" + str(yi) + ".fits"
                unbinned_file = params['file_name_base'] + "_unbinned_" + str(xi) + "_" + str(yi) + ".fits"
                
                generate_psf(params, xp, yp, psf_file, unbinned_file)
    except KeyError, e:
        raise Exception("params dictionary passed to load_psf_cache is missing needed key:\n" +
              str(e))
    
def initialise_moments_files(file_name_base):
    """Writes headers for the moments and moments_wings files.
       
       Requires: file_name_base (string for the root of the file name)
       
       Returns: moments_file, moments_wings_file (strings for the names of the files created)
       
       Side-effects: Overwrites moments_file and moments_wings_file and writes headers for them.
       
       Except: Will raise if the files cannot be accessed and written to
    """
    
    moments_file = file_name_base + "_moments.dat"
    with open(moments_file, "w") as fo:
        fo.write("# ID\tx_pix\ty_pix\tmonopole\tdipole_x\tdipole_y\t" + 
                 "quadrupole_xx\tquadrupole_xy\tmag\tstar_rad\tpsf_rad\t" +
                 "star_e1\tstar_e2\tpsf_e1\tpsf_e2\trad_diff\te1_diff\te2_diff\tstar_flux\n")
    moments_wings_file = file_name_base + "_moments_wings.dat"
    with open(moments_wings_file, "w") as fo:
        fo.write("# ID\tx_pix\ty_pix\tmonopole\tdipole_x\tdipole_y\t" + 
                 "quadrupole_xx\tquadrupole_xy\tmag\tstar_rad\tpsf_rad\t" +
                 "star_e1\tstar_e2\tpsf_e1\tpsf_e2\trad_diff\te1_diff\te2_diff\tstar_flux\n")
        
    return moments_file, moments_wings_file

def remove_psf_from_star(star, params):
    """Generates a postage stamp for a star, then removes the estimated psf from that stamp and
       generates a residual file. Also calculates moment data for the residual, and size and shape
       data for the star and psf.
       
       Requires: star <tuple> (tuple containing information on the star)
                 params <dictionary> (Must be properly loaded with needed keys:
                                          'file_name_base', 'quick_mode', 'image_file',
                                          'subsampling_factor', 'moments_file', 'moments_wings_file',
                                          'focus' (if params['quick_mode']==False))
                 
       Returns: stamp_file <string>, (name of postage stamp of star created)
                psf_file <string>, (name of undistorted psf file for this star)
                binned_file <string>, (name of rebinned distorted psf file for this star)
                residual_file <string> (name of residual file for this star, after psf subtracted off)
                
       Side-effects: All side-effects of generate_psf if params['quick_mode']==False
                     Creates and overwrites images for stamp_file, binned_file, and residual_file
                     All side-effects of cut_postage_stamp.cut_postage_stamp
                     All side-effects of subtract_psf.subtract_psf
                     
       Except: Will raise an exception if:
                   params is not set up with all needed keys
                   generate_psf raises an exception
                   cut_postage_stamp raises an exception
                   rebin_and_shift_psf raises an exception
                   subtract_psf raises an exception
    """
    
    try:    
        i = star[0]
        xp = star[1]
        yp = star[2]
        stamp_size = star[5]
        dx = star[6]
        dy = star[7]
        total_flux = star[8]
            
        
        binned_file = params['file_name_base'] + "_binned_" + str(i) + ".fits"
        stamp_file = params['file_name_base'] + "_stamp_" + str(i) + ".fits"
        residual_file = params['file_name_base'] + "_residual_" + str(i) + ".fits"
        
        if(params['quick_mode']):
            psf_file, unbinned_file = get_psf_file(params, xp, yp)
        else:
            psf_file = params['file_name_base'] + "_psf_" + str(i) + ".fits"
            unbinned_file = params['file_name_base'] + "_unbinned_" + str(i) + ".fits"
            
            generate_psf(params, xp, yp, psf_file, unbinned_file)
                
        # Cut a postage stamp around this star
        cut_postage_stamp(params['image_file'], (xp, yp), (stamp_size, stamp_size), stamp_file)
    
        # Rebin and cut the unbinned, distorted file for this star's offset
        rebin_and_shift_psf(unbinned_file, binned_file, params['subsampling_factor'], (dx, dy), stamp_size)
        
        # Subtract the psf and append the moments file with the moments of the residual
        try:
            subtract_psf(stamp_file, binned_file, residual_file, i, xp, yp, 
                        params['moments_file'], params['moments_wings_file'],total_flux,star[4])
        except KeyError, e:
            raise Exception("params dictionary passed to remove_psf_from_star is missing needed key:\n" +
                              str(e))
        except Exception, e:
            raise Exception("Cannot subtract psf for this stamp:\n" +
                        "Exception: " + str(e))
            
        return stamp_file, psf_file, binned_file, residual_file
    except KeyError, e:
        raise Exception("params dictionary passed to remove_psf_from_star is missing needed key:\n" +
                        str(e))
        
def get_psf_file(params, xp, yp):
    """Determines the name of the psf file to be used for this position in the image.
       This should be called only when params['quick_mode']==True
       
       Requires: params <dictionary> (Must be properly loaded with needed keys:
                                          'x_psfs', 'y_psfs', 'x_size', 'y_size',
                                          'quick_mode', 'file_name_base')
                 xp <int>
                 yp <int> (pixel coordinates in sextractor/tinytim order)
                 
       Returns: psf_file, unbinned_file (string names of undistorted psf file and distorted,
                    subsampled psf file)
                    
       Side-effects: (none)
       
       Except: Will raise an exception if:
                   params['quick_mode'] != True
                   params is not set up with all needed keys
    """
    
    try:
        if(not params['quick_mode']):
            raise Exception("get_psf_file must only be called when quick_mode is enabled.")
        
        xi = int(np.floor(xp * params['x_psfs'] / params['x_size']))
        yi = int(np.floor(yp * params['y_psfs'] / params['y_size']))
        
        # Check bounds
        if(xi < 0):
            xi = 0
        elif(xi >= params['x_psfs']):
            xi = params['x_psfs']-1
        if(yi < 0):
            yi = 0
        elif(yi >= params['y_psfs']):
            yi = params['y_psfs']-1
            
        psf_file = params['file_name_base'] + "_psf_" + str(xi) + "_" + str(yi) + ".fits"
        unbinned_file = params['file_name_base'] + "_unbinned_" + str(xi) + "_" + str(yi) + ".fits"
        
        return psf_file, unbinned_file
    except KeyError, e:
        raise Exception("params dictionary passed to get_psf_file is missing needed key:\n" +
                        str(e))
        

def remove_psf_no_fit(stars, params):
    """Removes the psf from a list of stars, using a fixed focus value, and calculates moments of the
       residuals and size and shape of the stars and psfs, printing and returning a chi-squared value.
       
       Requires: stars <list of tuples> (Each tuple in the list should contain information on a
                                         specific star)
                 params <dictionary> (Must be properly loaded with needed keys:
                                          'quick_mode', 'file_name_base', 'fit_focus',
                                          'x_size', 'y_size', 'x_psfs', 'y_psfs',
                                          'subsampling_factor', 'image_file',
                                          'focus' (if params['quick_mode']==False))
                 
       Returns: moments_file <string>, (name of file containing core moments)
                moments_wings_file <string>, (name of file containing wing moments)
                chi2 <float> (chi-squared value for the subtraction, assuming moments should be zero and
                    size/shape should be equal)
                    
       Side-effects: All side-effects of load_psf_cache
                     All side-effects of initialise_moments_files
                     All side-effects of remove_psf_from_star for each star in the list
                     
       Except: Will raise an exception if:
                   params is not set up with all needed keys
                   The moments files cannot be accessed and written to
                   remove_psf_from_star raises for all stars
                   The files needed for stack_residuals or draw_whisker_plot cannot be written to
       
    """
    
    try:
        # Set up psf cache if necessary
        if(params['quick_mode']):
            load_psf_cache(params)
                
        # Set up the moments residual files
        params['moments_file'], params['moments_wings_file'] = \
            initialise_moments_files(params['file_name_base'])
                
        processes = []
        # Now, go through each star, generate a PSF for it, and store the residual
        for star in stars:
            good_star_found = False
            try:
                # Call each thread
                p = mtp.Process(target=remove_psf_from_star, args=( star, params ) )
                p.start()
                processes.append(p)
                good_star_found = True
            except Exception, e:
                # raise # Uncomment this if you want to let any exception here result in a reraise
                removal_exception = e
                continue # Skip it. It's possible only certain files are write-protected
            
        # Join each thread back up so they can finish.
        for p in processes:
            p.join()
            
        if(not good_star_found):
            raise removal_exception
            # raise Exception("Cannot remove psfs from any stars in remove_psf_no_fit:\n" +
            #             str(removal_exception))
            
        
            
        # Set up ranges of various quantities within the moments files.
        # We don't actually need to make it this complicated, but it's here so that someone can more
        # easily to change to, for instance, only remove outliers on the monopole term. 
        mp_range = range(3,4)
        dp_range = range(4,6)
        qp_range = range(6,8)
        size_range = range(15,16)
        shape_range = range(16,18)
        
        dp_qp_range = dp_range+qp_range
        mp_dp_qp_range = mp_range+dp_qp_range
        size_shape_range = size_range+shape_range
        full_range = mp_dp_qp_range+size_shape_range
        
        remove_outliers_range = full_range
        
        # Dipole and quadrupole moments of core
        chi2_m_core, good_ids_and_fluxes = get_chi2(params['moments_file'],dp_qp_range,remove_outliers_range)
        
        # Dipole and quadrupole moments of wings +
        chi2_m_wings = get_chi2(params['moments_wings_file'],dp_qp_range,remove_outliers_range)[0]
        
        # Diff in radius and shape
        chi2_rad_shape = get_chi2(params['moments_file'],size_shape_range,remove_outliers_range)[0]
        
        chi2 = chi2_m_core+chi2_m_wings+chi2_rad_shape
        
        print("Chi^2 for core moments of this fit: " + str(chi2_m_core) + " (4 dof)")
        print("Chi^2 for wing moments of this fit: " + str(chi2_m_wings) + " (4 dof)")
        print("Chi^2 for size and shape: " + str(chi2_rad_shape) + " (3 dof)")
        print("Total Chi^2: " + str(chi2) + " (11 dof)")
        
        # Draw the whisker plot and stack residuals if we aren't fitting
        if(not params['fit_focus']):
            stack_residuals(params['file_name_base'], good_ids_and_fluxes)
            draw_whisker_plot(params['file_name_base'] + "_whisker.png",params['moments_file'],
                              good_ids_and_fluxes)
            draw_whisker_plot(params['file_name_base'] + "_whisker_wings.png",
                              params['moments_wings_file'],
                              good_ids_and_fluxes)
            
        return params['moments_file'], params['moments_wings_file'], chi2, good_ids_and_fluxes, \
            chi2_m_core, chi2_m_wings, chi2_rad_shape
    except KeyError, e:
        raise Exception("params dictionary passed to get_psf_file is missing needed key:\n" +
                        str(e))
        

def remove_psf_after_fit(stars, params):
    """Fits the best focus the set of stars, and uses it to generate the residual files and
       moments data files
       
       Requires: stars <list of tuples> (Each tuple in the list should contain information on a
                                         specific star)
                 params <dictionary> (Must be properly loaded with needed keys:
                                          'quick_mode', 'file_name_base', 'fit_focus',
                                          'x_size', 'y_size', 'x_psfs', 'y_psfs',
                                          'subsampling_factor', 'image_file',
                                          'focus' (if params['quick_mode']==False))
       
       Returns: (nothing)
       
       Except: Will raise an exception if:
                   params is not set up with all needed keys
                   The fitting procedure is unable to find any focii which give valid results,
                       possibly due to remove_psf_no_fit raising for all focii passed to it
    """

    try:
        # Using 'grid' fitting method
        
        # Magic numbers
        
        min_test_focus = -3.0
        max_test_focus = 3.0
        test_focus_points = 4
        target_precision = 0.25
        
        # Start by searching in a grid of evenly-spaced focus points
    
        best_chi2 = 1e99
        best_focus = -1e99
        
        for test_focus in np.linspace(min_test_focus,max_test_focus, num=test_focus_points):
            
            params['focus'] = test_focus
            try:
                # Call the removal method for this test focus
                unused_moments_file, unused_moments_wings_file, test_chi2, good_ids_and_fluxes, \
                    unused_chi2_core, unused_chi2_wings, unused_chi2_size_shape = \
                    remove_psf_no_fit(stars, params)
            except Exception, e:
                print(str(e))
            
            # Store this point if it's the best so far
            if(test_chi2 < best_chi2):
                best_chi2 = test_chi2
                best_focus = test_focus
                
        # Check if we found a suitable point to start from
        if(best_focus == -1e99):
            raise Exception("No suitable focus found.")
        
        # Now we go to a narrowing search
        
        # Initialise search_step and loop counter
        search_step = (max_test_focus-min_test_focus)/(test_focus_points-1) / 2.0
        loop_counter = 0
        
        # Go through a loop, finding the best of three points, then narrowing the search around it
        # until we reach the target precision
        while((search_step > target_precision/2) and (loop_counter < 100)):
            
            # Increment loop counter
            loop_counter += 1
            
            for test_focus in [best_focus-search_step,best_focus+search_step]:
                
                test_state = True
                
                params['focus'] = test_focus
                try:
                    # Call the removal method for this test focus
                    unused_moments_file, unused_moments_wings_file, test_chi2, good_ids_and_fluxes, \
                        best_chi2_core, best_chi2_wings, best_chi2_size_shape = \
                        remove_psf_no_fit(stars, params)
                
                    # Calculate the chi-squared for a null test of the residuals
                except Exception, e:
                    print(str(e))
                
                # Store this point if it's the best so far
                if(test_chi2 < best_chi2):
                    best_chi2 = test_chi2
                    best_focus = test_focus
                    test_state = False
    
            # Narrow the search step for the next loop
            search_step /= 4
            
        if(test_state):
            
            # In this case, the last results generated weren't the best, so regenerate the best results
            
            params['focus'] = best_focus
            try:
                # Call the removal method for this test focus
                unused_moments_file, unused_moments_wings_file, best_chi2, good_ids_and_fluxes, \
                        best_chi2_core, best_chi2_wings, best_chi2_size_shape = \
                    remove_psf_no_fit(stars, params)
            except Exception, e:
                print(str(e))
            
        # Print the results
        print("Best focus: " + str(best_focus))
        print("Best chi^2: " + str(best_chi2))
        
        # Output to summary file
        if(params['summary_file_name'] is not None):
            cmd = "echo '" + params['image_file'] + "\t" + str(best_focus) + "\t" + str(best_chi2) + \
                    "\t" + str(best_chi2_core) + "\t" + str(best_chi2_wings) + "\t" + str(best_chi2_size_shape) + \
                    "' >> " + params['summary_file_name']
            sbp.call(cmd,shell=True)
        
        # Draw the whisker plot now
        stack_residuals(params['file_name_base'], good_ids_and_fluxes)
        draw_whisker_plot(params['file_name_base'] + "_whisker.png",params['moments_file'],
                              good_ids_and_fluxes)
        draw_whisker_plot(params['file_name_base'] + "_whisker_wings.png",
                          params['moments_wings_file'],
                              good_ids_and_fluxes)
    except KeyError, e:
        raise Exception("params dictionary passed to get_psf_file is missing needed key:\n" +
                        str(e))
    
def flux_ratio(n):
    """Calculates the proper ratio for total_flux/flux_auto
    
       Requires: n <float> (Sersic index for profile)
    
       Returns: ratio <float> (total_flux/flux_auto)
       
       Side-effects: (none)
    """
    
    # Interpolate using Graham and Driver (2005)'s table
    
    n_list = [0.5,1,2,3,4,5,6,7,8,9,10]
    r_list = [.993,.96,.922,.908,.904,.905,.907,.910,.914,.919,.923]
    
    i = 0
    while((n_list[i] < n) and (i < len(n_list)-2)):
        i += 1
        
    ratio = r_list[i] + (n-n_list[i]) * ( r_list[i+1]-r_list[i] ) / ( n_list[i+1]-n_list[i] )
    
    if(ratio > 1):
        return 1.
    
    return 1./ratio

def stack_residuals(filename_base, ids_and_fluxes):
    """Generates a fits image which is a stack of the residuals of non-outlier stars.
    
       Requires: filename_base <string> (base for the residual files and output image)
                 ids_and_fluxes <list of tuples> (the first index of each should be the id of each
                                 star to use in the stacking of residuals, and the second should be
                                 the flux of that star (used for normalization))
                              
       Returns: (nothing)
       
       Side-effects: Writes/overwrites filename_base + '_residual_stack.fits' with a new stacked
                     residual file.
    
       Except: Will raise an exception if:
                      Any of the required residual files aren't accessible and readable
                      The ids_and_fluxes list isn't set up properly
                      Possibly if one of the residual files has an even size or isn't square
    """
    
    stacked_data = None
    
    for id_and_flux in ids_and_fluxes:
        
        star_id = id_and_flux[0]
        flux = id_and_flux[1]
        
        filename = filename_base + "_residual_" + str(int(star_id)) + ".fits"
        
        fits_struct = read_fits(filename)
        data = fits_struct[0].data
            
        if(stacked_data is None):
            stacked_data = data/flux
        else:
            shape_stacked = np.shape(stacked_data)
            shape_data = np.shape(data)
            try:
                if(shape_stacked[0]>shape_data[0]):
                    stacked_data = np.add(data/flux,stacked_data[(shape_stacked[0] - shape_data[0]) / 2:
                                                shape_stacked[0] - ((shape_stacked[0] - shape_data[0]) / 2),
                          (shape_stacked[1] - shape_data[1]) / 2:
                              shape_stacked[1] - ((shape_stacked[1] - shape_data[1]) / 2)
                          ])
                else:
                    stacked_data = np.add(stacked_data,data[(shape_data[0] - shape_stacked[0]) / 2:
                              shape_data[0] - ((shape_data[0] - shape_stacked[0]) / 2),
                          (shape_data[1] - shape_stacked[1]) / 2:
                              shape_data[1] - ((shape_data[1] - shape_stacked[1]) / 2)
                          ]/flux)
            except Exception, e:
                raise Exception("Likely due to size mismatch in stacking residuals.\n" +
                      "Ensure all residuals have odd sizes and are square:\n" +
                      str(e))

    
    out_fits = pyf.PrimaryHDU(stacked_data)
    out_fits.header = fits_struct[0].header
    out_fits.header.add_comment('Residual stack')
    out_fits.writeto(filename_base + "_residual_stack.fits", clobber=True)

def draw_whisker_plot(whisker_plot_filename, moments_file, good_ids_and_fluxes=None):
    """Draws a whisker plot for the ellipticity difference between a star and PSF.
       
       Requires: moments_file <string> (name of file containing moments information)
       
       Returns: (nothing)
       
       Side-effects: Overwrites whisker_plot_filename with new whisker plot
       
       Except: Will raise an exception if:
                      Possibly if the moments file isn't properly formatted
                      whisker_plot_filename cannot be written to
    """
    
    # Read in the moments file to get the g1, g2, x, y, and size lists
    objects = []
        
    with open(moments_file, 'r') as f:

        # Read in the file, except for comment lines
        for line in f:
            line.strip()
            if((line[0] != '#') and (len(line) > 0)):
                o = []
                for string in line.split():
                    o.append(float(string))
                if(good_ids_and_fluxes is None):
                    objects.append(o)
                else:
                    # Check if we find the id. If so, append
                    for id_and_flux in good_ids_and_fluxes:
                        if(float(id_and_flux[0])==o[0]):
                            objects.append(o)
                            break
                
    g1 = []
    g2 = []
    x = []
    y = []
    size = []
    
    for o in objects:
        g1.append(o[16])
        g2.append(o[17])
        x.append(o[1])
        y.append(o[2])
        size.append(o[15])
        
    draw_whisker_shear(g1, g2, x, y, size=size, title="Shear differentials",
                        filename=whisker_plot_filename, magnify=1.5, adjust_power=0.25)

if __name__ == "__main__":
    main(sys.argv)
