""" call_tinytim.py

    Created by Bryan Gillis 31 July, 2014
"""

import subprocess as sbp
from astropy.io import fits as pyf
from os.path import isfile

# Magic values
tinytim_path = "/home/user/Program_Files/tinytim-7.5"
psf_width = "3.0"
        
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

def get_psf(params, xp=-1, yp=-1, psf_file="temp_psf.fits", subsampled_file="test_subsampled.fits"):
    """ This function calls Tiny Tim to generate a PSF at a given pixel position (xp, yp),
        using the information in the passed params dictionary.
        
        Requires: params <dictionary> (parameter dictionary, properly set up with needed keys:
                                           'file_name_base', 'subsampling_factor', 'focus',
                                           'chip', 'detector', 'filter'. If xp and yp are not
                                           provided, then the keys 'x_size' and 'y_size' must
                                           be present in the dictionary.
        Optional: xp <int> (x pixel position for the PSF to be generated for)
                  yp <int> (y pixel '')
        
        Returns: Nothing
        
        Side-effects: Cleans and regenerates tinytim intermediate and output files
                      Creates and stores psf and subsampled files
        
        Except: The params dictionary is missing needed keys
                Tiny Tim fails to generate the expected psf files (possibly due to write-
                    protection or out of memory issues)
    """ 
    
    # Clean previous tinytim files if they exist
    cmd = "rm -f " + params['file_name_base'] + ".tt3 " + params['file_name_base'] + ".par " + \
          params['file_name_base'] + "00.fits" + params['file_name_base'] + "00_psf.fits"
    sbp.call(cmd, shell=True)
    
    # Check for default pixel values being passed in
    # If so, use the centre of the image
    if(xp==-1):
        xp = params['x_size']/2
    if(yp==-1):
        yp = params['y_size']/2
        
    # Set up the command to call tiny1
    cmd = "export TINYTIM=" + tinytim_path + "\n" + \
          tinytim_path + "/tiny1 " + params['file_name_base'] + ".par << EOF \n" + \
          str(params['detector']) + "\n" + \
          str(params['chip']) + "\n" + \
          str(xp) + " " + str(yp) + "\n" + \
          str(params['filter'])  + "\n" + \
          "2\n5000\n" + \
          psf_width + "\n" + \
          str(params['focus']) + "\n" + \
          params['file_name_base'] + "\nEOF"
    # Run the command to call tiny1
    sbp.call(cmd, shell=True)
    
    # Set up the command to call tiny2
    cmd = "export TINYTIM=" + tinytim_path + "\n" + \
          tinytim_path + "/tiny2 " + params['file_name_base'] + ".par"
    # Run the command to call tiny2
    sbp.call(cmd, shell=True)
    
    # Set up the command to call tiny3
    cmd = "export TINYTIM=" + tinytim_path + "\n" + \
          tinytim_path + "/tiny3 " + params['file_name_base'] + ".par SUB=" + \
        str(params['subsampling_factor'])
    # Run the command to call tiny3
    sbp.call(cmd, shell=True)
        
    # Check if the output from tinytim exists. If not, skip this star
    if(not isfile(params['file_name_base'] + "00_psf.fits")):
        raise Exception("Cannot create tinytim psfs.")
    if(not isfile(params['file_name_base'] + "00.fits")):
        raise Exception("Cannot create tinytim psfs.")
    
    # Move the psf file to store it
    cmd = "mv " + params['file_name_base'] + "00_psf.fits " + psf_file
    sbp.call(cmd, shell=True)
    
    # Move the unbinned distorted psf file to store it
    cmd = "mv " + params['file_name_base'] + "00.fits " + subsampled_file
    sbp.call(cmd, shell=True)
    
    return