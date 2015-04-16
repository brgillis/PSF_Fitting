""" call_tinytim.py

    Created by Bryan Gillis 31 July, 2014
"""

import subprocess as sbp

# Magic values
tinytim_path = "/home/brg/Program_Files/tinytim-7.5"
psf_width = "3.0"

#spec = (2, 5000) # Blackbody, 5000 K
spec = (1, 15) # K-type star
#spec = (1, 17)

def call_tinytim(params, xp=-1, yp=-1):
    """ This function calls Tiny Tim to generate a PSF at a given pixel position (xp, yp),
        using the information in the passed params dictionary.
        
        Requires: params <dictionary> (parameter dictionary, properly set up with needed keys:
                                           'file_name_base', 'subsampling_factor', 'focus',
                                           'chip', 'detector', 'filter'. If xp and yp are not
                                           provided, then the keys 'x_size' and 'y_size' must
                                           be present in the dictionary.
        Optional: xp <int> (x pixel position for the PSF to be generated for)
                  yp <int> (y pixel '')
        
        Returns: ???
        
        Side-effects: Cleans and regenerates tinytim intermediate and output files
                      Creates and stores psf and subsampled files
        
        Except: The params dictionary is missing needed keys
                Tiny Tim fails to generate the expected psf files (possibly due to write-
                    protection or out of memory issues)
    """ 
    
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
          str(spec[0]) + "\n" + \
          str(spec[1]) + "\n" + \
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