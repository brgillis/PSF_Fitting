#!/usr/bin/env python2
#
# Add a comment to a FITS header
#
# Created by Ole Marggraf 2013-12-10
# Edited by Bryan Gillis 2014

import pyfits as pyf

def read_fits(filename):
    """ Read the FITS file

        Requires: filename <string>

        Returns:  fits_struct <HDUList> (FITS data structure (Header & Data))
       
        Side-effects: (none)
       
        Except: Will raise an exception if the named file cannot be opened for some reason.
    """

    try:
        fits_struct = pyf.open(filename)
    except Exception, e:
        raise Exception("Could not open FITS file " + filename + ":\n" + str(e))

    return fits_struct


def add_comment(data_struct, text):
    """ Adds a comment line to a FITS header structure.

        Requires: data_struct <HDUList> FITS data structure
                  text <string>         Comment line to be added
                 
        Returns: (nothing)
        
        Side-effects: Adds a comment line to data_struct's header
    """

    data_struct[0].header.add_comment(text)

    return


def write_fits(fits_data, filename):
    """ Writes FITS data into file
    
        Requires: fits_data <HDUList> (FITS data structure (Header & Data))
                  filename <string> (File to write to)
                 
        Returns: (nothing)
       
        Side-effects: Overwrites filename with a new fits file
       
        Except: Will raise an exception if the named file cannot be accessed and written to.
    """

    try:
        fits_data.writeto(filename, clobber=True)
    except:
        raise Exception("Something went wrong writing the output file: " + filename)

    return
