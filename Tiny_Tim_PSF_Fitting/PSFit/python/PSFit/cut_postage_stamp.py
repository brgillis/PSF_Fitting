#!/usr/bin/env python2
#
# create psf postage stamps at a given position in the image frame
#
# Created by Ole Margraff 2014
# Edited by Bryan Gillis 2014

import numpy as np
import pyfits as pyf
from fits_functions import read_fits, write_fits


def extract_pstamp(fits, pos, size):
    """Extracts a postage stamp of given size and 
       at given position from the FITS file.

       Requires: fits <HDUList> (FITS data structure (Header & Data))
                 pos  (<int>,<int>) (Tuple of x,y coordinates for the postage stamp)
                 size (<int>,<int>) (Tuple with x,y size of the postage stamp)

       Returns:  stamp <HDUList> (Header & Data of <size>x<size>)
       
       Side-effects: (none)
       
       Except: Will raise an exception if the position is too close to an edge for a stamp to be
               extracted.
    """

    #print('Extracting postage stamp.')

    
    # Swap x/y due to numpy internals
    x1=pos[0]-np.floor_divide(size[0]-1, 2)
    y1=pos[1]-np.floor_divide(size[1]-1, 2)

    # Create the stamp data. This may fail if the stamp is too large and near and edge
    try:
        stamp_data = fits[0].data[y1-1:y1+size[1]-1, x1-1:x1+size[0]-1]
    except:
        print("ERROR: Cannot create postage stamp. Check that it isn't too large and/or close\n" +
              "to an edge.")
        raise

    stamp = pyf.PrimaryHDU(stamp_data)

    stamp.header.add_comment(
        "Postage stamp cut from file at image position (%.2f,%.2f)." % (pos[0],pos[1]))

    return stamp


def cut_postage_stamp(infile, pos, size, outfile):
    """Extracts and saves a postage stamp around a given position in a given file.
    
       Requires: infile <string> (name of file from which to extract the stamp)
                 pos  (<int>,<int>) (Tuple of x,y coordinates for the postage stamp)
                 size (<int>,<int>) (Tuple with x,y size of the postage stamp)
                 outfile <string> (name of the postage stamp file to create)
                 
       Returns: (nothing)
       
       Side-effects: Overwrites outfile with a new postage stamp
    
       Except: Will raise an exception if:
                   infile cannot be accessed
                   outfile cannot be accessed and written to
                   The position is too close to the edge of the image to extract a stamp
    """

    fitsfile = read_fits(infile)
    postage_stamp = extract_pstamp(fitsfile, pos, size)
    write_fits(postage_stamp, outfile)
