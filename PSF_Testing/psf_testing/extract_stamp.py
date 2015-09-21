#!/usr/bin/env python2
#
# create psf postage stamps at a given position in the image frame
#
# Created by Ole Margraff 2014
# Edited by Bryan Gillis 2014, 2015

import numpy as np

def extract_stamp_for_star(star, image):
    """Extracts a postage stamp of a star from the FITS image.

       Requires: star <star> (A star object, which contains the necessary info for the
                              position and size of the postage stamp)
                 image <HDU> (An astropy HDU containing the image data)

       Returns:  stamp <ndarray> 
       
       Side-effects: (none)
       
       Except: Will raise an exception if the position is too close to an edge for a stamp to be
               extracted.
    """
    
    # Swap x/y due to fits ordering
    x1=star.sky_object.y_pix-np.floor_divide(star.sky_object.stamp_size-1, 2)
    y1=star.sky_object.x_pix-np.floor_divide(star.sky_object.stamp_size-1, 2)

    # Create the stamp data. This may fail if the stamp is too large and near and edge
    try:
        stamp = image.data[x1-1:x1+star.sky_object.stamp_size-1,
                           y1-1:y1+star.sky_object.stamp_size-1]
    except:
        print("ERROR: Cannot create postage stamp. Check that it isn't too large and/or close\n" +
              "to an edge.")
        raise

    return stamp
