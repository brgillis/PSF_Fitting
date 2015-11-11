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
    
    stamp_size = 2*star.SkyObj.stamp_size + 1
    
    # Swap x/y due to fits ordering
    x1=star.SkyObj.y_pix-star.SkyObj.stamp_size
    y1=star.SkyObj.x_pix-star.SkyObj.stamp_size
    
    image_nx, image_ny = np.shape(image.data)
    
    # Check that the stamp isn't too close to an edge
    assert((x1>1) and (y1>1) and (x1+stamp_size-1<image_nx) and (y1+stamp_size-1<image_ny))

    # Create the stamp data
    stamp = image.data[x1-1:x1+stamp_size-1,
                       y1-1:y1+stamp_size-1]

    return stamp
