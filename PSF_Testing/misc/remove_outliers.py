""" @file remove_outliers.py

    Created 17 Sep 2015

    Function to remove outlier values from an array.

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

import numpy as np
from scipy.stats import norm

def remove_outliers(array,min_remaining_members=2):
    """ Removes outliers from a list of values, using Chauvenet's criterion.
    
        Requires: array <ndarray>
        Optional: min_remaining_membeters <int> (Minimum number of values that must remain. If we end
                                                 up with fewer than this, an exception will be raised.)
        
        Returns: <masked array> (array with outliers masked)
        
        Except: Less than min_remaining_members values remain in olist after outlier removal                            
    """
    
    assert(min_remaining_members>=2)
    
    marray = np.ma.masked_array(array) # Create a masked array from the array
    
    outliers_found = True
    num_tot = np.size(marray)
    
    while(outliers_found and (len(marray[~marray.mask])>2)):

        outliers_found = False
            
        # Find any outliers and remove them
        
        unmasked_elements = marray[~marray.mask]
        
        mean, sigma = np.mean(unmasked_elements), np.std(unmasked_elements)
        
        probs = 1-norm.cdf(abs(marray - mean) / sigma)
        
        outliers = probs*num_tot < 0.5
        
        # Do we have any new outliers?
        if(np.any(outliers[~marray.mask])):
            # Mark these
            marray.mask = np.logical_or(marray.mask,outliers)
            
            # Loop through again
            outliers_found = True
    
    # Check that we still have enough members left in the array
    if(len(marray[~marray.mask]) < min_remaining_members):
        raise Exception("WARNING: Too few arguments left after removing outliers. (Only " + str(len(list)) + ")")
    else:
        return marray