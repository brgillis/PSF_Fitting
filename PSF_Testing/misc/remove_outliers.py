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

def remove_outliers(olist,min_remaining_members=2):
    """ Removes outliers from a list of values, using Chauvenet's criterion.
    
        Requires: olist <list of floats>
        Optional: min_remaining_membeters <int> (Minimum number of values that must remain. If we end
                                                 up with fewer than this, an exception will be raised.)
        
        Returns: mean <float>, (mean of the list, after outlier removal)
                 sigma <float> (sigma of the list, after outlier removal)
                 
        Side-effects: Outliers removed from olist
        
        Except: Less than min_remaining_members values remain in olist after outlier removal                            
    """
    
    if(min_remaining_members<2):
        min_remaining_members = 2
    
    outliers_found = True
    num_tot = len(olist)
    
    while(outliers_found and (len(olist)>2)):

        outliers_found = False
            
        # Find any outliers and remove them
        
        mean, sigma = np.mean(olist), np.std(olist)
        
        for v in olist:
            
            p = 1-norm.cdf(abs(v - mean) / sigma)
            
            if( p*num_tot < 0.5 ):
                olist.remove(v)
                outliers_found = True
    
    # Check that we still have enough members left in the array
    if(len(olist) < min_remaining_members):
        raise Exception("WARNING: Too few arguments left after removing outliers. (Only " + str(len(list)) + ")")
    else:
        return mean, sigma