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
        #print("Removed " + str(num_tot-len(olist)) + "/" + str(num_tot) +  " outliers.")
        return mean, sigma
    

def remove_column_outliers(objects,ilist,min_remaining_members=2,p_threshold=0.5):
    """ Removes outliers from a set of objects, using Chauvenet's criterion.
    
        Requires: objects <list of lists of floats>
                  ilist <list of ints> (Indices of columns to check for outliers)
        Optional: p_threshold <float> (Minimum probability that we'd get a value this extreme for
                                       which we'll cut an object as an outlier)
                  min_remaining_members <int> (Minimum number of members that must remain after
                                               searching for outliers. If we don't have this many,
                                               an exception will be raised.)
                                
        Returns: (nothing)
                       
        Side-effects: Removes all outliers from objects list
        
        Except: Indices passed in ilist are out of range of the lists contained within objects
                Too few members remain after removing outliers
    """
    
    
    outliers_found = True
    num_tot = len(objects)
    
    ifactor = len(ilist)
    
    while(outliers_found and (len(objects)>2)):

        outliers_found = False
            
        # Find any outliers and remove them
        
        means = np.mean(objects,axis=0)
        stds = np.std(objects,axis=0)
        
        for index in ilist:
        
            mean = means[index]
            std = stds[index]
            
            for o in objects:
                
                p = 1-norm.cdf(abs(o[index] - mean) / std)
                
                if( p*num_tot*ifactor < p_threshold ):
                    objects.remove(o)
                    outliers_found = True
            
    
    # Check that we still have enough members left in the array
    if(len(objects) < min_remaining_members):
        raise Exception("WARNING: Too few arguments left after removing outliers. (Only " + str(len(objects)) + ")"
                        "\nReverting to non-trimmed dataset.\n" + 
                        "Chi-squared value will include presence of outliers.")
    else:
        #print("Removed " + str(num_tot-len(objects)) + "/" + str(num_tot) + " outliers.")
        return