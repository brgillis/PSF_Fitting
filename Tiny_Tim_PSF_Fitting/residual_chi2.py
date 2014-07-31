import numpy as np
import copy
from scipy.stats import norm


def remove_outliers(objects,ilist,min_remaining_members=2,p_threshold=0.5):
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
    
    while(outliers_found):

        outliers_found = False
            
        # Find any outliers and remove them
        
        for index in ilist:
        
            mean, sigma = get_mean_sigma(objects, index)
            
            for o in objects:
                
                p = 1-norm.cdf(abs(o[index] - mean) / sigma)
                
                if( p*num_tot*ifactor < p_threshold ):
                    objects.remove(o)
                    outliers_found = True
            
    
    # Check that we still have enough members left in the array
    if(len(objects) < min_remaining_members):
        raise Exception("WARNING: Too few arguments left after removing outliers. (Only " + str(len(objects)) + ")"
                        "\nReverting to non-trimmed dataset.\n" + 
                        "Chi-squared value will include presence of outliers.")
    else:
        print("Removed " + str(num_tot-len(objects)) + "/" + str(num_tot) + " outliers.")
    
def get_mean_sigma(objects, index):
    """ Gets mean and sigma (std deviation) for a column with a list of lists of floats.
    
        Requires: objects <list of lists of floats>
                  index <int> (Index to measure mean and sigma for)
                  
        Returns: mean <float>,
                 sigma <float>
                 
        Side-effects: (none)
    """
    
    mean = 0
    mean_sq = 0
    num = 0
    
    for o in objects:
        mean += o[index]
        mean_sq += np.square(o[index])
        num += 1
        
    if(num == 0):
        mean = 0
        mean_sq = 0
    else:
        mean /= num
        mean_sq /= num
        
    sigma = np.sqrt(mean_sq - np.square(mean))
    
    return mean, sigma

def get_chi2_for_index(objects, index, target):
    """ Calculates a Chi^2 value for a column within a list of lists of floats.
    
        Requires: objects <list of lists of floats>
                  index <int> (Index to measure mean and sigma for)
                  target <float> (Target value which would give Chi^2 == 0 if all values equal it)
                  
        Returns: chi2 <float> (Chi^2 value)
        
        Side-effects: (none)
        
    """
    mean, sigma = get_mean_sigma(objects, index)
    
    chi2 = np.square((mean-target)/(sigma/np.sqrt(len(objects))))
    
    return chi2

def get_chi2(file_name,ilist_chi2,ilist_outliers=None):
    """ Calculates a chi-squared value for the data contained within a file in certain columns,
        under the assumption that the values should ideally equal zero.
        
        Requires: file_name <string> (Must be standard ASCII catalog format with same number of values
                                      in each line.)
                  ilist_chi2 <list of ints> (indices of columns within the file to test)
                  ilist_outliers <list of ints> (indices of columns within the file to check for
                                                 outliers, which will be removed before chi-squared
                                                 calculation.)
                                                 
        Returns: chi2 <float>, (Resulting Chi-squared value)
                 good_ids_and_fluxes <list of (int, float)s> (list of tuples, each of which contains
                                                              an ID of a non-outlier galaxy in the
                                                              first index, and its flux in the second
                                                              index.)
                                                              
        Side-effects: (none)
        
        Except: file_name cannot be accessed or is improperly formatted
    
    """
    
    # Magic numbers
    target = 0 # Target value columns should equal for Chi^2 == 0.
    
    objects = np.loadtxt(file_name).tolist()
    
    if(ilist_outliers != None):
        
        # Trim outliers, but make sure we retain a quarter of the members
        # An exception will be thrown if we don't; if that's the case, restore backup
        min_remaining_members = len(objects)/4.
        
        # Create a deep copy of objects as a backup, in case something goes wrong
        objects_backup = copy.deepcopy(objects)
        
        try:
            remove_outliers(objects,ilist_outliers,min_remaining_members)
        except Exception, e:
            objects = objects_backup
            print(str(e))
        
    good_ids_and_fluxes = []    
    for o in objects:
        good_ids_and_fluxes.append((o[0],o[18]))
    
    # Sum up the chi2 for all dipole and quadrupole components
    chi2 = 0
    
    for index in ilist_chi2:
        chi2 += get_chi2_for_index(objects, index, target)
    
    return chi2, good_ids_and_fluxes
    
    