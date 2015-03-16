import numpy as np
import copy

from remove_outliers import remove_column_outliers
    
def get_mean_and_std_for_index(objects, index):
    """ Gets mean and sigma (std deviation) for a column with a list of lists of floats.
    
        Requires: objects <list of lists of floats>
                  index <int> (Index to measure mean and sigma for)
                  
        Returns: mean <float>,
                 sigma <float>
                 
        Side-effects: (none)
    """
    
    mean = np.mean(objects,axis=0)[index]
    std = np.std(objects,axis=0)[index]
    
#     mean = 0
#     mean_sq = 0
#     num = 0
#     
#     for o in objects:
#         mean += o[index]
#         mean_sq += np.square(o[index])
#         num += 1
#         
#     if(num == 0):
#         mean = 0
#         mean_sq = 0
#     else:
#         mean /= num
#         mean_sq /= num
#         
#     sigma = np.sqrt(mean_sq - np.square(mean))
    
    return mean, std

def get_chi2_and_mean_for_index(objects, index, target):
    """ Calculates a Chi^2 value for a column within a list of lists of floats.
    
        Requires: objects <list of lists of floats>
                  index <int> (Index to measure mean and sigma for)
                  target <float> (Target value which would give Chi^2 == 0 if all values equal it)
                  
        Returns: chi2 <float> (Chi^2 value)
        
        Side-effects: (none)
        
    """
    mean, std = get_mean_and_std_for_index(objects, index)
    
    chi2 = np.square((mean-target)/(std/np.sqrt(len(objects))))
    
    return chi2, mean

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
            remove_column_outliers(objects,ilist_outliers,min_remaining_members)
        except Exception, e:
            objects = objects_backup
            print(str(e))
        
    good_ids_and_fluxes = []    
    for o in objects:
        good_ids_and_fluxes.append((o[0],o[18]))
    
    # Sum up the chi2 for all dipole and quadrupole components
    chi2s = []
    means = []
    
    for index in ilist_chi2:
        chi2, mean = get_chi2_and_mean_for_index(objects, index, target)
        chi2s.append(chi2)
        means.append(mean)
    
    chi2 = np.sum(chi2s)
    
    return chi2, good_ids_and_fluxes, chi2s, means
    
    