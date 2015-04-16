import numpy as np
import copy

from fits_functions import write_fits
from remove_outliers import remove_outliers

def get_noise_array(counts, background_std=1.):
    return np.random.normal(size=np.shape(counts)) * np.sqrt(np.square(background_std) + counts)

def add_noise_from_stamp(stamp_data, psf_data, total_flux=1., gain=2.0):
    
    stamp_size = np.shape(stamp_data)
    
    # Go around the edge of the stamp, generating a list of values
    edge_data = []
    
    for v in stamp_data[:,0]:
        edge_data.append(v)
    for v in stamp_data[:,stamp_size[1]-1]:
        edge_data.append(v)
    for v in stamp_data[0,1:stamp_size[1]-2]:
        edge_data.append(v)
    for v in stamp_data[stamp_size[0]-1,1:stamp_size[1]-2]:
        edge_data.append(v)
        
    edge_data_copy = copy.deepcopy(edge_data)
    try:
        edge_std = remove_outliers(edge_data)[1]
    except:
        edge_data = edge_data_copy
        edge_std = np.std(edge_data)
        
    # Seed numpy from the first value in the array. This keeps it deterministic for fitting, but
    # random on a per-stamp case
    np.random.seed(int(abs(edge_data[0]*10000)))
    
    noise_data = get_noise_array(psf_data*total_flux/gain,edge_std/gain)*gain
     
    psf_data += noise_data/total_flux
    

def add_and_save_noise( stamp_struct, psf_struct, noisy_psf_file_name, total_flux=1., gain=2.0 ):
    """ Adds noise to a psf image, based on information from a star image.
    
        Requires: stamp_struct <HDUList> (stamp image's fits data structure)
                  psf_struct <HDUList> (psf image's fits data structure)
                  noise_psf_file_name <string> (file name in which to save the psf+noise)
        Optional: total_flux <float>  (Flux of the star, used to scale noise)
                  gain <float> (gain of the image)
        
        Returns: (nothing)
        
        Side-effects: psf_struct has noise added
                      Overwrites noisy_psf_file_name with new image (if not None)
        
        Except: noisy_psf_file_name cannot be written to
        
    """
    
    stamp_data = stamp_struct[0].data
    psf_data = psf_struct[0].data
    
    add_noise_from_stamp(stamp_data, psf_data, total_flux, gain)
    
    if(noisy_psf_file_name is not None):
        write_fits(psf_struct,noisy_psf_file_name)