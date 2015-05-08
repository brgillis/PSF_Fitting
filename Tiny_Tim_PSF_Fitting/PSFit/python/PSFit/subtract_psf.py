#!/usr/bin/env python2
#
# create psf postage stamps at a given position in the image frame
#
# OM 2014

import numpy as np
import pyfits as pyf
import copy
import os

import multipole_moments as mpm
from add_noise import add_and_save_noise
from fits_functions import read_fits, write_fits, add_comment
from remove_outliers import remove_outliers

# some magic numbers
aperture_size = 10
cut_off_scale = 5
wing_aperture_size = 20

def remove_psf(fits_data, psf_data, total_flux, calc_chi2=False, background_noise=0, gain=2.0):
    """ Subtracts the PSF from the data, scaling it by total_flux.
    
        Requires: fits_data <HDUList> (star's fits data structure)
                  psf_data <HDUList> (psf's fits data structure)
                  total_flux <float> (flux of the star, used here to scale the psf)
                  
        Optional: calc_chi2 <bool> (whether or not to calculate chi-squared for comparison of star
                                    with psf)
                  background_noise <float> (standard deviation in background level of star image)
                  gain <float> (gain of the image)
                  
        Returns: out_fits <HDUList> (residual's fits data structure),
                 chi2 <float> (chi-squared of comparison of star with PSF model)
        
        Side-effects: (none)
        
        Except: Either the star or PSF image is 1 pixel or less in size
    """
    
    data1 = fits_data[0].data
    data2 = psf_data[0].data

    size1 = np.shape(data1)
    size2 = np.shape(data2)
    
    # Check we have good data for both images
    if((size1[0] <= 1) or (size1[1] <= 1) or (size2[0] <= 1) or (size2[1] <= 1)):
        raise Exception("Bad size for psf or stamp image.")

    # we assume for now that data is square
    if size1[0] > size2[0]:
        data1 = data1[(size1[0] - size2[0]) / 2:
                          size1[0] - ((size1[0] - size2[0]) / 2),
                      (size1[1] - size2[1]) / 2:
                          size1[1] - ((size1[1] - size2[1]) / 2)
                      ]
        size1 = np.shape(data1)
        if(size1[0] > size2[0]):
            data1 = data1[1:,:]
            size1 = np.shape(data1)
        if(size1[1] > size2[1]):
            data1 = data1[:,1:]
            size1 = np.shape(data1)
    elif size1[0] < size2[0]:
        data2 = data2[(size2[0] - size1[0]) / 2:
                          size2[0] - ((size2[0] - size1[0]) / 2),
                      (size2[1] - size1[1]) / 2:
                          size2[1] - ((size2[1] - size1[1]) / 2)
                      ]
        size2 = np.shape(data2)
        if(size2[0] > size1[0]):
            data2 = data2[1:,:]
            size2 = np.shape(data2)
        if(size2[1] > size1[1]):
            data2 = data2[:,1:]
            size2 = np.shape(data2)
            
    # Check both images are now the same size
    if((size1[0] != size2[0]) or (size1[1] != size2[1])):
        # Should never occur, but just in case
        raise Exception("Could not cut images to same size.")

    scaled_psf_data = data2 * total_flux

    out = data1 - scaled_psf_data

    out_fits = pyf.HDUList(pyf.PrimaryHDU(out))
    #out_fits[0].header = fits_data[0].header
    add_comment(out_fits,"Postage stamp removed from central structure")
    
    if(calc_chi2):
        chi2 = calculate_chi2(data1,scaled_psf_data,background_noise,gain)
    else:
        chi2 = 0

    return out_fits, chi2
    

def remove_background( stamp_struct, background=None ):
    """ Removes mean background from a star's postage stamp, based on edge data
    
        Requires: stamp_struct <HDUList> (stamp image's fits data structure)
        
        Returns: (nothing)
        
        Side-effects: mean background removed from stamp_struct
        
    """
    
    stamp_data = stamp_struct[0].data
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
        edge_mean, edge_std = remove_outliers(edge_data)
    except:
        edge_data = edge_data_copy
        edge_mean, edge_std = np.mean(edge_data), np.std(edge_data)
        
    if(background is None):
        background = edge_mean
        
    stamp_struct[0].data = np.subtract(stamp_data,background)
    
    return edge_std


def subtract_psf(data_file_name, psf_file_name, residual_file_name, ID, x_pix, y_pix,
                 moments_file_name, moments_wings_file_name,total_flux,mag,calc_chi2=False,gain=2.0, background=None):
    """ Subtracts the psf from the star's image and calculates and stores multipole moment
        information - this is the function that should be called from external files.
        
        Requires: data_file_name <string> (Name of the star's postage stamp file)
                  psf_file_name <string> (Name of the psf's postage stamp file (after shifting and
                                          rebinning))
                  residual_file_name <string> (Name of the file in which to save the residual image)
                  ID <int> (ID of this star)
                  x_pix <int> (x-pixel location of this star)
                  y_pix <int> (y-pixel location of this star)
                  moments_file_name <string> (Name of the moments file we want to append this#
                                              residual's moment data to.)
                  moments_wings_file_name <string> (Ditto, but for the wing-weighted moments)
                  total_flux <float> (flux of the star, used for scaling the PSF)
                  mag <float> (magnitude of the star. Will be printed in the moments files)
        Optional: calc_chi2 <bool> (whether or not to calculate Chi-squared for comparison of this
                                    star with the PSF.)
                  gain <float> (gain of the image)
                  background <float> (local background level of the star)
                  
        Returns: chi2 <float> (0 if calc_chi2 is False, otherwise the Chi-squared value of the
                               star compared to the PSF.)
        
        Side-effects: Overwrites noisy_psf_file_name (see below) with a new noisy psf file
                      Overwrites residual_file_name with a new residual file
                      Appends moments data of this residual to moments_file_name and
                          moments_wings_file_name
        
        Except: residual_file_name cannot be written to
                moments_file_name or moments_wings_file_name cannot be written to
    """

    data_struct = read_fits(data_file_name)
    psf_struct = read_fits(psf_file_name)
    
    background_noise = remove_background(data_struct,background)
    
    scale = float(cut_off_scale)
    
    # Make sure we have a good scale. If the star's monopole term ends up being negative,
    # cut the scale in half and try again
    while(scale>float(cut_off_scale)/16):
        star_size, star_e1, star_e2, star_mp = mpm.get_size_and_shape(data_struct,-1,-1,
                                                                  scale,0,aperture_size)
        if(star_mp>0):
            break
        scale /= 2
        
    if(scale==float(cut_off_scale)/16):
        # If we get here, then the most likely scenario is that there's some issue with
        # the background estimation of this star. The safest is to undo removing the
        # background and try again
        
        # Restore the data structure and scale
        data_struct = read_fits(data_file_name)
        scale = float(cut_off_scale)
        
        while(scale>float(cut_off_scale)/16):
            star_size, star_e1, star_e2, star_mp = mpm.get_size_and_shape(data_struct,-1,-1,
                                                                      scale,0,aperture_size)
            if(star_mp>0):
                break
            scale /= 2
        
        if(scale==float(cut_off_scale)/16):
            # If we get here again, there's something seriously wrong with this star.
            # Throw an exception
            raise Exception("Star does not have positive monopole term; cannot scale PSF.")
        
    
    # Save the noisy psf. If this fails, just print a warning, since it's not needed
    try:
        noisy_psf_file_name = psf_file_name.replace("binned","noisy_psf")
        if(noisy_psf_file_name == psf_file_name):
            noisy_psf_file_name = os.path.splitext(psf_file_name)[0] + "_noisy.fits"
        add_and_save_noise(data_struct,psf_struct,noisy_psf_file_name,total_flux,gain)
    except Exception, _e:
        print("WARNING: Cannot output noisy PSF file to " + noisy_psf_file_name)
    
    
    psf_size, psf_e1, psf_e2, psf_mp = mpm.get_size_and_shape(psf_struct,-1,-1,
                                                              scale,0,aperture_size)
    
    scaling_factor = star_mp/psf_mp
    
    diff_struct, chi2 = remove_psf(data_struct, psf_struct, scaling_factor, calc_chi2, background_noise, gain)
    
    mm = mpm.get_2d_multipole_moments(diff_struct,-1,-1,scale,aperture_size,1./total_flux)
    mpm.append_moments(moments_file_name, ID, x_pix, y_pix, mag, mm,star_size,psf_size,
                 star_e1,star_e2,psf_e1,psf_e2,total_flux)
    
    mm = mpm.get_2d_multipole_moments(diff_struct,-1,-1,0,wing_aperture_size,1./total_flux)
    star_size, star_e1, star_e2, _mp = mpm.get_size_and_shape(data_struct,-1,-1,
                                                                    scale,0,wing_aperture_size)
    psf_size, psf_e1, psf_e2, _mp = mpm.get_size_and_shape(psf_struct,-1,-1,
                                                                 scale,0,wing_aperture_size)
    mpm.append_moments(moments_wings_file_name, ID, x_pix, y_pix, mag, mm,star_size,psf_size,
                 star_e1,star_e2,psf_e1,psf_e2,total_flux)

    write_fits(diff_struct, residual_file_name)
    
    return chi2

def calculate_chi2(data,model,background_noise=0,gain=2.0):
    """ Calculates the reduced Chi-squared comparison of a data ndarray with a model.
    
        Requires: data <ndarray> (Measured data)
                  model <ndarray> (Model to compare against)
                  
        Optional: background_noise <float> (Extra Gaussian noise which is expected to be present
                                            in data)
                  gain <float> (Ratio of counts to values, used to calculate Poisson noise)
        
        Returns: red_chi2 (reduced chi-squared value)
        
        Side-effects: (none)
        
        Except: data and model aren't ndarrays of the same shape
    """
    
    data_shape = np.shape(data)
    
    # Check data and model are the same shape
    if(data_shape != np.shape(model)):
        raise Exception("calculate_chi2 requires data and model to be ndarrays of the same shape.")
    
    num_points = np.product(data_shape)
    
    chi2 = 0
    
    for (data_value, model_value) in zip(data.ravel(),model.ravel()):
        
        if(model_value>=0):
            Poisson_noise_squared = (model_value/gain)
        else:
            Poisson_noise_squared = 0
        
        noise = np.sqrt(Poisson_noise_squared+np.square(background_noise))
        
        chi2 += np.square((data_value-model_value)/noise)
        
    
    red_chi2 = chi2/num_points
    
    return red_chi2
    