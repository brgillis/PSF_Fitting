#!/usr/bin/env python2
#
# create psf postage stamps at a given position in the image frame
#
# OM 2014

import numpy as np
from astropy.io import fits as pyf
import copy
from scipy.stats import norm

import multipole_moments as mpm
from fits_functions import read_fits, write_fits, add_comment

# some magic numbers
aperture_size = 10
cut_off_scale = 2
wing_aperture_size = 20

def center_intensity(data):
    """UNUSED - Calculates intensity at centre of image, taking the mean of the nine pixels
       nearest the centre.
    """

    size = np.shape(data)
    
    intensity = np.mean(data[size[0]/2-1:size[0]/2+2,size[1]/2-1:size[1]/2+2])

    return intensity

def max_intensity(data):
    """UNUSED - Get the maximum intensity on the image. 
    """

    maxpix = np.unravel_index(data.argmax(), data.shape)
    return data[maxpix[0],maxpix[1]]

def center_data(data):
    """UNUSED - Adjust the data to centre it on the barycentre.
    """
    
    centered = False
    
    # Iterate here - due to the constraints of the frame, we're unlikely to move more than
    # a pixel at a time, even if it's badly mis-centred
    while(not centered):
        cenpix = mean_pixel(data)
        size = copy.copy(np.shape(data))
        halfsize = min((cenpix[0], cenpix[1],
                        size[0] - cenpix[0] - 1, size[1] - cenpix[1] - 1))
    
        data = data[int(cenpix[0] - halfsize):int(cenpix[0] + halfsize+1),
                             int(cenpix[1] - halfsize):int(cenpix[1] + halfsize+1)]
        
        if(np.shape(data)[0] == size[0] and np.shape(data)[1] == size[1]):
            centered = True

    return data

def mean_pixel(data):
    """UNUSED - Determine the (unweighted) barycentre of the data.
    """
    
    # Construct two weight arrays, for x and y
    size = np.shape(data)
    x_w = np.fromfunction(lambda i,j : i, size)
    y_w = np.fromfunction(lambda i,j : j, size)
    
    try:
        xm = np.sum(np.multiply(data,x_w))/np.sum(data)
        xm = int(round(xm))
    except Exception:
        xm = int(round(size[0]/2))
    try:
        ym = np.sum(np.multiply(data,y_w))/np.sum(data)
        ym = int(round(ym))
    except Exception:
        ym = int(round(size[1]/2))
        
    return xm,ym

def remove_psf(fits_data, psf_data, total_flux):
    """ Subtracts the PSF from the data, scaling it by total_flux.
    
        Requires: fits_data <HDUList> (star's fits data structure)
                  psf_data <HDUList> (psf's fits data structure)
                  total_flux (flux of the star, used here to scale the psf)
                  
        Returns: out_fits <HDUList> (residual's fits data structure)
        
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

    out = data1 - data2 * total_flux

    out_fits = pyf.HDUList(pyf.PrimaryHDU(out))
    out_fits[0].header = fits_data[0].header
    add_comment(out_fits,"Postage stamp removed from central structure")

    return out_fits

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
    
    outliers_found = True
    num_tot = len(olist)
    
    while(outliers_found):

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
        raise Exception("WARNING: Too few arguments left after removing outliers. (Only " + str(len(list)) + ")"
                        "\nReverting to non-trimmed dataset.\n" + 
                        "Chi-squared value will include presence of outliers.")
    else:
        print("Removed " + str(num_tot-len(list)) + "/" + str(num_tot) +  " outliers.")
        return mean, sigma
    

def remove_background( stamp_struct ):
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
        edge_mean, unused_edge_std = remove_outliers(edge_data)
    except:
        edge_data = edge_data_copy
        edge_mean, unused_edge_std = np.mean(edge_data), np.std(edge_data)
        
    stamp_struct[0].data = np.subtract(stamp_data,edge_mean)
    

def add_noise( stamp_struct, psf_struct, noisy_psf_file_name, total_flux=1. ):
    """ Adds noise to a psf image, based on information from a star image.
    
        Requires: stamp_struct <HDUList> (stamp image's fits data structure)
                  psf_struct <HDUList> (psf image's fits data structure)
                  noise_psf_file_name <string> (file name in which to save the psf+noise)
        Optional: total_flux <float>  (Flux of the star, used to scale noise)
        
        Returns: (nothing)
        
        Side-effects: psf_struct has noise added
                      Overwrites noisy_psf_file_name with new image (if not None)
        
        Except: noisy_psf_file_name cannot be written to
        
    """
    
    stamp_data = stamp_struct[0].data
    psf_data = psf_struct[0].data
    
    stamp_size = np.shape(stamp_data)
    psf_size = np.shape(psf_data)
    
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
        
    # Seed numpy from the first value in the array. This keeps it deterministic for fitting, but
    # random on a per-stamp case
    np.random.seed(int(abs(edge_data[0]*10000)))
    
    noise_data = np.zeros(psf_size)
    
    for i in xrange(psf_size[0]):
        for j in xrange(psf_size[1]):
            noise_data[i,j] = (edge_std*np.random.randn() + edge_mean)/total_flux
    
    new_psf_data = np.add(psf_data,noise_data)
    
    psf_struct[0].data = new_psf_data
    
    if(noisy_psf_file_name is not None):
        write_fits(psf_struct,noisy_psf_file_name)


def subtract_psf(data_file_name, psf_file_name, residual_file_name, ID, x_pix, y_pix,
                 moments_file_name, moments_wings_file_name,total_flux,mag):
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
                  
        Returns: (nothing)
        
        Side-effects: Overwrites noisy_psf_file_name (see below) with a new noisy psf file
                      Overwrites residual_file_name with a new residual file
                      Appends moments data of this residual to moments_file_name and
                          moments_wings_file_name
        
        Except: residual_file_name cannot be written to
                moments_file_name or moments_wings_file_name cannot be written to
    """

    data_struct = read_fits(data_file_name)
    psf_struct = read_fits(psf_file_name)
    
    remove_background(data_struct)
    
    star_size, star_e1, star_e2, star_mp = mpm.get_size_and_shape(data_struct,-1,-1,
                                                                  cut_off_scale,0,aperture_size)
    psf_size, psf_e1, psf_e2, psf_mp = mpm.get_size_and_shape(psf_struct,-1,-1,
                                                              cut_off_scale,0,aperture_size)
    
    scaling_factor = star_mp/psf_mp
    
    diff_struct = remove_psf(data_struct, psf_struct, scaling_factor)
    
#     # Save the noisy psf. If this fails, just print a warning, since it's not needed
#     try:
#         noisy_psf_file_name = psf_file_name.replace("binned","noisy_psf")
#         if(noisy_psf_file_name == psf_file_name):
#             noisy_psf_file_name = os.path.splitext(psf_file_name)[0] + "_noisy.fits"
#         add_noise(data_struct,psf_struct,noisy_psf_file_name,total_flux)
#     except Exception:
#         print("WARNING: Cannot output noisy PSF file to " + noisy_psf_file_name)
    
    mm = mpm.get_2d_multipole_moments(diff_struct,-1,-1,cut_off_scale,aperture_size)
    mpm.append_moments(moments_file_name, ID, x_pix, y_pix, mag, mm,star_size,psf_size,
                 star_e1,star_e2,psf_e1,psf_e2,total_flux)
    
    mm = mpm.get_2d_multipole_moments(diff_struct,-1,-1,0,wing_aperture_size)
    star_size, star_e1, star_e2, unused_mp = mpm.get_size_and_shape(data_struct,-1,-1,
                                                                    cut_off_scale,0,wing_aperture_size)
    psf_size, psf_e1, psf_e2, unused_mp = mpm.get_size_and_shape(psf_struct,-1,-1,
                                                                 cut_off_scale,0,wing_aperture_size)
    mpm.append_moments(moments_wings_file_name, ID, x_pix, y_pix, mag, mm,star_size,psf_size,
                 star_e1,star_e2,psf_e1,psf_e2,total_flux)

    write_fits(diff_struct, residual_file_name)
