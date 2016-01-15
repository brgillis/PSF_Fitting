""" @file test_psf_for_focus.py

    Created 21 Sep 2015

    Function to test the PSF model for a single focus value.

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
from psf_testing import magic_values as mv
from psf_testing.get_model_psf import get_model_psf_for_star
from psf_testing.moments.get_Qs import get_m0_and_Qs
from psf_testing.psf_model_scheme import psf_model_scheme
from psf_testing.remove_outliers import remove_outliers
from psf_testing.parmap import parmap

def get_num_valid_stars(stars):
    
    num_valid_stars = 0
    
    for star in stars:
        if star.valid:
            num_valid_stars += 1
            
    return num_valid_stars

def get_star_arrays(stars):
    
    num_valid_stars = get_num_valid_stars(stars)
    
    properties = ("m0", "Qxy", "Qpcs", "x_pix", "y_pix")
    
    star_m0s = np.zeros((num_valid_stars,2),dtype=float)
    star_Qxys = np.zeros((num_valid_stars,2,2),dtype=float)
    star_Qpcss = np.zeros((num_valid_stars,3,2),dtype=float)
    star_x_pixs = np.zeros((num_valid_stars),dtype=float)
    star_y_pixs = np.zeros((num_valid_stars),dtype=float)
        
    i = 0
    for star in stars:
    
        if not star.valid:
            continue
    
        for prop in properties:
            exec("star_" + prop + "s[i] = star." + prop)
            
        i += 1
            
    result_exec_string = "res = ( "
    
    for prop in properties:
        result_exec_string += "star_" + prop + "s,"
        
    result_exec_string += ")"
        
    exec(result_exec_string)
    
    return res

def test_psf_for_focus(stars,

                       image_filename,
                       image=None,

                       test_focus=mv.default_test_focus,
                       num_grid_points=mv.default_num_grid_points,

                       prim_weight_func=mv.default_prim_weight_func,
                       sec_weight_func=mv.default_sec_weight_func,

                       tinytim_path=mv.default_tinytim_path,
                       tinytim_data_path=mv.default_tinytim_data_path,
                       
                       seed_factor=0,
                       gain=mv.gain,
                       fitted_params=0,
                       
                       save_models=False,
                       fitting_record=None,

                       outliers_mask=None,
                       files_to_cleanup=None,
                       
                       parallelize=False):

    if outliers_mask is None:
        outliers_mask = []

    if image is None:
        from astropy.io import fits
        image = fits.open(image_filename)[0]

    # Initialize arrays for results per star tested
    model_m0s = []
    model_m0_errs = []
    model_m0_covars = []
    model_Qxys = []
    model_Qxy_errs = []
    model_Qxy_covars = []
    model_Qpcss = []
    model_Qpcs_errs = []
    model_Qpcs_covars = []
    noisy_model_m0s = []
    noisy_model_m0_errs = []
    noisy_model_m0_covars = []
    noisy_model_Qxys = []
    noisy_model_Qxy_errs = []
    noisy_model_Qxy_covars = []
    noisy_model_Qpcss = []
    noisy_model_Qpcs_errs = []
    noisy_model_Qpcs_covars = []

    # Get the image shape, and reverse its ordering to x,y in fits ordering
    image_shape = np.shape(image)
    image_shape = image_shape[1], image_shape[0]

    # Set up the focus generating scheme
    model_scheme = psf_model_scheme(focus=test_focus,
                                    num_grid_points=num_grid_points,
                                    image_shape=image_shape)

    # Get results for each star
    num_stars = len(stars)
    assert get_num_valid_stars(stars) > 1
    
    def test_star_with_index(i):

        star = stars[i]

        if not star.valid:
            return star

        model_psf = get_model_psf_for_star(star=star,
                                           scheme=model_scheme,
                                           weight_func=prim_weight_func,
                                           tinytim_path=tinytim_path,
                                           tinytim_data_path=tinytim_data_path)

        (star.model_m0, star.model_m0_err, star.model_m0_covar,
         star.model_Qxy, star.model_Qxy_err, star.model_Qxy_covar,
         star.model_Qpcs, star.model_Qpcs_err, star.model_Qpcs_covar) = \
            get_m0_and_Qs(image=model_psf,
                          prim_weight_func=prim_weight_func,
                          sec_weight_func=sec_weight_func,
                          background_noise=star.background_noise,
                          gain=gain)

        # Now, get a noisy psf and test it (so we can test for noise bias)
        model_psf_noise = np.sqrt(model_psf / gain + np.square(star.background_noise))

        ny, nx = np.shape(model_psf)

        # Seed the random number generated with an arbitrary but stable integer
        np.random.seed(seed_factor*num_stars + i)
        noisy_model_psf = model_psf + model_psf_noise * np.random.randn(ny, nx)

        try:
            (star.noisy_model_m0, star.noisy_model_m0_err, star.noisy_model_m0_covar,
             star.noisy_model_Qxy, star.noisy_model_Qxy_err, star.noisy_model_Qxy_covar,
             star.noisy_model_Qpcs, star.noisy_model_Qpcs_err, star.noisy_model_Qpcs_covar) = \
                get_m0_and_Qs(image=noisy_model_psf,
                              prim_weight_func=prim_weight_func,
                              sec_weight_func=sec_weight_func,
                              background_noise=star.background_noise,
                              gain=gain)
        except AssertionError as _e:
            (star.noisy_model_m0, star.noisy_model_m0_err, star.noisy_model_m0_covar,
             star.noisy_model_Qxy, star.noisy_model_Qxy_err, star.noisy_model_Qxy_covar,
             star.noisy_model_Qpcs, star.noisy_model_Qpcs_err, star.noisy_model_Qpcs_covar) = \
                (star.model_m0, star.model_m0_err, star.model_m0_covar,
                 star.model_Qxy, star.model_Qxy_err, star.model_Qxy_covar,
                 star.model_Qpcs, star.model_Qpcs_err, star.model_Qpcs_covar)

        # Save the models if desired
        if save_models:
            star.model_psf = model_psf
            star.noisy_model_psf = noisy_model_psf
            
        return star
            
    if not parallelize:
        for i in range(num_stars):
            test_star_with_index(i)
    else:
        new_stars = parmap(test_star_with_index,range(num_stars))
        for i in range(num_stars):
            stars[i] = new_stars[i]
        del(new_stars)
        
    for i in range(num_stars):
            
        star = stars[i]
        
        if not star.valid:
            continue
        
        # Append m0 and Q data to the storage lists
        model_m0s.append(star.model_m0)
        model_m0_errs.append(star.model_m0_err)
        model_m0_covars.append(star.model_m0_covar)
        model_Qxys.append(star.model_Qxy)
        model_Qxy_errs.append(star.model_Qxy_err)
        model_Qxy_covars.append(star.model_Qxy_covar)
        model_Qpcss.append(star.model_Qpcs)
        model_Qpcs_errs.append(star.model_Qpcs_err)
        model_Qpcs_covars.append(star.model_Qpcs_covar)
        
        noisy_model_m0s.append(star.noisy_model_m0)
        noisy_model_m0_errs.append(star.noisy_model_m0_err)
        noisy_model_m0_covars.append(star.noisy_model_m0_covar)
        noisy_model_Qxys.append(star.noisy_model_Qxy)
        noisy_model_Qxy_errs.append(star.noisy_model_Qxy_err)
        noisy_model_Qxy_covars.append(star.noisy_model_Qxy_covar)
        noisy_model_Qpcss.append(star.noisy_model_Qpcs)
        noisy_model_Qpcs_errs.append(star.noisy_model_Qpcs_err)
        noisy_model_Qpcs_covars.append(star.noisy_model_Qpcs_covar)
        

    # Convert the psf m0 and Q storage lists to numpy arrays
    model_m0s = np.array(model_m0s)
    model_m0_errs = np.array(model_m0_errs)
    model_m0_covars = np.array(model_m0_covars)
    model_Qxys = np.array(model_Qxys)
    model_Qxy_errs = np.array(model_Qxy_errs)
    model_Qxy_covars = np.array(model_Qxy_covars)
    model_Qpcss = np.array(model_Qpcss)
    model_Qpcs_errs = np.array(model_Qpcs_errs)
    model_Qpcs_covars = np.array(model_Qpcs_covars)
    noisy_model_m0s = np.array(noisy_model_m0s)
    noisy_model_m0_errs = np.array(noisy_model_m0_errs)
    noisy_model_m0_covars = np.array(noisy_model_m0_covars)
    noisy_model_Qxys = np.array(noisy_model_Qxys)
    noisy_model_Qxy_errs = np.array(noisy_model_Qxy_errs)
    noisy_model_Qxy_covars = np.array(noisy_model_Qxy_covars)
    noisy_model_Qpcss = np.array(noisy_model_Qpcss)
    noisy_model_Qpcs_errs = np.array(noisy_model_Qpcs_errs)
    noisy_model_Qpcs_covars = np.array(noisy_model_Qpcs_covars)

    # Now get the comparison Z values for both the psf and noisy psf
    
    star_m0s, star_Qxys, star_Qpcss, star_x_pixs, star_y_pixs = get_star_arrays(stars)

    # Get the differences first
    m0_diffs = star_m0s - model_m0s
    Qxy_diffs = star_Qxys - model_Qxys
    Qpcs_diffs = star_Qpcss - model_Qpcss
    noisy_m0_diffs = star_m0s - noisy_model_m0s
    noisy_Qxy_diffs = star_Qxys - noisy_model_Qxys
    noisy_Qpcs_diffs = star_Qpcss - noisy_model_Qpcss

    # And then get the Z values. Use the noise-free psf error estimate instead of the stars'
    # error estimates, since it should be less biased
    
    m0_diff_diffs = (m0_diffs/model_m0_errs)[:,1] - (m0_diffs/model_m0_errs)[:,0]
    m0_diff_diff_errs = 2. - 2.*model_m0_covars[:,1,0]/np.product(model_m0_errs,axis=1)
    m0_diff_Zs = m0_diff_diffs/m0_diff_diff_errs
    
    Qxy_diff_diffs = (Qxy_diffs/model_Qxy_errs)[:,:,1] - (Qxy_diffs/model_Qxy_errs)[:,:,0]
    Qxy_diff_diff_errs = 2. - 2.*model_Qxy_covars[:,:,1,0]/np.product(model_Qxy_errs,axis=2)
    Qxy_diff_Zs = Qxy_diff_diffs/Qxy_diff_diff_errs
    
    Qpcs_diff_sums = (Qpcs_diffs/model_Qpcs_errs).sum(axis=2)
    Qpcs_diff_sum_errs = 2. + 2.*model_Qpcs_covars[:,:,1,0]/np.product(model_Qpcs_errs,axis=2)
    Qpcs_sum_Zs = Qpcs_diff_sums/Qpcs_diff_sum_errs
    
    Qpcs_diff_diffs = (Qpcs_diffs/model_Qpcs_errs)[:,:,1] - (Qpcs_diffs/model_Qpcs_errs)[:,:,0]
    Qpcs_diff_diff_errs = 2. - 2.*model_Qpcs_covars[:,:,1,0]/np.product(model_Qpcs_errs,axis=2)
    Qpcs_diff_Zs = Qpcs_diff_diffs/Qpcs_diff_diff_errs
    
    
    noisy_m0_diff_diffs = (noisy_m0_diffs/model_m0_errs/2.)[:,1] - (noisy_m0_diffs/model_m0_errs)[:,0]
    noisy_m0_diff_diff_errs = 2. - 2.*model_m0_covars[:,1,0]/np.product(model_m0_errs,axis=1)
    noisy_m0_diff_Zs = noisy_m0_diff_diffs/noisy_m0_diff_diff_errs
    
    noisy_Qxy_diff_diffs = (noisy_Qxy_diffs/model_Qxy_errs/2.)[:,:,1] - (noisy_Qxy_diffs/model_Qxy_errs)[:,:,0]
    noisy_Qxy_diff_diff_errs = 2. - 2.*model_Qxy_covars[:,:,1,0]/np.product(model_Qxy_errs,axis=2)
    noisy_Qxy_diff_Zs = noisy_Qxy_diff_diffs/noisy_Qxy_diff_diff_errs
    
    noisy_Qpcs_diff_sums = (noisy_Qpcs_diffs/model_Qpcs_errs/2.).sum(axis=2)
    noisy_Qpcs_diff_sum_errs = 2. + 2.*model_Qpcs_covars[:,:,1,0]/np.product(model_Qpcs_errs,axis=2)
    noisy_Qpcs_sum_Zs = noisy_Qpcs_diff_sums/noisy_Qpcs_diff_sum_errs
    
    noisy_Qpcs_diff_diffs = (noisy_Qpcs_diffs/model_Qpcs_errs/2.)[:,:,1] - (Qpcs_diffs/model_Qpcs_errs)[:,:,0]
    noisy_Qpcs_diff_diff_errs = 2. - 2.*model_Qpcs_covars[:,:,1,0]/np.product(model_Qpcs_errs,axis=2)
    noisy_Qpcs_diff_Zs = noisy_Qpcs_diff_diffs/noisy_Qpcs_diff_diff_errs

    if len(outliers_mask) == 0:

        # Look for outliers on each of the Q arrays
        Qx_mask = remove_outliers(Qxy_diff_diffs[:, 0]).mask
        Qy_mask = remove_outliers(Qxy_diff_diffs[:, 1]).mask
        Qplus_sum_mask = remove_outliers(Qpcs_diff_sums[:, 0]).mask
        Qcross_sum_mask = remove_outliers(Qpcs_diff_sums[:, 1]).mask
        Qsize_sum_mask = remove_outliers(Qpcs_diff_sums[:, 2]).mask
        Qplus_diff_mask = remove_outliers(Qpcs_diff_diffs[:, 0]).mask
        Qcross_diff_mask = remove_outliers(Qpcs_diff_diffs[:, 1]).mask
        Qsize_diff_mask = remove_outliers(Qpcs_diff_diffs[:, 2]).mask

        # If it's an outlier in any individual value, ignore it for all
        outliers_mask.append(np.logical_or.reduce((Qx_mask, Qy_mask,
                                                   Qplus_sum_mask, Qcross_sum_mask, Qsize_sum_mask,
                                                   Qplus_diff_mask, Qcross_diff_mask, Qsize_diff_mask)))

    omask = outliers_mask[0]
    
    assert not np.all(omask)

    # For the Q array, we need to duplicate the mask
    m0_omask = omask
    Qxy_omask = np.vstack((omask, omask)).transpose()
    Qpcs_omask = np.vstack((omask, omask, omask)).transpose()

    masked_m0_diff_diffs = np.ma.masked_array(m0_diff_diffs, mask=m0_omask)
    masked_Qxy_diff_diffs = np.ma.masked_array(Qxy_diff_diffs, mask=Qxy_omask)
    masked_Qpcs_diff_diffs = np.ma.masked_array(Qpcs_diff_diffs, mask=Qpcs_omask)
    masked_Qpcs_diff_sums = np.ma.masked_array(Qpcs_diff_sums, mask=Qpcs_omask)
    
    masked_m0_diff_Zs = np.ma.masked_array(m0_diff_Zs, mask=omask)
    masked_Qxy_diff_Zs = np.ma.masked_array(Qxy_diff_Zs, mask=Qxy_omask)
    masked_Qpcs_diff_Zs = np.ma.masked_array(Qpcs_diff_Zs, mask=Qpcs_omask)
    masked_Qpcs_sum_Zs = np.ma.masked_array(Qpcs_sum_Zs, mask=Qpcs_omask)
    

    masked_noisy_m0_diff_diffs = np.ma.masked_array(noisy_m0_diff_diffs, mask=m0_omask)
    masked_noisy_Qxy_diff_diffs = np.ma.masked_array(noisy_Qxy_diff_diffs, mask=Qxy_omask)
    masked_noisy_Qpcs_diff_diffs = np.ma.masked_array(noisy_Qpcs_diff_diffs, mask=Qpcs_omask)
    masked_noisy_Qpcs_diff_sums = np.ma.masked_array(noisy_Qpcs_diff_sums, mask=Qpcs_omask)
    
    masked_noisy_m0_diff_Zs = np.ma.masked_array(noisy_m0_diff_Zs, mask=omask)
    masked_noisy_Qxy_diff_Zs = np.ma.masked_array(noisy_Qxy_diff_Zs, mask=Qxy_omask)
    masked_noisy_Qpcs_diff_Zs = np.ma.masked_array(noisy_Qpcs_diff_Zs, mask=Qpcs_omask)
    masked_noisy_Qpcs_sum_Zs = np.ma.masked_array(noisy_Qpcs_sum_Zs, mask=Qpcs_omask)

    num_masked = np.sum(omask)
    num_unmasked = len(m0_diff_Zs) - num_masked

    # Get the mean Zs of the non-outliers
    m0_diff_diffs = masked_m0_diff_diffs[~m0_omask]
    mean_m0_diff_diff = np.mean(masked_m0_diff_diffs[~m0_omask], axis=0)
    
    mean_Qxy_diff_diffs = np.mean(masked_Qxy_diff_diffs[~Qxy_omask].reshape((num_unmasked, 2)), axis=0)
    mean_Qpcs_diff_sums = np.mean(masked_Qpcs_diff_sums[~Qpcs_omask].reshape((num_unmasked, 3)), axis=0)
    mean_Qpcs_diff_diffs = np.mean(masked_Qpcs_diff_diffs[~Qpcs_omask].reshape((num_unmasked, 3)), axis=0)
    
    noisy_m0_diff_diffs = masked_noisy_m0_diff_diffs[~m0_omask]
    mean_noisy_m0_diff_diff = np.mean(noisy_m0_diff_diffs, axis=0)
    
    mean_noisy_Qxy_diff_diffs = np.mean(masked_noisy_Qxy_diff_diffs[~Qxy_omask].reshape((num_unmasked, 2)), axis=0)
    mean_noisy_Qpcs_diff_sums = np.mean(masked_noisy_Qpcs_diff_sums[~Qpcs_omask].reshape((num_unmasked, 3)), axis=0)
    mean_noisy_Qpcs_diff_diffs = np.mean(masked_noisy_Qpcs_diff_diffs[~Qpcs_omask].reshape((num_unmasked, 3)), axis=0)
    
    
    m0_diff_diff_stddev = np.std(masked_m0_diff_diffs[~m0_omask], axis=0)
    Qxy_diff_diff_stddevs = np.std(masked_Qxy_diff_diffs[~Qxy_omask].reshape((num_unmasked, 2)), axis=0)
    Qpcs_diff_sum_stddevs = np.std(masked_Qpcs_diff_sums[~Qpcs_omask].reshape((num_unmasked, 3)), axis=0)
    Qpcs_diff_diff_stddevs = np.std(masked_Qpcs_diff_diffs[~Qpcs_omask].reshape((num_unmasked, 3)), axis=0)
    
    noisy_m0_diff_diff_stddev = np.std(masked_noisy_m0_diff_diffs[~m0_omask], axis=0)
    noisy_Qxy_diff_diff_stddevs = np.std(masked_noisy_Qxy_diff_diffs[~Qxy_omask].reshape((num_unmasked, 2)), axis=0)
    noisy_Qpcs_diff_sum_stddevs = np.std(masked_noisy_Qpcs_diff_sums[~Qpcs_omask].reshape((num_unmasked, 3)), axis=0)
    noisy_Qpcs_diff_diff_stddevs = np.std(masked_noisy_Qpcs_diff_diffs[~Qpcs_omask].reshape((num_unmasked, 3)), axis=0)

    m0_diff_emp_Z2 = np.sum(np.square(masked_m0_diff_diffs[~m0_omask]/m0_diff_diff_stddev), axis=0)
    Qxy_diff_emp_Z2s = np.sum(np.square(masked_Qxy_diff_diffs[~Qxy_omask].reshape((num_unmasked, 2))/Qxy_diff_diff_stddevs), axis=0)
    Qpcs_sum_emp_Z2s = np.sum(np.square(masked_Qpcs_diff_sums[~Qpcs_omask].reshape((num_unmasked, 3))/Qpcs_diff_sum_stddevs), axis=0)
    Qpcs_diff_emp_Z2s = np.sum(np.square(masked_Qpcs_diff_diffs[~Qpcs_omask].reshape((num_unmasked, 3))/Qpcs_diff_diff_stddevs), axis=0)

    noisy_m0_diff_emp_Z2 = np.sum(np.square(masked_noisy_m0_diff_diffs[~m0_omask]/noisy_m0_diff_diff_stddev)*(num_unmasked-1), axis=0)
    noisy_Qxy_diff_emp_Z2s = np.sum(np.square(masked_noisy_Qxy_diff_diffs[~Qxy_omask].reshape((num_unmasked, 2))/noisy_Qxy_diff_diff_stddevs), axis=0)
    noisy_Qpcs_sum_emp_Z2s = np.sum(np.square(masked_noisy_Qpcs_diff_sums[~Qpcs_omask].reshape((num_unmasked, 3))/noisy_Qpcs_diff_sum_stddevs), axis=0)
    noisy_Qpcs_diff_emp_Z2s = np.sum(np.square(masked_noisy_Qpcs_diff_diffs[~Qpcs_omask].reshape((num_unmasked, 3))/noisy_Qpcs_diff_diff_stddevs), axis=0)


    m0_diff_Z2 = np.sum(np.square(masked_m0_diff_Zs[~omask]), axis=0)
    noisy_m0_diff_Z2 = np.sum(np.square(masked_noisy_m0_diff_Zs[~omask]), axis=0)
    
    Qx_diff_Z2 = np.sum(np.square(masked_Qxy_diff_Zs[~Qxy_omask].reshape((num_unmasked, 2))[:, 0]), axis=0)
    Qy_diff_Z2 = np.sum(np.square(masked_Qxy_diff_Zs[~Qxy_omask].reshape((num_unmasked, 2))[:, 1]), axis=0)
    Qplus_sum_Z2 = np.sum(np.square(masked_Qpcs_sum_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 0]), axis=0)
    Qcross_sum_Z2 = np.sum(np.square(masked_Qpcs_sum_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 1]), axis=0)
    Qsize_sum_Z2 = np.sum(np.square(masked_Qpcs_sum_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 2]), axis=0)
    Qplus_diff_Z2 = np.sum(np.square(masked_Qpcs_diff_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 0]), axis=0)
    Qcross_diff_Z2 = np.sum(np.square(masked_Qpcs_diff_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 1]), axis=0)
    Qsize_diff_Z2 = np.sum(np.square(masked_Qpcs_diff_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 2]), axis=0)
    
    noisy_Qx_diff_Z2 = np.sum(np.square(masked_noisy_Qxy_diff_Zs[~Qxy_omask].reshape((num_unmasked, 2))[:, 0]), axis=0)
    noisy_Qy_diff_Z2 = np.sum(np.square(masked_noisy_Qxy_diff_Zs[~Qxy_omask].reshape((num_unmasked, 2))[:, 1]), axis=0)
    noisy_Qplus_sum_Z2 = np.sum(np.square(masked_noisy_Qpcs_sum_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 0]), axis=0)
    noisy_Qcross_sum_Z2 = np.sum(np.square(masked_noisy_Qpcs_sum_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 1]), axis=0)
    noisy_Qsize_sum_Z2 = np.sum(np.square(masked_noisy_Qpcs_sum_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 2]), axis=0)
    noisy_Qplus_diff_Z2 = np.sum(np.square(masked_noisy_Qpcs_diff_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 0]), axis=0)
    noisy_Qcross_diff_Z2 = np.sum(np.square(masked_noisy_Qpcs_diff_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 1]), axis=0)
    noisy_Qsize_diff_Z2 = np.sum(np.square(masked_noisy_Qpcs_diff_Zs[~Qpcs_omask].reshape((num_unmasked, 3))[:, 2]), axis=0)
  
    emp_chi2 = np.sum(Qxy_diff_emp_Z2s) + np.sum(Qpcs_diff_emp_Z2s[0:2]) + \
        np.sum(Qpcs_sum_emp_Z2s[0:2]) + noisy_Qpcs_sum_emp_Z2s[2] + noisy_Qpcs_diff_emp_Z2s[2]
    chi2 = Qx_diff_Z2 + Qy_diff_Z2 + Qplus_diff_Z2 + Qcross_diff_Z2 + Qplus_sum_Z2 + Qcross_sum_Z2 + \
        noisy_Qsize_diff_Z2 + noisy_Qsize_sum_Z2

    dof = 8 * num_unmasked - fitted_params
    emp_dof = 8 * num_unmasked - fitted_params
    
    test_results = (test_focus,
            ((chi2, emp_chi2), (dof, emp_dof)),
            ((mean_m0_diff_diff, (mean_Qxy_diff_diffs[0], mean_Qxy_diff_diffs[1],
                                  mean_Qpcs_diff_diffs[0], mean_Qpcs_diff_diffs[1], mean_Qpcs_diff_diffs[2],
                                  mean_Qpcs_diff_sums[0], mean_Qpcs_diff_sums[1], mean_Qpcs_diff_sums[2])),
             (mean_noisy_m0_diff_diff, (mean_noisy_Qxy_diff_diffs[0], mean_noisy_Qxy_diff_diffs[1],
                                  mean_noisy_Qpcs_diff_diffs[0], mean_noisy_Qpcs_diff_diffs[1], mean_noisy_Qpcs_diff_diffs[2],
                                  mean_noisy_Qpcs_diff_sums[0], mean_noisy_Qpcs_diff_sums[1], mean_noisy_Qpcs_diff_sums[2]))),
            ((m0_diff_Z2,
              (Qx_diff_Z2, Qy_diff_Z2, Qplus_sum_Z2, Qcross_sum_Z2, Qsize_sum_Z2,
               Qplus_diff_Z2, Qcross_diff_Z2, Qsize_diff_Z2)),
             (noisy_m0_diff_Z2,
              (noisy_Qx_diff_Z2, noisy_Qy_diff_Z2, noisy_Qplus_sum_Z2, noisy_Qcross_sum_Z2, noisy_Qsize_sum_Z2,
               noisy_Qplus_diff_Z2, noisy_Qcross_diff_Z2, noisy_Qsize_diff_Z2))),
            ((m0_diff_emp_Z2,
              (Qxy_diff_emp_Z2s[0], Qxy_diff_emp_Z2s[1],
               Qpcs_sum_emp_Z2s[0], Qpcs_sum_emp_Z2s[1], Qpcs_sum_emp_Z2s[2],
               Qpcs_diff_emp_Z2s[0], Qpcs_diff_emp_Z2s[1], Qpcs_diff_emp_Z2s[2])),
             (noisy_m0_diff_emp_Z2,
              (noisy_Qxy_diff_emp_Z2s[0], noisy_Qxy_diff_emp_Z2s[1],
               noisy_Qpcs_sum_emp_Z2s[0], noisy_Qpcs_sum_emp_Z2s[1], noisy_Qpcs_sum_emp_Z2s[2],
               noisy_Qpcs_diff_emp_Z2s[0], noisy_Qpcs_diff_emp_Z2s[1], noisy_Qpcs_diff_emp_Z2s[2]))),
            (star_m0s, star_Qxys, star_Qpcss),
            (model_m0s, model_Qxys, model_Qpcss),
            (noisy_model_m0s, noisy_model_Qxys, noisy_model_Qpcss),
            (model_m0_errs, model_Qxy_errs, model_Qpcs_errs),
            (omask, Qxy_omask, Qpcs_omask),
            (m0_diff_Zs, Qxy_diff_Zs, Qpcs_sum_Zs, Qpcs_diff_Zs),
            (star_x_pixs, star_y_pixs),
            )
    
    if fitting_record is not None:
        fitting_record.append(test_results)

    return test_results
