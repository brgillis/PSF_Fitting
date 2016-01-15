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

def get_star_arrays(stars,outliers_mask=None):

    if outliers_mask is None:
        outliers_mask = []
    
    num_valid_stars = get_num_valid_stars(stars)
    
    properties = ("m0", "Qxy", "Qpcs", "x_pix", "y_pix")
    
    star_props = {}
    
    star_props["m0"] = np.zeros((num_valid_stars,2),dtype=float)
    star_props["Qxy"] = np.zeros((num_valid_stars,2,2),dtype=float)
    star_props["Qpcs"] = np.zeros((num_valid_stars,3,2),dtype=float)
    star_props["x_pix"] = np.zeros((num_valid_stars),dtype=float)
    star_props["y_pix"] = np.zeros((num_valid_stars),dtype=float)
        
    i = 0
    for star in stars:
    
        if not star.valid:
            continue
    
        for prop in properties:
            exec("star_props[prop][i] = star." + prop)
            
        i += 1
            
    return star_props

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
    model_Qxys = []
    model_Qpcss = []
    noisy_model_m0s = []
    noisy_model_Qxys = []
    noisy_model_Qpcss = []

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

        (star.model_m0, star.model_Qxy, star.model_Qpcs) = \
            get_m0_and_Qs(image=model_psf,
                          prim_weight_func=prim_weight_func,
                          sec_weight_func=sec_weight_func)

        # Now, get a noisy psf and test it (so we can test for noise bias)
        model_psf_noise = np.sqrt(model_psf / gain + np.square(star.background_noise))

        ny, nx = np.shape(model_psf)

        # Seed the random number generated with an arbitrary but stable integer
        np.random.seed(seed_factor*num_stars + i)
        noisy_model_psf = model_psf + model_psf_noise * np.random.randn(ny, nx)

        try:
            (star.noisy_model_m0, star.noisy_model_Qxy, star.noisy_model_Qpcs ) = \
                get_m0_and_Qs(image=noisy_model_psf,
                              prim_weight_func=prim_weight_func,
                              sec_weight_func=sec_weight_func)
        except AssertionError as _e:
            (star.noisy_model_m0, star.noisy_model_Qxy, star.noisy_model_Qpcs, ) = \
                (star.model_m0, star.model_Qxy, star.model_Qpcs, )

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
        model_Qxys.append(star.model_Qxy)
        model_Qpcss.append(star.model_Qpcs)
        
        noisy_model_m0s.append(star.noisy_model_m0)
        noisy_model_Qxys.append(star.noisy_model_Qxy)
        noisy_model_Qpcss.append(star.noisy_model_Qpcs)
        

    # Convert the psf m0 and Q storage lists to numpy arrays
    model_m0s = np.array(model_m0s)
    model_Qxys = np.array(model_Qxys)
    model_Qpcss = np.array(model_Qpcss)
    noisy_model_m0s = np.array(noisy_model_m0s)
    noisy_model_Qxys = np.array(noisy_model_Qxys)
    noisy_model_Qpcss = np.array(noisy_model_Qpcss)

    # Now get the comparison Z values for both the psf and noisy psf
    
    star_props = get_star_arrays(stars)

    # Get the differences of the m0s and Qs between stars and models
    star_props["m0_diff"] = star_props["m0"] - model_m0s
    star_props["Qxy_diff"] = star_props["Qxy"] - model_Qxys
    star_props["Qpcs_diff"] = star_props["Qpcs"] - model_Qpcss
    star_props["noisy_m0_diff"] = star_props["m0"] - noisy_model_m0s
    star_props["noisy_Qxy_diff"] = star_props["Qxy"] - noisy_model_Qxys
    star_props["noisy_Qpcs_diff"] = star_props["Qpcs"] - noisy_model_Qpcss
        
    # Look for outliers in the Qxy and Qpcs diff values
    if len(outliers_mask) == 0:

        # Look for outliers on each of the Q arrays
        Qx1_mask = remove_outliers(star_props["Qxy_diff"][:, 0, 1]).mask
        Qy1_mask = remove_outliers(star_props["Qxy_diff"][:, 1, 1]).mask
        Qplus0_mask = remove_outliers(star_props["Qpcs_diff"][:, 0, 0]).mask
        Qplus1_mask = remove_outliers(star_props["Qpcs_diff"][:, 0, 1]).mask
        Qcross0_mask = remove_outliers(star_props["Qpcs_diff"][:, 1, 0]).mask
        Qcross1_mask = remove_outliers(star_props["Qpcs_diff"][:, 1, 1]).mask
        Qsize0_mask = remove_outliers(star_props["Qpcs_diff"][:, 2, 0]).mask
        Qsize1_mask = remove_outliers(star_props["Qpcs_diff"][:, 2, 1]).mask

        # If it's an outlier in any individual value, ignore it for all
        outliers_mask.append(np.logical_or.reduce((Qx1_mask, Qy1_mask,
                                                   Qplus0_mask, Qplus1_mask,
                                                   Qcross0_mask, Qcross1_mask,
                                                   Qsize0_mask, Qsize1_mask, )))
    omask = outliers_mask[0]
    
    assert not np.all(omask)

    # Set up multi-dimensional versions of the mask
    m0_omask = np.vstack((omask, omask)).transpose()
    Qxy_omask = np.vstack(np.vstack((omask, omask)).transpose(),
                          np.vstack((omask, omask)).transpose())
    Qpcs_omask = np.vstack(np.vstack((omask, omask, omask)).transpose(),
                           np.vstack((omask, omask, omask)).transpose())
    
    star_props["masked_m0_diff"] = np.ma.masked_array(star_props["m0_diff"], mask=m0_omask)
    star_props["masked_Qxy_diff"] = np.ma.masked_array(star_props["Qxy_diff"], mask=Qxy_omask)
    star_props["masked_Qpcs_diff"] = np.ma.masked_array(star_props["Qpcs_diff"], mask=Qpcs_omask)
    star_props["masked_noisy_m0_diff"] = np.ma.masked_array(star_props["m0_diff"], mask=m0_omask)
    star_props["masked_noisy_Qxy_diff"] = np.ma.masked_array(star_props["Qxy_diff"], mask=Qxy_omask)
    star_props["masked_noisy_Qpcs_diff"] = np.ma.masked_array(star_props["Qpcs_diff"], mask=Qpcs_omask)
    
    num_good_stars = len(star_props["masked_m0_diff"][~m0_omask])//2
    
    corr_factor = np.sqrt(num_good_stars/(num_good_stars-1))
    
    star_props["unmasked_m0_diff"] = np.reshape(star_props["masked_m0_diff"][~m0_omask],
                                                (num_good_stars,2))
    star_props["unmasked_Qxy_diff"] = np.reshape(star_props["masked_Qxy_diff"][~Qxy_omask],
                                                (num_good_stars,2,2))
    star_props["unmasked_Qpcs_diff"] = np.reshape(star_props["masked_Qpcs_diff"][~Qpcs_omask],
                                                (num_good_stars,3,2))
    
    star_props["unmasked_noisy_m0_diff"] = np.reshape(star_props["masked_noisy_m0_diff"][~m0_omask],
                                                (num_good_stars,2))
    star_props["unmasked_noisy_Qxy_diff"] = np.reshape(star_props["masked_noisy_Qxy_diff"][~Qxy_omask],
                                                (num_good_stars,2,2))
    star_props["unmasked_noisy_Qpcs_diff"] = np.reshape(star_props["masked_noisy_Qpcs_diff"][~Qpcs_omask],
                                                (num_good_stars,3,2))
    
    # Get the error for each property
    for prop in ("m0_diff","Qxy_diff","Qpcs_diff",
                 "noisy_m0_diff","noisy_Qxy_diff","noisy_Qpcs_diff"):
        star_props[prop + "_err"] = np.std(star_props["unmasked_" + prop],axis=0)*corr_factor

    # And now get the diffs, sums, and Z values
    
    for prop in ("m0_diff","noisy_m0_diff"):
        temp_array = star_props["unmasked_" + prop]/star_props[prop + "_err"]
        star_props[prop + "_sums"] = temp_array[:,1] + temp_array[:,0]
        star_props[prop + "_diffs"] = temp_array[:,1] - temp_array[:,0]
    
    for prop in ("Qxy_diff","noisy_Qxy_diff",
                 "Qpcs_diff","noisy_Qpcs_diff"):
        temp_array = star_props["unmasked_" + prop]/star_props[prop + "_err"]
        star_props[prop + "_sums"] = temp_array[:,:,1] + temp_array[:,:,0]
        star_props[prop + "_diffs"] = temp_array[:,:,1] - temp_array[:,:,0]
    
    for prop in ("m0_diff","noisy_m0_diff",
                 "Qxy_diff","noisy_Qxy_diff",
                 "Qpcs_diff","noisy_Qpcs_diff"):
        for comb in ("_sum","_diff"):
            star_props[prop + comb + "_mean"] = np.mean(star_props[prop + comb + "s"],axis=0)
            star_props[prop + comb + "_err"] = np.std(star_props[prop + comb + "s"],axis=0)*corr_factor
            star_props[prop + comb + "_Zs"] = (star_props[prop + comb + "s"]/
                                               star_props[prop + comb + "s_err"])
            star_props[prop + comb + "_Z2s"] = np.sum(np.square(star_props[prop + comb + "_Zs"],axis=0))
  
    chi2 = np.sum(star_props["Qxy_diff_diff_Z2s"]) + \
        np.sum(star_props["Qpcs_diff_diff_Z2s"][0:2]) + \
        np.sum(star_props["Qpcs_diff_sum_Z2s"][0:2]) + \
        np.sum(star_props["noisy_Qpcs_diff_diff_Z2s"][2]) + \
        np.sum(star_props["noisy_Qpcs_diff_sum_Z2s"][2])

    dof = 8 * num_good_stars - fitted_params
    
    test_results = (test_focus,
            (chi2, dof),
            ((star_props["m0_diff_diff_mean"], (star_props["Qxy_diff_diff_mean"][0],
                                                star_props["Qxy_diff_diff_mean"][1],
                                                star_props["Qpcs_diff_diff_mean"][0],
                                                star_props["Qpcs_diff_diff_mean"][1],
                                                star_props["Qpcs_diff_diff_mean"][2],
                                                star_props["Qpcs_diff_sum_mean"][0],
                                                star_props["Qpcs_diff_sum_mean"][1],
                                                star_props["Qpcs_diff_sum_mean"][2])),
             (star_props["noisy_m0_diff_diff_mean"], (star_props["noisy_Qxy_diff_diff_mean"][0],
                                                    star_props["noisy_Qxy_diff_diff_mean"][1],
                                                    star_props["noisy_Qpcs_diff_diff_mean"][0],
                                                    star_props["noisy_Qpcs_diff_diff_mean"][1],
                                                    star_props["noisy_Qpcs_diff_diff_mean"][2],
                                                    star_props["noisy_Qpcs_diff_sum_mean"][0],
                                                    star_props["noisy_Qpcs_diff_sum_mean"][1],
                                                    star_props["noisy_Qpcs_diff_sum_mean"][2]))),
            ((star_props["m0_diff_diff_Z2"], (star_props["Qxy_diff_diff_Z2"][0],
                                                star_props["Qxy_diff_diff_Z2"][1],
                                                star_props["Qpcs_diff_diff_Z2"][0],
                                                star_props["Qpcs_diff_diff_Z2"][1],
                                                star_props["Qpcs_diff_diff_Z2"][2],
                                                star_props["Qpcs_diff_sum_Z2"][0],
                                                star_props["Qpcs_diff_sum_Z2"][1],
                                                star_props["Qpcs_diff_sum_Z2"][2])),
             (star_props["noisy_m0_diff_diff_Z2"], (star_props["noisy_Qxy_diff_diff_Z2"][0],
                                                    star_props["noisy_Qxy_diff_diff_Z2"][1],
                                                    star_props["noisy_Qpcs_diff_diff_Z2"][0],
                                                    star_props["noisy_Qpcs_diff_diff_Z2"][1],
                                                    star_props["noisy_Qpcs_diff_diff_Z2"][2],
                                                    star_props["noisy_Qpcs_diff_sum_Z2"][0],
                                                    star_props["noisy_Qpcs_diff_sum_Z2"][1],
                                                    star_props["noisy_Qpcs_diff_sum_Z2"][2]))),
            (star_props["m0"], star_props["Qxy"], star_props["Qpcs"]),
            (model_m0s, model_Qxys, model_Qpcss),
            (noisy_model_m0s, noisy_model_Qxys, noisy_model_Qpcss),
            (star_props["m0_diff_err"], star_props["Qxy_diff_err"], star_props["Qpcs_diff_err"]),
            (omask, Qxy_omask, Qpcs_omask),
            (star_props["m0_diff_diff_Zs"],
             star_props["Qxy_diff_diff_Zs"],
             star_props["Qpcs_diff_diff_Zs"],
             star_props["Qpcs_diff_sum_Zs"]),
            (star_props["x_pix"], star_props["y_pix"]),
            )
    
    if fitting_record is not None:
        fitting_record.append(test_results)

    return test_results
