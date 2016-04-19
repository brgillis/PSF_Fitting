""" @file fit_params.py

    Created 15 Apr 2016

    Methods to fit a full set of params for the best PSF model

    ---------------------------------------------------------------------

    Copyright (C) 2016 Bryan Gillis

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
from scipy.optimize import minimize

from psf_testing import magic_values as mv
from psf_testing.test_psf_for_params import test_psf_for_params

# Magic values

# Offsets so that the solver can use relative errors without worrying that these params
# cross zero
focus_offset = 100
param_offset = 1


def get_X2_of_test_results(test_results):
    # Get X2
    return test_results[1][0]
    # Get chi2
    # return test_results[1][2]

def get_params_penalty(focus,
                       focus_penalty_sigma=mv.default_focus_penalty_sigma,
                       penalty_sigma=mv.default_penalty_sigma,
                       **params):
    
    # Use 0 as a flag to impose no penalty
    if focus_penalty_sigma==0:
        penalty = 0
    else:
        penalty = ((focus-mv.default_init_test_focus)/focus_penalty_sigma)**2
    
    if penalty_sigma > 0:
        for param in params:
            penalty += ((params[param]-mv.default_params[param])/penalty_sigma)**2
    
    return penalty

def fit_best_params_and_test_psf(stars,

                                image_filename,
                                image=None,

                                test_focus=mv.default_init_test_focus,
                                penalty_sigma=mv.default_penalty_sigma,
                                focus_penalty_sigma=mv.default_focus_penalty_sigma,
                                
                                num_grid_points=mv.default_num_grid_points,

                                prim_weight_func=mv.default_prim_weight_func,
                                sec_weight_func=mv.default_sec_weight_func,

                                tinytim_path=mv.default_tinytim_path,
                                tinytim_data_path=mv.default_tinytim_data_path,

                                gain=mv.gain,
                                save_models=True,

                                files_to_cleanup=None,
                                
                                parallelize=False,
                                
                                **params):

    # Use a set outliers mask for all tests?
    outliers_mask = []
    
    # Keep a record of the results for all tests
    fitting_record = []

    # Define a function that we can use for fitting the focus
    def get_X2_for_params(test_param_array):
        test_focus = test_param_array[0]-focus_offset
        test_results = test_psf_for_params(stars=stars,

                                            image_filename=image_filename,
                                            image=image,

                                            test_focus=test_focus,
                                            
                                            num_grid_points=num_grid_points,

                                            prim_weight_func=prim_weight_func,
                                            sec_weight_func=sec_weight_func,

                                            tinytim_path=tinytim_path,
                                            tinytim_data_path=tinytim_data_path,

                                            gain=gain,
                                            save_models=False,
                                            outliers_mask=outliers_mask,

                                            files_to_cleanup=files_to_cleanup,
                                            
                                            fitting_record=fitting_record,
                                            
                                            parallelize=parallelize,
                                            
                                            astigmatism_0 = test_param_array[1]-param_offset,
                                            astigmatism_45 = test_param_array[2]-param_offset,
                                            spherical_3rd = test_param_array[3]-param_offset,
                                            spherical_5th = test_param_array[4]-param_offset,)
        return get_X2_of_test_results(test_results) + \
            get_params_penalty(test_focus,focus_penalty_sigma=focus_penalty_sigma,
                               penalty_sigma=penalty_sigma,
                               astigmatism_0 = test_param_array[1]-param_offset,
                               astigmatism_45 = test_param_array[2]-param_offset,
                               spherical_3rd = test_param_array[3]-param_offset,
                               spherical_5th = test_param_array[4]-param_offset,)
    
    # Initialize the test
    param_array = np.empty(5)
    param_array[0] = test_focus + focus_offset
    param_array[1] = params["astigmatism_0"] + param_offset
    param_array[2] = params["astigmatism_45"] + param_offset
    param_array[3] = params["spherical_3rd"] + param_offset
    param_array[4] = params["spherical_5th"] + param_offset


    best_param_array = minimize(get_X2_for_params, param_array, method='Nelder-Mead',
                                options={'xtol':0.0001})
    best_focus = best_param_array.x[0] - focus_offset
    best_params = {}
    best_params["astigmatism_0"] = best_param_array.x[1] - param_offset
    best_params["astigmatism_45"] = best_param_array.x[2] - param_offset
    best_params["spherical_3rd"] = best_param_array.x[3] - param_offset
    best_params["spherical_5th"] = best_param_array.x[4] - param_offset

    test_results = test_psf_for_params(stars=stars,

                                            image_filename=image_filename,
                                            image=image,

                                            test_focus=best_focus,
                                            num_grid_points=num_grid_points,

                                            prim_weight_func=prim_weight_func,
                                            sec_weight_func=sec_weight_func,

                                            tinytim_path=tinytim_path,
                                            tinytim_data_path=tinytim_data_path,

                                            gain=gain,
                                            save_models=save_models,
                                            outliers_mask=outliers_mask,

                                            files_to_cleanup=files_to_cleanup,
                                            
                                            fitted_params=1,
                                            
                                            fitting_record=None,
                                            
                                            parallelize=parallelize)

    return test_results, fitting_record