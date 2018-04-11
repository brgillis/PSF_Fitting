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
focus_scale = 100
param_scale = 1


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
        penalty = ((focus-mv.default_focus)/focus_penalty_sigma)**2
    
    if penalty_sigma > 0:
        for param in params:
            penalty += ((params[param]-mv.default_params[param])/penalty_sigma)**2
    
    return penalty

def fit_best_params_and_test_psf(stars,

                                image_filename,
                                image=None,

                                focus=-1.0,
                                penalty_sigma=mv.default_penalty_sigma,
                                focus_penalty_sigma=mv.default_focus_penalty_sigma,
                                
                                num_grid_points=mv.default_num_grid_points,

                                prim_weight_func=mv.default_prim_weight_func,
                                sec_weight_func=mv.default_sec_weight_func,

                                tinytim_params=None,

                                gain=mv.gain,
                                save_models=True,

                                files_to_cleanup=None,
                                
                                parallelize=False,
                                          
                                norm_errors=False,
                                
                                seed=None,
                                
                                **params):

    # Use a set outliers mask for all tests?
    outliers_mask = []
    
    # Don't store fitting record just yet
    fitting_record = None

    # Define a function that we can use for fitting the parameters
    def get_X2_for_params(test_param_array):
        try:
            focus = focus_scale*(test_param_array[0]-1)
            test_results = test_psf_for_params(stars=stars,
    
                                                image_filename=image_filename,
                                                image=image,
    
                                                focus=focus,
                                                
                                                num_grid_points=num_grid_points,
    
                                                prim_weight_func=prim_weight_func,
                                                sec_weight_func=sec_weight_func,
    
                                                tinytim_params=tinytim_params,
    
                                                gain=gain,
                                                save_models=False,
                                                outliers_mask=outliers_mask,
    
                                                files_to_cleanup=files_to_cleanup,
                                                
                                                fitting_record=fitting_record,
                                                
                                                parallelize=parallelize,
                                                
                                                norm_errors=norm_errors,
                                                seed=seed,
                                                
                                                z2 = param_scale*(test_param_array[1]-1),
                                                z3 = param_scale*(test_param_array[2]-1),
                                                astigmatism_0 = param_scale*(test_param_array[3]-1),
                                                astigmatism_45 = param_scale*(test_param_array[4]-1),
                                                coma_x = param_scale*(test_param_array[5]-1),
                                                coma_y = param_scale*(test_param_array[6]-1),
                                                clover_x = param_scale*(test_param_array[7]-1),
                                                clover_y = param_scale*(test_param_array[8]-1),
                                                spherical_3rd = param_scale*(test_param_array[9]-1),
                                                z12 = param_scale*(test_param_array[10]-1),
                                                z13 = param_scale*(test_param_array[11]-1),
                                                z14 = param_scale*(test_param_array[12]-1),
                                                z15 = param_scale*(test_param_array[13]-1),
                                                z16 = param_scale*(test_param_array[14]-1),
                                                z17 = param_scale*(test_param_array[15]-1),
                                                z18 = param_scale*(test_param_array[16]-1),
                                                z19 = param_scale*(test_param_array[17]-1),
                                                z20 = param_scale*(test_param_array[18]-1),
                                                z21 = param_scale*(test_param_array[19]-1),
                                                spherical_5th = param_scale*(test_param_array[20]-1),
                                                kernel_adjustment = 1.5-param_scale*test_param_array[21])
            X2 = get_X2_of_test_results(test_results) + \
                get_params_penalty(focus,focus_penalty_sigma=focus_penalty_sigma,
                                   penalty_sigma=penalty_sigma,
                                    z2 = param_scale*(test_param_array[1]-1),
                                    z3 = param_scale*(test_param_array[2]-1),
                                    astigmatism_0 = param_scale*(test_param_array[3]-1),
                                    astigmatism_45 = param_scale*(test_param_array[4]-1),
                                    coma_x = param_scale*(test_param_array[5]-1),
                                    coma_y = param_scale*(test_param_array[6]-1),
                                    clover_x = param_scale*(test_param_array[7]-1),
                                    clover_y = param_scale*(test_param_array[8]-1),
                                    spherical_3rd = param_scale*(test_param_array[9]-1),
                                    z12 = param_scale*(test_param_array[10]-1),
                                    z13 = param_scale*(test_param_array[11]-1),
                                    z14 = param_scale*(test_param_array[12]-1),
                                    z15 = param_scale*(test_param_array[13]-1),
                                    z16 = param_scale*(test_param_array[14]-1),
                                    z17 = param_scale*(test_param_array[15]-1),
                                    z18 = param_scale*(test_param_array[16]-1),
                                    z19 = param_scale*(test_param_array[17]-1),
                                    z20 = param_scale*(test_param_array[18]-1),
                                    z21 = param_scale*(test_param_array[19]-1),
                                    spherical_5th = param_scale*(test_param_array[20]-1),
                                    kernel_adjustment = 1.5-param_scale*test_param_array[21])
        except Exception as _e:
            # Return a near-infinite value on exception
            X2 = 1.0e100
        
        print("X2 = " + str(X2))
        return X2
    
    # Initialize the test
    param_array = np.empty(23)
    param_array[0] = focus/focus_scale + 1
    param_array[1] = params["z2"]/param_scale + 1
    param_array[2] = params["z3"]/param_scale + 1
    param_array[3] = params["astigmatism_0"]/param_scale + 1
    param_array[4] = params["astigmatism_45"]/param_scale + 1
    param_array[5] = params["coma_x"]/param_scale + 1
    param_array[6] = params["coma_y"]/param_scale + 1
    param_array[7] = params["clover_x"]/param_scale + 1
    param_array[8] = params["clover_y"]/param_scale + 1
    param_array[9] = params["spherical_3rd"]/param_scale + 1
    param_array[10] = params["z12"]/param_scale + 1
    param_array[11] = params["z13"]/param_scale + 1
    param_array[12] = params["z14"]/param_scale + 1
    param_array[13] = params["z15"]/param_scale + 1
    param_array[14] = params["z16"]/param_scale + 1
    param_array[15] = params["z17"]/param_scale + 1
    param_array[16] = params["z18"]/param_scale + 1
    param_array[17] = params["z19"]/param_scale + 1
    param_array[18] = params["z20"]/param_scale + 1
    param_array[19] = params["z21"]/param_scale + 1
    param_array[20] = params["spherical_5th"]/param_scale + 1
    param_array[21] = (1.5 - params["kernel_adjustment"])/param_scale
    
#     # Calculate (and cache) the value for init params first, so we'll always use the
#     # outliers list for that
#     get_X2_for_params(param_array)

    # Store fitting record
    fitting_record = []
    
    best_param_array = minimize(get_X2_for_params, param_array, method='Nelder-Mead',
                                options={'xtol':0.001,'ftol':0.001}).x

    best_focus = focus_scale*(best_param_array[0] - 1)
    best_params = {}
    best_params["z2"] = param_scale*(best_param_array[1]-1)
    best_params["z3"] = param_scale*(best_param_array[2]-1)
    best_params["astigmatism_0"] = param_scale*(best_param_array[3]-1)
    best_params["astigmatism_45"] = param_scale*(best_param_array[4]-1)
    best_params["coma_x"] = param_scale*(best_param_array[5]-1)
    best_params["coma_y"] = param_scale*(best_param_array[6]-1)
    best_params["clover_x"] = param_scale*(best_param_array[7]-1)
    best_params["clover_y"] = param_scale*(best_param_array[8]-1)
    best_params["spherical_3rd"] = param_scale*(best_param_array[9]-1)
    best_params["z12"] = param_scale*(best_param_array[10]-1)
    best_params["z13"] = param_scale*(best_param_array[11]-1)
    best_params["z14"] = param_scale*(best_param_array[12]-1)
    best_params["z15"] = param_scale*(best_param_array[13]-1)
    best_params["z16"] = param_scale*(best_param_array[14]-1)
    best_params["z17"] = param_scale*(best_param_array[15]-1)
    best_params["z18"] = param_scale*(best_param_array[16]-1)
    best_params["z19"] = param_scale*(best_param_array[17]-1)
    best_params["z20"] = param_scale*(best_param_array[18]-1)
    best_params["z21"] = param_scale*(best_param_array[19]-1)
    best_params["spherical_5th"] = param_scale*(best_param_array[20]-1)
    best_params["kernel_adjustment"] = 1.5-param_scale*best_param_array[21]

    test_results = test_psf_for_params(stars=stars,

                                            image_filename=image_filename,
                                            image=image,

                                            focus=best_focus,
                                            num_grid_points=num_grid_points,

                                            prim_weight_func=prim_weight_func,
                                            sec_weight_func=sec_weight_func,

                                            tinytim_params=tinytim_params,

                                            gain=gain,
                                            save_models=save_models,
                                            outliers_mask=outliers_mask,

                                            files_to_cleanup=files_to_cleanup,
                                            
                                            fitted_params=1,
                                            
                                            fitting_record=None,
                                            
                                            parallelize=parallelize,
                                            
                                            norm_errors=norm_errors,
                                            
                                            z2 = param_scale*(best_param_array[1]-1),
                                            z3 = param_scale*(best_param_array[2]-1),
                                            astigmatism_0 = param_scale*(best_param_array[3]-1),
                                            astigmatism_45 = param_scale*(best_param_array[4]-1),
                                            coma_x = param_scale*(best_param_array[5]-1),
                                            coma_y = param_scale*(best_param_array[6]-1),
                                            clover_x = param_scale*(best_param_array[7]-1),
                                            clover_y = param_scale*(best_param_array[8]-1),
                                            spherical_3rd = param_scale*(best_param_array[9]-1),
                                            z12 = param_scale*(best_param_array[10]-1),
                                            z13 = param_scale*(best_param_array[11]-1),
                                            z14 = param_scale*(best_param_array[12]-1),
                                            z15 = param_scale*(best_param_array[13]-1),
                                            z16 = param_scale*(best_param_array[14]-1),
                                            z17 = param_scale*(best_param_array[15]-1),
                                            z18 = param_scale*(best_param_array[16]-1),
                                            z19 = param_scale*(best_param_array[17]-1),
                                            z20 = param_scale*(best_param_array[18]-1),
                                            z21 = param_scale*(best_param_array[19]-1),
                                            spherical_5th = param_scale*(best_param_array[20]-1),
                                            kernel_adjustment = 1.5-param_scale*best_param_array[21],)

    return test_results, fitting_record
