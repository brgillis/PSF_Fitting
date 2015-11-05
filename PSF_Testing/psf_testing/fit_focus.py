""" @file fit_focus.py

    Created 5 Nov 2015

    Functions used to fit the best focus parameter.

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

from psf_testing import magic_values as mv
from psf_testing.brute_force_minimize import bf_minimize
from psf_testing.test_psf_for_focus import test_psf_for_focus


def get_chi2_of_test_results(test_results):
    return test_results[1][0]

def fit_best_focus_and_test_psf(stars,

                                star_m0s,
                                star_Qs,

                                image_filename,
                                image=None,

                                init_test_focus=mv.default_init_test_focus,

                                test_focus=mv.default_test_focus,
                                num_grid_points=mv.default_num_grid_points,

                                prim_weight_func=mv.default_prim_weight_func,
                                sec_weight_func=mv.default_sec_weight_func,

                                tinytim_path=mv.default_tinytim_path,
                                tinytim_data_path=mv.default_tinytim_data_path,

                                gain=mv.gain,
                                save_models=True,

                                outliers_mask=None):

    # Define a function that we can use for fitting the focus
    def get_chi2_for_focus(test_focus):
        test_results = test_psf_for_focus(stars=stars,

                                            star_m0s=star_m0s,
                                            star_Qs=star_Qs,

                                            image_filename=image_filename,
                                            image=image,

                                            test_focus=test_focus,
                                            num_grid_points=num_grid_points,

                                            prim_weight_func=prim_weight_func,
                                            sec_weight_func=sec_weight_func,

                                            tinytim_path=tinytim_path,
                                            tinytim_data_path=tinytim_data_path,

                                            gain=gain,
                                            save_models=True)
        return get_chi2_of_test_results(test_results)

    # Initialize the test

    best_focus = bf_minimize(get_chi2_for_focus, init_test_focus)

    test_results = test_psf_for_focus(stars=stars,

                                            star_m0s=star_m0s,
                                            star_Qs=star_Qs,

                                            image_filename=image_filename,
                                            image=image,

                                            test_focus=best_focus,
                                            num_grid_points=num_grid_points,

                                            prim_weight_func=prim_weight_func,
                                            sec_weight_func=sec_weight_func,

                                            tinytim_path=tinytim_path,
                                            tinytim_data_path=tinytim_data_path,

                                            gain=gain,
                                            save_models=True,
                                            fitted_params=1)

    return test_results
