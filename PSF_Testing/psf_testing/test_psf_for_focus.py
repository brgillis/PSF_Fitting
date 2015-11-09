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


def test_psf_for_focus(stars,

                       star_m0s,
                       star_Qs,

                       image_filename,
                       image=None,

                       test_focus=mv.default_test_focus,
                       num_grid_points=mv.default_num_grid_points,

                       prim_weight_func=mv.default_prim_weight_func,
                       sec_weight_func=mv.default_sec_weight_func,

                       tinytim_path=mv.default_tinytim_path,
                       tinytim_data_path=mv.default_tinytim_data_path,

                       gain=mv.gain,
                       save_models=False,

                       outliers_mask=None,

                       fitted_params=0,

                       files_to_cleanup=None):

    if outliers_mask is None:
        outliers_mask = []

    if image is None:
        from astropy.io import fits
        image = fits.open(image_filename)[0]

    # Initialize arrays for results per star tested
    psf_m0s = []
    psf_m0_vars = []
    psf_Qs = []
    psf_Q_vars = []
    noisy_psf_m0s = []
    noisy_psf_m0_vars = []
    noisy_psf_Qs = []
    noisy_psf_Q_vars = []

    # Get the image shape, and reverse its ordering to x,y in fits ordering
    image_shape = np.shape(image)
    image_shape = image_shape[1], image_shape[0]

    # Set up the focus generating scheme
    model_scheme = psf_model_scheme(focus=test_focus,
                                    num_grid_points=num_grid_points,
                                    image_shape=image_shape)

    # Loop through stars and get results for each
    for star in stars:

        if not star.valid:
            continue

        model_psf = get_model_psf_for_star(star=star,
                                           scheme=model_scheme,
                                           weight_func=prim_weight_func,
                                           tinytim_path=tinytim_path,
                                           tinytim_data_path=tinytim_data_path)

        psf_m0, psf_m0_var, psf_Q, psf_Q_var = \
            get_m0_and_Qs(image=model_psf,
                          prim_weight_func=prim_weight_func,
                          sec_weight_func=sec_weight_func,
                          background_noise=star.background_noise,
                          gain=gain)

        # Append m0 and Q data to the storage lists
        psf_m0s.append(psf_m0)
        psf_m0_vars.append(psf_m0_var)
        psf_Qs.append(psf_Q)
        psf_Q_vars.append(psf_Q_var)

        # Now, get a noisy psf and test it (so we can test for noise bias)
        model_psf_noise = np.sqrt(model_psf / gain + np.square(star.background_noise))

        ny, nx = np.shape(model_psf)

        noisy_model_psf = model_psf + model_psf_noise * np.random.randn(ny, nx)

        noisy_psf_m0, noisy_psf_m0_var, noisy_psf_Q, noisy_psf_Q_var = \
            get_m0_and_Qs(image=noisy_model_psf,
                          prim_weight_func=prim_weight_func,
                          sec_weight_func=sec_weight_func,
                          background_noise=star.background_noise,
                          gain=gain)

        # Append m0 and Q data to the storage lists
        noisy_psf_m0s.append(noisy_psf_m0)
        noisy_psf_m0_vars.append(noisy_psf_m0_var)
        noisy_psf_Qs.append(noisy_psf_Q)
        noisy_psf_Q_vars.append(noisy_psf_Q_var)

        # Save the models if desired
        if save_models:
            star.model_psf = model_psf
            star.noisy_model_psf = noisy_model_psf

    # Convert the psf m0 and Q storage lists to numpy arrays
    psf_m0s = np.array(psf_m0s)
    psf_m0_vars = np.array(psf_m0_vars)
    psf_Qs = np.array(psf_Qs)
    psf_Q_vars = np.array(psf_Q_vars)
    noisy_psf_m0s = np.array(noisy_psf_m0s)
    noisy_psf_m0_vars = np.array(noisy_psf_m0_vars)
    noisy_psf_Qs = np.array(noisy_psf_Qs)
    noisy_psf_Q_vars = np.array(noisy_psf_Q_vars)

    # Now get the comparison Z values for both the psf and noisy psf

    # Get the differences first
    m0_diffs = star_m0s - psf_m0s
    Q_diffs = star_Qs - psf_Qs
    noisy_m0_diffs = star_m0s - noisy_psf_m0s
    noisy_Q_diffs = star_Qs - noisy_psf_Qs

    # And then get the Z values. Use the noise-free psf error estimate instead of the stars'
    # error estimates, since it should be less biased
    Q_signs = np.array([[-1, -1, 1, 1, 1]])
    comb_Q_diffs = (Q_signs * Q_diffs[:, :, 0] + Q_diffs[:, :, 1])
    comb_noisy_Q_diffs = (Q_signs * noisy_Q_diffs[:, :, 0] + noisy_Q_diffs[:, :, 1])

    m0_Zs = (m0_diffs[:, 0] + m0_diffs[:, 1]) / np.sqrt(psf_m0_vars[:, 0, 0] +
                                                          psf_m0_vars[:, 1, 1] +
                                                          2.0 * psf_m0_vars[:, 1, 0])
    Q_Zs = comb_Q_diffs / np.sqrt(psf_Q_vars[:, :, 0, 0] +
                                                          psf_Q_vars[:, :, 1, 1] +
                                                          2.0 * Q_signs * psf_Q_vars[:, :, 1, 0])
    noisy_m0_Zs = (noisy_m0_diffs[:, 0] + noisy_m0_diffs[:, 1]) / np.sqrt(psf_m0_vars[:, 0, 0] +
                                                          psf_m0_vars[:, 1, 1] +
                                                          2.0 * psf_m0_vars[:, 1, 0])
    noisy_Q_Zs = comb_noisy_Q_diffs / np.sqrt(psf_Q_vars[:, :, 0, 0] +
                                                          psf_Q_vars[:, :, 1, 1] +
                                                          2.0 * Q_signs * psf_Q_vars[:, :, 1, 0])

    if len(outliers_mask) == 0:

        # Look for outliers on each of the Q arrays
        Qx_mask = remove_outliers(Q_Zs[:, 0]).mask
        Qy_mask = remove_outliers(Q_Zs[:, 1]).mask
        Qplus_mask = remove_outliers(Q_Zs[:, 2]).mask
        Qcross_mask = remove_outliers(Q_Zs[:, 3]).mask
        Qsize_mask = remove_outliers(Q_Zs[:, 4]).mask

        # If it's an outlier in any individual value, ignore it for all
        outliers_mask.append(np.logical_or.reduce((Qx_mask, Qy_mask, Qplus_mask, Qcross_mask, Qsize_mask)))

    omask = outliers_mask[0]

    # For the Q array, we need to duplicate the mask
    m0_omask = np.vstack((omask,omask)).transpose()
    Q_omask = np.vstack((omask, omask, omask, omask, omask)).transpose()

    masked_m0_diffs = np.ma.masked_array(m0_diffs, mask=m0_omask)
    masked_Q_diffs = np.ma.masked_array(comb_Q_diffs, mask=Q_omask)
    masked_m0_Zs = np.ma.masked_array(m0_Zs, mask=omask)
    masked_Q_Zs = np.ma.masked_array(Q_Zs, mask=Q_omask)

    masked_noisy_m0_diffs = np.ma.masked_array(noisy_m0_diffs, mask=m0_omask)
    masked_noisy_Q_diffs = np.ma.masked_array(comb_noisy_Q_diffs, mask=Q_omask)
    masked_noisy_m0_Zs = np.ma.masked_array(noisy_m0_Zs, mask=omask)
    masked_noisy_Q_Zs = np.ma.masked_array(noisy_Q_Zs, mask=Q_omask)

    num_masked = np.sum(omask)
    num_unmasked = len(m0_Zs) - num_masked

    # Get the mean Zs of the non-outliers
    mean_m0_diff = np.mean(masked_m0_diffs[~m0_omask].reshape((num_unmasked, 2)), axis=0)
    mean_Q_diffs = np.mean(masked_Q_diffs[~Q_omask].reshape((num_unmasked, 5)), axis=0)
    mean_m0_Z = np.mean(masked_m0_Zs[~omask], axis=0)
    mean_Q_Zs = np.mean(masked_Q_Zs[~Q_omask].reshape((num_unmasked, 5)), axis=0)
    mean_noisy_m0_diff = np.mean(masked_noisy_m0_diffs[~m0_omask].reshape((num_unmasked, 2)), axis=0)
    mean_noisy_Q_diffs = np.mean(masked_noisy_Q_diffs[~Q_omask].reshape((num_unmasked, 5)), axis=0)
    mean_noisy_m0_Z = np.mean(masked_noisy_m0_Zs[~omask], axis=0)
    mean_noisy_Q_Zs = np.mean(masked_noisy_Q_Zs[~Q_omask].reshape((num_unmasked, 5)), axis=0)

    Qx_Z2 = np.sum(np.square(Q_Zs[~Q_omask].reshape((num_unmasked, 5))[:, 0]), axis=0)
    Qy_Z2 = np.sum(np.square(Q_Zs[~Q_omask].reshape((num_unmasked, 5))[:, 1]), axis=0)
    Qplus_Z2 = np.sum(np.square(Q_Zs[~Q_omask].reshape((num_unmasked, 5))[:, 2]), axis=0)
    Qcross_Z2 = np.sum(np.square(Q_Zs[~Q_omask].reshape((num_unmasked, 5))[:, 3]), axis=0)
    Qsize_Z2 = np.sum(np.square(noisy_Q_Zs[~Q_omask].reshape((num_unmasked, 5))[:, 4]), axis=0)

    chi2 = Qx_Z2 + Qy_Z2 + Qplus_Z2 + Qcross_Z2 + Qsize_Z2

    dof = 6 * num_unmasked - fitted_params

    return (test_focus,
            (chi2, dof),
            ((mean_m0_diff, mean_Q_diffs), (mean_noisy_m0_diff, mean_noisy_Q_diffs)),
            ((mean_m0_Z, mean_Q_Zs), (mean_noisy_m0_Z, mean_noisy_Q_Zs)),
            (star_m0s, star_Qs),
            (psf_m0s, psf_Qs),
            (noisy_psf_m0s, noisy_psf_Qs),
            (np.sqrt(psf_m0_vars), np.sqrt(psf_Q_vars)),
            (omask, Q_omask),
            (m0_Zs, Q_Zs),
            )
