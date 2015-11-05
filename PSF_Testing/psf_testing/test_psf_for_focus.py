""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing/test_psf_for_focus.py

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
                       star_m0_errs,
                       star_Qs,
                       star_Q_errs,

                       image_filename,
                       chip,
                       image=None,

                       test_focus=mv.default_test_focus,
                       num_grid_points=mv.default_num_grid_points,

                       prim_weight_func=mv.default_prim_weight_func,
                       sec_weight_func=mv.default_sec_weight_func,

                       tinytim_path=mv.default_tinytim_path,
                       tinytim_data_path=mv.default_tinytim_data_path,
                       cleanup_tinytim_files=False,

                       files_to_cleanup=None,

                       gain=mv.gain,
                       save_models=False):

    if(image is None):
        from astropy.io import fits
        image = fits.open(image_filename)[0]

    if(files_to_cleanup is None):
        files_to_cleanup = []

    # Initialize arrays for results per star tested
    psf_m0s = []
    psf_m0_errs = []
    psf_Qs = []
    psf_Q_errs = []
    noisy_psf_m0s = []
    noisy_psf_m0_errs = []
    noisy_psf_Qs = []
    noisy_psf_Q_errs = []

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
                                           weight_func=weight_func,
                                           tinytim_path=tinytim_path,
                                           tinytim_data_path=tinytim_data_path)

        psf_m0, psf_m0_err, psf_Q, psf_Q_err = \
            get_m0_and_Qs(image=model_psf,
                          weight_func=weight_func,
                          background_noise=star.background_noise,
                          gain=gain)

        # Append m0 and Q data to the storage lists
        psf_m0s.append(psf_m0)
        psf_m0_errs.append(psf_m0_err)
        psf_Qs.append(psf_Q)
        psf_Q_errs.append(psf_Q_err)

        # Now, get a noisy psf and test it (so we can test for noise bias)
        model_psf_noise = np.sqrt(model_psf / gain + np.square(star.background_noise))

        ny, nx = np.shape(model_psf)

        noisy_model_psf = model_psf + model_psf_noise * np.random.randn(ny, nx)

        noisy_psf_m0, noisy_psf_m0_err, noisy_psf_Q, noisy_psf_Q_err = \
            get_m0_and_Qs(image=noisy_model_psf,
                          weight_func=weight_func,
                          background_noise=star.background_noise,
                          gain=gain)

        if(not (noisy_psf_m0 > 0)):
            pass

        # Append m0 and Q data to the storage lists
        noisy_psf_m0s.append(noisy_psf_m0)
        noisy_psf_m0_errs.append(noisy_psf_m0_err)
        noisy_psf_Qs.append(noisy_psf_Q)
        noisy_psf_Q_errs.append(noisy_psf_Q_err)

        # Save the models if desired
        if(save_models):
            star.model_psf = model_psf
            star.noisy_model_psf = noisy_model_psf

    # Convert the psf m0 and Q storage lists to numpy arrays
    psf_m0s = np.array(psf_m0s)
    psf_m0_errs = np.array(psf_m0_errs)
    psf_Qs = np.array(psf_Qs)
    psf_Q_errs = np.array(psf_Q_errs)
    noisy_psf_m0s = np.array(noisy_psf_m0s)
    noisy_psf_m0_errs = np.array(noisy_psf_m0_errs)
    noisy_psf_Qs = np.array(noisy_psf_Qs)
    noisy_psf_Q_errs = np.array(noisy_psf_Q_errs)

    # Now get the comparison Z values for both the psf and noisy psf

    # Get the differences first
    m0_diffs = star_m0s - psf_m0s
    Q_diffs = star_Qs - psf_Qs
    noisy_m0_diffs = star_m0s - noisy_psf_m0s
    noisy_Q_diffs = star_Qs - noisy_psf_Qs

    # And then get the Z values. Use the noise-free psf error estimate instead of the stars'
    # error estimates, since it should be less biased
    m0_Zs = m0_diffs / psf_m0_errs
    Q_Zs = Q_diffs / psf_Q_errs
    noisy_m0_Zs = noisy_m0_diffs / (np.sqrt(2.) * psf_m0_errs)
    noisy_Q_Zs = noisy_Q_diffs / (np.sqrt(2.) * psf_Q_errs)

    # Look for outliers on each of the Q arrays
    Qx_mask = remove_outliers(Q_Zs[:, 0]).mask
    Qy_mask = remove_outliers(Q_Zs[:, 1]).mask
    Qplus_mask = remove_outliers(Q_Zs[:, 2]).mask
    Qcross_mask = remove_outliers(Q_Zs[:, 3]).mask
    Qsize_mask = remove_outliers(Q_Zs[:, 4]).mask

    # If it's an outlier in any individual value, ignore it for all
    outliers_mask = np.logical_or.reduce((Qx_mask, Qy_mask, Qplus_mask, Qcross_mask, Qsize_mask))

    # For the Q array, we need to duplicate the mask
    Q_outliers_mask = np.vstack((outliers_mask, outliers_mask, outliers_mask, outliers_mask, outliers_mask)).transpose()

    masked_m0_Zs = np.ma.masked_array(m0_Zs, mask=outliers_mask)
    masked_Q_Zs = np.ma.masked_array(Q_Zs, mask=Q_outliers_mask)
    masked_noisy_m0_Zs = np.ma.masked_array(noisy_m0_Zs, mask=outliers_mask)
    masked_noisy_Q_Zs = np.ma.masked_array(noisy_Q_Zs, mask=Q_outliers_mask)

    num_masked = np.sum(outliers_mask)
    num_unmasked = len(m0_Zs) - num_masked

    # Get the mean Zs of the non-outliers
    mean_m0_Z = np.mean(masked_m0_Zs[~outliers_mask])
    mean_Q_Zs = np.mean(masked_Q_Zs[~Q_outliers_mask].reshape((num_unmasked, 5)), axis=0)
    mean_noisy_m0_Z = np.mean(masked_noisy_m0_Zs[~outliers_mask])
    mean_noisy_Q_Zs = np.mean(masked_noisy_Q_Zs[~Q_outliers_mask].reshape((num_unmasked, 5)), axis=0)

    return mean_m0_Z, mean_Q_Zs, mean_noisy_m0_Z, mean_noisy_Q_Zs
