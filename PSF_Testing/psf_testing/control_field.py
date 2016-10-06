""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing/control_field.py

    Created 16 Sep 2016

    Function to generate a control field to test the PDF testing pipeline.

    ---------------------------------------------------------------------

    Copyright (C) 2016 Bryan R. Gillis

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

import multiprocessing

import galsim

import numpy as np
from psf_testing import magic_values as mv
from psf_testing.get_model_psf import get_cached_subsampled_psf
from psf_testing.get_model_psf import get_model_psf
from psf_testing.moments.centre_image import centre_image
from psf_testing.psf_model_scheme import psf_model_scheme


spec_types = np.linspace(0, 18, 19, endpoint=True)
spec_types_pdf = spec_types ** 3
spec_types_pdf[0] = spec_types_pdf[-1] = 0
spec_types_pdf /= spec_types_pdf.sum()

class get_psf_caller(object):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
    def __call__(self, x):
        return get_cached_subsampled_psf(psf_position=x[0], spec_type=x[1], *self.args, **self.kwargs)


def generate_star(image_shape, min_star_mag, max_star_mag, randomize_spectral_type, u, spec_f):

    star = {}
    star["xp"] = u() * image_shape[0]
    star["yp"] = u() * image_shape[1]
    star["mag"] = min_star_mag + u() * (max_star_mag - min_star_mag)

    # Choose a random spectral type if desired
    if randomize_spectral_type:
        st = int(round(spec_f()))
        st = max(1, min(17, st))
        star["spec_type"] = 1, st
    else:
        star["spec_type"] = mv.default_model_psf_spec_type
    return star

def make_control_field(image_filename,

                       random_seed=0,

                       image_shape=mv.default_image_shape,
                       exp_time=1298.,
                       scale=mv.pixel_scale,
                       gain=mv.gain,
                       zeropoint=mv.zeropoint,
                       read_noise=mv.read_noise,
                       sky_level=mv.sky_level,

                       num_stars=1000,
                       min_star_mag=mv.default_min_star_mag,
                       max_star_mag=mv.default_max_star_mag,
                       binary_fraction=0.,
                       binary_r_max=1.,

                       suppress_noise=False,

                       focus=mv.default_focus,
                       tinytim_params=None,
                       subsampling_factor=mv.default_subsampling_factor,
                       num_grid_points=(0, 0),
                       randomize_spectral_type=True,
                       use_cache=True,

                       parallelize=True):

    if tinytim_params is None:
        tinytim_params = mv.default_tinytim_params
    tinytim_params['subsampling_factor'] = subsampling_factor

    scheme = psf_model_scheme(focus=focus,
                              num_grid_points=num_grid_points,
                              image_shape=image_shape)

    image = galsim.Image(image_shape[0], image_shape[1], scale=scale)

    im_zeropoint = zeropoint + 2.5 * np.log10(exp_time)

    u = galsim.random.UniformDeviate(galsim.BaseDeviate(random_seed))
    spec_f = galsim.DistDeviate(galsim.BaseDeviate(random_seed + 1), function=galsim.LookupTable(spec_types, spec_types_pdf))

    # Set up star data first, to aid parallelization
    stars = []
    for _i in range(num_stars):
        star = generate_star(image_shape, min_star_mag, max_star_mag, randomize_spectral_type, u, spec_f)

        stars.append(star)

    # Generate binary stars now
    for i in range(num_stars):
        # Do a random check to see if we'll make this star a binary
        if u() > binary_fraction:
            continue

        # Generate another star
        extra_star = generate_star(image_shape, min_star_mag, max_star_mag, randomize_spectral_type, u, spec_f)

        # Override its position to be near the other star
        r_deviate = extra_star["xp"] / image_shape[0]
        theta_deviate = extra_star["yp"] / image_shape[1] # These will be uniform [0,1) now

        r = r_deviate ** 2 * binary_r_max
        theta = theta_deviate * 2 * np.pi

        extra_star["xp"] = star[i]["xp"] + r * np.cos(theta)
        extra_star["yp"] = star[i]["yp"] + r * np.sin(theta)

        stars.append(extra_star)

    # If we're parallelizing and using the cache, get possible psfs in parallel
    if parallelize and use_cache:
        points = set()
        for star in stars:
            new_point = (scheme.get_position_to_use(star["xp"], star["yp"]), star["spec_type"])
            points.add(new_point)

        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count(), maxtasksperchild=1)
        pool.map(get_psf_caller(tinytim_params_set=frozenset(tinytim_params.items()),
                                weight_func=mv.default_prim_weight_func,
                                focus=focus,
                                use_cache=use_cache), points, chunksize=1)
        pool.close()
        pool.join()
        pool.terminate()

    for i in range(num_stars):

        # Calculate the flux for this star given its magnitude
        flux = 10.0 ** (0.4 * (im_zeropoint - stars[i]["mag"]))

        psf_data = get_model_psf(stars[i]["xp"], stars[i]["yp"],
                                 scheme=scheme,
                                 tinytim_params=tinytim_params,
                                 use_cache=use_cache,
                                 spec_type=stars[i]["spec_type"])

        psf_data *= flux / psf_data.sum()

        # Centre the psf image
        psf_ny, psf_nx = np.shape(psf_data)
        psf_yc, psf_xc = centre_image(psf_data)[0:2]

        psf_x_arm = int(min([(psf_nx - 1.) / 2, psf_xc, psf_nx - psf_xc]))
        psf_y_arm = int(min([(psf_ny - 1.) / 2, psf_yc, psf_ny - psf_yc]))

        image_window = image.array
        psf_window = psf_data[int(psf_yc) - psf_y_arm : int(psf_yc) + psf_y_arm + 1 ,
                              int(psf_xc) - psf_x_arm : int(psf_xc) + psf_x_arm + 1 ]

        xp = int(stars[i]["xp"])
        yp = int(stars[i]["yp"])

        if yp < psf_y_arm:
            psf_window = psf_window[psf_y_arm - yp:, :]
            image_window = image_window[:yp + psf_y_arm + 1, :]
        elif yp > image_shape[1] - psf_y_arm - 1:
            psf_window = psf_window[:image_shape[1] - psf_y_arm - 1 - yp, :]
            image_window = image_window[yp - psf_y_arm:, :]
        else:
            image_window = image_window[yp - psf_y_arm:yp + psf_y_arm + 1, :]

        if xp < psf_x_arm:
            psf_window = psf_window[:, psf_x_arm - xp:]
            image_window = image_window[:, :xp + psf_x_arm + 1]
        elif xp > image_shape[0] - psf_x_arm - 1:
            psf_window = psf_window[:, :image_shape[0] - psf_x_arm - 1 - xp]
            image_window = image_window[:, xp - psf_x_arm:]
        else:
            image_window = image_window[:, xp - psf_x_arm:xp + psf_x_arm + 1]

        if not np.shape(image_window) == np.shape(psf_window):
            pass
        image_window += psf_window

        pass

    # Set up the image's header
    image.header = {}
    image.header[mv.header_chip_keyword] = tinytim_params["chip"]
    image.header[mv.header_exp_time_keyword] = exp_time
    image.header[mv.header_gain_keyword] = gain
    image.header[mv.header_obs_date_keyword] = 0.
    image.header[mv.header_obs_time_keyword] = 0.
    image.header[mv.header_ra_keyword] = 0.
    image.header[mv.header_dec_keyword] = 0.

    # Add noise to the image
    if not suppress_noise:
        image.addNoise(galsim.CCDNoise(galsim.BaseDeviate(random_seed + 2), gain=gain, read_noise=read_noise, sky_level=sky_level))

    # Output the image
    image.write(image_filename, clobber=True)

    return
