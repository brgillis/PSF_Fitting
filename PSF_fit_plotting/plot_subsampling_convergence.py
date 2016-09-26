#!/usr/bin/env python

""" @file plot_subsampling_convergence.py

    Created 26 Sep 2016

    This script is used to plot star/PSF/residual stacks which are
    generated as a result of running PSFit.

    ---------------------------------------------------------------------

    Copyright (C) 2016  Bryan R. Gillis

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

import os
import sys

from astropy.io import fits
import matplotlib

import matplotlib.pyplot as pyplot
import numpy as np


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

from psf_testing.magic_values import rel_weights

data_dir = "/disk2/brg/Data/HST_Fields/subsampling_convergence_testing"
results_file_root = "jb6v09shq_sci2_cor"

figsize = (6, 6)
fontsize = 12

def main(argv):
    """ @TODO main docstring
    """

    results_root = os.path.join(data_dir, results_file_root)

    ss_factors = np.linspace(1, 10, 10, endpoint=True).astype(int)

    X2s = np.zeros_like(ss_factors, dtype=float)
    Qp_sum_Z2s = np.zeros_like(ss_factors, dtype=float)
    Qp_diff_Z2s = np.zeros_like(ss_factors, dtype=float)
    Qc_sum_Z2s = np.zeros_like(ss_factors, dtype=float)
    Qc_diff_Z2s = np.zeros_like(ss_factors, dtype=float)
    Qs_sum_Z2s = np.zeros_like(ss_factors, dtype=float)
    Qs_diff_Z2s = np.zeros_like(ss_factors, dtype=float)

    for ssf in ss_factors:

        filename = results_root + "_ss" + str(ssf) + "_results.fits"

        header = fits.open(filename)[1].header

        X2s[ssf - 1] = header["X_SQR"]
        Qp_sum_Z2s[ssf - 1] = header["QPS_Z2"] / np.square(rel_weights["Qpcs_diff"][0])
        Qp_diff_Z2s[ssf - 1] = header["QPD_Z2"] / np.square(rel_weights["Qpcs_diff"][0])
        Qc_sum_Z2s[ssf - 1] = header["QPS_Z2"] / np.square(rel_weights["Qpcs_diff"][1])
        Qc_diff_Z2s[ssf - 1] = header["QPD_Z2"] / np.square(rel_weights["Qpcs_diff"][1])
        Qs_sum_Z2s[ssf - 1] = header["QSS_Z2"] / np.square(rel_weights["Qpcs_diff"][2])
        Qs_diff_Z2s[ssf - 1] = header["QSD_Z2"] / np.square(rel_weights["Qpcs_diff"][2])

    # Plot it up
    _fig = pyplot.figure(figsize=figsize)

    gs = matplotlib.gridspec.GridSpec(1, 1)
    gs.update(wspace=0.3, hspace=0.3, left=0.1, right=0.9, bottom=0.1, top=0.9)

    ax = pyplot.subplot(gs[0])

    ax.plot(ss_factors, X2s, label="X2")
    ax.plot(ss_factors, Qp_sum_Z2s, label="Qp sum Z2")
    ax.plot(ss_factors, Qp_diff_Z2s, label="Qp diff Z2")
    ax.plot(ss_factors, Qc_sum_Z2s, label="Qc sum Z2")
    ax.plot(ss_factors, Qc_diff_Z2s, label="Qc diff Z2")
    ax.plot(ss_factors, Qs_sum_Z2s, label="Qs sum Z2")
    ax.plot(ss_factors, Qs_diff_Z2s, label="Qs diff Z2")

    ax.set_xlim([0.5, 10.5])
    ax.set_ylim([0, 0.002])

    ax.legend(loc="upper right")

    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)
