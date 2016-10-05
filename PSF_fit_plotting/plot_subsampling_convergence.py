#!/usr/bin/env python

""" @file plot_subsampling_convergence.py

    Created 26 Sep 2016

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
from argparse import ArgumentParser

from astropy.io import fits
import matplotlib

import matplotlib.pyplot as pyplot
import numpy as np


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

data_dir = "/disk2/brg/Data/HST_Fields/subsampling_convergence_testing"
default_results_file_root = "control_image_n"

max_ssf = 15

figsize = (6, 6)
fontsize = 12

default_output_filename_root = "subsampling_convergence_test"
default_output_extension = "eps"

def main(argv):
    """ @TODO main docstring
    """
    
    parser = ArgumentParser()
    
    parser.add_argument("--results_file_root",default=default_results_file_root)
    parser.add_argument("--output_filename_root",default=default_output_filename_root)
    parser.add_argument("--output_extension",default=default_output_extension)
    
    args = parser.parse_args()

    results_root = os.path.join(data_dir, args.results_file_root)

    ss_factors = np.linspace(1, max_ssf, max_ssf, endpoint=True).astype(int)

    X2s = np.zeros_like(ss_factors, dtype=float)
    Qx_diff_Z2s = np.zeros_like(ss_factors, dtype=float)
    Qy_diff_Z2s = np.zeros_like(ss_factors, dtype=float)
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
        Qx_diff_Z2s[ssf - 1] = header["QXD_Z2"]
        Qy_diff_Z2s[ssf - 1] = header["QYD_Z2"]
        Qp_sum_Z2s[ssf - 1] = header["QPS_Z2"] 
        Qp_diff_Z2s[ssf - 1] = header["QPD_Z2"] 
        Qc_sum_Z2s[ssf - 1] = header["QPS_Z2"] 
        Qc_diff_Z2s[ssf - 1] = header["QPD_Z2"] 
        Qs_sum_Z2s[ssf - 1] = header["QSSNZ2"] 
        Qs_diff_Z2s[ssf - 1] = header["QSDNZ2"] 

    # Plot it up
    _fig = pyplot.figure(figsize=figsize)

    gs = matplotlib.gridspec.GridSpec(1, 1)
    gs.update(wspace=0.3, hspace=0.3, left=0.1, right=0.9, bottom=0.1, top=0.9)

    ax = pyplot.subplot(gs[0])

    ax.plot(ss_factors, X2s, label="$X^2$", color="black")
    ax.plot(ss_factors, Qx_diff_Z2s, label=r"$Q_x^{(-)} Z^2$", linestyle="dotted")
    ax.plot(ss_factors, Qy_diff_Z2s, label=r"$Q_y^{(-)} Z^2$", linestyle="dotted")
    ax.plot(ss_factors, Qp_sum_Z2s, label=r"$Q_+^{(+)} Z^2$", linestyle="dashed")
    ax.plot(ss_factors, Qp_diff_Z2s, label=r"$Q_+^{(-)} Z^2$", linestyle="dashed")
    ax.plot(ss_factors, Qc_sum_Z2s, label=r"$Q_{\times}^{(+)} Z^2$", linestyle="dashed")
    ax.plot(ss_factors, Qc_diff_Z2s, label=r"$Q_{\times}^{(-)} Z^2$", linestyle="dashed")
    ax.plot(ss_factors, Qs_sum_Z2s, label=r"$Q_s^{(+)} Z^2$", linestyle="dashdot")
    ax.plot(ss_factors, Qs_diff_Z2s, label=r"$Q_s^{(-)} Z^2$", linestyle="dashdot")

    ax.set_xlim([0.5, max_ssf+0.5])
    # ax.set_ylim([0, 0.02])
    ax.set_yscale("log", nonposy='clip')
    
    ax.set_xlabel("Subsampling Factor")

    ax.legend(loc="lower left",ncol=2)
    
    output_filename = args.output_filename_root + "." + args.output_extension

    pyplot.savefig(output_filename, format=args.output_extension, bbox_inches="tight", pad_inches=0.05)

    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)
