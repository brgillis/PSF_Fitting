#!/usr/bin/env python

""" @file plot_X_v_focus_v_chip.py

    Created 25 Apr 2016

    @TODO: File docstring

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

import argparse
import sys

from astropy.io import fits
import matplotlib
from scipy.interpolate import UnivariateSpline

import matplotlib.pyplot as pyplot
import numpy as np
import subprocess as sbp


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


default_summary_filename = "/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_all_params_summary.fits"
default_model_filename = "/home/brg/Data/HST_Fields/ModelFocus.fits"

default_plot_name = "model_focus_v_focus_v_chip_all_params"
default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"
default_file_type = "eps"

figsize = (6, 6)

fontsize = 12


def make_model_focus_v_focus_v_chip_plot(summary_filename=default_summary_filename,
                                         model_filename=default_model_filename,

                                         plot_name=default_plot_name,
                                         paper_location=default_paper_location,
                                         file_type=default_file_type,

                                         hide=False,
                                         ):

    summary_table = fits.open(summary_filename)[1].data
    model_table = fits.open(model_filename)[1].data

    fitted_focus = summary_table["focus"]

    # model_table.sort("t_sec")  # Make sure it's strictly increasing by time

    model_t_sec = model_table["t_sec"]
    model_focus = model_table["focus"]

    # We'll need to mask out any duplicate time values to make sure it's strictly increasing

    strictly_increasing = False

    while not strictly_increasing:

        t_diff = np.ones_like(model_t_sec)
        t_diff[0:-1] = model_t_sec[1:] - model_t_sec[:-1]
        t_mask = (t_diff <= 0)

        if t_mask.sum() == 0:
            strictly_increasing = True
        else:
            model_t_sec = np.ma.masked_array(model_t_sec, t_mask).compressed()
            model_focus = np.ma.masked_array(model_focus, t_mask).compressed()

    model_focus_spline = UnivariateSpline(model_t_sec, model_focus, s=0)

    interpolated_model_focus = np.zeros_like(fitted_focus)

    for i in range(len(fitted_focus)):
        interpolated_model_focus[i] = model_focus_spline(summary_table["obs_time"][i])

    _fig = pyplot.figure(figsize=figsize)

    gs = matplotlib.gridspec.GridSpec(1, 1)
    gs.update(wspace=0.3, hspace=0.3, left=0.1, right=0.9, bottom=0.1, top=0.9)

    ax = pyplot.subplot(gs[0])
    y_label = r"Model focus offset (microns)"

    chip1_mask = ~(summary_table["chip"] == 1)
    chip2_mask = ~(summary_table["chip"] == 2)

    for mask, label, color, marker in ((chip1_mask, "Chip 1", '#C00000', "o"),
                                       (chip2_mask, "Chip 2", '#4040FF', "^")):
        fitted_focii = np.ma.masked_array(fitted_focus, mask).compressed()
        model_focii = np.ma.masked_array(interpolated_model_focus, mask).compressed()

        pyplot.scatter(fitted_focii, model_focii, edgecolors=color,
                       label=label, alpha=1, marker=marker, facecolors='none')

    ax.set_xlim(ax.get_xlim())
    ax.set_ylim(ax.get_ylim())

    ax.set_xlabel("Best-fit focus offset (microns)", fontsize=fontsize)
    ax.set_ylabel(y_label, fontsize=fontsize, labelpad=-8)

    # Draw one-to-one line

    one_to_one_focus = np.linspace(-10, 10, 1000)

    pyplot.plot(one_to_one_focus, one_to_one_focus, linestyle="dashed", color="k")

    ax.legend(loc="upper left")

    # Save the figure
    outfile_name = plot_name + "." + file_type
    pyplot.savefig(outfile_name, format=file_type, bbox_inches="tight", pad_inches=0.05)

    # Copy it to the paper location if in eps format
    if file_type == "eps":
        cmd = "cp " + outfile_name + " " + paper_location
        sbp.call(cmd, shell=True)

    if not hide:
        pyplot.show()

    return


def main(argv):
    """ @TODO main docstring
    """

    parser = argparse.ArgumentParser()

    # Image filename
    parser.add_argument("--summary_filename", type=str, default=default_summary_filename)
    parser.add_argument("--model_filename", type=str, default=default_model_filename)

    parser.add_argument("--plot_name", type=str, default=default_plot_name)
    parser.add_argument("--paper_location", type=str, default=default_paper_location)
    parser.add_argument("--file_type", type=str, default=default_file_type)

    parser.add_argument("--hide", action="store_true")

    args = vars(parser.parse_args())

    make_model_focus_v_focus_v_chip_plot(**args)

    return


if __name__ == "__main__":
    main(sys.argv)
