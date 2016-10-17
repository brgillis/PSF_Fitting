#!/usr/bin/env python

""" @file plot_grid_convergence.py

    Created 27 Sep 2016

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
import subprocess as sbp
from argparse import ArgumentParser

from astropy.io import fits
import matplotlib

import matplotlib.pyplot as pyplot
import numpy as np


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

data_dir = "/home/brg/Data/HST_Fields/grid_convergence_testing"
default_results_file_root = "control_image_n"

grid_sizes = (2048, 1024, 512, 256, 128, 64, 32, 1)

figsize = (6, 6)
fontsize = 12

default_output_filename_root = "grid_convergence_test"
default_output_extension = "eps"

paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"

num_images = 100

def main(argv):
    """ @TODO main docstring
    """
    
    parser = ArgumentParser()
    
    parser.add_argument("--results_file_root",default=default_results_file_root)
    parser.add_argument("--output_filename_root",default=default_output_filename_root)
    parser.add_argument("--output_extension",default=default_output_extension)
    
    args = parser.parse_args()

    results_root = os.path.join(data_dir, args.results_file_root)
    
    indices = range(len(grid_sizes))
    
    labels = ("X2",
              "Qx_diff_Z2","Qy_diff_Z2",
              "Qp_sum_Z2","Qp_diff_Z2",
              "Qc_sum_Z2","Qc_diff_Z2",
              "Qs_sum_Z2","Qs_diff_Z2")
    
    vals = {}
    for label in labels:
        vals[label] = np.zeros((len(grid_sizes),num_images), dtype=float)
    
    for n in range(num_images):
        
        results_root_n = results_root.replace("_n","_"+str(n))
        
        for i in indices:
        
            filename = results_root_n + "_gp" + str(grid_sizes[i]) + "_results.fits"
        
            header = fits.open(filename)[1].header
        
            vals["X2"][i,n] = header["X_SQR"]
            vals["Qx_diff_Z2"][i,n] = header["QXD_Z2"]
            vals["Qy_diff_Z2"][i,n] = header["QYD_Z2"]
            vals["Qp_sum_Z2"][i,n] = header["QPS_Z2"] 
            vals["Qp_diff_Z2"][i,n] = header["QPD_Z2"] 
            vals["Qc_sum_Z2"][i,n] = header["QPS_Z2"] 
            vals["Qc_diff_Z2"][i,n] = header["QPD_Z2"] 
            vals["Qs_sum_Z2"][i,n] = header["QSS_Z2"] 
            vals["Qs_diff_Z2"][i,n] = header["QSD_Z2"]
            
    # Get the means
    val_means = {}
    for label in labels:
        val_means[label] = vals[label].mean(axis=1)

    # Plot it up
    _fig = pyplot.figure(figsize=figsize)

    gs = matplotlib.gridspec.GridSpec(1, 1)
    gs.update(wspace=0.3, hspace=0.3, left=0.1, right=0.9, bottom=0.1, top=0.9)

    ax = pyplot.subplot(gs[0])

    ax.plot(indices, val_means["X2"], label="$X^2$", color="black")
    ax.plot(indices, val_means["Qx_diff_Z2"], label=r"$Q_x^{(-)} Z^2$", linestyle="dotted")
    ax.plot(indices, val_means["Qy_diff_Z2"], label=r"$Q_y^{(-)} Z^2$", linestyle="dotted")
    ax.plot(indices, val_means["Qp_sum_Z2"], label=r"$Q_+^{(+)} Z^2$", linestyle="dashed")
    ax.plot(indices, val_means["Qp_diff_Z2"], label=r"$Q_+^{(-)} Z^2$", linestyle="dashed")
    ax.plot(indices, val_means["Qc_sum_Z2"], label=r"$Q_{\times}^{(+)} Z^2$", linestyle="dashed")
    ax.plot(indices, val_means["Qc_diff_Z2"], label=r"$Q_{\times}^{(-)} Z^2$", linestyle="dashed")
    ax.plot(indices, val_means["Qs_sum_Z2"], label=r"$Q_s^{(+)} Z^2$", linestyle="dashdot")
    ax.plot(indices, val_means["Qs_diff_Z2"], label=r"$Q_s^{(-)} Z^2$", linestyle="dashdot")

    ax.set_xlim([-0.5, len(indices)-0.5])
    # ax.set_ylim([0, 0.02])
    
    ax.set_yscale("log", nonposy='clip')
    
    ax.set_xlabel("Cell size",fontsize=fontsize)
    ax.set_xticks(indices)
    ax.set_xticklabels(grid_sizes)

    ax.legend(loc="lower left",ncol=2)
    
    output_filename = args.output_filename_root + "." + args.output_extension

    pyplot.savefig(output_filename, format=args.output_extension, bbox_inches="tight", pad_inches=0.05)
    
    if args.output_extension == "eps":
        cmd = "cp " + output_filename + " " + paper_location
        sbp.call(cmd,shell=True)

    pyplot.show()

if __name__ == "__main__":
    main(sys.argv)
