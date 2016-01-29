#!/usr/bin/env python

""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_fit_plotting/plot_fit_statistics.py

    Created 29 Jan 2016

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
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import subprocess as sbp
import sys

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

from astropy.io import fits

default_summary_filename = "/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_summary.fits"

default_plot_name = "fitting_paramater_hists"
default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"
default_file_type = "png"

default_red_X2_min = 0.1
default_red_X2_max = 100.0
default_nbins = 20

figsize = (8,8)

def make_fit_statistic_plots(summary_filename = default_summary_filename,
                             
                            plot_name = default_plot_name,
                            paper_location = default_paper_location,
                            file_type = default_file_type,
                            
                            red_X2_min = default_red_X2_min,
                            red_X2_max = default_red_X2_max,
                            nbins = default_nbins,
                            
                            fade_size = False,
                            
                            hide = False,
                            ):
    
    summary_table = fits.open(summary_filename)[1].data
    
    fig = pyplot.figure(figsize=figsize)
    gs = matplotlib.gridspec.GridSpec(3, 3)
    gs.update(wspace=0.3, hspace=0.3, left=0.1, right=0.9, bottom=0.1, top=0.9)
    
    lmin = np.log10(red_X2_min)
    lmax = np.log10(red_X2_max)
    dofs = summary_table["X2_dofs"]
    
    for param_name, label, color, i in (("X_squared", r"X^2", 'r', 0),
                                        ("Qx_diff_Z2", r"Z^2(Q_{x}^{-})", 'y', 1),
                                        ("Qy_diff_Z2", r"Z^2(Q_{y}^{-})", 'y', 2),
                                        ("Qplus_sum_Z2", r"Z^2(Q_{\rm +}^{+})", 'y', 3),
                                        ("Qcross_sum_Z2", r"Z^2(Q_{\times}^{+})", 'y', 4),
                                        ("Qsize_sum_Z2", r"Z^2(Q_{\rm s}^{+})", 'y', 5),
                                        ("Qplus_diff_Z2", r"Z^2(Q_{\rm +}^{-})", 'y', 6),
                                        ("Qcross_diff_Z2", r"Z^2(Q_{\times}^{-})", 'y', 7),
                                        ("Qsize_diff_Z2", r"Z^2(Q_{\rm s}^{-})", 'y', 8),):
        
        ax = pyplot.subplot(gs[i])
        
        lvals = np.log10(summary_table[param_name]/dofs)
        
        if fade_size and "size" in param_name:
            color = "w"
        
        pyplot.hist(lvals, bins=np.linspace(lmin,lmax,nbins+1), facecolor=color)
        
        ax.set_ylim([0,450])
        
        ax.set_xlabel(r"$\log_{10} " + label + r"/\nu$")
        ax.set_ylabel("\# of images")
        
    # Save the figure
    outfile_name = plot_name + "." + file_type
    pyplot.savefig(outfile_name, format=file_type, bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location if in eps format
    if file_type=="eps":
        cmd = "cp " + outfile_name + " " + paper_location
        sbp.call(cmd,shell=True)
    
    if not hide:
        pyplot.show()
     
    return

def main(argv):
    """ @TODO main docstring
    """
    
    parser = argparse.ArgumentParser()
    
    # Image filename
    parser.add_argument("--summary_filename",type=str, default=default_summary_filename)
    
    parser.add_argument("--plot_name",type=str, default=default_plot_name)
    parser.add_argument("--paper_location",type=str, default=default_paper_location)
    parser.add_argument("--file_type",type=str, default=default_file_type) 
    
    parser.add_argument("--red_X2_min",type=float, default=default_red_X2_min) 
    parser.add_argument("--red_X2_max",type=float, default=default_red_X2_max) 
    parser.add_argument("--nbins",type=int, default=default_nbins) 
    
    parser.add_argument("--fade_size", action="store_true")
    
    parser.add_argument("--hide", action="store_true")
    
    args = vars(parser.parse_args())
    
    make_fit_statistic_plots(**args)
    
    return

if __name__ == "__main__":
    main(sys.argv)
