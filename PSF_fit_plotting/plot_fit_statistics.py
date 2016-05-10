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

default_plot_name = "fitting_parameter_hists"
default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"
default_file_type = "png"

default_red_X2_min = 0.1
default_red_X2_max = 100.0
default_nbins = 20

figsize = (8,8)

base_fontsize = 12
base_tick_fontsize = 8

def make_fit_statistic_plots(summary_filename = default_summary_filename,
                             
                            plot_name = default_plot_name,
                            paper_location = default_paper_location,
                            file_type = default_file_type,
                            
                            red_X2_min = default_red_X2_min,
                            red_X2_max = default_red_X2_max,
                            nbins = default_nbins,
                            
                            fade_size = False,
                            
                            plot_chi2 = False,
                            
                            hide = False,
                            ):
    
    summary_table = fits.open(summary_filename)[1].data
    
    fig = pyplot.figure(figsize=figsize)
    
    if plot_chi2:
        dofs = summary_table["chi2_dofs"]+3
        plot_tuples = (("chi_squared", r"\chi^2", 'r', 0),)
        gs = matplotlib.gridspec.GridSpec(1, 1)
        ymax = 120
        yticks = [0,20,40,60,80,100,120]
        fontsize = 2*base_fontsize
        tick_fontsize = 2*base_tick_fontsize
        bins = np.logspace(np.log10(red_X2_min),np.log10(red_X2_max),
                                            nbins+1,base=10)
        xlog=True
    else:
        dofs = summary_table["X2_dofs"]
        gs = matplotlib.gridspec.GridSpec(3, 3)
        plot_tuples = (("X_squared", r"X^2", 'r', 0),
                        ("Qx_diff_Z2", r"Z^2(Q_{x}^{\perp})", 'y', 1),
                        ("Qy_diff_Z2", r"Z^2(Q_{y}^{\perp})", 'y', 2),
                        ("Qplus_sum_Z2", r"Z^2(Q_{\rm +}^{\parallel})", 'y', 3),
                        ("Qcross_sum_Z2", r"Z^2(Q_{\times}^{\parallel})", 'y', 4),
                        ("Qsize_sum_Z2", r"Z^2(Q_{\rm s}^{\parallel})", 'y', 5),
                        ("Qplus_diff_Z2", r"Z^2(Q_{\rm +}^{\perp})", 'y', 6),
                        ("Qcross_diff_Z2", r"Z^2(Q_{\times}^{\perp})", 'y', 7),
                        ("Qsize_diff_Z2", r"Z^2(Q_{\rm s}^{\perp})", 'y', 8),)
        ymax = 450
        yticks = [0,100,200,300,400]
        fontsize = base_fontsize
        tick_fontsize = base_tick_fontsize
        bins = np.logspace(np.log10(red_X2_min),np.log10(red_X2_max),
                                            nbins+1,base=10)
        xlog=True
        
    gs.update(wspace=0.3, hspace=0.3, left=0.1, right=0.9, bottom=0.1, top=0.9)
        
    
    for param_name, label, color, i in plot_tuples:
        
        ax = pyplot.subplot(gs[i])
        
        vals = summary_table[param_name]/dofs
        
        if fade_size and "size" in param_name:
            color = "w"
        
        pyplot.hist(vals, bins=bins, facecolor=color)
        
        ax.set_ylim([0,ymax])
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks,fontsize=tick_fontsize)
        
        ax.set_xticklabels(ax.get_xticklabels(),fontsize=tick_fontsize)
        
        ax.set_xlabel(r"$" + label + r"/\nu$",fontsize=fontsize)
        ax.set_ylabel("\# of images",fontsize=fontsize,labelpad=-1)
        
        if xlog:
            ax.set_xscale("log", nonposx='clip')
            
        # Get and print mean and sigma for this value
        mean = np.mean(vals)
        sigma = np.std(vals)
        
        print("For parameter " + label + "/nu:\n" +
              "mean = " + str(mean) + "\n" +
              "sigma = " + str(sigma))
        
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
    
    parser.add_argument("--plot_chi2", action="store_true")
    
    parser.add_argument("--hide", action="store_true")
    
    args = vars(parser.parse_args())
    
    make_fit_statistic_plots(**args)
    
    return

if __name__ == "__main__":
    main(sys.argv)
