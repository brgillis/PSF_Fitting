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

default_plot_name = "X_v_focus_v_chip"
default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"
default_file_type = "png"

default_red_X2_min = 0.1
default_red_X2_max = 10.0

figsize = (6,6)

fontsize = 12

def make_X_v_focus_v_chip_plot(summary_filename = default_summary_filename,
                             
                            plot_name = default_plot_name,
                            paper_location = default_paper_location,
                            file_type = default_file_type,
                            
                            red_X2_min = default_red_X2_min,
                            red_X2_max = default_red_X2_max,
                            
                            hide = False,
                            
                            plot_chi2 = False,
                            ):
    
    summary_table = fits.open(summary_filename)[1].data
    
    _fig = pyplot.figure(figsize=figsize)
        
    gs = matplotlib.gridspec.GridSpec(1, 1)
    gs.update(wspace=0.3, hspace=0.3, left=0.1, right=0.9, bottom=0.1, top=0.9)
        
    ax = pyplot.subplot(gs[0])
    
    if plot_chi2:
        y_label = r"Best $\chi^2_{\rm red}$"
        fitted_param = summary_table["chi_squared"]/(summary_table["chi2_dofs"]+3)
    else:
        y_label = r"Best $X^2$"
        fitted_param = summary_table["X_squared"]
    ax.set_yscale("log", nonposy='clip')

    chip1_mask = ~np.logical_and(summary_table["chip"]==1,fitted_param>0)
    chip2_mask = ~np.logical_and(summary_table["chip"]==2,fitted_param>0)

    for mask, label, color, marker in ((chip1_mask, "Chip 1", '#C00000', "o"),
                               (chip2_mask, "Chip 2", '#4040FF', "^")):
        focii = np.ma.masked_array(summary_table["focus"],mask).compressed()
        X2s = np.ma.masked_array(fitted_param,mask).compressed()

        pyplot.scatter(focii,X2s,edgecolors=color,label=label,alpha=1,marker=marker,facecolors='none')
    
    ax.set_ylim([red_X2_min,red_X2_max])
    
    ax.set_xlabel("Best-fit focus offset (microns)",fontsize=fontsize)
    ax.set_ylabel(y_label,fontsize=fontsize,labelpad=-8)
    
    ax.legend(loc="upper left")
        
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
    
    parser.add_argument("--plot_chi2", action="store_true")
    
    parser.add_argument("--hide", action="store_true")
    
    args = vars(parser.parse_args())
    
    make_X_v_focus_v_chip_plot(**args)
    
    return

if __name__ == "__main__":
    main(sys.argv)
