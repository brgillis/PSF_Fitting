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
from math import erfc, sqrt

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

from astropy.io import fits

default_summary_filename = "/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_all_params_summary.fits"

default_plot_name = "bestfit_param_hists"
default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"
default_file_type = "eps"

default_X2_max = 0.001

default_nx = 4
default_ny = 6

figsize = (8,12)

base_fontsize = 12
base_tick_fontsize = 8

param_colnames = (("Z 2",0.),
                  ("Z 3",0.),
                  ("0 degree astigmatism",0.031),
                  ("45 degree astigmatism",0.028),
                  ("X coma",0.003),
                  # ("Y coma",0.001),
                  ("X clover",0.008),
                  ("Y clover",0.018),
                  ("Spherical 3rd",-0.025),
                  ("Z 12",0.),
                  ("Z 13",0.),
                  ("Z 14",0.),
                  ("Z 15",0.),
                  ("Z 16",0.),
                  ("Z 17",0.),
                  ("Z 18",0.),
                  ("Z 19",0.),
                  ("Z 20",0.),
                  ("Z 21",0.),
                  ("Spherical 5th",0.009),
                  ("Kernel adjustment",1.0)
                  )

def make_bestfit_param_plots(summary_filename = default_summary_filename,
                             
                            plot_name = default_plot_name,
                            paper_location = default_paper_location,
                            file_type = default_file_type,
                            
                            X2_max = default_X2_max,
                            
                            nx = default_nx,
                            ny = default_ny,
                            
                            hide = False,
                            ):
    
    summary_table = fits.open(summary_filename)[1].data
    
    fig = pyplot.figure(figsize=figsize)
    
    gs = matplotlib.gridspec.GridSpec(ny, nx)
    fontsize = base_fontsize
    tick_fontsize = base_tick_fontsize
    xlog=True
        
    gs.update(wspace=0.05, hspace=0.6, left=0.1, right=0.9, bottom=0.1, top=0.9)
    
    i = 0
    
    for param_name, default_val in param_colnames:
        
        ax = pyplot.subplot(gs[i])
        i += 1
        
        vals = summary_table[param_name]
        
        # Draw the histogram
        pyplot.hist(vals,facecolor='y')
        
        ax.set_ylim([0,ax.get_ylim()[1]*1.1])
        ax.set_yticklabels([])
        
        # Draw the default
        ax.plot([default_val,default_val],ax.get_ylim(),color='k',linestyle='dashed')
        
        ax.set_xlabel(param_name,fontsize=fontsize)
            
        # Get and print mean and sigma for this value
        mean = np.mean(vals)
        sigma = np.std(vals)
        stderr = sigma / np.sqrt(len(vals)-1)
        
        Z = np.abs((mean-default_val)/stderr)
        
        p = erfc(Z/sqrt(2.))
            
        if p < 0.05/len(vals):
            ax.text(0.05,0.95,("$\\mathbf{ p =  %1.1e }}$" %  p).replace("e","\\times 10^{"),horizontalalignment='left',
                    verticalalignment='top',transform=ax.transAxes,
                    fontsize=fontsize)
        else:
            ax.text(0.05,0.95,("$ p =  %1.1e }$" %  p).replace("e","\\times 10^{"),horizontalalignment='left',
                    verticalalignment='top',transform=ax.transAxes,
                    fontsize=fontsize)
        
    if not hide:
        fig.show()
        
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
    
    parser.add_argument("--X2_max",type=float, default=default_X2_max) 
    
    parser.add_argument("--nx",type=int, default=default_nx) 
    parser.add_argument("--ny",type=int, default=default_ny) 
    
    parser.add_argument("--hide", action="store_true")
    
    args = vars(parser.parse_args())
    
    make_bestfit_param_plots(**args)
    
    return

if __name__ == "__main__":
    main(sys.argv)
