#!/usr/bin/env python

""" @file plot_focus_fitting_change.py

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
import subprocess as sbp
from argparse import ArgumentParser

from astropy.io import fits
import matplotlib
import matplotlib.pyplot as pyplot

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import numpy as np


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

data_dir = "/disk2/brg/Data/HST_Fields/"
default_results_file_root = "control_image_n"

default_output_filename_root = "focus_fitting_change"
default_output_extension = "eps"

paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"

num_images = 100

default_nx = 3
default_ny = 3

figsize = (10,12)

base_fontsize = 12
base_tick_fontsize = 8

def main(argv):
    """ @TODO main docstring
    """
    
    parser = ArgumentParser()
    
    parser.add_argument("--results_file_root",default=default_results_file_root)
    parser.add_argument("--output_filename_root",default=default_output_filename_root)
    parser.add_argument("--output_extension",default=default_output_extension)
    parser.add_argument("--nx",default=default_nx)
    parser.add_argument("--ny",default=default_ny)
    
    args = parser.parse_args()

    results_root = os.path.join(data_dir, args.results_file_root)
    
    control_types = (("Base",""),
                     ("Binaries", "_b3"),
                     ("Wide Binaries", "_b32"),
                     ("1D Guiding Error", "_blur"),
                     ("2D Guiding Error", "_blur2"),
                     ("Galaxy Background", "_gb"),
                     ("Varying Spec. Type", "_rs"),
                     ("Full", "_full"))
    
    vals = {}
    for control_type, _ in control_types:
        vals[control_type] = {"Focus": np.zeros(num_images, dtype=float),
                              "Fit Focus": np.zeros(num_images, dtype=float)}
    
    for n in range(num_images):
        
        results_root_n = results_root.replace("_n","_"+str(n))
        
        for results_tail, which_focus in (("_results.fits","Focus"),
                                          ("_fit_results.fits","Fit Focus")):
            for control_type, tag in control_types:
            
                filename = results_root_n + tag + results_tail
            
                header = fits.open(filename)[1].header
            
                vals[control_type][which_focus][n] = header["FOCUS"]
                
    # Set up the plot
    fig = pyplot.figure(figsize=figsize)
    
    gs = matplotlib.gridspec.GridSpec(args.ny, args.nx)
    fontsize = base_fontsize
    tick_fontsize = base_tick_fontsize
        
    gs.update(wspace=0.2, hspace=0.4, left=0.1, right=0.9, bottom=0.1, top=0.9)
    
    i = 0
    
    for control_type, _ in control_types:
        
        ax = pyplot.subplot(gs[i])
        i += 1
        
        pyplot.scatter(vals[control_type]["Focus"],vals[control_type]["Fit Focus"],color='r')
        pyplot.plot([-10,10],[-10,10],color='k',linestyle='dashed')
        
        ax.set_xlim([-6.5,4.5])
        ax.set_ylim([-8.5,6.5])
        
        ax.set_title(control_type)
        ax.set_xlabel("Actual Focus Offset",fontsize=fontsize)
        if i % args.nx == 1:
            ax.set_ylabel("Fit Focus Offset",fontsize=fontsize)
        
    # Save the figure
    outfile_name = args.output_filename_root + "." + args.output_extension
    pyplot.savefig(outfile_name, format=args.output_extension, bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location if in eps format
    if args.output_extension=="eps":
        cmd = "cp " + outfile_name + " " + paper_location
        sbp.call(cmd,shell=True)
    
    pyplot.show()
     
    return

if __name__ == "__main__":
    main(sys.argv)
