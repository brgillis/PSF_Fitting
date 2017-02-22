#!/usr/bin/env python

""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_fit_plotting/plot_bestfit_params_v_obs_time.py

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
from scipy.stats import linregress

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

from astropy.io import fits

default_summary_filename = "/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_all_params_summary.fits"

default_plot_name = "bestfit_param_v_focus"
default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"
default_file_type = "png"

default_X2_max = 0.001

default_nx = 5
default_ny = 5

figsize = (12,12)

base_fontsize = 10
base_tick_fontsize = 8

param_colnames = (("Z 2",0.,"$Z_{2}$ (Tip)","z2"),
                  ("Z 3",0.,"$Z_{3}$ (Tilt)","z3"),
                  ("0 degree astigmatism",0.031,"$Z_{5}$ (Oblique Astigmatism)","astigmatism_0"),
                  ("45 degree astigmatism",0.028,"$Z_{6}$ (Vertical Astigmatism)","astigmatism_45"),
                  ("X coma",0.003,"$Z_{7}$ (X coma)","coma_x"),
                  ("Y coma",0.001,"$Z_{8}$ (Y coma)","coma_y"),
                  ("X clover",0.008,"$Z_{9}$ (X clover)","clover_x"),
                  ("Y clover",0.018,"$Z_{10}$ (Y clover)","clover_y"),
                  ("Spherical 3rd",-0.025,"$Z_{11}$ (Primary Spherical)","spherical_3rd"),
                  ("Z 12",0.,"$Z_{12}$","z12"),
                  ("Z 13",0.,"$Z_{13}$","z13"),
                  ("Z 14",0.,"$Z_{14}$","z14"),
                  ("Z 15",0.,"$Z_{15}$","z15"),
                  ("Z 16",0.,"$Z_{16}$","z16"),
                  ("Z 17",0.,"$Z_{17}$","z17"),
                  ("Z 18",0.,"$Z_{18}$","z18"),
                  ("Z 19",0.,"$Z_{19}$","z19"),
                  ("Z 20",0.,"$Z_{20}$","z20"),
                  ("Z 21",0.,"$Z_{21}$","z21"),
                  ("Spherical 5th",0.009,"$Z_{22}$ (Secondary Spherical)","spherical_5th"),
                  ("Kernel adjustment",1.0,"Kernel Adjustment","kernel_adjustment")
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
        
    gs.update(wspace=0.05, hspace=0.4, left=0.1, right=0.9, bottom=0.1, top=0.9)
    
    i = 0
    bestfit_params_string = ""
    
    for param_name, default_val, col_name, param_key in param_colnames:
        
        ax = pyplot.subplot(gs[i])
        i += 1
        
        vals = summary_table[param_name]
        
        # Draw the plot
        
        ax.scatter(summary_table["focus"],vals,edgecolor='none',s=1.0)
        
        if param_name != "Kernel adjustment":
            # ax.set_ylim([-0.05,0.05])
            # ax.set_yticks(np.linspace(-0.04,0.04,5,True))
            # if i % nx == 1:
            #     ax.set_yticklabels(["-0.04","-0.02","0","0.02","0.04"],fontsize=tick_fontsize)
            # else:
            ax.set_yticklabels([])
        else:
            ax.set_ylim([.98,1.02])
            ax.set_yticks([0.98,0.99,1.00,1.01,1.02])
            if i % nx != 1:
                ax.set_yticklabels([])
                
        ax.set_xlim([-3,1])
        
        # Draw the default
        ax.plot(ax.get_xlim(),[default_val,default_val],color='k')
        
        # Get the linear regression and plot it
        regression = linregress(summary_table["focus"],vals)
        ax.plot(ax.get_xlim(),[regression[1]+regression[0]*ax.get_xlim()[0],regression[1]+regression[0]*ax.get_xlim()[1]],
                linestyle='dashed')
        
        ax.set_title(param_name,fontsize=fontsize)
            
        # Get and print mean and sigma for the intercept
        mean = np.mean(vals)
        sigma = np.std(vals)
        stderr = sigma / np.sqrt(len(vals)-1)
        
        Z = np.abs((mean-default_val)/stderr)
        
        p_int = erfc(Z/sqrt(2.))
            
        if p_int < 0.05/len(vals):
            p_label = ("$\\mathbf{ p_{\\rm m} =  %1.1e }}$" %  p_int).replace("e","\\times 10^{")
        else:
            p_label = ("$ p_{\\rm m} =  %1.1e }$" %  p_int).replace("e","\\times 10^{")
        ax.text(0.05,0.95,p_label,horizontalalignment='left',
                verticalalignment='top',transform=ax.transAxes,
                fontsize=fontsize)
        
        # Get and print p for the slope
        p_s = regression[3]
        if p_s < 0.05/len(vals):
            p_label = ("$\\mathbf{ p_{\\rm s} =  %1.1e }}$" %  p_s).replace("e","\\times 10^{")
        else:
            p_label = ("$ p_{\\rm s} =  %1.1e }$" %  p_s).replace("e","\\times 10^{")
        ax.text(0.05,0.85,p_label,horizontalalignment='left',
                verticalalignment='top',transform=ax.transAxes,
                fontsize=fontsize)
        
        # Append to the bestfit params string
        bestfit_params_string += "--" + param_key + " " + str(regression[1]) + " "
        bestfit_params_string += "--" + param_key + "_slope " + str(regression[0]) + " "
        
    if not hide:
        fig.show()
        
    # Save the figure
    outfile_name = plot_name + "." + file_type
    pyplot.savefig(outfile_name, format=file_type, bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location if in eps format
    if file_type=="png":
        cmd = "cp " + outfile_name + " " + paper_location
        sbp.call(cmd,shell=True)
        
    print(bestfit_params_string)
    
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
    
    parser.add_argument("--nx",type=int, default=default_nx) 
    parser.add_argument("--ny",type=int, default=default_ny) 
    
    parser.add_argument("--hide", action="store_true")
    
    args = vars(parser.parse_args())
    
    make_bestfit_param_plots(**args)
    
    return

if __name__ == "__main__":
    main(sys.argv)
