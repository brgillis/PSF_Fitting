#!/usr/bin/env python

""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_fit_plotting/plot_bestfit_params_FT.py

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
from scipy.stats import linregress

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

from astropy.io import fits

default_summary_filename = "/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_all_params_summary.fits"

default_plot_name = "bestfit_param_FT"
default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"
default_file_type = "png"

default_fmin = 1.0e-7
default_fmax = 1.0e-3
default_Nf = 100000
averaging = 100

default_nx = 5
default_ny = 5

figsize = (18,18)

base_fontsize = 15
base_tick_fontsize = 12

good_X2_limit = 5e-6

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

hst_period = 95.47 / 60 # in units of hours
hst_frequency = 1. / hst_period
# default_fmin = hst_frequency
# default_fmax = 6*hst_frequency
# default_Nf = 6
# averaging = 1

def moving_average(a, n) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def make_bestfit_param_plots(summary_filename = default_summary_filename,
                             
                            plot_name = default_plot_name,
                            paper_location = default_paper_location,
                            file_type = default_file_type,
                            
                            fmin = default_fmin,
                            fmax = default_fmax,
                            Nf = default_Nf,
                            
                            nx = default_nx,
                            ny = default_ny,
                            
                            hide = False,
                            ):
    
    summary_table = fits.open(summary_filename)[1].data
    
    good_X2 = summary_table["X_squared"] < good_X2_limit
    
    fig = pyplot.figure(figsize=figsize)
    
    gs = matplotlib.gridspec.GridSpec(ny, nx)
    fontsize = base_fontsize
    tick_fontsize = base_tick_fontsize
        
    gs.update(wspace=0.05, hspace=0.4, left=0.1, right=0.9, bottom=0.1, top=0.9)
    
    i = 0
    
    ws = np.linspace(2*np.pi*fmin,2*np.pi*fmax,Nf)
    cos_wave_basis = np.cos(np.outer(ws,summary_table["obs_time"][good_X2]))
    sin_wave_basis = np.sin(np.outer(ws,summary_table["obs_time"][good_X2]))
    
    ave_ws = moving_average(ws,averaging)[0::averaging]
    
    for param_name, _, param_title, _ in param_colnames:
        
        ax = pyplot.subplot(gs[i])
        i += 1
    
        vals = summary_table[param_name][good_X2]
    
        # Correct for linear relationship with focus
        regression = linregress(summary_table["focus"][good_X2],vals)
        zeroed_vals = vals - (regression[1]+regression[0]*summary_table["focus"][good_X2])
        
        cos_vals = moving_average(np.sum(cos_wave_basis*zeroed_vals,axis=1),averaging)[0::averaging]
        sin_vals = moving_average(np.sum(sin_wave_basis*zeroed_vals,axis=1),averaging)[0::averaging]
        
        amps = np.sqrt(cos_vals**2 + sin_vals**2)
        normed_amps = amps / np.mean(amps)
        
        phases = np.arctan2(sin_vals,cos_vals)
        
        ax.set_xlim(3600*fmin,3600*fmax)
        ax.set_ylim(0,5.0)
        
        # Draw multiples of HST's frequency
        for k in range(6):
            f = k*hst_frequency
            ax.plot([f,f],ax.get_ylim(),color=(0.5,0.5,0.5),linestyle='dashed',linewidth=0.5)
        
        # Draw the plot
        
        ax.scatter(3600*ave_ws/(2*np.pi),normed_amps,marker='.',color='b',s=16)
        # ax.scatter(3600*ave_ws/(2*np.pi),phases,marker='.',color='b',s=16)
        
        if i + nx > len(param_colnames):
            ax.set_xlabel(r"$f$ (Hr$^{-1}$)",fontsize=fontsize)
            
        if i % nx == 1:
            ax.set_ylabel(r"$A/\overline{A}$",fontsize=fontsize)
        else:
            ax.set_yticklabels([])
        
        # Draw zero
        ax.plot(ax.get_xlim(),[0.,0.],color='k',linestyle='solid')
        
        ax.set_title(param_title,fontsize=fontsize)
        
        
    if not hide:
        fig.show()
        
    # Save the figure
    outfile_name = plot_name + "." + file_type
    pyplot.savefig(outfile_name, format=file_type, bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location if in eps format
    if file_type=="png":
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
    
    parser.add_argument("--nx",type=int, default=default_nx) 
    parser.add_argument("--ny",type=int, default=default_ny) 
    
    parser.add_argument("--hide", action="store_true")
    
    args = vars(parser.parse_args())
    
    make_bestfit_param_plots(**args)
    
    return

if __name__ == "__main__":
    main(sys.argv)
