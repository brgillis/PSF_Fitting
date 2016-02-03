#!/usr/bin/env python

""" @file plot_stacks.py

    Created 16 Jul 2015

    This script is used to plot star/PSF/residual stacks which are
    generated as a result of running PSFit.

    ---------------------------------------------------------------------

    Copyright (C) 2015  Bryan R. Gillis

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
from os.path import join
import subprocess as sbp
import sys

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

from astropy.io import fits

default_image_location = "/disk2/brg/Data/HST_Fields"
default_image_name = "jb5d04qiq_sci1_cor"

default_plot_name_tail = "_stacks"
default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"
default_file_type = "png"

default_colormap = "bw"

default_i_min = 0.005
default_i_max = 8

default_res_i_min = -0.5
default_res_i_max = 0.5

figsize = (4,4)
labelsize = 8

def make_stacks(image_location = default_image_location,
                image_name = default_image_name,
                
                plot_name_tail = default_plot_name_tail,
                paper_location = default_paper_location,
                file_type = default_file_type,
                
                colormap = default_colormap,
                reverse_colormap = False,
                
                i_min = default_i_min,
                i_max = default_i_max,
                
                res_i_min = default_res_i_min,
                res_i_max = default_res_i_max,
                
                hide = False,
                ):
    
    stack_filename_base = join(image_location,image_name)
    
    stacks = {}
    
    for stack_name in ("star", "model", "noisy_model", "residual"):
        
        stack_filename = stack_filename_base + "_" + stack_name + "_stack.fits"
        
        stack = fits.open(stack_filename)[0]
    
        stacks[stack_name] = np.flipud(stack.data)
            
    # Do the plotting now
    
    fig = pyplot.figure(figsize=figsize)
    gs = matplotlib.gridspec.GridSpec(2, 2)
    gs.update(wspace=0.0, hspace=0, left=0.1, right=0.9, bottom=0.1, top=0.9) 
    
    for stack_name, label, i in (("star","Star",0),
                          ("model","Model",1),
                          ("noisy_model","Noisy model",2)):
        ax = pyplot.subplot(gs[i])
        clipped_data = np.clip(stacks[stack_name],i_min,i_max)
        ax.matshow(np.log(clipped_data),
                   cmap=pyplot.cm.gray,
                   interpolation='nearest',
                   vmin=np.log(i_min),
                   vmax=np.log(i_max),
                   )
        fig.gca().xaxis.set_major_locator(pyplot.NullLocator())
        fig.gca().yaxis.set_major_locator(pyplot.NullLocator())
        
        
        if(i<2):
            ax2 = ax.twiny()
            ax2.set_xlabel(label)
            fig.gca().xaxis.set_major_locator(pyplot.NullLocator())
            fig.gca().yaxis.set_major_locator(pyplot.NullLocator())
        else:
            ax.set_xlabel(label)

        
    ax = pyplot.subplot(gs[3])
    clipped_data = np.clip(stacks["residual"],res_i_min,res_i_max)
    ax.matshow(clipped_data,
               cmap=pyplot.cm.gray,
               interpolation='nearest',
               vmin=res_i_min,
               vmax=res_i_max,
               )
    fig.gca().xaxis.set_major_locator(pyplot.NullLocator())
    fig.gca().yaxis.set_major_locator(pyplot.NullLocator())
        
    ax.set_xlabel("Residual")
        
    # Save the figure
    outfile_name = image_name + plot_name_tail + "." + file_type
    pyplot.savefig(outfile_name, format=file_type, bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location if in eps format
    if file_type=="eps":
        cmd = "cp " + outfile_name + " " + paper_location
        sbp.call(cmd,shell=True)
    
    if not hide:
        fig.show()
    
    return

def main(argv):
    """ @TODO main docstring
    """
    
    parser = argparse.ArgumentParser()
    
    # Image filename
    parser.add_argument("--image_location",type=str, default=default_image_location)
    parser.add_argument("--image_name",type=str, default=[], action='append')
    
    parser.add_argument("--plot_name_tail",type=str, default=default_plot_name_tail)
    parser.add_argument("--paper_location",type=str, default=default_paper_location)
    parser.add_argument("--file_type",type=str, default=default_file_type) 
    
    parser.add_argument("--colormap",type=str, default=default_colormap) 
    parser.add_argument("--reverse_colormap", action="store_true")
    
    parser.add_argument("--hide", action="store_true")
    
    args = vars(parser.parse_args())
    
    images = args['image_name']
    
    if(len(images)==0):
        images.append(default_image_name)
    
    for image in images:
        args['image_name'] = image
        make_stacks(**args)
    
    if not args['hide']:
        pyplot.show()
    
    return

if __name__ == "__main__":
    main(sys.argv)
