#!/usr/bin/env python

""" @file plot_Qsize.py

    Created 22 Jan 2016

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

import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import os
from scipy import signal
import subprocess as sbp
import sys

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

from psf_testing import magic_values as mv
from psf_testing.moments.get_Qs import get_moments, get_Qsize

plot_name = "Qsize_comparison"
noisy_plot_name = "Qsize_noise_comparison"

file_type = "png"
    
paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"

nx = ny = 41
max_image_sigma = 5.0
num_test_points = 50

gain = 2.0
flux = 10000.0
background_noise = 10.0
bg_level = 1

test_shear = 0.2

seed = 151357

def main(argv):
    """ @TODO main docstring
    """
    
    test_sigmas = np.linspace(max_image_sigma/num_test_points, max_image_sigma, num_test_points)
    
    Qsizes = np.zeros_like(test_sigmas)
    quad_moments = np.zeros_like(test_sigmas)
    quad_moment_with_dets = np.zeros_like(test_sigmas)
    noisy_Qsizes = np.zeros_like(test_sigmas)
    noisy_quad_moments = np.zeros_like(test_sigmas)
    noisy_quad_moment_with_dets = np.zeros_like(test_sigmas)
    
    np.random.seed(seed)
    
    for i in range(len(test_sigmas)):
        
        # Make a Gaussian image with the proper sigma
        sigma = test_sigmas[i]
        unnormed_image = np.outer(signal.gaussian(nx, sigma), signal.gaussian(ny, sigma))
        image = flux*unnormed_image/unnormed_image.sum()
    
        noise_per_pixel = np.sqrt(image/gain + np.square(background_noise))
        
        noisy_image = image + noise_per_pixel * np.random.randn(nx,ny)
        
        Qsizes[i] = get_Qsize(image, xc=(nx-1.)/2, yc=(nx-1.)/2)[1]
        noisy_Qsizes[i] = get_Qsize(noisy_image, xc=(nx-1.)/2, yc=(nx-1.)/2)[1]
        
        moments = get_moments(image, xc=(nx-1.)/2, yc=(nx-1.)/2)
        noisy_moments = get_moments(noisy_image, xc=(nx-1.)/2, yc=(nx-1.)/2)
        
        Mxx = moments[2][0][1]
        Myy = moments[2][1][1]
        Mxy = moments[2][2][1]
        noisy_Mxx = noisy_moments[2][0][1]
        noisy_Myy = noisy_moments[2][1][1]
        noisy_Mxy = noisy_moments[2][2][1]
        
        quad_moments[i] = np.sqrt(np.abs(Mxx + Myy))
        quad_moment_with_dets[i] = np.sqrt(np.abs(Mxx + Myy + 2 * np.sqrt( Mxx*Myy - np.square(Mxy) )))
        noisy_quad_moments[i] = np.sqrt(np.abs(noisy_Mxx + noisy_Myy))
        noisy_quad_moment_with_dets[i] = np.sqrt(np.abs(noisy_Mxx + noisy_Myy + \
            2 * np.sqrt( np.abs(noisy_Mxx*noisy_Myy - np.square(noisy_Mxy) ) ) ) )
        
    # Do the plotting now
    
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r"$\sigma/r_{\mathrm max}$$",labelpad=10,fontsize=16)
    ax.set_ylabel("$r/r_{\mathrm{max}}$",labelpad=10,fontsize=16)
    
    ax.plot(test_sigmas/mv.default_weight_rmax,
            quad_moments/mv.default_weight_rmax,
            label="Quad. moments w/o det")
    ax.plot(test_sigmas/mv.default_weight_rmax,
            quad_moment_with_dets/mv.default_weight_rmax,
            label="Quad. moments with det")
    ax.plot(test_sigmas/mv.default_weight_rmax,
            Qsizes/mv.pixel_scale/mv.default_weight_rmax,
            label=r"$Q_{\rm size}$")
    
    ax.legend(loc='upper left')
    
    # Save the figure
    outfile_name = os.path.splitext(plot_name)[0] + "." + file_type
    pyplot.savefig(outfile_name, format=file_type, bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location if in eps format
    if file_type=="eps":
        cmd = "cp " + outfile_name + " " + paper_location
        sbp.call(cmd,shell=True)
    
    fig.show()
    
    fig = pyplot.figure()
    fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r"$\sigma/r_{\mathrm max}$$",labelpad=10,fontsize=16)
    ax.set_ylabel("$r/r_{\mathrm{max}}$",labelpad=10,fontsize=16)
    
    ax.plot(test_sigmas/mv.default_weight_rmax,
            noisy_quad_moments/mv.default_weight_rmax,
            label="Quad. moments w/o det")
    ax.plot(test_sigmas/mv.default_weight_rmax,
            noisy_quad_moment_with_dets/mv.default_weight_rmax,
            label="Quad. moments with det")
    ax.plot(test_sigmas/mv.default_weight_rmax,
            noisy_Qsizes/mv.pixel_scale/mv.default_weight_rmax,
            label=r"$Q_{\rm size}$")
    
    ax.legend(loc='upper left')
    
    # Save the figure
    outfile_name = os.path.splitext(noisy_plot_name)[0] + "." + file_type
    pyplot.savefig(outfile_name, format=file_type, bbox_inches="tight", pad_inches=0.05)
    
    # Copy it to the paper location if in eps format
    if file_type=="eps":
        cmd = "cp " + outfile_name + " " + paper_location
        sbp.call(cmd,shell=True)
    
    pyplot.show()
    

if __name__ == "__main__":
    main(sys.argv)
