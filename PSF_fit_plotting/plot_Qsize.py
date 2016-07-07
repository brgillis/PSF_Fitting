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

import argparse
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import os
import subprocess as sbp
import sys
from scipy.interpolate import InterpolatedUnivariateSpline

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

from psf_testing import magic_values as mv
from psf_testing.moments.get_Qs import get_moments, get_Qsize
from psf_testing.moments.coords import get_x_and_y_of_array

default_plot_name = "Qsize_comparison"
default_noisy_plot_name = "Qsize_noise_comparison"

default_file_type = "png"
    
default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"

default_nx = default_ny = 41
default_max_image_sigma = 5.0
default_num_test_points = 100
n = 4
k = 7.669

default_gain = 2.0
default_flux = 10000.0
default_background_noise = 10.0
default_bg_level = 0

default_seed = 151357

colors = ("#600060","#00C0C0","#FFFF40","k")
point_size = 4

fontsize = 24
tick_fontsize = 16

Qsize_scale = {"core":5, "wings":3}

star_Qsizes_pix = {"core":0.055599, "wings":0.097552}
model_Qsizes_pix = {"core":0.048102, "wings":0.094873}

def make_Qsize_plots(plot_name = default_plot_name,
                     noisy_plot_name = default_noisy_plot_name,
                     
                     file_type = default_file_type,
                     
                     paper_location = default_paper_location,
                     
                     nx = default_nx,
                     ny = default_ny,
                     max_image_sigma = default_max_image_sigma,
                     num_test_points = default_num_test_points,
                     fit_centre = False,
                     
                     gain = default_gain,
                     flux = default_flux,
                     background_noise = default_background_noise,
                     bg_level = default_bg_level,
                     
                     seed = default_seed,
                     
                     hide = False,
                     ):
    """ @TODO make_Qsize_plots docstring
    """
    
    if(fit_centre):
        xc = None
        yc = None
    else:
        xc = (nx-1)/2.
        yc = (ny-1)/2.
    
    for wf_i, wf_name, wf_scale, wf_scale_name in ((0,"core",mv.default_weight_sigma,r"\sigma_{\mathrm w}"),
                                                   (1,"wings",mv.default_weight_rmax,r"r_{\mathrm max}")):
    
        test_sigmas = np.linspace(max_image_sigma/num_test_points, max_image_sigma, num_test_points)
        
        Qsizes = np.zeros_like(test_sigmas)
        quad_moments = np.zeros_like(test_sigmas)
        quad_moment_with_dets = np.zeros_like(test_sigmas)
        noisy_Qsizes = np.zeros_like(test_sigmas)
        noisy_quad_moments = np.zeros_like(test_sigmas)
        noisy_quad_moment_with_dets = np.zeros_like(test_sigmas)
        
        np.random.seed(seed)
            
        x_array, y_array = get_x_and_y_of_array(nx, ny, xc, yc)
        
        x2_array = np.square(x_array)
        y2_array = np.square(y_array)
        
        r_array = np.sqrt(x2_array+y2_array)
        r_array = np.clip(r_array,0.25,nx+ny)
        
        for i in range(len(test_sigmas)):
            
            # Make a Gaussian image with the proper sigma
            sigma = test_sigmas[i]
            
            unnormed_image = np.exp(-k*np.power(r_array/sigma,1./n))
            image = flux*unnormed_image/unnormed_image.sum() + bg_level
        
            noise_per_pixel = np.sqrt(image/gain + np.square(background_noise))
            
            noisy_image = image + noise_per_pixel * np.random.randn(nx,ny)
            
            Qsizes[i] = get_Qsize(image, xc=xc, yc=yc)[wf_i]/wf_scale/mv.pixel_scale
            noisy_Qsizes[i] = get_Qsize(noisy_image, xc=xc, yc=yc)[wf_i]/wf_scale/mv.pixel_scale
            
            moments = get_moments(image, xc=xc, yc=yc)
            noisy_moments = get_moments(noisy_image, xc=xc, yc=yc)
            
            Mxx = moments[2][0][wf_i]
            Myy = moments[2][1][wf_i]
            Mxy = moments[2][2][wf_i]
            noisy_Mxx = noisy_moments[2][0][wf_i]
            noisy_Myy = noisy_moments[2][1][wf_i]
            noisy_Mxy = noisy_moments[2][2][wf_i]
            
            quad_moments[i] = np.abs(Mxx + Myy)/wf_scale**2
            quad_moment_with_dets[i] = np.abs(Mxx + Myy + 2 * np.sqrt( Mxx*Myy - np.square(Mxy) ))/wf_scale**2
            noisy_quad_moments[i] = np.abs(noisy_Mxx + noisy_Myy)/wf_scale**2
            noisy_quad_moment_with_dets[i] = np.abs(noisy_Mxx + noisy_Myy + \
                2 * np.sqrt( np.abs(noisy_Mxx*noisy_Myy - np.square(noisy_Mxy) ) ) ) /wf_scale**2
        
        normed_sigmas = test_sigmas/wf_scale
        
        # Calculate some equivalent values needed for the paper
        
        star_Qsize = star_Qsizes_pix[wf_name]/wf_scale/mv.pixel_scale
        model_Qsize = model_Qsizes_pix[wf_name]/wf_scale/mv.pixel_scale
        
        quad_moments_spline = InterpolatedUnivariateSpline(Qsizes,quad_moments)
        
        star_quad_moments = quad_moments_spline(star_Qsize)
        model_quad_moments = quad_moments_spline(model_Qsize)
        
        print("Equivalent moment size for Qs of " + str(star_Qsize*wf_scale*mv.pixel_scale)
              + " = " + str(star_quad_moments*(wf_scale*mv.pixel_scale)**2))
        print("Equivalent moment size for Qs of " + str(model_Qsize*wf_scale*mv.pixel_scale)
              + " = " + str(model_quad_moments*(wf_scale*mv.pixel_scale)**2))
        
        quad_moment_with_dets_spline = InterpolatedUnivariateSpline(Qsizes,quad_moment_with_dets)
        
        star_quad_moment_with_dets = quad_moment_with_dets_spline(star_Qsize)
        model_quad_moment_with_dets = quad_moment_with_dets_spline(model_Qsize)
        
        print("Equivalent moment with det size for Qs of " + str(star_Qsize*wf_scale*mv.pixel_scale)
              + " = " + str(star_quad_moment_with_dets*(wf_scale*mv.pixel_scale)**2))
        print("Equivalent moment with det size for Qs of " + str(model_Qsize*wf_scale*mv.pixel_scale)
              + " = " + str(model_quad_moment_with_dets*(wf_scale*mv.pixel_scale)**2))
            
        # Do the plotting now
        
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(r"$r_{\mathrm{s}}/" + wf_scale_name + r"$",labelpad=10,fontsize=fontsize)
        ax.set_ylabel(r"$(r_{\mathrm{" + wf_name + r"}}/" + wf_scale_name + ")^2$",
                      labelpad=10,fontsize=fontsize)
        
        ax.plot(normed_sigmas,
                (Qsize_scale[wf_name]*Qsizes)**2,
                label=r"$(" + str(Qsize_scale[wf_name]) + r"Q_{\rm s})^2$",
                color=colors[3],
                linewidth=2)
        ax.plot(normed_sigmas,
                (Qsize_scale[wf_name]*noisy_Qsizes)**2,
                label=r"$(" + str(Qsize_scale[wf_name]) + r"Q_{\rm s})^2$ + noise",
                color=colors[3],
                linewidth=1)
        ax.plot(normed_sigmas,
                quad_moment_with_dets,
                label="Moments with det",
                color=colors[1],
                linewidth=2)
        ax.plot(normed_sigmas,
                noisy_quad_moment_with_dets,
                label="Moments with det + noise",
                color=colors[1],
                linewidth=1)
        ax.plot(normed_sigmas,
                quad_moments,
                label="Moments w/o det",
                color=colors[0],
                linewidth=2)
        ax.plot(normed_sigmas,
                noisy_quad_moments,
                label="Moments w/o det + noise",
                color=colors[0],
                linewidth=1)
        
        ax.legend(loc='upper left')
        
        ax.set_xlim(np.min(normed_sigmas),np.max(normed_sigmas))
        ax.set_ylim(0,1.1*np.max(noisy_quad_moment_with_dets))
        
        ax.set_xticklabels(ax.get_xticks(),fontsize=tick_fontsize)
        ax.set_yticklabels(ax.get_yticks(),fontsize=tick_fontsize)
        
        # Save the figure
        outfile_name = os.path.splitext(plot_name)[0] + "_" + wf_name + "." + file_type
        pyplot.savefig(outfile_name, format=file_type, bbox_inches="tight", pad_inches=0.05)
        
        # Copy it to the paper location if in eps format
        if file_type=="eps":
            cmd = "cp " + outfile_name + " " + paper_location
            sbp.call(cmd,shell=True)
            
        if not hide:
            fig.show()
            
        # Check for a possible bias      
        if not hide:
            Qsize_noise_ratio = np.mean(noisy_Qsizes/Qsizes)
            print("Mean bias ratio: " + str(Qsize_noise_ratio))
            
        # Now olot the relation between the estimators
            
        # Do the plotting now for the slopes
        
        fig = pyplot.figure()
        fig.subplots_adjust(wspace=0.5, hspace=0, bottom=0.1, right=0.95, top=0.95, left=0.12)
        
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel(r"$Q_{\mathrm s}/" + wf_scale_name + "$",labelpad=10,fontsize=fontsize)
        ax.set_ylabel(r"$(r_{\mathrm{" + wf_name + r"}}/" + wf_scale_name + ")^2$",
                      labelpad=10,fontsize=fontsize)
        
        normed_sigmas = test_sigmas/wf_scale
        
        ax.plot(Qsizes,
                quad_moment_with_dets,
                label="Moments with det",
                color=colors[1],
                linewidth=2)
        ax.plot(Qsizes,
                quad_moments,
                label="Moments w/o det",
                color=colors[0],
                linewidth=2)
        
        ax.legend(loc='upper left')
        
        ax.set_xlim(np.min(Qsizes),np.max(Qsizes))
        ax.set_ylim(np.min(quad_moments),
                    1.1*np.max(quad_moment_with_dets))
        
        ax.set_xticklabels(ax.get_xticks(),fontsize=tick_fontsize)
        ax.set_yticklabels(ax.get_yticks(),fontsize=tick_fontsize)
        
        # Save the figure
        outfile_name = os.path.splitext(plot_name)[0] + "_" + wf_name + "_comps." + file_type
        pyplot.savefig(outfile_name, format=file_type, bbox_inches="tight", pad_inches=0.05)
        
        # Copy it to the paper location if in eps format
        if file_type=="eps":
            cmd = "cp " + outfile_name + " " + paper_location
            sbp.call(cmd,shell=True)
            
        if not hide:
            fig.show()
        
    if not hide:
        pyplot.show()
    
def main(argv):
    
    parser = argparse.ArgumentParser()
    
    # Image filename
    parser.add_argument("--plot_name",type=str, default=default_plot_name)
    parser.add_argument("--noisy_plot_name",type=str, default=default_noisy_plot_name)
    
    parser.add_argument("--file_type",type=str, default=default_file_type)
    
    parser.add_argument("--paper_location",type=str, default=default_paper_location)
    
    parser.add_argument("--nx",type=int, default=default_nx)
    parser.add_argument("--ny",type=int, default=default_ny)
    parser.add_argument("--max_image_sigma",type=float, default=default_max_image_sigma)
    parser.add_argument("--num_test_points",type=int, default=default_num_test_points)
    parser.add_argument("--fit_centre", action="store_true")
    
    parser.add_argument("--gain",type=float, default=default_gain)
    parser.add_argument("--flux",type=float, default=default_flux)
    parser.add_argument("--background_noise",type=float, default=default_background_noise)
    parser.add_argument("--bg_level",type=float, default=default_bg_level)
    
    parser.add_argument("--seed",type=int, default=default_seed)
    
    parser.add_argument("--hide", action="store_true")
    
    args = vars(parser.parse_args())
    
    return make_Qsize_plots(**args)

if __name__ == "__main__":
    main(sys.argv)
