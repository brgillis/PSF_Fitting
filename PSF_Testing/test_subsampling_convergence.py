#!/usr/bin/env python

""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/test_subsampling_convergence.py

    Created 22 Sep 2016

    Runs a series of focus fittings to see how subsampling converges.

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

import os
import sys
import numpy as np

from psf_testing.arg_parser import get_arg_parser, parse_and_get_kwargs
from psf_testing.test_psf import test_psf
from utility.smart_logging import get_default_logger

class test_psf_caller(object):
    def __init__(self,image_filename,extra_tag,*args,**kwargs):
        self.image_filename = image_filename
        self.extra_tag = extra_tag
        self.args = args
        self.kwargs = kwargs
    def __call__(self,x):
        try:
            test_psf(self.image_filename,subsampling_factor=x,results_tag=self.extra_tag+"ss"+str(x),
                     *self.args,**self.kwargs)
        except Exception as e:
            logger = get_default_logger()
            logger.error("Exception when trying subsampling factor " + str(x) + ": " + str(e))
            raise

default_image_dir = "/home/brg/Data/HST_Fields/"
default_image_filename = "control_image_n.fits"
default_results_dir = "/home/brg/Data/HST_Fields/subsampling_convergence_testing"

default_num_images = 1

def main(argv):
    """ @TODO main docstring
    """
    
    subsampling_factors = np.linspace(1,10,10,endpoint=True).astype(int)
    
    # Check if we're debugging
    try:
        import pydevd as _
        debugging = True
    except ImportError:
        debugging = False
        
    # Execute command-line parsing
    parser = get_arg_parser()
    
    parser.add_argument("--num_images",default=default_num_images,type=int)
    
    kwargs, special_kwargs = parse_and_get_kwargs(parser,
                                                  special_keys=("image_filename","image_list_filename",
                                                                "image_dir","logging_level","disable_parallelization",
                                                                "subsampling_factor","results_dir","results_tag",
                                                                "num_images"))
        
    if special_kwargs["image_dir"] is None:
        image_dir = default_image_dir
    else:
        image_dir = special_kwargs["image_dir"]
        
    if special_kwargs["results_dir"] is None:
        results_dir = default_results_dir
    else:
        results_dir = special_kwargs["results_dir"]
        
    if kwargs["galsim_rebin"]:
        extra_tag = "grb_"
    else:
        extra_tag = ""
        
    if special_kwargs["num_images"] is None:
        num_images = default_num_images
    else:
        num_images = special_kwargs["num_images"]
        
    parallelize = not (special_kwargs["disable_parallelization"] or debugging)

    for i in range(num_images):
    
        if special_kwargs["image_filename"] is None:
            image_filename = default_image_filename
        else:
            image_filename = special_kwargs["image_filename"]
        if num_images != 1:
            image_filename = image_filename.replace("_n","_"+str(i))
        
        caller = test_psf_caller(os.path.join(image_dir, image_filename),
                                 extra_tag,
                                 results_dir=results_dir,
                                 parallelize=parallelize,
                                 **kwargs)
            
        for subsampling_factor in subsampling_factors:
            caller(subsampling_factor)
        
        print("Finished subsampling convergence testing on " + image_filename + " with extra results tag " + extra_tag + ".")

if __name__ == "__main__":
    main(sys.argv)
