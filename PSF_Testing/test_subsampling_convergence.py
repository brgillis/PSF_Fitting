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
from multiprocessing import Pool, cpu_count
from copy import deepcopy

from psf_testing import magic_values as mv
from psf_testing.test_psf import test_psf
from psf_testing.smart_logging import get_default_logger

class test_psf_caller(object):
    def __init__(self,image_filename,*args,**kwargs):
        self.image_filename = image_filename
        self.args = args
        self.kwargs = kwargs
    def __call__(self,x):
        try:
            test_psf(self.image_filename,subsampling_factor=x,results_tag="ss"+str(x),*self.args,**self.kwargs)
        except Exception as e:
            logger = get_default_logger()
            logger.error("Exception when trying subsampling factor " + str(x) + ": " + str(e))
            raise

image_dir = "/disk2/brg/Data/HST_Fields/"
image_filename = "control_image_n_rs.fits"
results_dir = "/disk2/brg/Data/HST_Fields/subsampling_convergence_testing"
num_grid_points = (32,16)
parallelize = True

def main(argv):
    """ @TODO main docstring
    """
    
    subsampling_factors = np.linspace(1,20,20,endpoint=True).astype(int)
    
    # Check if we're debugging
    try:
        import pydevd as unused_import
        debugging = True
    except ImportError:
        debugging = False
        
    caller = test_psf_caller(os.path.join(image_dir,image_filename),
                             num_grid_points=num_grid_points,
                             min_lowest_separation=1.0,
                             min_class_star=0.01,
                             min_star_mag=22.,
                             max_star_mag=25.,
                             focus=mv.default_focus,
                             norm_errors=True,
                             min_star_snr=50.,
                             results_dir=results_dir,
                             parallelize=not debugging)
        
#     if debugging or not parallelize or True:
    for subsampling_factor in subsampling_factors:
        caller(subsampling_factor)
#     else:
#         nproc = max((cpu_count()-1,1))
#         pool = Pool(processes=nproc,maxtasksperchild=1)
#         pool.map(caller, subsampling_factors,chunksize=1)
#         pool.close()
#         pool.join()
#         pool.terminate()
    

if __name__ == "__main__":
    main(sys.argv)
