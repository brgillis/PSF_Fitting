#!/usr/bin/env python

""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/Run_PSF_Testing.py

    Created 16 Sep 2015

    Main module for running PSF testing.

    ---------------------------------------------------------------------

    Copyright (C) 2015 Bryan R. Gillis

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

import sys
from os.path import join
from multiprocessing import Pool, cpu_count

from psf_testing.arg_parser import get_arg_parser, parse_and_get_kwargs
from psf_testing.test_psf import test_psf
from utility.smart_logging import get_default_logger

class test_psf_caller(object):
    def __init__(self,*args,**kwargs):
        self.args = args,
        self.kwargs = kwargs
    def __call__(self,x):
        try:
            test_psf(x,*self.args,**self.kwargs)
        except AssertionError as e:
            logger = get_default_logger()
            logger.error("Exception when processing file " + str(x) + ": " + str(e))
        except Exception as e:
            logger = get_default_logger()
            logger.error("Exception when processing file " + str(x) + ": " + str(e))
            if "usable stars" not in str(e):
                raise
            
def main(argv):
    """ @TODO main docstring
    """
    
    # Check if we're debugging
    try:
        import pydevd as _
        debugging = True
    except ImportError:
        debugging = False
    debugging = False # FIXME
        
    # Execute command-line parsing
    parser = get_arg_parser()
    kwargs, special_kwargs = parse_and_get_kwargs(parser,
                                                  special_keys=("image_filename","image_list_filename",
                                                                "image_dir","logging_level","disable_parallelization"))
        
    # Set up logger
    logger = get_default_logger()
    logger.setLevel(special_kwargs["logging_level"].upper())
    
    # Pass the cline-args to the test_psf function, which carries out the testing
    
    if special_kwargs["image_filename"] is not None:
        test_psf_caller(
                 parallelize = (not debugging) and (not special_kwargs["disable_parallelization"]),
                 **kwargs )(special_kwargs["image_filename"])
        image_filenames = [special_kwargs["image_filename"]]
    else:
        # We'll test a list of images
        image_filenames = []
        with open(special_kwargs["image_list_filename"]) as fi:
            for line in fi:
                for word in line.split():
                    image_filename = join(special_kwargs["image_dir"], word)
                    image_filenames.append(image_filename)
        
        if kwargs["refresh_only"]:
            nproc = 1
            parallelize_each_image = (not debugging) and (not special_kwargs["disable_parallelization"])
        else:
            nproc = max((cpu_count()-1,1))
            parallelize_each_image = False
        
        if (nproc == 1) or (special_kwargs["disable_parallelization"]):
            for image_filename in image_filenames:
                test_psf_caller(parallelize = parallelize_each_image,
                                **kwargs )(image_filename)
        else:
        
            pool = Pool(processes=nproc,maxtasksperchild=1)
            _ = pool.map(test_psf_caller(parallelize = False,
                                         **kwargs ),
                     image_filenames,chunksize=1)
            pool.close()
            pool.join()
    
    logger.info("Execution complete.")

if __name__ == "__main__":
    main(sys.argv)
