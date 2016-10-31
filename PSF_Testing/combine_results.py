#!/usr/bin/env python

""" @file combine_results.py

    Created 11 Nov 2015

    @TODO: File docstring

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
import argparse

from utility.smart_logging import get_default_logger
from psf_testing.summarize_results import make_results_summary, make_stack_stacks
from psf_testing import magic_values as mv

def main(argv):
    """ @TODO main docstring
    """
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--image_list_filename",type=str, default=None,
                        help="The filename which contains a list of images to test.")
    parser.add_argument("--image_dir",type=str, default=None,
                        help="The directory where images can be found.")
    parser.add_argument("--out",type=str, default="psf_testing_results_summary.fits",
                        help="The desired output filename, containing the results summary.")
    parser.add_argument("--stack_root",type=str, default="psf_testing",
                        help="The desired root of stacks of stacks.")
    
    args = parser.parse_args()
        
    logger = get_default_logger()
    
    # Check that we have both needed arguments
    if (args.image_list_filename is None) or (args.image_dir is None):
        logger.error("Both the '--image_list_filename' and '--image_dir' arguments must be passed " +
                     "at the command-line. eg. ' python combine_results.py --image_list_filename " +
                     " /path/to/my_list.dat --image_dir /path/to/images/'.")
        
    results_filename_roots = []
    with open(args.image_list_filename) as fi:
        for line in fi:
            for word in line.split():
                image_filename = join(args.image_dir, word).replace(mv.image_extension,"")
                results_filename_roots.append(image_filename)
                
    make_stack_stacks(results_filename_roots=results_filename_roots,
                         stack_stack_filename_root=args.stack_root)
                
    make_results_summary(results_filename_roots=results_filename_roots,
                         summary_filename=args.out)
    
    logger.info("Execution complete.")

if __name__ == "__main__":
    main(sys.argv)
