#!/usr/bin/env python

""" @file plot_subsampling_convergence.py

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

import numpy as np


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

data_dir = "/disk2/brg/Data/HST_Fields/"
default_results_file_root = "control_image_n"

default_output_filename_root = "control_field_tests"
default_output_extension = "tex"

paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"

num_images = 10

def main(argv):
    """ @TODO main docstring
    """
    
    parser = ArgumentParser()
    
    parser.add_argument("--results_file_root",default=default_results_file_root)
    parser.add_argument("--output_filename_root",default=default_output_filename_root)
    parser.add_argument("--output_extension",default=default_output_extension)
    
    args = parser.parse_args()

    results_root = os.path.join(data_dir, args.results_file_root)
    
    control_types = (("Base",""),
                     ("Binaries", "_b3"),
                     ("Wide Binaries", "_b32"),
                     ("Guiding Error", "_blur"),
                     ("Double Guiding Error", "_blur2"),
                     ("Galaxy Background", "_gb"),
                     ("Varying Spectral Type", "_rs"),
                     ("Full", "_full"))
    
    labels_and_colnames = (("X2",r"$X^2$"),
                           ("Qx_diff_Z2",r"$Z^{2(-)}_{x}$"),
                           ("Qy_diff_Z2",r"$Z^{2(-)}_{y}$"),
                           ("Qp_sum_Z2",r"$Z^{2(+)}_{+}$"),
                           ("Qp_diff_Z2",r"$Z^{2(-)}_{+}$"),
                           ("Qc_sum_Z2",r"$Z^{2(+)}_{\times}$"),
                           ("Qc_diff_Z2",r"$Z^{2(-)}_{\times}$"),
                           ("Qs_sum_Z2",r"$Z^{2(+)}_{s}$"),
                           ("Qs_diff_Z2",r"$Z^{2(-)}_{s}$"))
    
    vals = {}
    for control_type, _ in control_types:
        vals[control_type] = {}
        for label, _ in labels_and_colnames:
            vals[control_type][label] = np.zeros(num_images, dtype=float)
    
    for n in range(num_images):
        
        results_root_n = results_root.replace("_n","_"+str(n))
        
        for control_type, tag in control_types:
        
            filename = results_root_n + tag + "_results.fits"
        
            header = fits.open(filename)[1].header
        
            vals[control_type]["X2"][n] = header["X_SQR"]
            vals[control_type]["Qx_diff_Z2"][n] = header["QXD_Z2"]
            vals[control_type]["Qy_diff_Z2"][n] = header["QYD_Z2"]
            vals[control_type]["Qp_sum_Z2"][n] = header["QPS_Z2"] 
            vals[control_type]["Qp_diff_Z2"][n] = header["QPD_Z2"] 
            vals[control_type]["Qc_sum_Z2"][n] = header["QPS_Z2"] 
            vals[control_type]["Qc_diff_Z2"][n] = header["QPD_Z2"] 
            vals[control_type]["Qs_sum_Z2"][n] = header["QSS_Z2"] 
            vals[control_type]["Qs_diff_Z2"][n] = header["QSD_Z2"]
            
    # Get the means
    val_means = {}
    for control_type, _ in control_types:
        val_means[control_type] = {}
        for label, _ in labels_and_colnames:
            val_means[control_type][label] = vals[control_type][label].mean()

    # Print it in table format
    table_filename = args.output_filename_root + "." + args.output_extension
    
    with open(table_filename,'w') as fo:
        
        # Write the header portion
        fo.write("\\begin{center}\n")
        fo.write("\t\\begin{tabular}{llllllllll}\n")
        
        line = "\t\tType"
        for _, colname in labels_and_colnames:
            line += " & " + colname
        line += " \\\\ \\hline \n"
        fo.write(line)
        
        # Print a line for each control type
        for control_type, _ in control_types:
            line = "\t\t" + control_type
            for label, _ in labels_and_colnames:
                line += (" & $%1.1e}$" % val_means[control_type][label]).replace("e","\\times 10^{").replace("{+0","{").replace("{-0","{-")
            line += " \\\\\n"
            
            fo.write(line)
        
        # Write the tail portion
        fo.write("\t\\end{tabular}\n")
        fo.write("\\end{center}\n")
    
    # Copy it to the paper location if in tex format
    if args.output_extension=="tex":
        cmd = "cp " + table_filename + " " + paper_location
        sbp.call(cmd,shell=True)

if __name__ == "__main__":
    main(sys.argv)
