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

data_dir = "/disk2/brg/Data/HST_Fields/control_results"
default_results_file_root = "control_image_n"

default_output_filename_root = "control_field_tests"
default_output_extension = "tex"

paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"

num_images = 100

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
                     ("1D Guiding Error", "_blur"),
                     ("2D Guiding Error", "_blur2"),
                     ("Galaxy Background", "_gb"),
                     ("Varying Spec. Type", "_rs"),
                     ("Full", "_full"))
    
    labels_and_colnames = (("X2",r"$X^2$"),
                           ("Qx_diff_Z2",r"$Z^{2(-)}_{x}$"),
                           ("Qy_diff_Z2",r"$Z^{2(-)}_{y}$"),
                           ("Qp_sum_Z2",r"$Z^{2(+)}_{+}$"),
                           ("Qp_diff_Z2",r"$Z^{2(-)}_{+}$"),
                           ("Qc_sum_Z2",r"$Z^{2(+)}_{\times}$"),
                           ("Qc_diff_Z2",r"$Z^{2(-)}_{\times}$"),
                           ("Qs_sum_Z2",r"$Z^{2(+)}_{\rm s}$"),
                           ("Qs_diff_Z2",r"$Z^{2(-)}_{\rm s}$"))
    
    for results_tail in "_results.fits", "_fit_results.fits", "_fit_all_params_results.fits":
    
        vals = {}
        for control_type, _ in control_types:
            vals[control_type] = {}
            vals[control_type]["NSTAR"] = np.zeros(num_images, dtype=int)
            for label, _ in labels_and_colnames:
                vals[control_type][label] = np.zeros(num_images, dtype=float)
                vals[control_type]["Focus"] = np.zeros(num_images, dtype=float)
        
        for n in range(num_images):
            
            results_root_n = results_root.replace("_n","_"+str(n))
            
            for control_type, tag in control_types:
                
                if results_tail=="_fit_all_params_results.fits":
                    if control_type not in ("Base","Full"):
                        continue
            
                filename = results_root_n + tag + results_tail
                
                if not os.path.exists(filename):
                    raise Exception("Expected file doesn't exist: " + filename)
            
                try:
                    header = fits.open(filename)[1].header
                except Exception:
                    raise
            
                vals[control_type]["NSTAR"][n] = header["NSTAR"]
                vals[control_type]["X2"][n] = header["X_SQR"]
                vals[control_type]["Qx_diff_Z2"][n] = header["QXD_Z2"]
                vals[control_type]["Qy_diff_Z2"][n] = header["QYD_Z2"]
                vals[control_type]["Qp_sum_Z2"][n] = header["QPS_Z2"] 
                vals[control_type]["Qp_diff_Z2"][n] = header["QPD_Z2"] 
                vals[control_type]["Qc_sum_Z2"][n] = header["QCS_Z2"] 
                vals[control_type]["Qc_diff_Z2"][n] = header["QCD_Z2"] 
                vals[control_type]["Qs_sum_Z2"][n] = header["QSS_Z2"] 
                vals[control_type]["Qs_diff_Z2"][n] = header["QSD_Z2"]
                vals[control_type]["Focus"][n] = header["FOCUS"]
                
        # Get the means and stderrs
        val_means = {}
        val_stderrs = {}
        for control_type, _ in control_types:
            val_means[control_type] = {}
            val_stderrs[control_type] = {}
            for label, _ in labels_and_colnames:
                val_means[control_type][label] = vals[control_type][label].mean()
                val_stderrs[control_type][label] = vals[control_type][label].std()
            val_means[control_type]["Focus"] = vals[control_type]["Focus"].mean()
            val_stderrs[control_type]["Focus"] = vals[control_type]["Focus"].std()
                
        if results_tail == "_results.fits":
            output_root = args.output_filename_root
            table_label = "Known Focus Offset"
        elif results_tail == "_fit_results.fits": 
            output_root = args.output_filename_root + "_fit"
            table_label = "Fit Focus Offset"
        elif results_tail == "_fit_all_params_results.fits": 
            output_root = args.output_filename_root + "_fit_all_params"
            table_label = "Fit All Params"
        else:
            raise ValueError("Unknown results_tail value: " + str(results_tail))
    
        # Print the means in table format
        table_filename = output_root + "." + args.output_extension
        
        with open(table_filename,'w') as fo:
            
            # Write the header portion
            fo.write("\\begin{center}\n")
            fo.write("\t\\begin{tabular}{llllllllll}\n")
            fo.write("\t\t\\multicolumn{10}{c}{\\textbf{"+table_label+"}} \\\\")
            
            line = "\t\tType"
            for _, colname in labels_and_colnames:
                line += " & " + colname
            line += " \\\\ \\hline \n"
            fo.write(line)
            
            # Print a line for each control type
            for control_type, _ in control_types:
                
                if results_tail=="_fit_all_params_results.fits":
                    if control_type not in ("Base","Full"):
                        continue
                    
                line = "\t\t" + control_type
                for label, _ in labels_and_colnames:
                    line += (" & $%1.1e}$" % val_means[control_type][label]).replace("e","\\times 10^{").replace("{+0","{").replace("{-0","{-")
                line += " \\\\\n"
                
                fo.write(line)
            
            # Write the tail portion
            fo.write("\t\\end{tabular}\n")
            fo.write("\\end{center}\n")
    
        # Print the standard errors in table format
        stderr_table_filename = output_root + "_stderr." + args.output_extension
        
        with open(stderr_table_filename,'w') as fo:
            
            # Write the header portion
            fo.write("\\begin{center}\n")
            fo.write("\t\\begin{tabular}{llllllllll}\n")
            fo.write("\t\t\\multicolumn{10}{c}{\\textbf{"+table_label+" Errors}} \\\\")
            
            line = "\t\tType"
            for _, colname in labels_and_colnames:
                line += " & " + colname
            line += " \\\\ \\hline \n"
            fo.write(line)
            
            # Print a line for each control type
            for control_type, _ in control_types:
                
                if results_tail=="_fit_all_params_results.fits":
                    if control_type not in ("Base","Full"):
                        continue
                    
                line = "\t\t" + control_type
                for label, _ in labels_and_colnames:
                    line += (" & $%1.1e}$" % val_stderrs[control_type][label]).replace("e","\\times 10^{").replace("{+0","{").replace("{-0","{-")
                line += " \\\\\n"
                
                fo.write(line)
            
            # Write the tail portion
            fo.write("\t\\end{tabular}\n")
            fo.write("\\end{center}\n")
        
        # Copy them to the paper location if in tex format
        if args.output_extension=="tex":
            cmd = "cp " + table_filename + " " + paper_location
            sbp.call(cmd,shell=True)
            cmd = "cp " + stderr_table_filename + " " + paper_location
            sbp.call(cmd,shell=True)

if __name__ == "__main__":
    main(sys.argv)
