#!/usr/bin/env python

""" @file /disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_fit_plotting/make_all_paper_plots.py

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

import sys

from plot_fit_statistics import make_fit_statistic_plots
from plot_Qsize import make_Qsize_plots
from plot_stacks import make_stacks
    
default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"

def make_all_plots(paper_location = default_paper_location):
    
    make_fit_statistic_plots(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_inc_size_summary.fits",
                             plot_name="fitting_paramater_inc_size_hists",
                             red_X2_max=100.0,
                             fade_size=False,
                             paper_location=paper_location,
                             file_type="eps",
                             hide=True)
    make_fit_statistic_plots(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_noinc_size_summary.fits",
                             plot_name="fitting_paramater_noinc_size_hists",
                             red_X2_max=10.0,
                             fade_size=True,
                             file_type="eps",
                             paper_location=paper_location,
                             hide=True)
    
    make_Qsize_plots(file_type="eps",
                     paper_location=paper_location,
                     hide=True)
    
    # Inc size "good"
    make_stacks(image_location="/disk2/brg/Data/HST_Fields/results_inc_size/",
                image_name="jb6v21f7q_sci1_cor",
                plot_name_tail="_stacks_inc_size",
                paper_location=paper_location,
                file_type="eps",
                hide=True)
    
    # Inc size "bad"
    make_stacks(image_location="/disk2/brg/Data/HST_Fields/results_inc_size/",
                image_name="jb6v09shq_sci2_cor",
                plot_name_tail="_stacks_inc_size",
                paper_location=paper_location,
                file_type="eps",
                hide=True)
    
    # Noinc size "good"
    make_stacks(image_location="/disk2/brg/Data/HST_Fields/results_noinc_size/",
                image_name="jb6v25azq_sci2_cor",
                plot_name_tail="_stacks_noinc_size",
                paper_location=paper_location,
                file_type="eps",
                hide=True)
    
    # Noinc size "bad"
    make_stacks(image_location="/disk2/brg/Data/HST_Fields/results_noinc_size/",
                image_name="jb6v24ugq_sci1_cor",
                plot_name_tail="_stacks_noinc_size",
                paper_location=paper_location,
                file_type="eps",
                hide=True)

def main(argv):
    """ @TODO main docstring
    """
    
    make_all_plots()
    
    return

if __name__ == "__main__":
    main(sys.argv)
