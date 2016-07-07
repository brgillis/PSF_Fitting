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
from plot_X_v_focus_v_chip import make_X_v_focus_v_chip_plot
from plot_X_v_obs_time_v_nstar import make_X_v_obs_time_v_nstar_plot
    
default_paper_location = "/disk2/brg/Dropbox/gillis-comp-shared/Papers/PSF_Model_Testing/"

def make_all_plots(paper_location = default_paper_location):
    
    make_fit_statistic_plots(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_inc_size_summary.fits",
                             plot_name="fitting_paramater_inc_size_hists",
                             red_X2_min=0.1,
                             red_X2_max=10.0,
                             fade_size=False,
                             paper_location=paper_location,
                             file_type="eps",
                             hide=True)
    make_fit_statistic_plots(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_noinc_size_summary.fits",
                             plot_name="fitting_paramater_noinc_size_hists",
                             red_X2_min=0.1,
                             red_X2_max=10.0,
                             fade_size=True,
                             file_type="eps",
                             paper_location=paper_location,
                             hide=True)
    
    make_fit_statistic_plots(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_inc_size_summary.fits",
                             plot_name="fitting_paramater_chi_inc_size_hists",
                             red_X2_min=0.1,
                             red_X2_max=1000.0,
                             paper_location=paper_location,
                             plot_chi2=True,
                             file_type="eps",
                             hide=True)
    make_fit_statistic_plots(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_noinc_size_summary.fits",
                             plot_name="fitting_paramater_chi_noinc_size_hists",
                             red_X2_min=0.1,
                             red_X2_max=1000.0,
                             paper_location=paper_location,
                             plot_chi2=True,
                             file_type="eps",
                             hide=True)
    
    make_Qsize_plots(file_type="eps",
                     seed=1,
                     paper_location=paper_location,
                     hide=True)
    
    # Stacks inc size "good"
    make_stacks(image_location="/disk2/brg/Data/HST_Fields/results_inc_size/",
                image_name="jb6v21f7q_sci1_cor",
                plot_name_tail="_stacks_inc_size",
                paper_location=paper_location,
                file_type="eps",
                hide=True)
    
    # Stacks inc size "bad"
    make_stacks(image_location="/disk2/brg/Data/HST_Fields/results_inc_size/",
                image_name="jb6v09shq_sci2_cor",
                plot_name_tail="_stacks_inc_size",
                paper_location=paper_location,
                file_type="eps",
                hide=True)
    
    # Stacks noinc size "good"
    make_stacks(image_location="/disk2/brg/Data/HST_Fields/results_noinc_size/",
                image_name="jb6v21f7q_sci1_cor",
                plot_name_tail="_stacks_noinc_size",
                paper_location=paper_location,
                file_type="eps",
                hide=True)
    
    # Stacks noinc size "bad"
    make_stacks(image_location="/disk2/brg/Data/HST_Fields/results_noinc_size/",
                image_name="jb6v09shq_sci2_cor",
                plot_name_tail="_stacks_noinc_size",
                paper_location=paper_location,
                file_type="eps",
                hide=True)
    
    # X v focus v chip inc size
    make_X_v_focus_v_chip_plot(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_inc_size_summary.fits",
                               plot_name="X_v_focus_v_chip_inc_size",
                               red_X2_max=10.0,
                               red_X2_min=1,
                               file_type="eps",
                               hide=True)
    
    # X v focus v chip noinc size
    make_X_v_focus_v_chip_plot(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_noinc_size_summary.fits",
                               plot_name="X_v_focus_v_chip_noinc_size",
                               red_X2_max=10.0,
                               red_X2_min=1,
                               file_type="eps",
                               hide=True)
    
    # chi v focus v chip inc size
    make_X_v_focus_v_chip_plot(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_inc_size_summary.fits",
                               plot_name="chi_v_focus_v_chip_inc_size",
                               plot_chi2=True,
                               red_X2_max=1000.0,
                               red_X2_min=0.1,
                               file_type="eps",
                               hide=True)
    
    # chi v focus v chip noinc size
    make_X_v_focus_v_chip_plot(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_noinc_size_summary.fits",
                               plot_name="chi_v_focus_v_chip_noinc_size",
                               plot_chi2=True,
                               red_X2_max=1000.0,
                               red_X2_min=0.1,
                               file_type="eps",
                               hide=True)
    
    # X v obs_time v nstar inc size
    make_X_v_obs_time_v_nstar_plot(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_inc_size_summary.fits",
                               plot_name="X_v_obs_time_v_nstar_inc_size",
                               red_X2_max=10.0,
                               red_X2_min=1,
                               file_type="eps",
                               hide=True)
    
    # X v obs_time v nstar noinc size
    make_X_v_obs_time_v_nstar_plot(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_noinc_size_summary.fits",
                               plot_name="X_v_obs_time_v_nstar_noinc_size",
                               red_X2_max=10.0,
                               red_X2_min=1,
                               file_type="eps",
                               hide=True)
    
    # chi v obs_time v nstar inc size
    make_X_v_obs_time_v_nstar_plot(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_inc_size_summary.fits",
                               plot_name="chi_v_obs_time_v_nstar_inc_size",
                               plot_chi2=True,
                               red_X2_max=1000.0,
                               red_X2_min=0.1,
                               file_type="eps",
                               hide=True)
    
    # chi v obs_time v nstar noinc size
    make_X_v_obs_time_v_nstar_plot(summary_filename="/disk2/brg/git/Tiny_Tim_PSF_Fitting/PSF_Testing/psf_testing_results_noinc_size_summary.fits",
                               plot_name="chi_v_obs_time_v_nstar_noinc_size",
                               plot_chi2=True,
                               red_X2_max=1000.0,
                               red_X2_min=0.1,
                               file_type="eps",
                               hide=True)

def main(argv):
    """ @TODO main docstring
    """
    
    make_all_plots()
    
    return

if __name__ == "__main__":
    main(sys.argv)
