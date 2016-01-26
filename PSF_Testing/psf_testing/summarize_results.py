""" @file summarize_results.py

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

from astropy.io import fits
from psf_testing import smart_logging
from psf_testing import magic_values as mv

def make_results_summary(results_filenames,
                         summary_filename):
    
    image_filenames = []
    chips = []
    obs_times = []
    focii = []
    
    X_squareds = []
    chi_squareds = []
    
    X2_dofs = []
    chi2_dofs = []
    
    m0_diff_diffs = []
    m0_Zs = []

    m0_noisy_diff_diffs = []
    m0_noisy_Zs = []
    
    Qs = {}
    
    logger = smart_logging.get_default_logger()

    for Q_label in ("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"):
        Qs[Q_label + "_DIF"] = []
        Qs[Q_label + "_Z2"] = []
        Qs[Q_label + "NDIF"] = []
        Qs[Q_label + "NZ2"] = []
    
    for results_filename in results_filenames:
        try:
            results_file = fits.open(results_filename)
        except IOError as _e:
            logger.warn("File " + results_filename + " cannot be opened and will be skipped.")
            continue
            
        header = results_file[1].header
        
        image_filenames.append(results_filename.split('/')[-1].replace(mv.results_tail,
                                                                       mv.image_extension))
        
        chips.append(header['CCDCHIP'])
        obs_times.append(header['OBS_TIME'])

        focii.append(header["FOCUS"])
        
        X_squareds.append(header["X_SQR"])
        chi_squareds.append(header["CHI_SQR"])
        
        X2_dofs.append(header["XDOF"])
        chi2_dofs.append(header["CDOF"])
        
        m0_diff_diffs.append(header["M0D_DIF"])
        m0_Zs.append(header["M0D_Z2"])
    
        m0_noisy_diff_diffs.append(header["M0DNDIF"])
        m0_noisy_Zs.append(header["M0DNZ2"])
    
        for Q_label in ("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"):
            Qs[Q_label + "_DIF"].append(header[Q_label + "_DIF"])
            Qs[Q_label + "_Z2"].append(header[Q_label + "_Z2"])
            Qs[Q_label + "NDIF"].append(header[Q_label + "NDIF"])
            Qs[Q_label + "NZ2"].append(header[Q_label + "NZ2"])
    
    columns = [fits.Column(name="filename", format='30A', array=image_filenames),
               fits.Column(name="chip", format='B', array=chips),
               fits.Column(name="obs_time", format='E', array=obs_times),
               fits.Column(name="focus", format='E', array=focii),
               fits.Column(name="X_squared", format='E', array=X_squareds),
               fits.Column(name="chi_squared", format='E', array=chi_squareds),
               fits.Column(name="X2_dofs", format='E', array=X2_dofs),
               fits.Column(name="chi2_dofs", format='E', array=chi2_dofs),
               fits.Column(name="m0_diff_diff", format='E', array=m0_diff_diffs),
               fits.Column(name="m0_Z2", format='E', array=m0_Zs)]
    
    for Q_label, colname in zip(("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"), 
                                ("Qx_diff", "Qy_diff",
                                 "Qplus_sum", "Qcross_sum", "Qsize_sum",
                                 "Qplus_diff", "Qcross_diff", "Qsize_diff")):
        columns.append(fits.Column(name=colname + "_diff", format='E', array=Qs[Q_label + "_DIF"]))
        columns.append(fits.Column(name=colname + "_Z2", format='E', array=Qs[Q_label + "_Z2"]))
        
    columns += [fits.Column(name="m0_noisy_diff_diff", format='E', array=m0_noisy_diff_diffs),
                fits.Column(name="m0_noisy_Z2", format='E', array=m0_noisy_Zs)]
    
    for Q_label, colname in zip(("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"), 
                                ("Qx_diff", "Qy_diff",
                                 "Qplus_sum", "Qcross_sum", "Qsize_sum",
                                 "Qplus_diff", "Qcross_diff", "Qsize_diff")):
        columns.append(fits.Column(name=colname + "_noisy_diff", format='E', array=Qs[Q_label + "NDIF"]))
        columns.append(fits.Column(name=colname + "_noisy_Z2", format='E', array=Qs[Q_label + "NZ2"]))

    tbhdu = fits.BinTableHDU.from_columns(columns)

    tbhdu.writeto(summary_filename,clobber=True)
    
    logger.info("Results summary output to " + summary_filename + ".")
    
    return
        
        