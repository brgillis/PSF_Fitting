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

def make_results_summary(results_filenames,
                         summary_filename):
    
    image_filenames = []
    chips = []
    focii = []
    chi_squareds = []
    emp_chi_squareds = []
    
    dofs = []
    edofs = []
    
    m0_core_diffs = []
    m0_wings_diffs = []
    m0_Zs = []
    m0_emp_Zs = []

    m0_noisy_core_diffs = []
    m0_noisy_wings_diffs = []
    m0_noisy_Zs = []
    m0_noisy_emp_Zs = []
    
    Qs = {}
    
    logger = smart_logging.get_default_logger()

    for Q_label in ("QX", "QY", "QP", "QC", "QS"):
        Qs[Q_label + "_DIFF"] = []
        Qs[Q_label + "_Z2"] = []
        Qs[Q_label + "_EZ2"] = []
        Qs[Q_label + "NDIFF"] = []
        Qs[Q_label + "NZ2"] = []
        Qs[Q_label + "NEZ2"] = []
    
    for results_filename in results_filename:
        try:
            results_file = fits.open(results_filename)
        except IOError as _e:
            logger.warn("File " + results_filename + " cannot be opened and will be skipped.")
            
        header = results_file[1].header
        
        image_filenames.append(results_filename.split('/')[-1])
        
        chips.append(header['CCDCHIP'])

        focii.append(header["FOCUS"])
        
        chi_squareds.append(header["CHI_SQR"])
        emp_chi_squareds.append(header["ECHI_SQR"])
        
        dofs.append(header["DOF"])
        edofs.append(header["EDOF"])
    
        m0_core_diffs.append(header["M0_CDIFF"])
        m0_wings_diffs.append(header["M0_WDIFF"])
        m0_Zs.append(header["M0_Z2"])
        m0_emp_Zs.append(header["M0_EZ2"])
    
        m0_noisy_core_diffs.append(header["M0NCDIFF"])
        m0_noisy_wings_diffs.append(header["M0NWDIFF"])
        m0_noisy_Zs.append(header["M0NZ2"])
        m0_noisy_emp_Zs.append(header["M0NEZ2"])
    
        for Q_label in ("QX", "QY", "QP", "QC", "QS"):
            Qs[Q_label + "_DIFF"].append(header[Q_label + "_DIFF"])
            Qs[Q_label + "_Z2"].append(header[Q_label + "_Z2"])
            Qs[Q_label + "_EZ2"].append(header[Q_label + "_EZ2"])
            Qs[Q_label + "NDIFF"].append(header[Q_label + "NDIFF"])
            Qs[Q_label + "NZ2"].append(header[Q_label + "NZ2"])
            Qs[Q_label + "NEZ2"].append(header[Q_label + "NEZ2"])
    
    columns = [fits.Column(name="filename", format='30A', array=image_filenames),
               fits.Column(name="CHIP", format='B', array=chips),
               fits.Column(name="focus", format='E', array=focii),
               fits.Column(name="chi_squared", format='E', array=chi_squareds),
               fits.Column(name="dofs", format='E', array=dofs),
               fits.Column(name="emp_chi_squared", format='E', array=emp_chi_squareds),
               fits.Column(name="emp_dofs", format='E', array=edofs),
               fits.Column(name="m0_core_diff", format='E', array=m0_core_diffs),
               fits.Column(name="m0_wings_diff", format='E', array=m0_wings_diffs),
               fits.Column(name="m0_Z2", format='E', array=m0_Zs),
               fits.Column(name="m0_emp_Z2", format='E', array=m0_emp_Zs)]
    
    for Q_label, colname in zip(("QX", "QY", "QP", "QC", "QS"), 
                                ("Qx", "Qy", "Qplus", "Qcross", "Qsize")):
        columns.append(fits.Column(name=colname + "_diff", format='E', array=Qs[Q_label + "_DIFF"]))
        columns.append(fits.Column(name=colname + "_Z2", format='E', array=Qs[Q_label + "_Z2"]))
        columns.append(fits.Column(name=colname + "_emp_Z2", format='E', array=Qs[Q_label + "_EZ2"]))
        
    columns += [fits.Column(name="m0_noisy_core_diff", format='E', array=m0_noisy_core_diffs),
                fits.Column(name="m0_noisy_wings_diff", format='E', array=m0_noisy_wings_diffs),
                fits.Column(name="m0_noisy_Z2", format='E', array=m0_noisy_Zs),
                fits.Column(name="m0_noisy_emp_Z2", format='E', array=m0_noisy_emp_Zs)]
    
    for Q_label, colname in zip(("QX", "QY", "QP", "QC", "QS"), 
                                ("Qx", "Qy", "Qplus", "Qcross", "Qsize")):
        columns.append(fits.Column(name=colname + "_noisy_diff", format='E', array=Qs[Q_label + "NDIFF"]))
        columns.append(fits.Column(name=colname + "_noisy_Z2", format='E', array=Qs[Q_label + "NZ2"]))
        columns.append(fits.Column(name=colname + "_noisy_emp_Z2", format='E', array=Qs[Q_label + "NEZ2"]))

    tbhdu = fits.BinTableHDU.from_columns(columns)

    tbhdu.writeto(summary_filename,clobber=True)
    
    logger.info("Results summary output to " + summary_filename + ".")
    
    return
        
        