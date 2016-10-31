""" @file report_results.py

    Created 5 Nov 2015

    Reports the results of a PSF test, printing basic results in the
    console and detailed results in a fits table.

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
import numpy as np

from psf_testing import magic_values as mv
from utility.smart_logging import get_default_logger

def save_fitting_record(fitting_record,
                        filename_root):
    params = {}
    for param in mv.default_params:
        params[param] = []
    params["focus"] = []
    X_squareds = []
    chi_squareds = []
    
    m0_diffs = []
    m0_Zs = []

    m0_noisy_diffs = []
    m0_noisy_Zs = []
    
    Qs = {}

    for Q_label in ("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"):
        Qs[Q_label + "_DIF"] = []
        Qs[Q_label + "_Z2"] = []
        Qs[Q_label + "NDIFF"] = []
        Qs[Q_label + "NZ2"] = []

    for Q_label in ("QXC", "QYC", "QPC", "QCC", "QSC", "QXW", "QYW", "QPW", "QCW", "QSW"):
        Qs[Q_label + "_DIF"] = []
        
    for test_results in fitting_record:
        
        for param in mv.default_params:
            params[param].append(test_results[0][param])
        params["focus"].append(test_results[0]["focus"])
        
        X_squareds.append(test_results[1][0])
        chi_squareds.append(test_results[1][2])
    
        m0_diffs.append(test_results[2][0][0])
        m0_Zs.append(test_results[3][0][0])
    
        m0_noisy_diffs.append(test_results[2][1][0])
        m0_noisy_Zs.append(test_results[3][1][0])
    
        for Q_label, j in zip(("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"), range(8)):
            Qs[Q_label + "_DIF"].append(test_results[2][0][1][j])
            Qs[Q_label + "_Z2"].append(test_results[3][0][1][j])
            Qs[Q_label + "NDIFF"].append(test_results[2][1][1][j])
            Qs[Q_label + "NZ2"].append(test_results[3][1][1][j])
    
        for Q_label, j in zip(("QXC", "QYC", "QPC", "QCC", "QSC", "QXW", "QYW", "QPW", "QCW", "QSW",),
                              range(10)):
            Qs[Q_label + "_DIF"].append(np.mean(test_results[11][j]))
    
    columns = []
    
    columns.append(fits.Column(name="focus", format='E', array=params["focus"]))
    for param in mv.default_params:
        columns.append(fits.Column(name=param, format='E', array=params[param]))
        
    columns.append(fits.Column(name="X_squared", format='E', array=X_squareds))
    columns.append(fits.Column(name="chi_squared", format='E', array=chi_squareds))
    columns.append(fits.Column(name="m0_diff", format='E', array=m0_diffs))
    columns.append(fits.Column(name="m0_Z2", format='E', array=m0_Zs))
    
    for Q_label, colname in zip(("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"), 
                                ("Qx_diff", "Qy_diff", "Qplus_sum", "Qcross_sum", "Qsize_sum",
                                 "Qplus_diff", "Qcross_diff", "Qsize_diff")):
        columns.append(fits.Column(name=colname + "_diff", format='E', array=Qs[Q_label + "_DIF"]))
        columns.append(fits.Column(name=colname + "_Z2", format='E', array=Qs[Q_label + "_Z2"]))
    
    for Q_label, colname in zip(("QXC", "QYC", "QPC", "QCC", "QSC", "QXW", "QYW", "QPW", "QCW", "QSW",), 
                                ("Qx_core", "Qy_core", "Qplus_core", "Qcross_core", "Qsize_core",
                                 "Qx_wings", "Qy_wings", "Qplus_wings", "Qcross_wings", "Qsize_wings")):
        columns.append(fits.Column(name=colname + "_diff", format='E', array=Qs[Q_label + "_DIF"]))
        
    columns += [fits.Column(name="m0_noisy_diff", format='E', array=m0_noisy_diffs),
                fits.Column(name="m0_noisy_Z2", format='E', array=m0_noisy_Zs)]
    
    for Q_label, colname in zip(("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"), 
                                ("Qx_diff", "Qy_diff", "Qplus_sum", "Qcross_sum", "Qsize_sum",
                                 "Qplus_diff", "Qcross_diff", "Qsize_diff")):
        columns.append(fits.Column(name=colname + "_noisy_diff", format='E', array=Qs[Q_label + "NDIFF"]))
        columns.append(fits.Column(name=colname + "_noisy_Z2", format='E', array=Qs[Q_label + "NZ2"]))

    tbhdu = fits.BinTableHDU.from_columns(columns)
    
    tbhdu.header["NSTAR"] = test_results[1][1]
    tbhdu.header["DOF"] = test_results[1][3]

    fitting_record_filename = filename_root + "_fitting_record" + mv.table_extension

    tbhdu.writeto(fitting_record_filename,clobber=True)

def report_results(test_results,
                   filename_root,
                   chip,
                   obs_time,
                   exp_time,
                   ra,
                   dec,
                   fitting_record=None):

    # Set up the columns first
    columns = []
        
    Q_err_index = 7

    for weight_label, k in zip(("core_", "wings_"), (0,1)):

        for label, i, err_i in zip(("star_", "model_", "noisy_model_"),
                                   (4, 5, 6),
                                   (Q_err_index, None, Q_err_index)):
    
            if err_i is None:
                get_errs = np.zeros_like
                err_i = Q_err_index
            else:
                get_errs = lambda x : x
    
            columns.append(fits.Column(name=weight_label + label + "m0", format='E',
                                       array=test_results[i][0][:,k]))
            columns.append(fits.Column(name=weight_label + label + "m0" + "_err", format='E',
                                       array=get_errs(test_results[err_i][0][:,k])))
    
            for Q_label, j in zip(("Qx", "Qy"), range(2)):
    
                columns.append(fits.Column(name=weight_label + label + Q_label, format='E',
                                           array=test_results[i][1][:,j,k]))
                columns.append(fits.Column(name=weight_label + label + Q_label + "_err", format='E',
                                           array=get_errs(test_results[err_i][1][:,j,k])))
    
            for Q_label, j in zip(("Qplus", "Qcross", "Qsize"), range(3)):
    
                columns.append(fits.Column(name=weight_label + label + Q_label, format='E',
                                           array=test_results[i][2][:,j,k]))
                columns.append(fits.Column(name=weight_label + label + Q_label + "_err", format='E',
                                           array=get_errs(test_results[err_i][2][:,j,k])))
        
    columns.append(fits.Column(name="star_x_pix", format='E', array=test_results[10][0]))
    columns.append(fits.Column(name="star_y_pix", format='E', array=test_results[10][1]))

    columns.append(fits.Column(name="is_not_outlier", format='L',
                               array=np.logical_not(test_results[8][0])))

    tbhdu = fits.BinTableHDU.from_columns(columns)

    # Add metadata to the header

    tbhdu.header["CCDCHIP"] = chip
    tbhdu.header["OBS_TIME"] = obs_time
    tbhdu.header["EXP_TIME"] = exp_time
    tbhdu.header["RA_TARG"] = ra
    tbhdu.header["DEC_TARG"] = dec
    tbhdu.header["FOCUS"] = test_results[0]["focus"]
    tbhdu.header["Z2"] = test_results[0]["z2"]
    tbhdu.header["Z3"] = test_results[0]["z3"]
    tbhdu.header["ASTIG_0"] = test_results[0]["astigmatism_0"]
    tbhdu.header["ASTIG_45"] = test_results[0]["astigmatism_45"]
    tbhdu.header["COMA_X"] = test_results[0]["coma_x"]
    tbhdu.header["COMA_X"] = test_results[0]["coma_x"]
    tbhdu.header["CLOVER_X"] = test_results[0]["clover_x"]
    tbhdu.header["CLOVER_Y"] = test_results[0]["clover_y"]
    tbhdu.header["SPHERE_3"] = test_results[0]["spherical_3rd"]
    tbhdu.header["Z12"] = test_results[0]["z12"]
    tbhdu.header["Z13"] = test_results[0]["z13"]
    tbhdu.header["Z14"] = test_results[0]["z14"]
    tbhdu.header["Z15"] = test_results[0]["z15"]
    tbhdu.header["Z16"] = test_results[0]["z16"]
    tbhdu.header["Z17"] = test_results[0]["z17"]
    tbhdu.header["Z18"] = test_results[0]["z18"]
    tbhdu.header["Z19"] = test_results[0]["z19"]
    tbhdu.header["Z20"] = test_results[0]["z20"]
    tbhdu.header["Z21"] = test_results[0]["z21"]
    tbhdu.header["SPHERE_5"] = test_results[0]["spherical_5th"]
    tbhdu.header["KERN_ADJ"] = test_results[0]["kernel_adjustment"]
    tbhdu.header["X_SQR"] = test_results[1][0]
    tbhdu.header["CHI_SQR"] = test_results[1][2]
    tbhdu.header["DOF"] = test_results[1][3]
    tbhdu.header["NSTAR"] = test_results[1][1]

    tbhdu.header["M0D_DIF"] = test_results[2][0][0]
    tbhdu.header["M0D_Z2"] = test_results[3][0][0]

    for Q_label, j in zip(("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"), range(8)):
        tbhdu.header[Q_label + "_DIF"] = test_results[2][0][1][j]
        tbhdu.header[Q_label + "_Z2"] = test_results[3][0][1][j]

    for Q_label, j in zip(("QXC", "QYC", "QPC", "QCC", "QSC", "QXW", "QYW", "QPW", "QCW", "QSW",),
                          range(10)):
        tbhdu.header[Q_label + "_DIF"] = np.mean(test_results[11][j])

    tbhdu.header["M0DNDIF"] = test_results[2][1][0]
    tbhdu.header["M0DNZ2"] = test_results[3][1][0]

    for Q_label, j in zip(("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"), range(8)):
        tbhdu.header[Q_label + "NDIF"] = test_results[2][1][1][j]
        tbhdu.header[Q_label + "NZ2"] = test_results[3][1][1][j]


    results_filename = filename_root + mv.results_tail

    tbhdu.writeto(results_filename,clobber=True)
    
    # Print summary
    logger = get_default_logger()
    logger.info(str(tbhdu.header["NSTAR"]) + " stars tested.")
    log_string = ("X^2 = " + str(tbhdu.header["X_SQR"]) + " for " + str(tbhdu.header["DOF"]) +
          " degrees of freedom.")
    for param in test_results[0]:
        log_string += "\n" + param + " = " + str(test_results[0][param])
    logger.info(log_string)
    logger.info("chi^2  = " + str(tbhdu.header["CHI_SQR"]) + ", for " + str(tbhdu.header["DOF"]) +
          " degrees of freedom.")
    
    if fitting_record is not None:
        if len(fitting_record)>0:
            save_fitting_record(fitting_record,filename_root)

    return
